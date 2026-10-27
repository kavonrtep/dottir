//! Spaced seeds over 2-bit packed DNA.
//!
//! A seed pattern is a string of `1` (position must match) and `0` (don't
//! care). A window of `span` bases is packed into a `u64` with the leftmost
//! base in the highest bits; the key is the gather of the `1` positions,
//! leftmost first, so keys of the same pattern compare equal iff the masked
//! bases are equal.

use crate::error::DottirError;

/// Maximum weight supported by the direct-address [`super::index::SeedIndex`]
/// (`4^13` buckets × 4 bytes = 256 MiB for the offset table).
pub const MAX_WEIGHT: u8 = 13;
/// Minimum weight; below this nearly every bucket is hit by chance.
pub const MIN_WEIGHT: u8 = 6;
/// Maximum span, bounded by the 64-bit packed window.
pub const MAX_SPAN: u8 = 32;

/// LASTZ's default 12-of-19 pattern (Harris 2007).
pub const DEFAULT_PATTERN: &str = "1110100110010101111";

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpacedSeed {
    pattern: String,
    span: u8,
    weight: u8,
    /// For each `1` in the pattern (leftmost first), the right-shift that
    /// brings that base to the low 2 bits of the packed window.
    shifts: Vec<u8>,
}

impl SpacedSeed {
    pub fn parse(pattern: &str) -> Result<Self, DottirError> {
        let bytes = pattern.as_bytes();
        if bytes.is_empty() || bytes.len() > MAX_SPAN as usize {
            return Err(DottirError::InvalidConfig(format!(
                "seed pattern must be 1..={MAX_SPAN} characters, got {}",
                bytes.len()
            )));
        }
        if let Some(bad) = bytes.iter().find(|&&b| b != b'0' && b != b'1') {
            return Err(DottirError::InvalidConfig(format!(
                "seed pattern may only contain '0' and '1', found {:?}",
                *bad as char
            )));
        }
        if bytes[0] != b'1' || bytes[bytes.len() - 1] != b'1' {
            return Err(DottirError::InvalidConfig(
                "seed pattern must start and end with '1'".into(),
            ));
        }
        let span = bytes.len() as u8;
        let shifts: Vec<u8> = bytes
            .iter()
            .enumerate()
            .filter(|(_, &b)| b == b'1')
            .map(|(j, _)| 2 * (span - 1 - j as u8))
            .collect();
        let weight = shifts.len() as u8;
        if !(MIN_WEIGHT..=MAX_WEIGHT).contains(&weight) {
            return Err(DottirError::InvalidConfig(format!(
                "seed weight must be {MIN_WEIGHT}..={MAX_WEIGHT}, pattern {pattern:?} has {weight}"
            )));
        }
        Ok(Self {
            pattern: pattern.to_string(),
            span,
            weight,
            shifts,
        })
    }

    /// A contiguous seed of length `k` (`"111…1"`).
    pub fn contiguous(k: u8) -> Result<Self, DottirError> {
        Self::parse(&"1".repeat(k as usize))
    }

    /// Contiguous k-mer "seed" without the index weight limits; used by
    /// the depth counter, which never builds a direct-address table.
    pub(crate) fn contiguous_unchecked(k: u8) -> Self {
        let span = k.clamp(1, MAX_SPAN);
        Self {
            pattern: "1".repeat(span as usize),
            span,
            weight: span,
            shifts: (0..span).map(|j| 2 * (span - 1 - j)).collect(),
        }
    }

    pub fn pattern(&self) -> &str {
        &self.pattern
    }
    pub fn span(&self) -> usize {
        self.span as usize
    }
    pub fn weight(&self) -> u8 {
        self.weight
    }
    /// Number of distinct keys, `4^weight`.
    pub fn n_keys(&self) -> usize {
        1usize << (2 * self.weight as usize)
    }

    /// Mask selecting the low `2 * span` bits of a packed window.
    #[inline]
    pub(crate) fn window_mask(&self) -> u64 {
        if self.span == 32 {
            u64::MAX
        } else {
            (1u64 << (2 * self.span as u32)) - 1
        }
    }

    /// Gather the `1` positions of a packed window into a key.
    #[inline]
    pub fn key(&self, packed: u64) -> u32 {
        let mut k = 0u32;
        for &sh in &self.shifts {
            k = (k << 2) | ((packed >> sh) & 3) as u32;
        }
        k
    }
}

impl Default for SpacedSeed {
    fn default() -> Self {
        Self::parse(DEFAULT_PATTERN).expect("default pattern is valid")
    }
}

/// Rolling packer over an encoded DNA stream (codes 0..=3 = ACGT; anything
/// else invalidates the window). Yields forward and reverse-complement
/// packed windows; the reverse window is the reverse complement of the
/// forward one, packed in the same leftmost-high convention.
#[derive(Debug, Clone)]
pub(crate) struct Roller {
    span: usize,
    mask: u64,
    rc_top: u32,
    fwd: u64,
    rc: u64,
    valid: usize,
}

impl Roller {
    pub(crate) fn new(seed: &SpacedSeed) -> Self {
        Self {
            span: seed.span(),
            mask: seed.window_mask(),
            rc_top: 2 * (seed.span() as u32 - 1),
            fwd: 0,
            rc: 0,
            valid: 0,
        }
    }

    #[inline]
    pub(crate) fn reset(&mut self) {
        self.valid = 0;
    }

    /// Push one encoded base. Returns `true` when the last `span` bases
    /// are all valid, i.e. [`Self::fwd`] / [`Self::rc`] hold a full window.
    #[inline]
    pub(crate) fn push(&mut self, code: u8) -> bool {
        if code > 3 {
            self.valid = 0;
            return false;
        }
        let c = code as u64;
        self.fwd = ((self.fwd << 2) | c) & self.mask;
        self.rc = (self.rc >> 2) | ((3 - c) << self.rc_top);
        self.valid += 1;
        self.valid >= self.span
    }

    #[inline]
    pub(crate) fn fwd(&self) -> u64 {
        self.fwd
    }
    #[inline]
    pub(crate) fn rc(&self) -> u64 {
        self.rc
    }
}

/// Call `f(start, fwd_packed, rc_packed)` for every window of `span`
/// valid bases whose start lies in `starts` and that does not cross a
/// record boundary. `breaks` is the sorted list of record start offsets
/// (excluding 0); windows spanning one are skipped.
pub(crate) fn for_each_window<F: FnMut(usize, u64, u64)>(
    seed: &SpacedSeed,
    encoded: &[u8],
    breaks: &[usize],
    starts: std::ops::Range<usize>,
    mut f: F,
) {
    let span = seed.span();
    let end = (starts.end + span - 1).min(encoded.len());
    if starts.start >= end {
        return;
    }
    let mut roller = Roller::new(seed);
    let mut next_break_idx = breaks.partition_point(|&b| b <= starts.start);
    for (e, &code) in encoded.iter().enumerate().take(end).skip(starts.start) {
        if next_break_idx < breaks.len() && breaks[next_break_idx] == e {
            roller.reset();
            next_break_idx += 1;
        }
        if roller.push(code) {
            f(e + 1 - span, roller.fwd(), roller.rc());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::{encode, reverse_complement_dna, AlphabetKind};

    fn naive_key(seed: &SpacedSeed, win: &[u8]) -> u32 {
        let enc = encode(win, AlphabetKind::Dna);
        let mut k = 0u32;
        for (j, b) in seed.pattern().bytes().enumerate() {
            if b == b'1' {
                k = (k << 2) | enc[j] as u32;
            }
        }
        k
    }

    #[test]
    fn default_is_12_of_19() {
        let s = SpacedSeed::default();
        assert_eq!(s.span(), 19);
        assert_eq!(s.weight(), 12);
    }

    #[test]
    fn rejects_bad_patterns() {
        assert!(SpacedSeed::parse("").is_err());
        assert!(SpacedSeed::parse("1102").is_err());
        assert!(SpacedSeed::parse("0111111111").is_err());
        assert!(SpacedSeed::parse("11111").is_err()); // weight 5
        assert!(SpacedSeed::parse(&"1".repeat(14)).is_err());
    }

    #[test]
    fn rolling_keys_match_naive_and_rc() {
        let seed = SpacedSeed::default();
        let seq = b"ACGTTGCAAGGCTTACGATCGATCGGATCCATGCATGCAAACCCGGGTTT";
        let enc = encode(seq, AlphabetKind::Dna);
        let mut n = 0;
        for_each_window(&seed, &enc, &[], 0..seq.len(), |p, fwd, rc| {
            let win = &seq[p..p + seed.span()];
            assert_eq!(seed.key(fwd), naive_key(&seed, win), "fwd at {p}");
            let rcw = reverse_complement_dna(win);
            assert_eq!(seed.key(rc), naive_key(&seed, &rcw), "rc at {p}");
            n += 1;
        });
        assert_eq!(n, seq.len() - seed.span() + 1);
    }

    #[test]
    fn windows_skip_n_and_breaks() {
        let seed = SpacedSeed::contiguous(6).unwrap();
        let enc = encode(b"ACGTACNACGTACGTACGTA", AlphabetKind::Dna);
        let mut starts = vec![];
        for_each_window(&seed, &enc, &[12], 0..enc.len(), |p, _, _| starts.push(p));
        // N at 6 kills windows 1..=6; break at 12 kills windows 7..=11.
        assert_eq!(starts, vec![0, 12, 13, 14]);
        // A sub-range only yields windows starting inside it.
        let mut sub = vec![];
        for_each_window(&seed, &enc, &[12], 5..13, |p, _, _| sub.push(p));
        assert_eq!(sub, vec![12]);
    }
}
