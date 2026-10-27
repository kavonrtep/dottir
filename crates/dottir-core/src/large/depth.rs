//! Per-base repeat depth from canonical k-mer copy number, and the repeat
//! classes derived from it.
//!
//! A base's depth is the maximum copy number (within its own sequence) of
//! any canonical k-mer covering it. Counting is exact and bounded in memory:
//! k-mers are hash-partitioned and each partition is sorted and run-length
//! counted on its own.

use super::seed::SpacedSeed;
use crate::error::DottirError;

/// Repeat class of a base or a dot, ordered from least to most repetitive.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RepeatClass {
    Unique = 0,
    Low = 1,
    Mid = 2,
    High = 3,
}

impl RepeatClass {
    pub const ALL: [RepeatClass; 4] = [Self::Unique, Self::Low, Self::Mid, Self::High];
    pub const COUNT: usize = 4;

    pub fn index(self) -> usize {
        self as usize
    }
    pub fn from_index(i: usize) -> Option<Self> {
        Self::ALL.get(i).copied()
    }
    pub fn name(self) -> &'static str {
        match self {
            Self::Unique => "unique",
            Self::Low => "low",
            Self::Mid => "mid",
            Self::High => "high",
        }
    }
    pub fn from_name(s: &str) -> Option<Self> {
        Self::ALL
            .into_iter()
            .find(|c| c.name().eq_ignore_ascii_case(s))
    }
}

/// Copy-number thresholds: `depth < low` → unique, `< mid` → low,
/// `< high` → mid, otherwise high.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ClassThresholds {
    pub low: u16,
    pub mid: u16,
    pub high: u16,
}

impl Default for ClassThresholds {
    fn default() -> Self {
        Self {
            low: 2,
            mid: 5,
            high: 51,
        }
    }
}

impl ClassThresholds {
    pub fn validate(&self) -> Result<(), DottirError> {
        if 1 <= self.low && self.low < self.mid && self.mid < self.high {
            Ok(())
        } else {
            Err(DottirError::InvalidConfig(format!(
                "class thresholds must satisfy 1 <= low < mid < high, got {}/{}/{}",
                self.low, self.mid, self.high
            )))
        }
    }

    #[inline]
    pub fn classify(&self, depth: u16) -> RepeatClass {
        if depth < self.low {
            RepeatClass::Unique
        } else if depth < self.mid {
            RepeatClass::Low
        } else if depth < self.high {
            RepeatClass::Mid
        } else {
            RepeatClass::High
        }
    }

    /// Smallest depth that classifies as `class`.
    pub fn floor_of(&self, class: RepeatClass) -> u16 {
        match class {
            RepeatClass::Unique => 0,
            RepeatClass::Low => self.low,
            RepeatClass::Mid => self.mid,
            RepeatClass::High => self.high,
        }
    }
}

/// How per-base copy number is estimated.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DepthMethod {
    /// No estimate: every base is unique unless a floor applies.
    None,
    /// Copy number of exact canonical k-mers. Fast; blind to repeats
    /// diverged beyond roughly `100/k` %.
    Kmer { k: u8 },
    /// Number of two-hit diagonal runs covering the base in a
    /// self-comparison, plus one; seeds above `max_occ` count as high-copy.
    /// As divergence-tolerant as the plot; costs one extra seeding pass
    /// per sequence.
    Seed,
}

impl DepthMethod {
    pub fn describe(self) -> String {
        match self {
            Self::None => "none".into(),
            Self::Kmer { k } => format!("kmer:{k}"),
            Self::Seed => "seed".into(),
        }
    }
}

/// Per-base depth of one (concatenated, multi-record) sequence.
#[derive(Debug, Clone)]
pub struct DepthTrack {
    pub method: DepthMethod,
    /// Saturating copy number per base; 0 or 1 where nothing repeats.
    pub depth: Vec<u16>,
}

impl DepthTrack {
    /// A track of all-zero depth (everything classifies as unique). Used
    /// when depth classification is disabled.
    pub fn flat(len: usize) -> Self {
        Self {
            method: DepthMethod::None,
            depth: vec![0; len],
        }
    }

    /// Raise every base to at least `floor[i]` (e.g. soft-mask or
    /// annotation-derived repeat floors).
    pub fn apply_floor(&mut self, floor: &[u16]) -> Result<(), DottirError> {
        if floor.len() != self.depth.len() {
            return Err(DottirError::InvalidConfig(format!(
                "depth floor has {} entries, sequence has {}",
                floor.len(),
                self.depth.len()
            )));
        }
        for (d, &f) in self.depth.iter_mut().zip(floor) {
            *d = (*d).max(f);
        }
        Ok(())
    }
}

/// Roughly how many `(kmer, pos)` pairs one counting partition holds.
const PARTITION_TARGET: usize = 1 << 24;

pub const MIN_DEPTH_K: u8 = 11;
pub const MAX_DEPTH_K: u8 = 31;

/// Compute the depth track of an encoded sequence. `breaks` are record
/// start offsets (excluding 0); k-mers never span a record boundary.
pub fn compute_depth(encoded: &[u8], breaks: &[usize], k: u8) -> Result<DepthTrack, DottirError> {
    if !(MIN_DEPTH_K..=MAX_DEPTH_K).contains(&k) {
        return Err(DottirError::InvalidConfig(format!(
            "depth k must be {MIN_DEPTH_K}..={MAX_DEPTH_K}, got {k}"
        )));
    }
    if encoded.len() >= u32::MAX as usize {
        return Err(DottirError::InvalidConfig(
            "sequence exceeds the 4 Gb large-mode limit".into(),
        ));
    }
    let n = encoded.len();
    let n_parts = n.div_ceil(PARTITION_TARGET).next_power_of_two().max(1);
    let part_bits = n_parts.trailing_zeros();

    // `counts[i]` = copy number of the k-mer starting at i (0 = none).
    let mut counts = vec![0u16; n];
    for part in 0..n_parts as u64 {
        let mut pairs = collect_partition(encoded, breaks, k, part_bits, part);
        sort_pairs(&mut pairs);
        let mut i = 0;
        while i < pairs.len() {
            let key = pairs[i].0;
            let mut j = i + 1;
            while j < pairs.len() && pairs[j].0 == key {
                j += 1;
            }
            let c = (j - i).min(u16::MAX as usize) as u16;
            for &(_, pos) in &pairs[i..j] {
                counts[pos as usize] = c;
            }
            i = j;
        }
    }

    // Sliding max over the k k-mers covering each base, in place: the
    // deque holds (start, count) so overwriting `counts[i]` is safe.
    let k_us = k as usize;
    let mut deque: std::collections::VecDeque<(usize, u16)> = std::collections::VecDeque::new();
    for (i, slot) in counts.iter_mut().enumerate() {
        let c = *slot;
        if c > 0 {
            while deque.back().is_some_and(|&(_, v)| v <= c) {
                deque.pop_back();
            }
            deque.push_back((i, c));
        }
        while deque.front().is_some_and(|&(s, _)| s + k_us <= i) {
            deque.pop_front();
        }
        *slot = deque.front().map_or(0, |&(_, v)| v);
    }
    Ok(DepthTrack {
        method: DepthMethod::Kmer { k },
        depth: counts,
    })
}

#[inline]
fn partition_of(kmer: u64, part_bits: u32) -> u64 {
    if part_bits == 0 {
        0
    } else {
        kmer.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - part_bits)
    }
}

fn collect_partition(
    encoded: &[u8],
    breaks: &[usize],
    k: u8,
    part_bits: u32,
    part: u64,
) -> Vec<(u64, u32)> {
    let seed = SpacedSeed::contiguous_unchecked(k);
    let scan = |range: std::ops::Range<usize>| {
        let mut out = Vec::new();
        super::seed::for_each_window(&seed, encoded, breaks, range, |p, fwd, rc| {
            let canon = fwd.min(rc);
            if partition_of(canon, part_bits) == part {
                out.push((canon, p as u32));
            }
        });
        out
    };
    #[cfg(feature = "rayon")]
    {
        use rayon::prelude::*;
        const CHUNK: usize = 1 << 20;
        let chunks: Vec<_> = (0..encoded.len())
            .step_by(CHUNK)
            .map(|lo| lo..(lo + CHUNK).min(encoded.len()))
            .collect();
        let parts: Vec<Vec<(u64, u32)>> = chunks.into_par_iter().map(scan).collect();
        parts.concat()
    }
    #[cfg(not(feature = "rayon"))]
    {
        scan(0..encoded.len())
    }
}

fn sort_pairs(pairs: &mut [(u64, u32)]) {
    #[cfg(feature = "rayon")]
    {
        use rayon::prelude::*;
        pairs.par_sort_unstable_by_key(|p| p.0);
    }
    #[cfg(not(feature = "rayon"))]
    pairs.sort_unstable_by_key(|p| p.0);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::{encode, reverse_complement_dna, AlphabetKind};
    use std::collections::BTreeMap;

    fn naive_depth(seq: &[u8], k: usize) -> Vec<u16> {
        let canon = |w: &[u8]| {
            let rc = reverse_complement_dna(w);
            if rc.as_slice() < w {
                rc
            } else {
                w.to_vec()
            }
        };
        let valid = |w: &[u8]| w.iter().all(|b| b"ACGT".contains(b));
        let mut counts: BTreeMap<Vec<u8>, u16> = BTreeMap::new();
        for w in seq.windows(k).filter(|w| valid(w)) {
            *counts.entry(canon(w)).or_default() += 1;
        }
        let mut d = vec![0u16; seq.len()];
        for (i, w) in seq.windows(k).enumerate() {
            if !valid(w) {
                continue;
            }
            let c = counts[&canon(w)];
            for x in &mut d[i..i + k] {
                *x = (*x).max(c);
            }
        }
        d
    }

    fn lcg_dna(n: usize, mut state: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(state >> 62) as usize]
            })
            .collect()
    }

    #[test]
    fn matches_naive_counter_with_planted_repeats() {
        let unit = lcg_dna(300, 7);
        let mut seq = lcg_dna(2000, 1);
        // three forward copies and one reverse-complement copy
        seq.splice(100..100, unit.iter().copied());
        seq.splice(900..900, unit.iter().copied());
        seq.splice(1500..1500, reverse_complement_dna(&unit));
        seq.splice(1550..1550, b"NNNN".iter().copied());
        let enc = encode(&seq, AlphabetKind::Dna);
        let got = compute_depth(&enc, &[], 15).unwrap();
        assert_eq!(got.depth, naive_depth(&seq, 15));
        assert!(got.depth[150] >= 3);
        assert_eq!(got.depth[50], 1);
    }

    #[test]
    fn kmers_do_not_cross_records() {
        let a = lcg_dna(40, 3);
        let mut seq = a.clone();
        seq.extend_from_slice(&a);
        let enc = encode(&seq, AlphabetKind::Dna);
        let d = compute_depth(&enc, &[40], 11).unwrap();
        assert!(d.depth.iter().all(|&x| x == 2));
    }

    #[test]
    fn classes() {
        let t = ClassThresholds::default();
        assert_eq!(t.classify(0), RepeatClass::Unique);
        assert_eq!(t.classify(1), RepeatClass::Unique);
        assert_eq!(t.classify(2), RepeatClass::Low);
        assert_eq!(t.classify(50), RepeatClass::Mid);
        assert_eq!(t.classify(51), RepeatClass::High);
        assert!(ClassThresholds {
            low: 3,
            mid: 3,
            high: 9
        }
        .validate()
        .is_err());
    }
}
