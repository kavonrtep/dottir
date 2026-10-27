//! Direct-address spaced-seed index of the subject sequence.
//!
//! `offsets[key]..offsets[key + 1]` delimits the bucket of subject window
//! starts with that key. Buckets are filled in ascending subject order, so
//! each bucket is sorted and can be restricted to a subject range by
//! binary search (viewport recompute).

use std::ops::Range;

use super::seed::{for_each_window, SpacedSeed};
use crate::error::DottirError;

#[derive(Debug, Clone)]
pub struct SeedIndex {
    seed: SpacedSeed,
    offsets: Vec<u32>,
    positions: Vec<u32>,
    subject_len: usize,
}

impl SeedIndex {
    /// Bytes [`Self::build`] will allocate for a subject of `len` bases.
    /// An upper bound: every position is assumed to start a valid window.
    pub fn estimate_bytes(seed: &SpacedSeed, len: usize) -> u64 {
        (seed.n_keys() as u64 + 1) * 4 + len as u64 * 4
    }

    /// Index an encoded subject (codes from [`crate::alphabet::encode_dna`]).
    /// `breaks` are record start offsets (excluding 0), sorted.
    pub fn build(
        seed: &SpacedSeed,
        subject_encoded: &[u8],
        breaks: &[usize],
    ) -> Result<Self, DottirError> {
        if subject_encoded.len() >= u32::MAX as usize {
            return Err(DottirError::InvalidConfig(format!(
                "subject of {} bases exceeds the 4 Gb large-mode limit",
                subject_encoded.len()
            )));
        }
        let n_keys = seed.n_keys();
        let mut offsets = vec![0u32; n_keys + 1];
        let all = 0..subject_encoded.len();
        for_each_window(seed, subject_encoded, breaks, all.clone(), |_, fwd, _| {
            offsets[seed.key(fwd) as usize + 1] += 1;
        });
        for i in 1..=n_keys {
            offsets[i] += offsets[i - 1];
        }
        let total = offsets[n_keys] as usize;
        let mut positions = vec![0u32; total];
        // `cursor[key]` starts at the bucket start; reuse a copy of offsets.
        let mut cursor = offsets[..n_keys].to_vec();
        for_each_window(seed, subject_encoded, breaks, all, |p, fwd, _| {
            let k = seed.key(fwd) as usize;
            positions[cursor[k] as usize] = p as u32;
            cursor[k] += 1;
        });
        Ok(Self {
            seed: seed.clone(),
            offsets,
            positions,
            subject_len: subject_encoded.len(),
        })
    }

    pub fn seed(&self) -> &SpacedSeed {
        &self.seed
    }
    pub fn subject_len(&self) -> usize {
        self.subject_len
    }
    /// Number of indexed windows.
    pub fn n_positions(&self) -> usize {
        self.positions.len()
    }
    pub fn bytes(&self) -> u64 {
        (self.offsets.len() as u64 + self.positions.len() as u64) * 4
    }

    /// Total occurrences of `key` in the whole subject.
    #[inline]
    pub fn count(&self, key: u32) -> u32 {
        let k = key as usize;
        self.offsets[k + 1] - self.offsets[k]
    }

    /// Ascending subject starts with `key`.
    #[inline]
    pub fn bucket(&self, key: u32) -> &[u32] {
        let k = key as usize;
        &self.positions[self.offsets[k] as usize..self.offsets[k + 1] as usize]
    }

    /// Subject starts with `key` that lie in `range`.
    #[inline]
    pub fn bucket_in(&self, key: u32, range: &Range<usize>) -> &[u32] {
        let b = self.bucket(key);
        if range.start == 0 && range.end >= self.subject_len {
            return b;
        }
        let lo = b.partition_point(|&p| (p as usize) < range.start);
        let hi = b.partition_point(|&p| (p as usize) < range.end);
        &b[lo..hi]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::{encode, AlphabetKind};

    #[test]
    fn buckets_list_every_occurrence_in_order() {
        let seed = SpacedSeed::contiguous(6).unwrap();
        let seq = b"ACGTACGTACGTNNACGTAC";
        let enc = encode(seq, AlphabetKind::Dna);
        let idx = SeedIndex::build(&seed, &enc, &[]).unwrap();
        let mut n = 0;
        for_each_window(&seed, &enc, &[], 0..enc.len(), |p, fwd, _| {
            let b = idx.bucket(seed.key(fwd));
            assert!(b.contains(&(p as u32)));
            assert!(b.windows(2).all(|w| w[0] < w[1]));
            n += 1;
        });
        assert_eq!(idx.n_positions(), n);
        let k = {
            let mut k = 0;
            for_each_window(&seed, &enc, &[], 0..1, |_, fwd, _| k = seed.key(fwd));
            k
        };
        // ACGTAC occurs at 0, 4 and 14.
        assert_eq!(idx.bucket(k), &[0, 4, 14]);
        assert_eq!(idx.bucket_in(k, &(1..15)), &[4, 14]);
        assert_eq!(idx.bucket_in(k, &(5..14)), &[] as &[u32]);
        assert_eq!(idx.count(k), 3);
    }
}
