//! Direct-address spaced-seed index of the subject sequence.
//!
//! `offsets[key]..offsets[key + 1]` delimits the bucket of subject window
//! starts with that key. Buckets are filled in ascending subject order, so
//! each bucket is sorted and can be restricted to a subject range by
//! binary search (viewport recompute).

use std::ops::Range;
use std::sync::atomic::{AtomicU32, Ordering::Relaxed};

use super::seed::{for_each_window, SpacedSeed};
use crate::error::DottirError;

/// Window starts walked per task by [`SeedIndex::build`].
const BUILD_CHUNK: usize = if cfg!(test) { 1 << 10 } else { 1 << 20 };

/// Sort every bucket `positions[offsets[k]..offsets[k + 1]]`. Parallel
/// tasks take contiguous key ranges, so they own disjoint slices.
fn sort_buckets(offsets: &[u32], positions: &mut [u32]) {
    let n_keys = offsets.len() - 1;
    let sort_range = |keys: Range<usize>, slice: &mut [u32]| {
        let base = offsets[keys.start] as usize;
        for k in keys {
            let b = &mut slice[offsets[k] as usize - base..offsets[k + 1] as usize - base];
            if !b.is_sorted() {
                b.sort_unstable();
            }
        }
    };
    #[cfg(feature = "rayon")]
    {
        use rayon::prelude::*;
        const KEYS_PER_TASK: usize = 1 << 14;
        let mut parts = Vec::with_capacity(n_keys.div_ceil(KEYS_PER_TASK));
        let mut rest = positions;
        for k0 in (0..n_keys).step_by(KEYS_PER_TASK) {
            let k1 = (k0 + KEYS_PER_TASK).min(n_keys);
            let (head, tail) = rest.split_at_mut((offsets[k1] - offsets[k0]) as usize);
            parts.push((k0..k1, head));
            rest = tail;
        }
        parts
            .into_par_iter()
            .for_each(|(keys, slice)| sort_range(keys, slice));
    }
    #[cfg(not(feature = "rayon"))]
    sort_range(0..n_keys, positions);
}

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
        let n = subject_encoded.len();
        let n_keys = seed.n_keys();
        // Window starts are split into chunks walked independently; the
        // bucket counters are atomic so chunks can run on any thread.
        let chunks: Vec<Range<usize>> = (0..n.div_ceil(BUILD_CHUNK))
            .map(|c| c * BUILD_CHUNK..((c + 1) * BUILD_CHUNK).min(n))
            .collect();
        let walk = |f: &(dyn Fn(usize, usize) + Sync)| {
            let per_chunk = |r: &Range<usize>| {
                for_each_window(seed, subject_encoded, breaks, r.clone(), |p, fwd, _| {
                    f(p, seed.key(fwd) as usize)
                })
            };
            #[cfg(feature = "rayon")]
            {
                use rayon::prelude::*;
                chunks.par_iter().for_each(per_chunk);
            }
            #[cfg(not(feature = "rayon"))]
            chunks.iter().for_each(per_chunk);
        };

        let counts: Vec<AtomicU32> = (0..=n_keys).map(|_| AtomicU32::new(0)).collect();
        walk(&|_, k| {
            counts[k + 1].fetch_add(1, Relaxed);
        });
        let mut offsets: Vec<u32> = counts.into_iter().map(AtomicU32::into_inner).collect();
        for i in 1..=n_keys {
            offsets[i] += offsets[i - 1];
        }
        let total = offsets[n_keys] as usize;
        let cursor: Vec<AtomicU32> = offsets[..n_keys]
            .iter()
            .map(|&o| AtomicU32::new(o))
            .collect();
        let slots: Vec<AtomicU32> = (0..total).map(|_| AtomicU32::new(0)).collect();
        walk(&|p, k| {
            slots[cursor[k].fetch_add(1, Relaxed) as usize].store(p as u32, Relaxed);
        });
        drop(cursor);
        let mut positions: Vec<u32> = slots.into_iter().map(AtomicU32::into_inner).collect();
        // Chunks interleave within a bucket; sorting restores ascending
        // subject order, so the index is identical whatever the scheduling.
        sort_buckets(&offsets, &mut positions);
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

    /// The original serial build: buckets filled in subject order.
    fn serial_build(seed: &SpacedSeed, enc: &[u8], breaks: &[usize]) -> (Vec<u32>, Vec<u32>) {
        let n_keys = seed.n_keys();
        let mut offsets = vec![0u32; n_keys + 1];
        for_each_window(seed, enc, breaks, 0..enc.len(), |_, fwd, _| {
            offsets[seed.key(fwd) as usize + 1] += 1;
        });
        for i in 1..=n_keys {
            offsets[i] += offsets[i - 1];
        }
        let mut positions = vec![0u32; offsets[n_keys] as usize];
        let mut cursor = offsets[..n_keys].to_vec();
        for_each_window(seed, enc, breaks, 0..enc.len(), |p, fwd, _| {
            let k = seed.key(fwd) as usize;
            positions[cursor[k] as usize] = p as u32;
            cursor[k] += 1;
        });
        (offsets, positions)
    }

    #[test]
    fn chunked_build_matches_the_serial_one() {
        // Low-entropy sequence with Ns and record breaks, long enough for
        // many chunks, so buckets collect windows from several of them.
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut seq = Vec::new();
        for i in 0..40 * BUILD_CHUNK {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let b = if i % 997 < 3 {
                b'N'
            } else if state % 4 == 0 {
                b"ACGTTGCA"[i % 8]
            } else {
                b"ACGT"[(state >> 8) as usize % 4]
            };
            seq.push(b);
        }
        let enc = encode(&seq, AlphabetKind::Dna);
        let breaks = [5, BUILD_CHUNK - 3, BUILD_CHUNK, 7 * BUILD_CHUNK + 11];
        for seed in [SpacedSeed::contiguous(6).unwrap(), SpacedSeed::default()] {
            let idx = SeedIndex::build(&seed, &enc, &breaks).unwrap();
            let (offsets, positions) = serial_build(&seed, &enc, &breaks);
            assert_eq!(idx.offsets, offsets);
            assert_eq!(idx.positions, positions);
        }
    }
}
