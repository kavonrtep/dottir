//! Pixel max-merge and scaling helpers (spec §4.1.5).
//!
//! The pixelmap is a `width × height` row-major buffer; each pixel
//! covers a `zoom × zoom` block of the underlying score matrix.
//! Multiple diagonals fall into the same pixel and the maximum value
//! wins.
//!
//! ## Concurrency model (Phase A1)
//!
//! `PixelMap` stores its data as `Vec<AtomicU8>` so [`max_merge`] takes
//! `&self`. The parallel driver in [`crate::plot`] shares one
//! [`PixelMap`] across all rayon workers — there is no per-chunk
//! pixelmap allocation, so peak memory equals just the output map plus
//! each worker's small ping-pong sum buffers. This honours the
//! `memory_limit_bytes` cap (spec §4.5.2) regardless of thread count.
//!
//! Determinism is preserved: integer `max` is associative and
//! commutative, and the `compare_exchange_weak` loop guarantees every
//! `max_merge` ultimately leaves the cell at `max(previous, value)`
//! independent of interleaving (spec §4.1.11).

use std::sync::atomic::{AtomicU8, Ordering};

use crate::error::DottirError;

/// A row-major byte buffer storing pixelmap values.
///
/// Internally `Vec<AtomicU8>` so that [`Self::max_merge`] can be called
/// concurrently through `&self`. After parallel compute completes,
/// `&mut self`-taking methods view the same storage as plain `&mut [u8]`
/// via [`AtomicU8::get_mut_slice`].
#[derive(Debug)]
pub struct PixelMap {
    width: u32,
    height: u32,
    data: Vec<AtomicU8>,
}

impl PixelMap {
    /// Allocate a new zero-filled pixelmap. Checked against
    /// `memory_limit_bytes` — refuses to allocate larger and returns
    /// [`DottirError::OutOfMemory`] (spec §4.5.2).
    ///
    /// After Phase A1 the cap is honest: the parallel driver shares one
    /// pixelmap across all workers, so the actual peak is this size plus
    /// `n_threads × (2 × qlen × 4)` for the ping-pong sum buffers, which
    /// is negligible.
    pub fn new_checked(
        width: u32,
        height: u32,
        memory_limit_bytes: u64,
    ) -> Result<Self, DottirError> {
        let n = width as u64 * height as u64;
        if n > memory_limit_bytes {
            return Err(DottirError::OutOfMemory {
                requested: n,
                limit: memory_limit_bytes,
            });
        }
        let mut data = Vec::with_capacity(n as usize);
        data.resize_with(n as usize, || AtomicU8::new(0));
        Ok(Self { width, height, data })
    }

    #[inline]
    pub fn width(&self) -> usize {
        self.width as usize
    }
    #[inline]
    pub fn height(&self) -> usize {
        self.height as usize
    }
    #[inline]
    pub fn width_u32(&self) -> u32 {
        self.width
    }
    #[inline]
    pub fn height_u32(&self) -> u32 {
        self.height
    }

    /// Consume the pixelmap and return its raw bytes. Useful at the end
    /// of compute when ownership moves into [`crate::plot::DotPlot`]
    /// — all post-processing (self-comparison mirror, forward/reverse
    /// combine) runs on `Vec<u8>` rather than through the atomic API.
    pub fn into_vec(self) -> Vec<u8> {
        self.data.into_iter().map(|a| a.into_inner()).collect()
    }

    /// Atomic max-merge: `pixel[y * width + x] = max(pixel[y * width + x],
    /// value)`. Safe to call concurrently through `&self`. Out-of-bounds
    /// coordinates are silently dropped; the caller is expected to have
    /// gated on bounds already.
    #[inline]
    pub fn max_merge(&self, x: usize, y: usize, value: u8) {
        if x >= self.width() || y >= self.height() {
            return;
        }
        let idx = y * self.width() + x;
        let cell = &self.data[idx];
        // `Relaxed` is sufficient: the kernel doesn't need to observe any
        // memory ordering between max_merge calls — every pixel is
        // independent of every other pixel, and the post-compute reader
        // synchronises via Arc::try_unwrap / &mut self anyway.
        let mut current = cell.load(Ordering::Relaxed);
        while value > current {
            match cell.compare_exchange_weak(
                current,
                value,
                Ordering::Relaxed,
                Ordering::Relaxed,
            ) {
                Ok(_) => break,
                Err(observed) => current = observed,
            }
        }
    }

    /// Snapshot the buffer to an owned `Vec<u8>` without consuming the
    /// pixelmap. Used by tests for inspection.
    pub fn to_vec(&self) -> Vec<u8> {
        self.data.iter().map(|a| a.load(Ordering::Relaxed)).collect()
    }
}

/// Element-wise max-merge `src` into `dst` on raw byte slices. Used by
/// [`crate::plot`] to combine forward and reverse strand passes after
/// they've been consumed out of their [`PixelMap`]s — see spec §4.1.7
/// and the Phase A1 notes in this module.
pub fn merge_max_into(dst: &mut [u8], src: &[u8]) {
    debug_assert_eq!(dst.len(), src.len());
    for (a, &b) in dst.iter_mut().zip(src.iter()) {
        if b > *a {
            *a = b;
        }
    }
}

/// Compute the output image dimension for a sequence of length `seq_len`
/// at zoom factor `zoom`. Matches `getImageDimension` (dotplot.c:724):
/// `ceil(seq_len / zoom)`.
pub fn image_dimension(seq_len: usize, zoom: u32) -> u32 {
    let zoom = zoom.max(1) as usize;
    ((seq_len + zoom - 1) / zoom) as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn max_merge_keeps_max() {
        let p = PixelMap::new_checked(3, 2, 1024).unwrap();
        p.max_merge(0, 0, 100);
        p.max_merge(0, 0, 50); // should be ignored
        p.max_merge(0, 0, 200);
        assert_eq!(p.to_vec()[0], 200);
    }

    #[test]
    fn out_of_bounds_max_merge_is_noop() {
        let p = PixelMap::new_checked(2, 2, 1024).unwrap();
        p.max_merge(5, 5, 100); // should not panic
        assert!(p.to_vec().iter().all(|&v| v == 0));
    }

    #[test]
    fn memory_limit_refuses_oversize() {
        let err = PixelMap::new_checked(1000, 1000, 999_999).unwrap_err();
        match err {
            DottirError::OutOfMemory { requested, limit } => {
                assert_eq!(requested, 1_000_000);
                assert_eq!(limit, 999_999);
            }
            _ => panic!("unexpected error: {err:?}"),
        }
    }

    #[test]
    fn merge_max_into_slice() {
        let mut a = vec![10, 20, 30, 0];
        let b = [80, 0, 40, 25];
        merge_max_into(&mut a, &b);
        assert_eq!(a, vec![80, 20, 40, 25]);
    }

    /// Concurrent writers calling max_merge through `&self` reach the
    /// same final state as a serial run.
    #[test]
    fn concurrent_max_merge_is_deterministic() {
        use std::sync::Arc;
        let pm = Arc::new(PixelMap::new_checked(64, 64, 1 << 20).unwrap());
        let n_threads = 8;
        let n_writes_per_thread = 5000;
        let mut handles = Vec::new();
        for t in 0..n_threads {
            let pm = Arc::clone(&pm);
            handles.push(std::thread::spawn(move || {
                // Each thread writes a different pattern; the final
                // value at each pixel must equal the max across all
                // writers' contributions.
                for i in 0..n_writes_per_thread {
                    let x = (i * 7 + t * 3) % 64;
                    let y = (i * 11 + t) % 64;
                    let v = ((i ^ t) as u8).wrapping_mul(3);
                    pm.max_merge(x, y, v);
                }
            }));
        }
        for h in handles {
            h.join().unwrap();
        }
        // Re-derive the expected value sequentially.
        let mut expected = vec![0_u8; 64 * 64];
        for t in 0..n_threads {
            for i in 0..n_writes_per_thread {
                let x = (i * 7 + t * 3) % 64;
                let y = (i * 11 + t) % 64;
                let v = ((i ^ t) as u8).wrapping_mul(3);
                let idx = y * 64 + x;
                if v > expected[idx] {
                    expected[idx] = v;
                }
            }
        }
        let pm = Arc::try_unwrap(pm).expect("only one Arc remaining");
        assert_eq!(pm.to_vec(), expected);
    }

    #[test]
    fn image_dimension_ceil_div() {
        assert_eq!(image_dimension(100, 1), 100);
        assert_eq!(image_dimension(100, 10), 10);
        assert_eq!(image_dimension(101, 10), 11);
        assert_eq!(image_dimension(0, 10), 0);
        assert_eq!(image_dimension(1, 1), 1);
        // Defensive: zoom=0 is treated as zoom=1.
        assert_eq!(image_dimension(5, 0), 5);
    }
}
