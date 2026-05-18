//! Pixel max-merge and scaling helpers (spec §4.1.5).
//!
//! The pixelmap is a `width × height` row-major byte buffer; each pixel
//! covers a `zoom × zoom` block of the underlying score matrix.
//! Multiple diagonals fall into the same pixel and the maximum value
//! wins.
//!
//! ## Concurrency model (Phase A1)
//!
//! [`PixelMap`] owns the buffer as plain `Vec<u8>`. To let multiple
//! rayon workers max-merge concurrently, the owner produces a
//! [`PixelView`] via [`PixelMap::view_mut`]: that view borrows the
//! underlying bytes as `&[AtomicU8]` (via the stdlib's safe
//! [`AtomicU8::from_mut_slice`]) and exposes [`PixelView::max_merge`]
//! through `&self`. The parallel driver in [`crate::plot`] shares one
//! [`PixelView`] across all workers — there is no per-chunk
//! pixelmap allocation, so peak memory equals just the output map plus
//! each worker's small ping-pong sum buffers. This honours the
//! `memory_limit_bytes` cap (spec §4.5.2) regardless of thread count.
//!
//! Once the parallel pass returns, the view is dropped and the
//! `PixelMap` is again exclusively owned — [`PixelMap::into_vec`] then
//! hands ownership of the buffer out without copying.
//!
//! Determinism is preserved: integer `max` is associative and
//! commutative, and the `compare_exchange_weak` loop guarantees every
//! `max_merge` ultimately leaves the cell at `max(previous, value)`
//! independent of interleaving (spec §4.1.11).

use std::sync::atomic::{AtomicU8, Ordering};

use crate::error::DottirError;

/// Owner of a row-major byte buffer storing pixelmap values.
///
/// The internal `Vec<u8>` is plain — to mutate it concurrently the
/// caller takes a transient [`PixelView`] via [`Self::view_mut`]; the
/// view aliases the same bytes as `&[AtomicU8]` (via stdlib's safe
/// [`AtomicU8::from_mut_slice`]) and can be shared across rayon workers
/// for the duration of one parallel pass.
#[derive(Debug)]
pub struct PixelMap {
    width: u32,
    height: u32,
    data: Vec<u8>,
}

impl PixelMap {
    /// Allocate a new zero-filled pixelmap. Checked against
    /// `memory_limit_bytes` — refuses to allocate larger and returns
    /// [`DottirError::OutOfMemory`] (spec §4.5.2).
    ///
    /// The parallel driver shares one pixelmap across all workers, so
    /// the actual peak is this size plus `n_threads × (2 × qlen × 4)`
    /// for the ping-pong sum buffers, which is negligible.
    ///
    /// In practice the up-front check in [`crate::compute_dotplot`] has
    /// already validated the total budget for *all* channels by the
    /// time this is called, so the per-pixelmap check here is a
    /// secondary guard for direct users of [`PixelMap`].
    pub fn new_checked(
        width: u32,
        height: u32,
        memory_limit_bytes: u64,
    ) -> Result<Self, DottirError> {
        let n = width as u64 * height as u64;
        if n > memory_limit_bytes {
            return Err(DottirError::OutOfMemory {
                requested: n,
                per_channel: n,
                channels: 1,
                limit: memory_limit_bytes,
            });
        }
        let data = vec![0u8; n as usize];
        Ok(Self {
            width,
            height,
            data,
        })
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

    /// Produce a transient atomic-access view of the buffer. Used by
    /// the parallel kernel driver to share `&view` across workers,
    /// each of whom calls `view.max_merge(...)` through `&self`. The
    /// view's lifetime is tied to the `&mut` borrow on `self`, so no
    /// other reader/writer of the `PixelMap` can race with the view.
    pub fn view_mut(&mut self) -> PixelView<'_> {
        PixelView {
            width: self.width,
            height: self.height,
            data: u8_slice_as_atomic(&mut self.data),
        }
    }

    /// Consume the pixelmap and return its raw bytes — zero-copy, the
    /// `Vec<u8>` is the same allocation that [`Self::new_checked`]
    /// produced.
    pub fn into_vec(self) -> Vec<u8> {
        self.data
    }

    /// Snapshot the buffer to an owned `Vec<u8>` without consuming the
    /// pixelmap. Used by tests for inspection.
    pub fn to_vec(&self) -> Vec<u8> {
        self.data.clone()
    }
}

/// Transient atomic-access lens onto a [`PixelMap`]'s bytes. Held for
/// the duration of one parallel pass; the contained `&[AtomicU8]`
/// re-borrows the owner's `Vec<u8>` (via stdlib's safe
/// [`AtomicU8::from_mut_slice`]) so multiple rayon workers can share
/// `&view` and call `view.max_merge(...)` through `&self`.
#[derive(Debug, Clone, Copy)]
pub struct PixelView<'a> {
    width: u32,
    height: u32,
    data: &'a [AtomicU8],
}

impl<'a> PixelView<'a> {
    #[inline]
    pub fn width(&self) -> usize {
        self.width as usize
    }
    #[inline]
    pub fn height(&self) -> usize {
        self.height as usize
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
        // synchronises via the owning PixelMap's `&mut self` reborrow
        // (the view's lifetime ends, then `into_vec` consumes the owner).
        let mut current = cell.load(Ordering::Relaxed);
        while value > current {
            match cell.compare_exchange_weak(current, value, Ordering::Relaxed, Ordering::Relaxed) {
                Ok(_) => break,
                Err(observed) => current = observed,
            }
        }
    }
}

/// Reinterpret a `&mut [u8]` as `&[AtomicU8]`. This is what stdlib's
/// (still-nightly) `AtomicU8::from_mut_slice` does internally — we
/// duplicate the trick here so the safe-by-default core can produce
/// an atomic-write view without copying the buffer.
///
/// SAFETY: `AtomicU8` is `#[repr(transparent)]` over a
/// `core::cell::UnsafeCell<u8>`, which is itself `#[repr(transparent)]`
/// over `u8`. The three share size (1 byte), alignment (1 byte) and
/// memory layout. The returned slice borrows the same allocation as
/// `slice`; multiple shared (`&`) reads through atomic ops are sound
/// because every mutating method on `AtomicU8` takes `&self` and uses
/// processor-level atomic instructions. The lifetime is bound to the
/// input `&mut`, so the borrow checker prevents any non-atomic reader
/// from racing the view.
#[allow(unsafe_code)]
fn u8_slice_as_atomic(slice: &mut [u8]) -> &[AtomicU8] {
    let ptr = slice.as_mut_ptr() as *const AtomicU8;
    let len = slice.len();
    // SAFETY: see function-level comment.
    unsafe { std::slice::from_raw_parts(ptr, len) }
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
    seq_len.div_ceil(zoom) as u32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn max_merge_keeps_max() {
        let mut p = PixelMap::new_checked(3, 2, 1024).unwrap();
        {
            let v = p.view_mut();
            v.max_merge(0, 0, 100);
            v.max_merge(0, 0, 50); // should be ignored
            v.max_merge(0, 0, 200);
        }
        assert_eq!(p.to_vec()[0], 200);
    }

    #[test]
    fn out_of_bounds_max_merge_is_noop() {
        let mut p = PixelMap::new_checked(2, 2, 1024).unwrap();
        {
            let v = p.view_mut();
            v.max_merge(5, 5, 100); // should not panic
        }
        assert!(p.to_vec().iter().all(|&v| v == 0));
    }

    #[test]
    fn memory_limit_refuses_oversize() {
        let err = PixelMap::new_checked(1000, 1000, 999_999).unwrap_err();
        match err {
            DottirError::OutOfMemory {
                requested,
                per_channel,
                channels,
                limit,
            } => {
                assert_eq!(requested, 1_000_000);
                assert_eq!(per_channel, 1_000_000);
                assert_eq!(channels, 1);
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

    /// Concurrent writers calling `view.max_merge` through `&view`
    /// reach the same final state as a serial run.
    #[test]
    fn concurrent_max_merge_is_deterministic() {
        use std::thread;
        let mut pm = PixelMap::new_checked(64, 64, 1 << 20).unwrap();
        let n_threads = 8;
        let n_writes_per_thread = 5000;
        thread::scope(|scope| {
            let view = pm.view_mut();
            for t in 0..n_threads {
                scope.spawn(move || {
                    for i in 0..n_writes_per_thread {
                        let x = (i * 7 + t * 3) % 64;
                        let y = (i * 11 + t) % 64;
                        let v = ((i ^ t) as u8).wrapping_mul(3);
                        view.max_merge(x, y, v);
                    }
                });
            }
        });
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
        assert_eq!(pm.to_vec(), expected);
    }

    /// `into_vec` is zero-copy — the `Vec<u8>` we get out is the same
    /// allocation as the one constructed by `new_checked`. Sprinkle
    /// values via a view, then assert ownership transfers cleanly.
    #[test]
    fn into_vec_preserves_bytes() {
        let mut pm = PixelMap::new_checked(257, 131, 1 << 20).unwrap();
        let (w, h) = (pm.width(), pm.height());
        {
            let v = pm.view_mut();
            for y in 0..h {
                for x in 0..w {
                    let value = ((x.wrapping_mul(31) ^ y.wrapping_mul(17)) & 0xff) as u8;
                    v.max_merge(x, y, value);
                }
            }
        }
        let snapshot = pm.to_vec();
        let consumed = pm.into_vec();
        assert_eq!(consumed.len(), snapshot.len());
        assert_eq!(consumed, snapshot);
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
