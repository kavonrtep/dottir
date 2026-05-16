//! Pixel max-merge and scaling helpers (spec §4.1.5).
//!
//! The pixelmap is a `width × height` `Vec<u8>` row-major; each pixel
//! covers a `zoom × zoom` block of the underlying score matrix. Multiple
//! diagonals fall into the same pixel and the maximum value wins.

use crate::error::DottirError;

/// A row-major `Vec<u8>` pixel buffer with explicit width/height. Keeps
/// the per-call indexing in [`crate::sliding`] readable.
#[derive(Debug, Clone)]
pub struct PixelMap {
    width: u32,
    height: u32,
    data: Vec<u8>,
}

impl PixelMap {
    /// Allocate a new zero-filled pixelmap. Checked against
    /// `memory_limit_bytes` — refuses to allocate larger and returns
    /// [`DottirError::OutOfMemory`] (spec §4.5.2).
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
        Ok(Self {
            width,
            height,
            data: vec![0_u8; n as usize],
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

    /// Consume the pixelmap and return its raw bytes.
    pub fn into_vec(self) -> Vec<u8> {
        self.data
    }

    /// Borrow the raw bytes.
    pub fn as_slice(&self) -> &[u8] {
        &self.data
    }

    /// `pixel[y * width + x] = max(pixel[y * width + x], value)`. Out-of-
    /// bounds coordinates are silently dropped; the caller is expected to
    /// have gated on bounds already.
    #[inline]
    pub fn max_merge(&mut self, x: usize, y: usize, value: u8) {
        if x >= self.width() || y >= self.height() {
            return;
        }
        let idx = y * self.width() + x;
        if value > self.data[idx] {
            self.data[idx] = value;
        }
    }

    /// Element-wise max-merge another pixelmap of the same dimensions
    /// into this one. Used to combine forward and reverse strand passes
    /// (spec §4.1.7) and to merge parallel-chunked outputs (Phase 3).
    pub fn merge_from(&mut self, other: &PixelMap) {
        assert_eq!(self.width, other.width);
        assert_eq!(self.height, other.height);
        for (a, b) in self.data.iter_mut().zip(other.data.iter()) {
            if *b > *a {
                *a = *b;
            }
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
        let mut p = PixelMap::new_checked(3, 2, 1024).unwrap();
        p.max_merge(0, 0, 100);
        p.max_merge(0, 0, 50); // should be ignored
        p.max_merge(0, 0, 200);
        assert_eq!(p.as_slice()[0], 200);
    }

    #[test]
    fn out_of_bounds_max_merge_is_noop() {
        let mut p = PixelMap::new_checked(2, 2, 1024).unwrap();
        p.max_merge(5, 5, 100); // should not panic
        assert!(p.as_slice().iter().all(|&v| v == 0));
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
    fn merge_from_takes_max_elementwise() {
        let mut a = PixelMap::new_checked(2, 2, 1024).unwrap();
        let mut b = PixelMap::new_checked(2, 2, 1024).unwrap();
        a.max_merge(0, 0, 10);
        a.max_merge(1, 1, 50);
        b.max_merge(0, 0, 80);
        b.max_merge(1, 0, 20);
        a.merge_from(&b);
        assert_eq!(a.as_slice(), &[80, 20, 0, 50]);
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
