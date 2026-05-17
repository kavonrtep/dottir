//! Read and write the C-dotter `.dot` binary format (D2).
//!
//! ## File layout
//!
//! Bytes are little-endian (the original C is host-byte-order, with a
//! reverse-bytes shim for big-endian Alpha hosts that nobody runs any
//! more; little-endian is the modern reality).
//!
//! | Offset | Bytes | Field |
//! |--------|-------|-------|
//! | 0      | 1     | format version (1, 2, or 3) |
//! | 1      | 4 or 8| zoom: i32 (formats 1/2) or f64 (format 3) |
//! |        | 4     | image_width  (i32) |
//! |        | 4     | image_height (i32) |
//! | ↓ formats 2 & 3 only: |
//! |        | 4     | pixel_fac (i32) |
//! |        | 4     | sliding_window_size (i32) |
//! |        | 4     | matrix_name_len (i32) |
//! |        | N     | matrix_name (UTF-8 bytes, length matrix_name_len) |
//! |        | 4×576 | matrix (24×24 i32, row-major) |
//! | ↓ pixelmap: |
//! |        | w×h   | pixels (u8, row-major) |
//!
//! Format 1 is missing the matrix and parameters block; dottir-io
//! tolerates it by substituting C-dotter's default sliding_window_size = 25
//! and pixel_fac = 50 (matching the C reader's behaviour).
//!
//! Reference: `third_party/seqtools/dotterApp/dotplot.c:1610` (loadPlot)
//! and `:1812` (savePlot).

use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use thiserror::Error;

const MAX_MATRIX_NAME_LENGTH: usize = 256;

#[derive(Debug, Error)]
pub enum DotError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("unknown .dot format version: {0} (expected 1, 2, or 3)")]
    UnknownFormat(u8),
    #[error("matrix name length {0} exceeds limit {MAX_MATRIX_NAME_LENGTH}")]
    MatrixNameTooLong(u32),
    #[error("expected {expected} bytes of pixelmap, found {got}")]
    PixelmapSizeMismatch { expected: u64, got: u64 },
    #[error("invalid UTF-8 in matrix name")]
    InvalidMatrixName,
}

/// In-memory representation of a `.dot` file.
#[derive(Debug, Clone)]
pub struct DotFile {
    pub format: u8,
    pub zoom: f64,
    pub width: u32,
    pub height: u32,
    /// Pixel factor. `None` for format 1 (caller substitutes the
    /// default 50).
    pub pixel_fac: Option<i32>,
    /// Sliding window size. `None` for format 1 (caller substitutes
    /// the default 25).
    pub sliding_window_size: Option<i32>,
    /// Score matrix name, e.g. "BLOSUM62" or "DNA+5/-4". `None` for
    /// format 1.
    pub matrix_name: Option<String>,
    /// 24×24 score matrix in row-major order. `None` for format 1.
    pub matrix: Option<Vec<i32>>,
    /// Pixelmap: row-major u8 array of length `width * height`.
    pub pixels: Vec<u8>,
}

impl DotFile {
    /// Read a `.dot` file from disk. Auto-detects format 1/2/3.
    pub fn read<P: AsRef<Path>>(path: P) -> Result<Self, DotError> {
        let file = File::open(path)?;
        Self::from_reader(BufReader::new(file))
    }

    /// Read from any `Read` source.
    pub fn from_reader<R: Read>(mut r: R) -> Result<Self, DotError> {
        let format = read_u8(&mut r)?;
        if !(format == 1 || format == 2 || format == 3) {
            return Err(DotError::UnknownFormat(format));
        }
        let zoom = if format == 3 {
            read_f64(&mut r)?
        } else {
            read_i32(&mut r)? as f64
        };
        let width = read_i32(&mut r)? as u32;
        let height = read_i32(&mut r)? as u32;

        let (pixel_fac, sliding_window_size, matrix_name, matrix) = if format == 1 {
            (None, None, None, None)
        } else {
            let pf = read_i32(&mut r)?;
            let w = read_i32(&mut r)?;
            let name_len = read_i32(&mut r)? as u32;
            if (name_len as usize) > MAX_MATRIX_NAME_LENGTH {
                return Err(DotError::MatrixNameTooLong(name_len));
            }
            let mut name_buf = vec![0_u8; name_len as usize];
            r.read_exact(&mut name_buf)?;
            let name = String::from_utf8(name_buf).map_err(|_| DotError::InvalidMatrixName)?;
            let mut matrix = Vec::with_capacity(24 * 24);
            for _ in 0..(24 * 24) {
                matrix.push(read_i32(&mut r)?);
            }
            (Some(pf), Some(w), Some(name), Some(matrix))
        };

        // Trailing pixelmap.
        let mut pixels = vec![0_u8; (width as usize) * (height as usize)];
        r.read_exact(&mut pixels)?;

        Ok(DotFile {
            format,
            zoom,
            width,
            height,
            pixel_fac,
            sliding_window_size,
            matrix_name,
            matrix,
            pixels,
        })
    }

    /// Write as format-3 `.dot` (the current C-dotter format).
    /// `matrix_name` / `matrix` are required for format 3.
    pub fn write_format3<P: AsRef<Path>>(&self, path: P) -> Result<(), DotError> {
        let file = File::create(path)?;
        let mut w = BufWriter::new(file);
        self.write_format3_to(&mut w)
    }

    pub fn write_format3_to<W: Write>(&self, w: &mut W) -> Result<(), DotError> {
        if (self.width as usize) * (self.height as usize) != self.pixels.len() {
            return Err(DotError::PixelmapSizeMismatch {
                expected: (self.width as u64) * (self.height as u64),
                got: self.pixels.len() as u64,
            });
        }
        w.write_all(&[3_u8])?; // format
        w.write_all(&self.zoom.to_le_bytes())?;
        w.write_all(&(self.width as i32).to_le_bytes())?;
        w.write_all(&(self.height as i32).to_le_bytes())?;
        w.write_all(&self.pixel_fac.unwrap_or(50).to_le_bytes())?;
        w.write_all(&self.sliding_window_size.unwrap_or(25).to_le_bytes())?;
        let name = self.matrix_name.as_deref().unwrap_or("UNKNOWN");
        let name_bytes = name.as_bytes();
        if name_bytes.len() > MAX_MATRIX_NAME_LENGTH {
            return Err(DotError::MatrixNameTooLong(name_bytes.len() as u32));
        }
        w.write_all(&(name_bytes.len() as i32).to_le_bytes())?;
        w.write_all(name_bytes)?;
        let matrix = self.matrix.clone().unwrap_or_else(|| vec![0_i32; 24 * 24]);
        for v in matrix.iter().take(24 * 24) {
            w.write_all(&v.to_le_bytes())?;
        }
        // Pad if caller supplied a smaller matrix.
        if matrix.len() < 24 * 24 {
            for _ in matrix.len()..(24 * 24) {
                w.write_all(&0_i32.to_le_bytes())?;
            }
        }
        w.write_all(&self.pixels)?;
        Ok(())
    }
}

#[inline]
fn read_u8<R: Read>(r: &mut R) -> Result<u8, DotError> {
    let mut b = [0_u8; 1];
    r.read_exact(&mut b)?;
    Ok(b[0])
}

#[inline]
fn read_i32<R: Read>(r: &mut R) -> Result<i32, DotError> {
    let mut b = [0_u8; 4];
    r.read_exact(&mut b)?;
    Ok(i32::from_le_bytes(b))
}

#[inline]
fn read_f64<R: Read>(r: &mut R) -> Result<f64, DotError> {
    let mut b = [0_u8; 8];
    r.read_exact(&mut b)?;
    Ok(f64::from_le_bytes(b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn fixture() -> DotFile {
        DotFile {
            format: 3,
            zoom: 1.5,
            width: 4,
            height: 3,
            pixel_fac: Some(60),
            sliding_window_size: Some(20),
            matrix_name: Some("TOY".into()),
            matrix: Some((0..(24 * 24)).collect()),
            pixels: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        }
    }

    #[test]
    fn format3_round_trip() {
        let f = fixture();
        let mut buf = Vec::new();
        f.write_format3_to(&mut buf).unwrap();
        let g = DotFile::from_reader(Cursor::new(&buf)).unwrap();
        assert_eq!(g.format, 3);
        assert_eq!(g.zoom, 1.5);
        assert_eq!(g.width, 4);
        assert_eq!(g.height, 3);
        assert_eq!(g.pixel_fac, Some(60));
        assert_eq!(g.sliding_window_size, Some(20));
        assert_eq!(g.matrix_name.as_deref(), Some("TOY"));
        assert_eq!(g.pixels, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]);
    }

    /// Read a hand-built format-1 stream (just version + i32 zoom +
    /// dims + pixels). Should parse with matrix fields = None.
    #[test]
    fn reads_format_1() {
        let mut buf = Vec::new();
        buf.push(1_u8); // format
        buf.extend_from_slice(&2_i32.to_le_bytes()); // zoom = 2 (int)
        buf.extend_from_slice(&3_i32.to_le_bytes()); // width
        buf.extend_from_slice(&2_i32.to_le_bytes()); // height
        buf.extend_from_slice(&[10, 20, 30, 40, 50, 60]); // pixels
        let f = DotFile::from_reader(Cursor::new(&buf)).unwrap();
        assert_eq!(f.format, 1);
        assert_eq!(f.zoom, 2.0);
        assert_eq!(f.width, 3);
        assert_eq!(f.height, 2);
        assert!(f.pixel_fac.is_none());
        assert!(f.matrix.is_none());
        assert_eq!(f.pixels, vec![10, 20, 30, 40, 50, 60]);
    }

    #[test]
    fn rejects_unknown_format() {
        let buf = [99_u8, 0, 0, 0, 0];
        let err = DotFile::from_reader(Cursor::new(&buf[..])).unwrap_err();
        assert!(matches!(err, DotError::UnknownFormat(99)));
    }

    #[test]
    fn rejects_oversized_matrix_name() {
        let mut buf = Vec::new();
        buf.push(2_u8); // format
        buf.extend_from_slice(&1_i32.to_le_bytes());
        buf.extend_from_slice(&1_i32.to_le_bytes());
        buf.extend_from_slice(&1_i32.to_le_bytes());
        buf.extend_from_slice(&50_i32.to_le_bytes()); // pixel_fac
        buf.extend_from_slice(&25_i32.to_le_bytes()); // win size
        buf.extend_from_slice(&(MAX_MATRIX_NAME_LENGTH as i32 + 1).to_le_bytes());
        let err = DotFile::from_reader(Cursor::new(&buf)).unwrap_err();
        assert!(matches!(err, DotError::MatrixNameTooLong(_)));
    }

    #[test]
    fn write_default_matrix_pads_zero() {
        // Caller omits the matrix; writer pads with zeros.
        let mut f = fixture();
        f.matrix = None;
        let mut buf = Vec::new();
        f.write_format3_to(&mut buf).unwrap();
        let g = DotFile::from_reader(Cursor::new(&buf)).unwrap();
        assert_eq!(g.matrix.unwrap(), vec![0_i32; 24 * 24]);
    }
}
