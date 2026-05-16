//! PNG writer for a dottir pixelmap.
//!
//! Greyscale 8-bit PNG. Optionally embeds a `tEXt`/`zTXt` chunk per spec
//! §4.4.4 capturing the dotplot parameters (matrix, window size, zoom,
//! pixel factor, input file SHA-256s, dottir version) so the image is
//! self-describing.

use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum PngError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("PNG encoder error: {0}")]
    Encoding(#[from] png::EncodingError),
    #[error("pixelmap dimensions {width}×{height} don't match data length {len}")]
    DimensionMismatch { width: u32, height: u32, len: usize },
}

/// Write a greyscale 8-bit PNG. `pixels` must be exactly `width * height`
/// bytes, row-major. Optional `text_chunks` are added as `tEXt` entries
/// for parameter provenance.
pub fn write_grayscale_png<P: AsRef<Path>>(
    path: P,
    width: u32,
    height: u32,
    pixels: &[u8],
    text_chunks: &[(&str, &str)],
) -> Result<(), PngError> {
    if pixels.len() != (width as usize).saturating_mul(height as usize) {
        return Err(PngError::DimensionMismatch {
            width,
            height,
            len: pixels.len(),
        });
    }
    let file = File::create(path)?;
    let w = BufWriter::new(file);
    let mut encoder = png::Encoder::new(w, width, height);
    encoder.set_color(png::ColorType::Grayscale);
    encoder.set_depth(png::BitDepth::Eight);
    for (k, v) in text_chunks {
        // Keep tEXt keys ASCII-clean and within PNG's 79-byte limit; the
        // png crate enforces this internally and returns an error if not.
        encoder.add_text_chunk((*k).to_string(), (*v).to_string())?;
    }
    let mut writer = encoder.write_header()?;
    writer.write_image_data(pixels)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_greyscale() {
        let dir = std::env::temp_dir();
        let path = dir.join("dottir_png_test.png");
        let pixels: Vec<u8> = (0..255_u8).collect();
        write_grayscale_png(&path, 17, 15, &pixels, &[("dottir", "test")]).unwrap();
        // Decode it back and check the data.
        let file = File::open(&path).unwrap();
        let decoder = png::Decoder::new(file);
        let mut reader = decoder.read_info().unwrap();
        let mut out = vec![0_u8; reader.output_buffer_size()];
        let info = reader.next_frame(&mut out).unwrap();
        assert_eq!(info.width, 17);
        assert_eq!(info.height, 15);
        assert_eq!(&out[..info.buffer_size()], pixels.as_slice());
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn dimension_mismatch_errors() {
        let dir = std::env::temp_dir();
        let path = dir.join("dottir_png_test_mismatch.png");
        let err = write_grayscale_png(&path, 10, 10, &[0_u8; 50], &[]).unwrap_err();
        assert!(matches!(err, PngError::DimensionMismatch { .. }));
        std::fs::remove_file(&path).ok();
    }
}
