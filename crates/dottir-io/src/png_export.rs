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

/// Write a greyscale 8-bit PNG with a margin, axis ticks, and
/// numeric labels around the pixelmap. The input `pixels` is
/// inverted on the way out (raw kernel output is "0 = no hit", we
/// emit "no hit = white") so the result is the analysis-friendly
/// light-background plot.
///
/// `coord_w` and `coord_h` are the coordinate-space dimensions
/// shown along the axes (typically the input sequence lengths in
/// residues). They can differ from the pixel dimensions `width` /
/// `height` when the caller has nearest-neighbour upscaled the
/// pixelmap. For the no-rescale case pass `coord_w = width,
/// coord_h = height`.
pub fn write_grayscale_png_with_axes<P: AsRef<Path>>(
    path: P,
    width: u32,
    height: u32,
    coord_w: u32,
    coord_h: u32,
    pixels: &[u8],
    margin: u32,
    text_chunks: &[(&str, &str)],
) -> Result<(), PngError> {
    if pixels.len() != (width as usize).saturating_mul(height as usize) {
        return Err(PngError::DimensionMismatch {
            width,
            height,
            len: pixels.len(),
        });
    }
    let inverted = crate::text_overlay::inverted(pixels);
    let (canvas, total_w, total_h) = crate::text_overlay::compose_image_with_axes(
        &inverted, width, height, coord_w, coord_h, margin,
    );
    write_grayscale_png_raw(path, total_w, total_h, &canvas, text_chunks)
}

/// Lower-level: write a raw greyscale PNG with whatever pixels you
/// give it, no inversion, no axes. Use this when you've already
/// composed your own canvas (e.g. the GUI's "Save PNG" applies the
/// current greyramp LUT before calling this).
pub fn write_grayscale_png<P: AsRef<Path>>(
    path: P,
    width: u32,
    height: u32,
    pixels: &[u8],
    text_chunks: &[(&str, &str)],
) -> Result<(), PngError> {
    write_grayscale_png_raw(path, width, height, pixels, text_chunks)
}

fn write_grayscale_png_raw<P: AsRef<Path>>(
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
