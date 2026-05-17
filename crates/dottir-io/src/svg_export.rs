//! SVG exporter (D1).
//!
//! Writes an `<svg>` document that embeds the pixelmap as a base64
//! `data:image/png` `<image>` element, plus axis lines and tick labels
//! as SVG primitives. The result opens cleanly in any browser and is
//! resolution-independent for the axes (the pixelmap itself is still
//! a bitmap).
//!
//! Metadata is embedded in a `<metadata>` block so the SVG is
//! self-describing in the same way as the PNG `tEXt` chunks.

use std::fs::File;
use std::io::{BufWriter, Cursor, Write};
use std::path::Path;

use base64::Engine as _;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum SvgError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("PNG encoding error: {0}")]
    Png(#[from] png::EncodingError),
    #[error("pixelmap dimensions {width}×{height} don't match data length {len}")]
    DimensionMismatch { width: u32, height: u32, len: usize },
}

/// Render the pixelmap to an SVG file at `path`. The image is
/// embedded as a base64 PNG; axis ticks are drawn as SVG lines.
///
/// `margin` is the screen-pixel margin around the plot reserved for
/// axis labels (default 50 reads well at typical web zoom).
///
/// `metadata` is an optional list of (key, value) pairs included as
/// `<dottir:KEY>VALUE</dottir:KEY>` elements inside the SVG
/// `<metadata>` block. Mirror what the PNG `tEXt` writer accepts.
pub fn write_svg<P: AsRef<Path>>(
    path: P,
    width: u32,
    height: u32,
    pixels: &[u8],
    margin: u32,
    metadata: &[(&str, &str)],
) -> Result<(), SvgError> {
    if pixels.len() != (width as usize) * (height as usize) {
        return Err(SvgError::DimensionMismatch {
            width,
            height,
            len: pixels.len(),
        });
    }

    // Invert before encoding: raw kernel output uses 0 = no hit
    // (black), but SVG viewers expect the analysis-conventional
    // "white background, dark hits" look.
    let inverted = crate::text_overlay::inverted(pixels);
    let mut png_buf = Vec::with_capacity(inverted.len() + 1024);
    {
        let cursor = Cursor::new(&mut png_buf);
        let mut encoder = png::Encoder::new(cursor, width, height);
        encoder.set_color(png::ColorType::Grayscale);
        encoder.set_depth(png::BitDepth::Eight);
        let mut writer = encoder.write_header()?;
        writer.write_image_data(&inverted)?;
    }
    let png_b64 = base64::engine::general_purpose::STANDARD.encode(&png_buf);

    let total_w = width + 2 * margin;
    let total_h = height + 2 * margin;
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);

    // SVG attributes use double-quoted strings inside the format
    // string, so we use `write!` with `{}` placeholders rather than
    // raw strings (raw-string delimiters interact badly with the
    // literal `"#` in `fill="#222"`).
    writeln!(w, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>")?;
    writeln!(
        w,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" \
         xmlns:xlink=\"http://www.w3.org/1999/xlink\" \
         xmlns:dottir=\"https://github.com/petr/dottir/ns/1\" \
         width=\"{total_w}\" height=\"{total_h}\" \
         viewBox=\"0 0 {total_w} {total_h}\">"
    )?;

    if !metadata.is_empty() {
        writeln!(w, "  <metadata>")?;
        for (k, v) in metadata {
            writeln!(
                w,
                "    <dottir:{k}>{}</dottir:{k}>",
                xml_escape(v)
            )?;
        }
        writeln!(w, "  </metadata>")?;
    }
    writeln!(
        w,
        "  <rect x=\"0\" y=\"0\" width=\"{total_w}\" height=\"{total_h}\" fill=\"white\"/>"
    )?;
    writeln!(
        w,
        "  <image x=\"{margin}\" y=\"{margin}\" width=\"{width}\" height=\"{height}\" \
         image-rendering=\"pixelated\" \
         xlink:href=\"data:image/png;base64,{png_b64}\"/>"
    )?;
    // Axis box.
    writeln!(
        w,
        "  <rect x=\"{margin}\" y=\"{margin}\" width=\"{width}\" height=\"{height}\" \
         fill=\"none\" stroke=\"#444\" stroke-width=\"1\"/>"
    )?;

    // Top-axis ticks.
    let step_x = nice_step(width as u64);
    let mut t = 0u64;
    while t <= width as u64 {
        let x = margin + t as u32;
        writeln!(
            w,
            "  <line x1=\"{x}\" y1=\"{}\" x2=\"{x}\" y2=\"{}\" stroke=\"#444\" stroke-width=\"1\"/>",
            margin,
            margin.saturating_sub(6),
        )?;
        writeln!(
            w,
            "  <text x=\"{x}\" y=\"{}\" text-anchor=\"middle\" \
             font-family=\"monospace\" font-size=\"11\" fill=\"#222\">{}</text>",
            margin.saturating_sub(10),
            format_kb(t),
        )?;
        t += step_x;
    }
    // Left-axis ticks.
    let step_y = nice_step(height as u64);
    let mut t = 0u64;
    while t <= height as u64 {
        let y = margin + t as u32;
        writeln!(
            w,
            "  <line x1=\"{}\" y1=\"{y}\" x2=\"{}\" y2=\"{y}\" stroke=\"#444\" stroke-width=\"1\"/>",
            margin,
            margin.saturating_sub(6),
        )?;
        writeln!(
            w,
            "  <text x=\"{}\" y=\"{}\" text-anchor=\"end\" \
             font-family=\"monospace\" font-size=\"11\" fill=\"#222\">{}</text>",
            margin.saturating_sub(10),
            y + 4,
            format_kb(t),
        )?;
        t += step_y;
    }

    writeln!(w, "</svg>")?;
    w.flush()?;
    Ok(())
}

fn xml_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '&' => out.push_str("&amp;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            _ => out.push(c),
        }
    }
    out
}

fn nice_step(span: u64) -> u64 {
    if span == 0 {
        return 1;
    }
    let target = span as f64 / 8.0;
    let exp = target.log10().floor();
    let base = 10f64.powf(exp) as u64;
    for &m in &[1, 2, 5, 10] {
        let s = m * base;
        if (s as f64) >= target {
            return s.max(1);
        }
    }
    base.max(1)
}

fn format_kb(n: u64) -> String {
    if n >= 1_000_000 {
        format!("{:.1}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{}k", n / 1_000)
    } else {
        format!("{n}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_valid_svg_with_metadata() {
        let dir = std::env::temp_dir();
        let path = dir.join("dottir_svg_test.svg");
        let pixels: Vec<u8> = (0..255_u8).collect();
        write_svg(
            &path,
            17,
            15,
            &pixels,
            40,
            &[("matrix", "BLOSUM62"), ("window", "25")],
        )
        .unwrap();
        let s = std::fs::read_to_string(&path).unwrap();
        assert!(s.contains("<svg"));
        assert!(s.contains("data:image/png;base64,"));
        assert!(s.contains("BLOSUM62"));
        assert!(s.contains("dottir:window"));
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn dimension_mismatch_errors() {
        let dir = std::env::temp_dir();
        let path = dir.join("dottir_svg_test_mismatch.svg");
        let err = write_svg(&path, 10, 10, &[0_u8; 50], 30, &[]).unwrap_err();
        assert!(matches!(err, SvgError::DimensionMismatch { .. }));
    }

    #[test]
    fn xml_escape_basic() {
        assert_eq!(xml_escape("a&b<c>d"), "a&amp;b&lt;c&gt;d");
    }

    #[test]
    fn nice_step_sane_choices() {
        assert!(nice_step(100) > 0);
        assert!(nice_step(1_000_000) >= 100_000);
        assert_eq!(nice_step(0), 1);
    }
}
