//! PDF writer (spec §4.4.2): one page holding the composed plot image
//! (axes included) as a Flate-compressed image, sized to fit an A4-ish
//! page while keeping every pixel. The `text_chunks` that PNG exports
//! carry go into the document information dictionary.
//!
//! A raster page needs only an image XObject, so the PDF is written
//! directly rather than through a PDF library.

use std::io::Write;
use std::path::Path;

use flate2::write::ZlibEncoder;
use flate2::Compression;

/// Longest page side, in points (A4 is 595 × 842).
const MAX_PAGE_PT: f64 = 842.0;

/// Pixels of a composed image: greyscale or RGB, row-major.
#[derive(Debug, Clone, Copy)]
pub enum PdfPixels<'a> {
    Grey(&'a [u8]),
    Rgb(&'a [u8]),
}

/// A PDF literal string: parentheses and backslashes escaped, non-ASCII
/// replaced (the information dictionary uses PDFDocEncoding).
fn pdf_string(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    out.push('(');
    for c in s.chars() {
        match c {
            '(' | ')' | '\\' => {
                out.push('\\');
                out.push(c);
            }
            '\n' => out.push_str("\\n"),
            c if c.is_ascii() && !c.is_ascii_control() => out.push(c),
            _ => out.push('?'),
        }
    }
    out.push(')');
    out
}

/// A PDF name from a metadata key: characters outside `[A-Za-z0-9-_.]`
/// written as `#xx`.
fn pdf_name(s: &str) -> String {
    let mut out = String::from("/");
    for b in s.bytes() {
        if b.is_ascii_alphanumeric() || b"-_.".contains(&b) {
            out.push(b as char);
        } else {
            out.push_str(&format!("#{b:02X}"));
        }
    }
    out
}

/// The PDF as bytes.
pub fn encode_pdf(
    width: u32,
    height: u32,
    pixels: PdfPixels<'_>,
    text_chunks: &[(&str, &str)],
) -> std::io::Result<Vec<u8>> {
    let (data, channels, space) = match pixels {
        PdfPixels::Grey(p) => (p, 1, "/DeviceGray"),
        PdfPixels::Rgb(p) => (p, 3, "/DeviceRGB"),
    };
    if data.len() != width as usize * height as usize * channels {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("{width}×{height}×{channels} image, {} bytes", data.len()),
        ));
    }
    let mut z = ZlibEncoder::new(Vec::new(), Compression::default());
    z.write_all(data)?;
    let image = z.finish()?;

    let scale = (MAX_PAGE_PT / width.max(height).max(1) as f64).min(1.0);
    let (pw, ph) = (width as f64 * scale, height as f64 * scale);
    let content = format!("q {pw:.3} 0 0 {ph:.3} 0 0 cm /Im0 Do Q\n");
    let mut info = String::from("<< /Producer (dottir)");
    for (k, v) in text_chunks {
        info.push(' ');
        info.push_str(&pdf_name(k));
        info.push(' ');
        info.push_str(&pdf_string(v));
    }
    info.push_str(" >>");

    let mut out: Vec<u8> = b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n".to_vec();
    let mut offsets = Vec::new();
    let mut object = |out: &mut Vec<u8>, body: &[u8]| {
        offsets.push(out.len());
        out.extend_from_slice(format!("{} 0 obj\n", offsets.len()).as_bytes());
        out.extend_from_slice(body);
        out.extend_from_slice(b"\nendobj\n");
    };
    object(&mut out, b"<< /Type /Catalog /Pages 2 0 R >>");
    object(&mut out, b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>");
    object(
        &mut out,
        format!(
            "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {pw:.3} {ph:.3}] \
             /Resources << /XObject << /Im0 4 0 R >> >> /Contents 5 0 R >>"
        )
        .as_bytes(),
    );
    let mut img = format!(
        "<< /Type /XObject /Subtype /Image /Width {width} /Height {height} \
         /ColorSpace {space} /BitsPerComponent 8 /Filter /FlateDecode \
         /Length {} >>\nstream\n",
        image.len()
    )
    .into_bytes();
    img.extend_from_slice(&image);
    img.extend_from_slice(b"\nendstream");
    object(&mut out, &img);
    let mut c = format!("<< /Length {} >>\nstream\n", content.len()).into_bytes();
    c.extend_from_slice(content.as_bytes());
    c.extend_from_slice(b"endstream");
    object(&mut out, &c);
    object(&mut out, info.as_bytes());

    let xref = out.len();
    out.extend_from_slice(
        format!("xref\n0 {}\n0000000000 65535 f \n", offsets.len() + 1).as_bytes(),
    );
    for o in &offsets {
        out.extend_from_slice(format!("{o:010} 00000 n \n").as_bytes());
    }
    out.extend_from_slice(
        format!(
            "trailer\n<< /Size {} /Root 1 0 R /Info 6 0 R >>\nstartxref\n{xref}\n%%EOF\n",
            offsets.len() + 1
        )
        .as_bytes(),
    );
    Ok(out)
}

/// Write a composed image (see [`crate::text_overlay::compose_image_with_axes`])
/// as a one-page PDF.
pub fn write_pdf<P: AsRef<Path>>(
    path: P,
    width: u32,
    height: u32,
    pixels: PdfPixels<'_>,
    text_chunks: &[(&str, &str)],
) -> std::io::Result<()> {
    std::fs::write(path, encode_pdf(width, height, pixels, text_chunks)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    #[test]
    fn writes_a_well_formed_single_page_pdf() {
        let rgb: Vec<u8> = (0..6 * 4 * 3).map(|i| i as u8).collect();
        let pdf = encode_pdf(6, 4, PdfPixels::Rgb(&rgb), &[("dottir-version", "0.6 (x)")]).unwrap();
        let text = String::from_utf8_lossy(&pdf);
        assert!(text.starts_with("%PDF-1.4"));
        assert!(text.ends_with("%%EOF\n"));
        assert!(text.contains("/dottir-version (0.6 \\(x\\))"));
        assert!(text.contains("/Width 6 /Height 4 /ColorSpace /DeviceRGB"));
        // Every xref offset points at its object header.
        let rfind = |pat: &[u8]| pdf.windows(pat.len()).rposition(|w| w == pat).unwrap();
        let xref = rfind(b"\nxref\n") + 1;
        let tail = std::str::from_utf8(&pdf[xref..]).unwrap();
        let start: usize = tail[tail.find("startxref\n").unwrap() + 10..]
            .lines()
            .next()
            .unwrap()
            .parse()
            .unwrap();
        assert_eq!(start, xref);
        for (i, line) in tail.lines().skip(3).take(6).enumerate() {
            let off: usize = line[..10].parse().unwrap();
            assert!(pdf[off..].starts_with(format!("{} 0 obj", i + 1).as_bytes()));
        }
        // The image stream inflates to the pixels.
        let s = pdf.windows(7).position(|w| w == b"stream\n").unwrap() + 7;
        let mut back = Vec::new();
        flate2::read::ZlibDecoder::new(&pdf[s..])
            .read_to_end(&mut back)
            .unwrap();
        assert_eq!(back, rgb);
        assert!(encode_pdf(6, 4, PdfPixels::Grey(&rgb), &[]).is_err());
    }
}
