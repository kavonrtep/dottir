//! Genome-mode outputs (ADR 0008): the PNG and the SVG (identity-coloured
//! alignments, record separators, shaded mix blocks, record names along
//! the top (target) and left (query) edges), the association table and
//! FASTA exports.

use std::io::Write;
use std::path::Path;

use dottir_core::genome::order::Placement;
use dottir_core::genome::raster::identity_class;
use dottir_core::genome::{AxisLayout, GenomeRaster, GenomeRecord};

use crate::faidx::{FaidxError, FastaIndex};
use crate::png_export::{write_rgb_png, PngError};
use crate::text_overlay::{draw_text, text_width, FONT_H};

/// Colours of the four D-GENIES identity classes (< 25, < 50, < 75,
/// ≥ 75 %), light to dark (viridis stops), readable on white.
pub const IDENTITY_COLOURS: [[u8; 3]; 4] = [
    [0xb5, 0xb0, 0x1e],
    [0x29, 0xaf, 0x7f],
    [0x2e, 0x6f, 0x8e],
    [0x44, 0x01, 0x54],
];

const MARGIN: usize = 24;
const MAX_LABEL_W: usize = 160;
const SEPARATOR: u8 = 205;
const MIX_SHADE: u8 = 238;
const FRAME: u8 = 80;
const INK: u8 = 30;

fn label(layout: &AxisLayout, recs: &[GenomeRecord], seg: usize) -> String {
    let s = &layout.segments[seg];
    let name = &recs[s.rec as usize].name;
    if s.reverse {
        format!("{name}-")
    } else {
        name.clone()
    }
}

/// Compose the plot with margins and labels; returns RGB pixels and size.
pub fn compose_genome_image(
    raster: &GenomeRaster,
    x: &AxisLayout,
    x_recs: &[GenomeRecord],
    y: &AxisLayout,
    y_recs: &[GenomeRecord],
) -> (Vec<u8>, u32, u32) {
    let (pw, ph) = (raster.width as usize, raster.height as usize);
    let bp = raster.bp_per_px;
    let px = |v: u64, max: usize| ((v / bp) as usize).min(max);
    let y_label_w = (0..y.segments.len())
        .map(|i| text_width(&label(y, y_recs, i)))
        .max()
        .unwrap_or(0)
        .min(MAX_LABEL_W);
    let (left, top) = (MARGIN + y_label_w + 6, MARGIN + FONT_H + 6);
    let (w, h) = (left + pw + MARGIN, top + ph + MARGIN);
    let mut grey = vec![255u8; w * h];
    let set = |g: &mut [u8], xx: usize, yy: usize, v: u8| {
        let p = &mut g[yy * w + xx];
        *p = (*p).min(v);
    };

    for m in &x.mixes {
        for xx in px(m.start, pw)..px(m.end, pw) {
            for yy in 0..ph {
                set(&mut grey, left + xx, top + yy, MIX_SHADE);
            }
        }
    }
    for m in &y.mixes {
        for yy in px(m.start, ph)..px(m.end, ph) {
            for xx in 0..pw {
                set(&mut grey, left + xx, top + yy, MIX_SHADE);
            }
        }
    }
    let boundaries = |l: &AxisLayout| -> Vec<u64> {
        l.segments
            .windows(2)
            .filter(|s| s[0].mix.is_none() || s[0].mix != s[1].mix)
            .map(|s| s[1].start)
            .collect()
    };
    for b in boundaries(x) {
        let xx = px(b, pw - 1);
        for yy in 0..ph {
            set(&mut grey, left + xx, top + yy, SEPARATOR);
        }
    }
    for b in boundaries(y) {
        let yy = px(b, ph - 1);
        for xx in 0..pw {
            set(&mut grey, left + xx, top + yy, SEPARATOR);
        }
    }
    for xx in left - 1..=left + pw {
        set(&mut grey, xx, top - 1, FRAME);
        set(&mut grey, xx, top + ph, FRAME);
    }
    for yy in top - 1..=top + ph {
        set(&mut grey, left - 1, yy, FRAME);
        set(&mut grey, left + pw, yy, FRAME);
    }

    // Labels: one per record, or one per mix block; skipped when they
    // do not fit their span.
    let spans = |l: &AxisLayout, recs: &[GenomeRecord]| -> Vec<(u64, u64, String)> {
        let mut out: Vec<(u64, u64, String)> = Vec::new();
        for (i, s) in l.segments.iter().enumerate() {
            match s.mix {
                None => out.push((s.start, s.start + s.len, label(l, recs, i))),
                Some(m) if i == 0 || l.segments[i - 1].mix != Some(m) => {
                    let r = &l.mixes[m as usize];
                    let n = l.segments.iter().filter(|t| t.mix == Some(m)).count();
                    out.push((r.start, r.end, format!("{n} small")));
                }
                Some(_) => {}
            }
        }
        out
    };
    for (a, b, text) in spans(x, x_recs) {
        let (xa, xb) = (px(a, pw), px(b, pw));
        let tw = text_width(&text);
        if tw + 4 <= xb - xa {
            let lx = left + (xa + xb) / 2 - tw / 2;
            draw_text(&mut grey, w, h, lx, top - 4 - FONT_H, &text, INK);
        }
    }
    for (a, b, text) in spans(y, y_recs) {
        let (ya, yb) = (px(a, ph), px(b, ph));
        let tw = text_width(&text);
        if FONT_H + 2 <= yb - ya && tw <= MAX_LABEL_W {
            let ly = top + (ya + yb) / 2 - FONT_H / 2;
            draw_text(&mut grey, w, h, left - 4 - tw, ly, &text, INK);
        }
    }

    let mut rgb: Vec<u8> = grey.iter().flat_map(|&g| [g, g, g]).collect();
    for yy in 0..ph {
        for xx in 0..pw {
            if let Some(c) = identity_class(raster.pixels[yy * pw + xx]) {
                let o = ((top + yy) * w + left + xx) * 3;
                rgb[o..o + 3].copy_from_slice(&IDENTITY_COLOURS[c]);
            }
        }
    }
    (rgb, w as u32, h as u32)
}

/// Compose and write the plot as an RGB PNG.
pub fn write_genome_png(
    path: &Path,
    raster: &GenomeRaster,
    x: &AxisLayout,
    x_recs: &[GenomeRecord],
    y: &AxisLayout,
    y_recs: &[GenomeRecord],
    text_chunks: &[(&str, &str)],
) -> Result<(), PngError> {
    let (rgb, w, h) = compose_genome_image(raster, x, x_recs, y, y_recs);
    write_rgb_png(path, w, h, &rgb, text_chunks)
}

/// Longest side of the SVG plot area, in user units.
const SVG_PLOT: f64 = 1000.0;
const SVG_FONT: f64 = 10.0;

/// The plot as SVG: the same elements as the PNG, with alignments as
/// vector lines (`lines` from [`dottir_core::genome::plot_lines`]),
/// `line_width` units wide, coloured by identity class.
#[allow(clippy::too_many_arguments)]
pub fn write_genome_svg<W: Write>(
    mut w: W,
    lines: &[dottir_core::genome::PlotLine],
    line_width: f64,
    x: &AxisLayout,
    x_recs: &[GenomeRecord],
    y: &AxisLayout,
    y_recs: &[GenomeRecord],
    metadata: &[(&str, &str)],
) -> std::io::Result<()> {
    use crate::svg_export::xml_escape;
    let scale = SVG_PLOT / x.total.max(y.total).max(1) as f64;
    let (pw, ph) = (x.total as f64 * scale, y.total as f64 * scale);
    // Text width estimate: monospace-ish 0.6 em per character.
    let tw = |t: &str| t.chars().count() as f64 * SVG_FONT * 0.6;
    let y_label_w = (0..y.segments.len())
        .map(|i| tw(&label(y, y_recs, i)))
        .fold(0.0, f64::max)
        .min(MAX_LABEL_W as f64);
    let (left, top) = (
        MARGIN as f64 + y_label_w + 6.0,
        MARGIN as f64 + SVG_FONT + 6.0,
    );
    let (total_w, total_h) = (left + pw + MARGIN as f64, top + ph + MARGIN as f64);
    let grey = |g: u8| format!("#{g:02x}{g:02x}{g:02x}");

    writeln!(
        w,
        "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>"
    )?;
    writeln!(
        w,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" \
         xmlns:dottir=\"https://github.com/petr/dottir/ns/1\" \
         width=\"{total_w:.0}\" height=\"{total_h:.0}\" viewBox=\"0 0 {total_w:.1} {total_h:.1}\">"
    )?;
    if !metadata.is_empty() {
        writeln!(w, "  <metadata>")?;
        for (k, v) in metadata {
            writeln!(w, "    <dottir:{k}>{}</dottir:{k}>", xml_escape(v))?;
        }
        writeln!(w, "  </metadata>")?;
    }
    writeln!(w, "  <rect width=\"100%\" height=\"100%\" fill=\"white\"/>")?;
    writeln!(w, "  <g transform=\"translate({left:.2} {top:.2})\">")?;
    // Mix blocks.
    writeln!(w, "    <g fill=\"{}\">", grey(MIX_SHADE))?;
    for m in &x.mixes {
        let (a, b) = (m.start as f64 * scale, m.end as f64 * scale);
        writeln!(
            w,
            "      <rect x=\"{a:.2}\" y=\"0\" width=\"{:.2}\" height=\"{ph:.2}\"/>",
            b - a
        )?;
    }
    for m in &y.mixes {
        let (a, b) = (m.start as f64 * scale, m.end as f64 * scale);
        writeln!(
            w,
            "      <rect x=\"0\" y=\"{a:.2}\" width=\"{pw:.2}\" height=\"{:.2}\"/>",
            b - a
        )?;
    }
    writeln!(w, "    </g>")?;
    // Record separators.
    let boundaries = |l: &AxisLayout| -> Vec<f64> {
        l.segments
            .windows(2)
            .filter(|s| s[0].mix.is_none() || s[0].mix != s[1].mix)
            .map(|s| s[1].start as f64 * scale)
            .collect()
    };
    writeln!(
        w,
        "    <g stroke=\"{}\" stroke-width=\"0.5\">",
        grey(SEPARATOR)
    )?;
    for b in boundaries(x) {
        writeln!(
            w,
            "      <line x1=\"{b:.2}\" y1=\"0\" x2=\"{b:.2}\" y2=\"{ph:.2}\"/>"
        )?;
    }
    for b in boundaries(y) {
        writeln!(
            w,
            "      <line x1=\"0\" y1=\"{b:.2}\" x2=\"{pw:.2}\" y2=\"{b:.2}\"/>"
        )?;
    }
    writeln!(w, "    </g>")?;
    // Alignments, one group per identity class, low identity first so
    // high identity is drawn on top (the PNG keeps the highest).
    let mut by_class: [Vec<&dottir_core::genome::PlotLine>; 4] = Default::default();
    for l in lines {
        let c = identity_class(dottir_core::genome::raster::identity_code(l.2)).unwrap_or(0);
        by_class[c].push(l);
    }
    for (c, ls) in by_class.iter().enumerate() {
        if ls.is_empty() {
            continue;
        }
        let [r, g, b] = IDENTITY_COLOURS[c];
        writeln!(
            w,
            "    <g stroke=\"#{r:02x}{g:02x}{b:02x}\" stroke-width=\"{line_width}\" \
             stroke-linecap=\"round\">"
        )?;
        for ((x0, y0), (x1, y1), _) in ls {
            writeln!(
                w,
                "      <line x1=\"{:.2}\" y1=\"{:.2}\" x2=\"{:.2}\" y2=\"{:.2}\"/>",
                *x0 as f64 * scale,
                *y0 as f64 * scale,
                *x1 as f64 * scale,
                *y1 as f64 * scale
            )?;
        }
        writeln!(w, "    </g>")?;
    }
    writeln!(
        w,
        "    <rect width=\"{pw:.2}\" height=\"{ph:.2}\" fill=\"none\" stroke=\"{}\"/>",
        grey(FRAME)
    )?;
    writeln!(w, "  </g>")?;
    // Labels, skipped when they do not fit their span.
    let spans = |l: &AxisLayout, recs: &[GenomeRecord]| -> Vec<(u64, u64, String)> {
        let mut out: Vec<(u64, u64, String)> = Vec::new();
        for (i, s) in l.segments.iter().enumerate() {
            match s.mix {
                None => out.push((s.start, s.start + s.len, label(l, recs, i))),
                Some(m) if i == 0 || l.segments[i - 1].mix != Some(m) => {
                    let r = &l.mixes[m as usize];
                    let n = l.segments.iter().filter(|t| t.mix == Some(m)).count();
                    out.push((r.start, r.end, format!("{n} small")));
                }
                Some(_) => {}
            }
        }
        out
    };
    writeln!(
        w,
        "  <g font-family=\"monospace\" font-size=\"{SVG_FONT}\" fill=\"{}\">",
        grey(INK)
    )?;
    for (a, b, text) in spans(x, x_recs) {
        let (xa, xb) = (a as f64 * scale, b as f64 * scale);
        if tw(&text) + 4.0 <= xb - xa {
            writeln!(
                w,
                "    <text x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"middle\">{}</text>",
                left + (xa + xb) / 2.0,
                top - 4.0,
                xml_escape(&text)
            )?;
        }
    }
    for (a, b, text) in spans(y, y_recs) {
        let (ya, yb) = (a as f64 * scale, b as f64 * scale);
        if SVG_FONT + 2.0 <= yb - ya && tw(&text) <= MAX_LABEL_W as f64 {
            writeln!(
                w,
                "    <text x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"end\" \
                 dominant-baseline=\"middle\">{}</text>",
                left - 4.0,
                top + (ya + yb) / 2.0,
                xml_escape(&text)
            )?;
        }
    }
    writeln!(w, "  </g>")?;
    writeln!(w, "</svg>")?;
    Ok(())
}

/// Complement of an IUPAC nucleotide code, preserving case.
pub(crate) fn complement(b: u8) -> u8 {
    let c = match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        x => x,
    };
    if b.is_ascii_lowercase() {
        c.to_ascii_lowercase()
    } else {
        c
    }
}

/// Association table (D-GENIES columns): one row per vertical-axis
/// (query) record in plot order, with its best target record and the
/// extent of the alignments to it. Coordinates are 0-based, forward
/// strand, end-exclusive; `na` without alignments.
pub fn write_association<W: Write>(
    mut w: W,
    y: &AxisLayout,
    y_recs: &[GenomeRecord],
    x_recs: &[GenomeRecord],
    placements: &[Option<Placement>],
) -> std::io::Result<()> {
    writeln!(
        w,
        "Query\tTarget\tStrand\tQ-len\tQ-start\tQ-stop\tT-len\tT-start\tT-stop"
    )?;
    for s in &y.segments {
        let q = &y_recs[s.rec as usize];
        let strand = if s.reverse { '-' } else { '+' };
        match &placements[s.rec as usize] {
            Some(p) => {
                let t = &x_recs[p.t_rec as usize];
                writeln!(
                    w,
                    "{}\t{}\t{strand}\t{}\t{}\t{}\t{}\t{}\t{}",
                    q.name, t.name, q.len, p.q_span.0, p.q_span.1, t.len, p.t_span.0, p.t_span.1
                )?;
            }
            None => writeln!(
                w,
                "{}\tNone\t{strand}\t{}\tna\tna\tna\tna\tna",
                q.name, q.len
            )?,
        }
    }
    w.flush()
}

/// Write `records` of an indexed FASTA in `order`, reverse-complementing
/// those with `flips[rec]` (their header gets ` dottir:revcomp`).
pub fn write_fasta_records<W: Write>(
    mut w: W,
    index: &FastaIndex,
    order: &[u32],
    flips: &[bool],
) -> Result<(), FaidxError> {
    let io = |e: std::io::Error| FaidxError::Io {
        path: index.path.clone(),
        source: e,
    };
    for &rec in order {
        let e = &index.entries[rec as usize];
        let mut seq = index.fetch_raw(rec as usize, 0..e.len)?;
        let rc = flips.get(rec as usize).copied().unwrap_or(false);
        if rc {
            seq.reverse();
            seq.iter_mut().for_each(|b| *b = complement(*b));
        }
        writeln!(w, ">{}{}", e.name, if rc { " dottir:revcomp" } else { "" }).map_err(io)?;
        for line in seq.chunks(80) {
            w.write_all(line).map_err(io)?;
            w.write_all(b"\n").map_err(io)?;
        }
    }
    w.flush().map_err(io)
}

#[cfg(test)]
mod tests {
    use super::*;
    use dottir_core::genome::raster::identity_code;

    #[test]
    fn plot_pixels_are_coloured_by_class_inside_the_frame() {
        let recs = vec![
            GenomeRecord {
                name: "chr1".into(),
                len: 1000,
            },
            GenomeRecord {
                name: "chr2".into(),
                len: 1000,
            },
        ];
        let l = AxisLayout::new(&recs, None, &[], 0.0);
        let mut raster = GenomeRaster {
            width: 200,
            height: 200,
            bp_per_px: 10,
            pixels: vec![0; 200 * 200],
        };
        raster.pixels[0] = identity_code(0.9);
        let (rgb, w, _) = compose_genome_image(&raster, &l, &recs, &l, &recs);
        let top = MARGIN + FONT_H + 6;
        let left = MARGIN + text_width("chr1") + 6;
        let o = (top * w as usize + left) * 3;
        assert_eq!(&rgb[o..o + 3], &IDENTITY_COLOURS[3]);
        // Record separator at x = 100 px.
        let o = ((top + 50) * w as usize + left + 100) * 3;
        assert_eq!(rgb[o], SEPARATOR);
    }

    #[test]
    fn svg_draws_alignments_as_lines_by_identity_class() {
        let recs = vec![
            GenomeRecord {
                name: "chr<1>".into(),
                len: 1000,
            },
            GenomeRecord {
                name: "chr2".into(),
                len: 1000,
            },
        ];
        let l = AxisLayout::new(&recs, None, &[], 0.0);
        let aln = dottir_core::genome::Alignment {
            q_rec: 1,
            q_start: 0,
            q_end: 500,
            t_rec: 0,
            t_start: 100,
            t_end: 600,
            reverse: true,
            matches: 450,
            block_len: 500,
            mapq: 60,
            primary: true,
            divergence: None,
        };
        let lines = dottir_core::genome::plot_lines([&aln], &l, &l, false);
        let mut out = Vec::new();
        write_genome_svg(&mut out, &lines, 1.5, &l, &recs, &l, &recs, &[("k", "v")]).unwrap();
        let svg = String::from_utf8(out).unwrap();
        // 2000 bp on 1000 units: reverse, so target 600 at query 0.
        assert!(svg.contains(r#"<line x1="300.00" y1="500.00" x2="50.00" y2="750.00"/>"#));
        let [r, g, b] = IDENTITY_COLOURS[3];
        assert!(svg.contains(&format!("stroke=\"#{r:02x}{g:02x}{b:02x}\"")));
        assert!(svg.contains("chr&lt;1&gt;"));
        assert!(svg.contains("<dottir:k>v</dottir:k>"));
        assert!(svg.trim_end().ends_with("</svg>"));
    }

    #[test]
    fn fasta_export_orders_and_reverse_complements_keeping_case() {
        let d = std::env::temp_dir().join(format!("dottir-gexp-{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let p = d.join("a.fa");
        std::fs::write(&p, ">a\nACGTN\n>b\nAAcR\n").unwrap();
        let idx = FastaIndex::open(&p, &d).unwrap();
        let mut out = Vec::new();
        write_fasta_records(&mut out, &idx, &[1, 0], &[true, false]).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            ">b\nAAcR\n>a dottir:revcomp\nNACGT\n"
        );
    }
}
