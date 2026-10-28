//! Genome-mode outputs (ADR 0008): the PNG (identity-coloured
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
