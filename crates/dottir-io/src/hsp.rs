//! HSP overlays (spec §4.2.3, §4.4.2): alignments from a PAF file drawn
//! on a dotplot. PAF positions are in source records; sorted, flipped or
//! region records of the plotted sequences are handled through their
//! [`RegionOrigin`](crate::sequence::RegionOrigin), so lines land where
//! the plot shows those bases.

use std::collections::BTreeMap;
use std::path::Path;

use dottir_core::genome::{Alignment, GenomeRecord};

use crate::paf::{read_paf_subset_oriented, AlignmentPath, PafError};
use crate::sequence::{RecordSpan, Sequence};

/// How HSPs are shown, after C dotter's `DOTTER_HSPS_*` modes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum HspMode {
    #[default]
    Off,
    /// One colour.
    Line,
    /// Coloured by identity (the genome-mode bins).
    Score,
    /// Instead of the dot image: grey by identity on white.
    Greyscale,
}

impl HspMode {
    pub const ALL: [HspMode; 4] = [
        HspMode::Off,
        HspMode::Line,
        HspMode::Score,
        HspMode::Greyscale,
    ];

    pub fn name(self) -> &'static str {
        match self {
            HspMode::Off => "off",
            HspMode::Line => "line",
            HspMode::Score => "score",
            HspMode::Greyscale => "greyscale",
        }
    }

    pub fn parse(s: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|m| m.name() == s)
    }
}

/// One HSP as a segment in plot coordinates (concatenated residues;
/// x = query, horizontal; y = subject, vertical).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hsp {
    pub x0: f64,
    pub y0: f64,
    pub x1: f64,
    pub y1: f64,
    pub identity: f64,
    pub primary: bool,
}

/// HSPs of a PAF placed on a plot of `x` (horizontal) against `y`.
#[derive(Debug, Clone, Default)]
pub struct HspSet {
    pub path: std::path::PathBuf,
    /// Segments: one per alignment, or one per gapless piece when the
    /// PAF has CIGARs (`cg:Z`).
    pub hsps: Vec<Hsp>,
    /// PAF lines read, and those naming records not in the plot.
    pub alignments: usize,
    pub skipped: u64,
    /// Alignments with at least one segment on the plot.
    pub placed: usize,
    /// Alignments drawn along their CIGAR.
    pub with_cigar: usize,
}

/// Source records of `seq`: name and full length, first appearance first.
fn source_records(seq: &Sequence) -> Vec<GenomeRecord> {
    let mut seen = BTreeMap::new();
    let mut out = Vec::new();
    for r in &seq.records {
        let (name, len) = match &r.origin {
            Some(o) => (o.record.clone(), o.record_len),
            None => (r.id.clone(), r.len() as u64),
        };
        if seen.insert(name.clone(), ()).is_none() {
            out.push(GenomeRecord { name, len });
        }
    }
    out
}

/// Plotted records showing each source record, by source name.
fn shown_by_source(seq: &Sequence) -> BTreeMap<&str, Vec<&RecordSpan>> {
    let mut m: BTreeMap<&str, Vec<&RecordSpan>> = BTreeMap::new();
    for r in &seq.records {
        let name = r
            .origin
            .as_ref()
            .map_or(r.id.as_str(), |o| o.record.as_str());
        m.entry(name).or_default().push(r);
    }
    m
}

/// Source interval a plotted record shows, and its orientation.
fn shown_interval(r: &RecordSpan) -> (f64, f64, bool) {
    match &r.origin {
        Some(o) => (o.range.start as f64, o.range.end as f64, o.reverse),
        None => (0.0, r.len() as f64, false),
    }
}

/// Plot coordinate of source position `p` in plotted record `r`.
fn to_plot(r: &RecordSpan, p: f64) -> f64 {
    let (lo, hi, rev) = shown_interval(r);
    r.range.start as f64 + if rev { hi - p } else { p - lo }
}

/// Clip the segment `(x0, y0) → (x1, y1)` to `x ∈ [xl, xh]`,
/// `y ∈ [yl, yh]` (Liang–Barsky on the parameter).
fn clip(
    (x0, y0, x1, y1): (f64, f64, f64, f64),
    (xl, xh): (f64, f64),
    (yl, yh): (f64, f64),
) -> Option<(f64, f64, f64, f64)> {
    let (mut t0, mut t1) = (0.0_f64, 1.0_f64);
    for (p, q) in [
        (-(x1 - x0), x0 - xl),
        (x1 - x0, xh - x0),
        (-(y1 - y0), y0 - yl),
        (y1 - y0, yh - y0),
    ] {
        if p == 0.0 {
            if q < 0.0 {
                return None;
            }
        } else {
            let t = q / p;
            if p < 0.0 {
                t0 = t0.max(t);
            } else {
                t1 = t1.min(t);
            }
        }
    }
    (t0 < t1).then(|| {
        let (dx, dy) = (x1 - x0, y1 - y0);
        (x0 + t0 * dx, y0 + t0 * dy, x0 + t1 * dx, y0 + t1 * dy)
    })
}

/// Place `alns` (query = vertical records `yr`, target = horizontal `xr`)
/// on the plot of `x` against `y`. An alignment with a path (`paths[i]`,
/// from its CIGAR) is drawn as its gapless pieces, so a long alignment
/// with indels follows its dots; without one it is a straight line from
/// start to end.
pub fn place(
    alns: &[Alignment],
    paths: &[Option<AlignmentPath>],
    xr: &[GenomeRecord],
    yr: &[GenomeRecord],
    x: &Sequence,
    y: &Sequence,
    n_placed: &mut usize,
) -> Vec<Hsp> {
    let (xs, ys) = (shown_by_source(x), shown_by_source(y));
    let mut out = Vec::new();
    for (i, a) in alns.iter().enumerate() {
        let (Some(xrecs), Some(yrecs)) = (
            xs.get(xr[a.t_rec as usize].name.as_str()),
            ys.get(yr[a.q_rec as usize].name.as_str()),
        ) else {
            continue;
        };
        // In source coordinates: target along x, query along y.
        let segs: Vec<(f64, f64, f64, f64)> = match paths.get(i).and_then(Option::as_ref) {
            Some(p) => p
                .iter()
                .map(|&[q0, t0, q1, t1]| (t0 as f64, q0 as f64, t1 as f64, q1 as f64))
                .collect(),
            None if a.reverse => vec![(
                a.t_end as f64,
                a.q_start as f64,
                a.t_start as f64,
                a.q_end as f64,
            )],
            None => vec![(
                a.t_start as f64,
                a.q_start as f64,
                a.t_end as f64,
                a.q_end as f64,
            )],
        };
        let mut placed = false;
        for rx in xrecs {
            let (xl, xh, _) = shown_interval(rx);
            for ry in yrecs {
                let (yl, yh, _) = shown_interval(ry);
                for &seg in &segs {
                    if let Some((x0, y0, x1, y1)) = clip(seg, (xl, xh), (yl, yh)) {
                        out.push(Hsp {
                            x0: to_plot(rx, x0),
                            y0: to_plot(ry, y0),
                            x1: to_plot(rx, x1),
                            y1: to_plot(ry, y1),
                            identity: a.identity(),
                            primary: a.primary,
                        });
                        placed = true;
                    }
                }
            }
        }
        *n_placed += usize::from(placed);
    }
    out
}

/// Read a PAF (either orientation; records not plotted are skipped) and
/// place it on the plot of `x` (horizontal) against `y` (vertical).
pub fn load(path: &Path, x: &Sequence, y: &Sequence) -> Result<HspSet, PafError> {
    let (xr, yr) = (source_records(x), source_records(y));
    let (alns, skipped, _, paths) = read_paf_subset_oriented(path, &yr, &xr)?;
    let mut placed = 0;
    let hsps = place(&alns, &paths, &xr, &yr, x, y, &mut placed);
    Ok(HspSet {
        path: path.to_path_buf(),
        hsps,
        alignments: alns.len(),
        skipped,
        placed,
        with_cigar: paths.iter().filter(|p| p.is_some()).count(),
    })
}

/// Colour of an HSP line in `mode`.
pub fn colour(mode: HspMode, identity: f64) -> [u8; 3] {
    match mode {
        HspMode::Off | HspMode::Line => [220, 30, 30],
        // Genome-mode bins (D-GENIES): < 25, 25–50, 50–75, ≥ 75 %.
        HspMode::Score => {
            crate::genome_export::IDENTITY_COLOURS[((identity * 4.0) as usize).min(3)]
        }
        HspMode::Greyscale => {
            // Darker for higher identity; 50 % identity and below light.
            let g = 230.0 - 230.0 * ((identity - 0.5) / 0.5).clamp(0.0, 1.0);
            [g as u8; 3]
        }
    }
}

/// Draw `hsps` into an RGBA pixelmap of `w × h` showing residues from
/// `(q_off, s_off)` at `zoom` residues per pixel. Greyscale replaces the
/// image with white first.
#[allow(clippy::too_many_arguments)]
pub fn draw_rgba(
    rgba: &mut [u8],
    w: usize,
    h: usize,
    hsps: &[Hsp],
    mode: HspMode,
    q_off: usize,
    s_off: usize,
    zoom: usize,
) {
    if mode == HspMode::Off {
        return;
    }
    if mode == HspMode::Greyscale {
        rgba.fill(255);
    }
    let z = zoom.max(1) as f64;
    let px = |v: f64, off: usize| (v - off as f64) / z;
    for hsp in hsps {
        let c = colour(mode, hsp.identity);
        let (x0, y0) = (px(hsp.x0, q_off), px(hsp.y0, s_off));
        let (x1, y1) = (px(hsp.x1, q_off), px(hsp.y1, s_off));
        let steps = (x1 - x0).abs().max((y1 - y0).abs()).ceil().max(1.0) as usize;
        for i in 0..=steps {
            let t = i as f64 / steps as f64;
            let (x, y) = (x0 + t * (x1 - x0), y0 + t * (y1 - y0));
            if x < 0.0 || y < 0.0 {
                continue;
            }
            let (xi, yi) = (x as usize, y as usize);
            if xi < w && yi < h {
                let k = (yi * w + xi) * 4;
                rgba[k..k + 3].copy_from_slice(&c);
                rgba[k + 3] = 255;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::parse_fasta;

    fn seq(text: &str) -> Sequence {
        Sequence::from_records(parse_fasta(text).unwrap(), None)
    }

    fn write(name: &str, text: &str) -> std::path::PathBuf {
        let d = std::env::temp_dir().join(format!("dottir_hsp_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let p = d.join(name);
        std::fs::write(&p, text).unwrap();
        p
    }

    #[test]
    fn places_hsps_in_either_paf_orientation_and_skips_other_records() {
        let x = seq(">a\nACGTACGTAC\n>b\nACGTACGTACGTACGTACGT\n");
        let y = seq(">s\nACGTACGTACGTACGTACGTACGTACGTAC\n");
        // Query s 0..10 on target b 5..15, forward; one line for a record
        // not in the plot.
        let paf = "s\t30\t0\t10\t+\tb\t20\t5\t15\t9\t10\t60\n\
                   s\t30\t0\t10\t+\tzz\t50\t0\t10\t9\t10\t60\n";
        let set = load(&write("a.paf", paf), &x, &y).unwrap();
        assert_eq!((set.alignments, set.skipped), (1, 1));
        // b starts at plot x 10.
        assert_eq!(
            set.hsps,
            [Hsp {
                x0: 15.0,
                y0: 0.0,
                x1: 25.0,
                y1: 10.0,
                identity: 0.9,
                primary: true
            }]
        );
        // The same alignment written the other way round.
        let swapped = "b\t20\t5\t15\t+\ts\t30\t0\t10\t9\t10\t60\n";
        let set2 = load(&write("b.paf", swapped), &x, &y).unwrap();
        assert_eq!(set2.hsps, set.hsps);
    }

    #[test]
    fn hsps_follow_reversed_and_region_records() {
        // x shows source b reversed (whole record), y a region s:10-20.
        let x = seq(">b\nACGTACGTACGTACGTACGT\n").reordered(&[0], &[true]);
        let mut y = seq(">s:11-20\nACGTACGTAC\n");
        y.records[0].origin = Some(crate::sequence::RegionOrigin {
            record: "s".into(),
            record_len: 30,
            range: 10..20,
            reverse: false,
        });
        // Forward alignment s 0..20 on b 0..20: clipped to s 10..20
        // (b 10..20), then b mirrored: x 10..0.
        let paf = "s\t30\t0\t20\t+\tb\t20\t0\t20\t20\t20\t60\n";
        let set = load(&write("c.paf", paf), &x, &y).unwrap();
        assert_eq!(set.hsps.len(), 1);
        let h = set.hsps[0];
        assert_eq!((h.x0, h.y0, h.x1, h.y1), (10.0, 0.0, 0.0, 10.0));
    }

    #[test]
    fn hsps_with_a_cigar_follow_it_across_indels() {
        let x = seq(&format!(">b\n{}\n", "ACGT".repeat(10)));
        let y = seq(&format!(">s\n{}\n", "ACG".repeat(10)));
        let seg = |h: &Hsp| (h.x0, h.y0, h.x1, h.y1);
        // A 20 bp deletion splits the line into its two gapless pieces.
        let fwd = "s\t30\t0\t10\t+\tb\t40\t0\t30\t10\t30\t60\tcg:Z:5M20D5M\n";
        let set = load(&write("cg1.paf", fwd), &x, &y).unwrap();
        assert_eq!((set.placed, set.with_cigar), (1, 1));
        let segs: Vec<_> = set.hsps.iter().map(seg).collect();
        assert_eq!(segs, [(0.0, 0.0, 5.0, 5.0), (25.0, 5.0, 30.0, 10.0)]);
        // The same PAF line written the other way round (I for D).
        let swapped = "b\t40\t0\t30\t+\ts\t30\t0\t10\t10\t30\t60\tcg:Z:5M20I5M\n";
        let set2 = load(&write("cg2.paf", swapped), &x, &y).unwrap();
        assert_eq!(set2.hsps, set.hsps);
        // Reverse strand: the query runs down from its end.
        let rev = "s\t30\t0\t10\t-\tb\t40\t0\t30\t10\t30\t60\tcg:Z:5M20D5M\n";
        let set3 = load(&write("cg3.paf", rev), &x, &y).unwrap();
        let segs: Vec<_> = set3.hsps.iter().map(seg).collect();
        assert_eq!(segs, [(0.0, 10.0, 5.0, 5.0), (25.0, 5.0, 30.0, 0.0)]);
        // Small indels stay within one line; a CIGAR not spanning the
        // alignment is ignored (straight line).
        let small = "s\t30\t0\t10\t+\tb\t40\t0\t13\t10\t13\t60\tcg:Z:5M3D5M\n\
                     s\t30\t0\t10\t+\tb\t40\t20\t30\t10\t10\t60\tcg:Z:9M\n";
        let set4 = load(&write("cg4.paf", small), &x, &y).unwrap();
        let segs: Vec<_> = set4.hsps.iter().map(seg).collect();
        assert_eq!(segs, [(0.0, 0.0, 13.0, 10.0), (20.0, 0.0, 30.0, 10.0)]);
        assert_eq!(set4.with_cigar, 1);
    }

    #[test]
    fn rasterises_lines_and_greyscale_replaces_the_image() {
        let hsps = [Hsp {
            x0: 0.0,
            y0: 0.0,
            x1: 8.0,
            y1: 8.0,
            identity: 1.0,
            primary: true,
        }];
        let mut rgba = vec![0_u8; 4 * 4 * 4];
        draw_rgba(&mut rgba, 4, 4, &hsps, HspMode::Line, 0, 0, 2);
        for i in 0..4 {
            assert_eq!(&rgba[(i * 4 + i) * 4..(i * 4 + i) * 4 + 3], &[220, 30, 30]);
        }
        assert_eq!(&rgba[4..8], &[0, 0, 0, 0], "off-diagonal untouched");
        draw_rgba(&mut rgba, 4, 4, &hsps, HspMode::Greyscale, 0, 0, 2);
        assert_eq!(&rgba[4..8], &[255, 255, 255, 255], "background cleared");
        assert_eq!(&rgba[0..3], &[0, 0, 0], "100 % identity is black");
    }
}
