//! I/O glue for large-sequence mode (ADR 0007): depth floors from
//! soft-masking and repeat annotations, chain export as PAF, and depth
//! tracks as bedGraph.

use std::io::{BufWriter, Write};
use std::path::Path;

use dottir_core::large::{Chain, DepthTrack};

use crate::annotation::AnnotSet;
use crate::sequence::Sequence;
use crate::svg_export::SvgLine;

/// Per-base depth floor: `value` on soft-masked bases, 0 elsewhere.
/// `None` when the sequence has no soft-masked bases.
pub fn soft_mask_floor(seq: &Sequence, value: u16) -> Option<Vec<u16>> {
    if seq.soft_mask.is_empty() {
        return None;
    }
    let mut floor = vec![0u16; seq.len()];
    for r in &seq.soft_mask {
        floor[r.start.min(seq.len())..r.end.min(seq.len())].fill(value);
    }
    Some(floor)
}

/// Per-base depth floor: `value` inside features whose type is in `types`
/// (case-insensitive; empty = every feature). Features on records not in
/// `seq` are skipped; the count of skipped features is returned too.
pub fn annotation_floor(
    seq: &Sequence,
    annots: &AnnotSet,
    types: &[String],
    value: u16,
) -> (Vec<u16>, usize) {
    let locator = seq.feature_locator();
    let mut floor = vec![0u16; seq.len()];
    let mut skipped = 0;
    for f in &annots.features {
        if !types.is_empty()
            && !types
                .iter()
                .any(|t| t.eq_ignore_ascii_case(&f.feature_type))
        {
            continue;
        }
        let at = locator.locate(&f.record, f.range.clone());
        if at.is_empty() {
            skipped += 1;
        }
        for r in at {
            floor[r].fill(value);
        }
    }
    (floor, skipped)
}

/// Element-wise max of two optional floors.
pub fn merge_floors(a: Option<Vec<u16>>, b: Option<Vec<u16>>) -> Option<Vec<u16>> {
    match (a, b) {
        (Some(mut a), Some(b)) => {
            for (x, y) in a.iter_mut().zip(b) {
                *x = (*x).max(y);
            }
            Some(a)
        }
        (a, None) => a,
        (None, b) => b,
    }
}

/// Write chains as PAF (columns 1–12 plus `tp:A:P`, `cm:i`, `s1:i`; no
/// CIGAR — chains are not alignments). Coordinates are record-local.
/// Column 10 (matches) is unknown and set to the anchor-covered estimate
/// `min(qlen, slen)`; column 12 (MAPQ) is 255 (unavailable).
pub fn write_chains_paf<P: AsRef<Path>>(
    path: P,
    chains: &[Chain],
    query: &Sequence,
    subject: &Sequence,
) -> std::io::Result<()> {
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    for c in chains {
        let (Some(qr), Some(sr)) = (
            query.records.get(c.q_rec as usize),
            subject.records.get(c.s_rec as usize),
        ) else {
            continue;
        };
        // Sorted or drill-down records are reported in their source
        // records; a reversed record flips the strand.
        let (qn, ql, qs, qe, qrev) = qr.source_interval(
            c.q_start as usize - qr.range.start,
            c.q_end as usize - qr.range.start,
        );
        let (sn, sl, ss, se, srev) = sr.source_interval(
            c.s_start as usize - sr.range.start,
            c.s_end as usize - sr.range.start,
        );
        let block = (qe - qs).max(se - ss);
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\ttp:A:P\tcm:i:{}\ts1:i:{}",
            qn,
            ql,
            qs,
            qe,
            if c.reverse ^ qrev ^ srev { '-' } else { '+' },
            sn,
            sl,
            ss,
            se,
            (qe - qs).min(se - ss),
            block,
            c.n_anchors,
            c.score.round() as i64,
        )?;
    }
    w.flush()
}

/// Chain colours: forward dark green, reverse magenta (as the canvas).
pub const CHAIN_FORWARD: [u8; 3] = [0, 140, 60];
pub const CHAIN_REVERSE: [u8; 3] = [190, 0, 190];

/// Chains as vector lines in pixelmap pixels at `zoom` residues per
/// pixel; `mirror` also draws them across the diagonal (self-comparison).
pub fn chain_svg_lines(chains: &[Chain], zoom: usize, mirror: bool) -> Vec<SvgLine> {
    let z = zoom.max(1) as f64;
    let mut out = Vec::new();
    for c in chains {
        let (q0, q1) = (c.q_start as f64 / z, c.q_end as f64 / z);
        let (s0, s1) = (c.s_start as f64 / z, c.s_end as f64 / z);
        let (from, to, rgb) = if c.reverse {
            ((q0, s1), (q1, s0), CHAIN_REVERSE)
        } else {
            ((q0, s0), (q1, s1), CHAIN_FORWARD)
        };
        out.push(SvgLine {
            from,
            to,
            rgb,
            width: 1.5,
        });
        if mirror {
            out.push(SvgLine {
                from: (from.1, from.0),
                to: (to.1, to.0),
                rgb,
                width: 1.5,
            });
        }
    }
    out
}

/// Write chains as PAF, or as BEDPE when `path` ends in `.bedpe`.
pub fn write_chains<P: AsRef<Path>>(
    path: P,
    chains: &[Chain],
    query: &Sequence,
    subject: &Sequence,
) -> std::io::Result<()> {
    let bedpe = path
        .as_ref()
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("bedpe"));
    if bedpe {
        write_chains_bedpe(path, chains, query, subject)
    } else {
        write_chains_paf(path, chains, query, subject)
    }
}

/// Write chains as BEDPE: query interval, subject interval (0-based,
/// end-exclusive, source records), name `chainN`, score, and strands
/// (`+` for the query, the chain's relative strand for the subject).
pub fn write_chains_bedpe<P: AsRef<Path>>(
    path: P,
    chains: &[Chain],
    query: &Sequence,
    subject: &Sequence,
) -> std::io::Result<()> {
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    for (i, c) in chains.iter().enumerate() {
        let (Some(qr), Some(sr)) = (
            query.records.get(c.q_rec as usize),
            subject.records.get(c.s_rec as usize),
        ) else {
            continue;
        };
        let (qn, _, qs, qe, qrev) = qr.source_interval(
            c.q_start as usize - qr.range.start,
            c.q_end as usize - qr.range.start,
        );
        let (sn, _, ss, se, srev) = sr.source_interval(
            c.s_start as usize - sr.range.start,
            c.s_end as usize - sr.range.start,
        );
        writeln!(
            w,
            "{qn}\t{qs}\t{qe}\t{sn}\t{ss}\t{se}\tchain{}\t{}\t+\t{}",
            i + 1,
            c.score.round() as i64,
            if c.reverse ^ qrev ^ srev { '-' } else { '+' },
        )?;
    }
    w.flush()
}

/// Write a depth track as bedGraph, averaged over `bin`-base bins (source-
/// record coordinates; bins never span records).
pub fn write_depth_bedgraph<P: AsRef<Path>>(
    path: P,
    depth: &DepthTrack,
    seq: &Sequence,
    bin: usize,
) -> std::io::Result<()> {
    let bin = bin.max(1);
    let mut w = BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        w,
        "track type=bedGraph name=\"dottir depth {}\"",
        depth.method.describe()
    )?;
    for r in &seq.records {
        let mut lines = Vec::new();
        let mut lo = r.range.start;
        while lo < r.range.end {
            let hi = (lo + bin).min(r.range.end);
            let sum: u64 = depth.depth[lo..hi].iter().map(|&d| d as u64).sum();
            let mean = sum as f64 / (hi - lo) as f64;
            let (name, _, s, e, _) = r.source_interval(lo - r.range.start, hi - r.range.start);
            lines.push(format!("{name}\t{s}\t{e}\t{mean:.2}"));
            lo = hi;
        }
        // Reversed records run backwards in the source; bedGraph is sorted.
        if r.origin.as_ref().is_some_and(|o| o.reverse) {
            lines.reverse();
        }
        for l in lines {
            writeln!(w, "{l}")?;
        }
    }
    w.flush()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::parse_fasta;

    fn seq(text: &str) -> Sequence {
        Sequence::from_records(parse_fasta(text).unwrap(), None)
    }

    #[test]
    fn exports_report_source_records_after_sorting() {
        use crate::sequence::RegionOrigin;
        use dottir_core::large::DepthMethod;
        // Query: one 10-base record, reversed by --sort-records.
        let q = seq(">q\nAACCGGTTAC\n").reordered(&[0], &[true]);
        // Subject: a drill-down region, bases 100..108 of a 1000-base "s".
        let mut s = seq(">s:101-108\nACGTACGT\n");
        s.records[0].origin = Some(RegionOrigin {
            record: "s".into(),
            record_len: 1000,
            range: 100..108,
            reverse: false,
        });
        let chain = Chain {
            reverse: false,
            q_rec: 0,
            s_rec: 0,
            q_start: 0,
            q_end: 4,
            s_start: 2,
            s_end: 6,
            score: 400.0,
            n_segments: 1,
            n_anchors: 4,
        };
        let dir = std::env::temp_dir().join(format!("dottir_sorted_export_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let paf = dir.join("c.paf");
        write_chains_paf(&paf, std::slice::from_ref(&chain), &q, &s).unwrap();
        let line = std::fs::read_to_string(&paf).unwrap();
        let f: Vec<&str> = line.trim().split('\t').collect();
        // Displayed 0..4 of the reversed query is source 6..10, minus strand.
        assert_eq!(
            &f[..9],
            ["q", "10", "6", "10", "-", "s", "1000", "102", "106"]
        );
        let bedpe = dir.join("c.bedpe");
        write_chains(&bedpe, std::slice::from_ref(&chain), &q, &s).unwrap();
        assert_eq!(
            std::fs::read_to_string(&bedpe).unwrap(),
            "q\t6\t10\ts\t102\t106\tchain1\t400\t+\t-\n"
        );
        let lines = chain_svg_lines(std::slice::from_ref(&chain), 2, true);
        assert_eq!(lines.len(), 2);
        assert_eq!((lines[0].from, lines[0].to), ((0.0, 1.0), (2.0, 3.0)));
        assert_eq!((lines[1].from, lines[1].to), ((1.0, 0.0), (3.0, 2.0)));
        assert_eq!(lines[0].rgb, CHAIN_FORWARD);

        let depth = DepthTrack {
            method: DepthMethod::None,
            depth: vec![1, 1, 1, 1, 1, 1, 1, 1, 1, 3],
        };
        let bg = dir.join("d.bedgraph");
        write_depth_bedgraph(&bg, &depth, &q, 5).unwrap();
        let text = std::fs::read_to_string(&bg).unwrap();
        let rows: Vec<&str> = text.lines().skip(1).collect();
        // Displayed bins 0..5 and 5..10 are source 5..10 and 0..5; sorted.
        assert_eq!(rows, ["q\t0\t5\t1.40", "q\t5\t10\t1.00"]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn floors_from_soft_mask_and_annotation() {
        let s = seq(">a\nACgtAC\n>b\nacGT\n");
        assert_eq!(
            soft_mask_floor(&s, 9).unwrap(),
            vec![0, 0, 9, 9, 0, 0, 9, 9, 0, 0]
        );
        let gff = "##gff-version 3\nb\t.\tLTR\t2\t3\t.\t+\t.\tName=x\nzz\t.\tLTR\t1\t2\t.\t+\t.\t.\na\t.\tgene\t1\t2\t.\t+\t.\t.\n";
        let annots = AnnotSet {
            source_path: None,
            source: crate::AnnotSource::Gff3,
            features: crate::annotation::parse_gff3_bytes(gff.as_bytes()).unwrap(),
        };
        let (f, skipped) = annotation_floor(&s, &annots, &["ltr".to_string()], 7);
        assert_eq!(f, vec![0, 0, 0, 0, 0, 0, 0, 7, 7, 0]);
        assert_eq!(skipped, 1);
        let m = merge_floors(soft_mask_floor(&s, 9), Some(f)).unwrap();
        assert_eq!(m, vec![0, 0, 9, 9, 0, 0, 9, 9, 7, 0]);
    }

    #[test]
    fn chain_paf_is_record_local() {
        let q = seq(">q1\nACGTACGTAC\n>q2\nACGTACGTACGTACGTACGT\n");
        let s = seq(">s1\nACGTACGTACGTACGTACGTACGTACGTAC\n");
        let chains = vec![Chain {
            reverse: true,
            q_rec: 1,
            s_rec: 0,
            q_start: 12,
            q_end: 25,
            s_start: 3,
            s_end: 17,
            score: 123.4,
            n_segments: 2,
            n_anchors: 9,
        }];
        let dir = std::env::temp_dir().join(format!("dottir-paf-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("c.paf");
        write_chains_paf(&p, &chains, &q, &s).unwrap();
        let text = std::fs::read_to_string(&p).unwrap();
        assert_eq!(
            text,
            "q2\t20\t2\t15\t-\ts1\t30\t3\t17\t13\t14\t255\ttp:A:P\tcm:i:9\ts1:i:123\n"
        );
        std::fs::remove_dir_all(dir).ok();
    }
}
