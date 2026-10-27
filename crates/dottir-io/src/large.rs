//! I/O glue for large-sequence mode (ADR 0007): depth floors from
//! soft-masking and repeat annotations, chain export as PAF, and depth
//! tracks as bedGraph.

use std::collections::BTreeMap;
use std::io::{BufWriter, Write};
use std::path::Path;

use dottir_core::large::{Chain, DepthTrack};

use crate::annotation::AnnotSet;
use crate::sequence::Sequence;

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
    let by_id: BTreeMap<&str, std::ops::Range<usize>> = seq
        .records
        .iter()
        .map(|r| (r.id.as_str(), r.range.clone()))
        .collect();
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
        let Some(rec) = by_id.get(f.record.as_str()) else {
            skipped += 1;
            continue;
        };
        let lo = (rec.start + f.range.start).min(rec.end);
        let hi = (rec.start + f.range.end).min(rec.end);
        floor[lo..hi].fill(value);
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
        let qs = c.q_start as usize - qr.range.start;
        let qe = c.q_end as usize - qr.range.start;
        let ss = c.s_start as usize - sr.range.start;
        let se = c.s_end as usize - sr.range.start;
        let block = (qe - qs).max(se - ss);
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t255\ttp:A:P\tcm:i:{}\ts1:i:{}",
            qr.id,
            qr.len(),
            qs,
            qe,
            if c.reverse { '-' } else { '+' },
            sr.id,
            sr.len(),
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

/// Write a depth track as bedGraph, averaged over `bin`-base bins (record-
/// local coordinates; bins never span records).
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
        let mut lo = r.range.start;
        while lo < r.range.end {
            let hi = (lo + bin).min(r.range.end);
            let sum: u64 = depth.depth[lo..hi].iter().map(|&d| d as u64).sum();
            let mean = sum as f64 / (hi - lo) as f64;
            writeln!(
                w,
                "{}\t{}\t{}\t{:.2}",
                r.id,
                lo - r.range.start,
                hi - r.range.start,
                mean
            )?;
            lo = hi;
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
