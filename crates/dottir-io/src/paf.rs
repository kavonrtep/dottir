//! PAF (minimap2 pairwise mapping format) reader and writer for genome
//! mode (ADR 0008). Plain or gzipped input. Only the 12 mandatory columns
//! and the `tp:A` / `de:f` tags are interpreted (plus `cg:Z` for HSP
//! overlays, [`read_paf_subset_oriented`]); other tags are ignored.
//!
//! There is no noodles PAF crate; the format is 12 tab-separated columns.

use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use dottir_core::genome::alignment::{assign_primary, MASK_LEVEL};
use dottir_core::genome::{Alignment, GenomeRecord};
use flate2::read::MultiGzDecoder;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum PafError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("PAF line {line}: {msg}")]
    Malformed { line: u64, msg: String },
    #[error(
        "PAF line {line}: {side} record {name} has length {paf_len} in the PAF \
         but {fasta_len} in the FASTA — the PAF was made from different sequences"
    )]
    LengthMismatch {
        line: u64,
        side: &'static str,
        name: String,
        paf_len: u64,
        fasta_len: u64,
    },
    #[error(
        "{count} PAF line(s) name {side} records not in the FASTA, e.g. {}",
        examples.join(", ")
    )]
    UnknownRecords {
        side: &'static str,
        count: u64,
        examples: Vec<String>,
    },
}

/// Open a file, transparently gunzipping it when it starts with the gzip
/// magic bytes.
pub fn open_maybe_gz(path: &Path) -> std::io::Result<Box<dyn BufRead>> {
    let mut r = BufReader::new(std::fs::File::open(path)?);
    let gz = r.fill_buf()?.starts_with(&[0x1f, 0x8b]);
    Ok(if gz {
        Box::new(BufReader::new(MultiGzDecoder::new(r)))
    } else {
        Box::new(r)
    })
}

/// Read a PAF file against the query and target record lists.
///
/// Alignments shorter than `min_len` (the longer of the two spans, as
/// the display filter measures) are not kept, so a huge PAF (e.g.
/// minimap2 `-P` on large genomes) need not fit in memory; their record
/// names are still checked.
pub fn read_paf_file(
    path: &Path,
    query: &[GenomeRecord],
    target: &[GenomeRecord],
    min_len: u64,
) -> Result<Vec<Alignment>, PafError> {
    read_paf_impl(open_maybe_gz(path)?, query, target, min_len)
}

/// Read a PAF for a plot with `vertical` records on the query side and
/// `horizontal` ones on the target side. A PAF made the other way round
/// (its query names unknown, but valid when swapped) is read swapped; the
/// flag says so. `min_len` as in [`read_paf_file`].
pub fn read_paf_oriented(
    path: &Path,
    vertical: &[GenomeRecord],
    horizontal: &[GenomeRecord],
    min_len: u64,
) -> Result<(Vec<Alignment>, bool), PafError> {
    match read_paf_file(path, vertical, horizontal, min_len) {
        Ok(a) => Ok((a, false)),
        // Shared names (chr1 in both) fail on lengths, not names.
        Err(first @ (PafError::UnknownRecords { .. } | PafError::LengthMismatch { .. })) => {
            match read_paf_file(path, horizontal, vertical, min_len) {
                Ok(a) => Ok((a.iter().map(Alignment::swapped).collect(), true)),
                Err(_) => Err(first),
            }
        }
        Err(e) => Err(e),
    }
}

struct Names<'a> {
    side: &'static str,
    recs: &'a [GenomeRecord],
    by_name: BTreeMap<&'a str, u32>,
    unknown: u64,
    examples: Vec<String>,
}

impl<'a> Names<'a> {
    fn new(side: &'static str, recs: &'a [GenomeRecord]) -> Self {
        let by_name = recs
            .iter()
            .enumerate()
            .map(|(i, r)| (r.name.as_str(), i as u32))
            .collect();
        Self {
            side,
            recs,
            by_name,
            unknown: 0,
            examples: Vec::new(),
        }
    }

    fn resolve(&mut self, line: u64, name: &str, paf_len: u64) -> Result<Option<u32>, PafError> {
        let Some(&i) = self.by_name.get(name) else {
            self.unknown += 1;
            if self.examples.len() < 5 && !self.examples.iter().any(|e| e == name) {
                self.examples.push(name.to_string());
            }
            return Ok(None);
        };
        let fasta_len = self.recs[i as usize].len;
        if fasta_len != paf_len {
            return Err(PafError::LengthMismatch {
                line,
                side: self.side,
                name: name.to_string(),
                paf_len,
                fasta_len,
            });
        }
        Ok(Some(i))
    }

    fn finish(self) -> Result<(), PafError> {
        if self.unknown == 0 {
            return Ok(());
        }
        Err(PafError::UnknownRecords {
            side: self.side,
            count: self.unknown,
            examples: self.examples,
        })
    }
}

/// Read PAF from `reader`. Record names are resolved against `query` and
/// `target`; unknown names, length disagreements and out-of-range
/// coordinates are errors. Empty lines and lines starting with `#` are
/// skipped.
pub fn read_paf<R: BufRead>(
    reader: R,
    query: &[GenomeRecord],
    target: &[GenomeRecord],
) -> Result<Vec<Alignment>, PafError> {
    read_paf_impl(reader, query, target, 0)
}

/// Like [`read_paf_oriented`], but lines naming records outside the
/// lists are skipped instead of failing (overlays on a dotplot of some of
/// the records). Returns the alignments, the number of lines skipped, and
/// whether the PAF was read swapped: the orientation that places more
/// alignments wins. Also each alignment's [`AlignmentPath`] when it has a
/// valid `cg:Z` CIGAR, in the returned orientation.
pub fn read_paf_subset_oriented(
    path: &Path,
    vertical: &[GenomeRecord],
    horizontal: &[GenomeRecord],
) -> Result<SubsetRead, PafError> {
    // One pass resolving each line both ways round.
    let mut a = Orientation::new(vertical, horizontal);
    let mut b = Orientation::new(horizontal, vertical);
    for_each_line(open_maybe_gz(path)?, |l, line| {
        let p = parse_line(l, line, true)?;
        let path = p.cigar.and_then(|cg| cigar_path(cg, &p.aln));
        a.take(line, &p, &path);
        b.take(line, &p, &path);
        Ok(())
    })?;
    let swap = |b: Orientation<'_>| {
        let (alns, skipped, paths) = b.finish();
        let paths = paths
            .into_iter()
            .map(|p| {
                p.map(|p| {
                    p.into_iter()
                        .map(|[q0, t0, q1, t1]| [t0, q0, t1, q1])
                        .collect()
                })
            })
            .collect();
        Ok((
            alns.iter().map(Alignment::swapped).collect(),
            skipped,
            true,
            paths,
        ))
    };
    match (a.error.take(), b.error.is_none()) {
        (None, true) if b.out.len() > a.out.len() => swap(b),
        (None, _) => {
            let (alns, skipped, paths) = a.finish();
            Ok((alns, skipped, false, paths))
        }
        (Some(_), true) => swap(b),
        (Some(e), false) => Err(e),
    }
}

/// One way of reading a PAF onto the plotted records, for
/// [`read_paf_subset_oriented`]: lines naming other records are skipped,
/// the first length disagreement rules this orientation out.
struct Orientation<'a> {
    qn: Names<'a>,
    tn: Names<'a>,
    out: Vec<Alignment>,
    paths: Vec<Option<AlignmentPath>>,
    skipped: u64,
    error: Option<PafError>,
}

impl<'a> Orientation<'a> {
    fn new(query: &'a [GenomeRecord], target: &'a [GenomeRecord]) -> Self {
        Self {
            qn: Names::new("query", query),
            tn: Names::new("target", target),
            out: Vec::new(),
            paths: Vec::new(),
            skipped: 0,
            error: None,
        }
    }

    fn take(&mut self, line: u64, p: &PafLine<'_>, path: &Option<AlignmentPath>) {
        if self.error.is_some() {
            return;
        }
        let resolved = self
            .qn
            .resolve(line, p.q_name, p.q_len)
            .and_then(|q| Ok((q, self.tn.resolve(line, p.t_name, p.t_len)?)));
        match resolved {
            Err(e) => self.error = Some(e),
            Ok((Some(q_rec), Some(t_rec))) => {
                self.out.push(Alignment {
                    q_rec,
                    t_rec,
                    ..p.aln.clone()
                });
                self.paths.push(path.clone());
            }
            Ok(_) => self.skipped += 1,
        }
    }

    fn finish(mut self) -> PafRead {
        derive_primaries(&mut self.out);
        (self.out, self.skipped, self.paths)
    }
}

/// minimap2 -P marks every chain secondary; derive the primaries it
/// would have chosen so "primary only" still means something.
fn derive_primaries(alns: &mut [Alignment]) {
    if !alns.is_empty() && alns.iter().all(|a| !a.primary) {
        assign_primary(alns, MASK_LEVEL);
    }
}

/// The fields of one PAF line; `aln` has record indices 0 until the
/// names are resolved.
struct PafLine<'l> {
    q_name: &'l str,
    q_len: u64,
    t_name: &'l str,
    t_len: u64,
    aln: Alignment,
    cigar: Option<&'l str>,
}

/// Parse one data line (no trailing newline). `cg:Z` is kept only when
/// `keep_cigar`.
fn parse_line(l: &str, line: u64, keep_cigar: bool) -> Result<PafLine<'_>, PafError> {
    let bad = |msg: String| PafError::Malformed { line, msg };
    let mut fields = l.split('\t');
    let mut col = [""; 12];
    for (i, c) in col.iter_mut().enumerate() {
        *c = fields
            .next()
            .ok_or_else(|| bad(format!("{i} columns, expected at least 12")))?;
    }
    let num = |i: usize| -> Result<u64, PafError> {
        col[i]
            .parse()
            .map_err(|_| bad(format!("column {} is not a number: {:?}", i + 1, col[i])))
    };
    let (q_len, q_start, q_end) = (num(1)?, num(2)?, num(3)?);
    let reverse = match col[4] {
        "+" => false,
        "-" => true,
        s => return Err(bad(format!("strand must be + or -, got {s:?}"))),
    };
    let (t_len, t_start, t_end) = (num(6)?, num(7)?, num(8)?);
    let (matches, block_len) = (num(9)?, num(10)?);
    let mapq = num(11)?.min(255) as u8;
    if q_start > q_end || q_end > q_len || t_start > t_end || t_end > t_len {
        return Err(bad("coordinates out of range".into()));
    }
    let mut primary = true;
    let mut divergence = None;
    let mut cigar = None;
    for tag in fields {
        if let Some(v) = tag.strip_prefix("tp:A:") {
            primary = matches!(v, "P" | "I");
        } else if let Some(v) = tag.strip_prefix("de:f:") {
            divergence = Some(v.parse().map_err(|_| bad(format!("bad de:f tag {v:?}")))?);
        } else if let Some(v) = tag.strip_prefix("cg:Z:").filter(|_| keep_cigar) {
            cigar = Some(v);
        }
    }
    Ok(PafLine {
        q_name: col[0],
        q_len,
        t_name: col[5],
        t_len,
        aln: Alignment {
            q_rec: 0,
            q_start,
            q_end,
            t_rec: 0,
            t_start,
            t_end,
            reverse,
            matches,
            block_len,
            mapq,
            primary,
            divergence,
        },
        cigar,
    })
}

/// Call `f` with each data line of `reader` (empty lines and `#`
/// comments skipped) and its 1-based line number.
fn for_each_line<R: BufRead>(
    mut reader: R,
    mut f: impl FnMut(&str, u64) -> Result<(), PafError>,
) -> Result<(), PafError> {
    let mut buf = String::new();
    let mut line = 0u64;
    loop {
        buf.clear();
        if reader.read_line(&mut buf)? == 0 {
            return Ok(());
        }
        line += 1;
        let l = buf.trim_end_matches(['\n', '\r']);
        if !l.is_empty() && !l.starts_with('#') {
            f(l, line)?;
        }
    }
}

/// Alignments, lines skipped and their paths.
type PafRead = (Vec<Alignment>, u64, Vec<Option<AlignmentPath>>);

/// [`read_paf_subset_oriented`]: alignments, lines skipped, read
/// swapped, paths.
pub type SubsetRead = (Vec<Alignment>, u64, bool, Vec<Option<AlignmentPath>>);

/// Gapless pieces of an alignment from its `cg:Z` CIGAR, as
/// `[q0, t0, q1, t1]` source positions (query descending on the reverse
/// strand). Pieces separated by indels are merged while the diagonal
/// drifts by at most [`PATH_TOLERANCE`] bases.
pub type AlignmentPath = Vec<[u64; 4]>;

/// Diagonal drift (bp) allowed within one piece of an [`AlignmentPath`].
pub const PATH_TOLERANCE: u64 = 10;

/// The [`AlignmentPath`] of `a` from CIGAR `cg`; `None` when the CIGAR is
/// malformed or does not span the alignment.
fn cigar_path(cg: &str, a: &Alignment) -> Option<AlignmentPath> {
    let (mut qp, mut tp) = (0_u64, 0_u64);
    let mut out: AlignmentPath = Vec::new();
    // Current piece, in alignment-relative positions.
    let mut piece: Option<(u64, u64, u64, u64)> = None;
    let at = |qp: u64, tp: u64| {
        let q = if a.reverse {
            a.q_end - qp
        } else {
            a.q_start + qp
        };
        (q, a.t_start + tp)
    };
    let close = |out: &mut AlignmentPath, (q0, t0, q1, t1): (u64, u64, u64, u64)| {
        let (a0, b0) = at(q0, t0);
        let (a1, b1) = at(q1, t1);
        out.push([a0, b0, a1, b1]);
    };
    let mut n = 0_u64;
    for c in cg.bytes() {
        if c.is_ascii_digit() {
            n = n.checked_mul(10)?.checked_add(u64::from(c - b'0'))?;
            continue;
        }
        match c {
            b'M' | b'=' | b'X' => {
                piece = match piece {
                    Some((q0, t0, _, _)) if (qp - q0).abs_diff(tp - t0) <= PATH_TOLERANCE => {
                        Some((q0, t0, qp + n, tp + n))
                    }
                    Some(p) => {
                        close(&mut out, p);
                        Some((qp, tp, qp + n, tp + n))
                    }
                    None => Some((qp, tp, qp + n, tp + n)),
                };
                qp += n;
                tp += n;
            }
            b'I' => qp += n,
            b'D' | b'N' => tp += n,
            _ => return None,
        }
        n = 0;
    }
    if n != 0 || qp != a.q_end - a.q_start || tp != a.t_end - a.t_start {
        return None;
    }
    if let Some(p) = piece {
        close(&mut out, p);
    }
    Some(out)
}

fn read_paf_impl<R: BufRead>(
    reader: R,
    query: &[GenomeRecord],
    target: &[GenomeRecord],
    min_len: u64,
) -> Result<Vec<Alignment>, PafError> {
    let mut parser = PafParser::new(query, target);
    parser.min_len = min_len;
    for_each_line(reader, |l, line| parser.data_line(l, line))?;
    parser.finish()
}

/// Strict PAF reading ([`read_paf`]) one line at a time, for output that
/// arrives as a stream (an aligner's): only the alignments are kept.
pub struct PafParser<'a> {
    qn: Names<'a>,
    tn: Names<'a>,
    out: Vec<Alignment>,
    line: u64,
    min_len: u64,
}

impl<'a> PafParser<'a> {
    pub fn new(query: &'a [GenomeRecord], target: &'a [GenomeRecord]) -> Self {
        Self {
            qn: Names::new("query", query),
            tn: Names::new("target", target),
            out: Vec::new(),
            line: 0,
            min_len: 0,
        }
    }

    /// Read the next line (a trailing newline is ignored; empty lines and
    /// `#` comments are skipped).
    pub fn push(&mut self, l: &str) -> Result<(), PafError> {
        self.line += 1;
        let l = l.trim_end_matches(['\n', '\r']);
        if l.is_empty() || l.starts_with('#') {
            return Ok(());
        }
        self.data_line(l, self.line)
    }

    fn data_line(&mut self, l: &str, line: u64) -> Result<(), PafError> {
        let p = parse_line(l, line, false)?;
        let q = self.qn.resolve(line, p.q_name, p.q_len)?;
        let t = self.tn.resolve(line, p.t_name, p.t_len)?;
        if let (Some(q_rec), Some(t_rec)) = (q, t) {
            let a = Alignment {
                q_rec,
                t_rec,
                ..p.aln
            };
            if a.len() >= self.min_len {
                self.out.push(a);
            }
        }
        Ok(())
    }

    /// The alignments; an error if any line named an unknown record.
    pub fn finish(self) -> Result<Vec<Alignment>, PafError> {
        self.qn.finish()?;
        self.tn.finish()?;
        let mut out = self.out;
        derive_primaries(&mut out);
        Ok(out)
    }
}

/// Write alignments as PAF: 12 columns, `tp:A:P|S`, and `de:f` when known.
pub fn write_paf<W: Write>(
    mut w: W,
    alns: &[Alignment],
    query: &[GenomeRecord],
    target: &[GenomeRecord],
) -> std::io::Result<()> {
    for a in alns {
        let q = &query[a.q_rec as usize];
        let t = &target[a.t_rec as usize];
        write!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\ttp:A:{}",
            q.name,
            q.len,
            a.q_start,
            a.q_end,
            if a.reverse { '-' } else { '+' },
            t.name,
            t.len,
            a.t_start,
            a.t_end,
            a.matches,
            a.block_len,
            a.mapq,
            if a.primary { 'P' } else { 'S' }
        )?;
        if let Some(de) = a.divergence {
            write!(w, "\tde:f:{de:.4}")?;
        }
        writeln!(w)?;
    }
    w.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn recs(v: &[(&str, u64)]) -> Vec<GenomeRecord> {
        v.iter()
            .map(|&(n, len)| GenomeRecord {
                name: n.to_string(),
                len,
            })
            .collect()
    }

    const PAF: &str = "\
# comment
q1\t1000\t10\t500\t+\tt1\t2000\t100\t590\t450\t490\t60\ttp:A:P\tcm:i:40
q2\t500\t0\t500\t-\tt1\t2000\t1000\t1500\t480\t500\t0\ttp:A:S\tde:f:0.0123

";

    #[test]
    fn orientation_detected_with_shared_record_names() {
        // Review 2026-10-28 #6: query chr1 is the 100-base horizontal one.
        let d = std::env::temp_dir().join(format!("dottir_paf_orient_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let p = d.join("a.paf");
        std::fs::write(&p, "chr1\t100\t0\t80\t+\tchr1\t80\t0\t80\t78\t80\t60\n").unwrap();
        let vertical = recs(&[("chr1", 80)]);
        let horizontal = recs(&[("chr1", 100)]);
        let (a, swapped) = read_paf_oriented(&p, &vertical, &horizontal, 0).unwrap();
        assert!(swapped);
        assert_eq!((a[0].q_end, a[0].t_end), (80, 80));
        assert_eq!(
            (a[0].q_start, a[0].t_start, a[0].q_rec, a[0].t_rec),
            (0, 0, 0, 0)
        );
        // Valid neither way: the first error is kept.
        let other = recs(&[("chr1", 90)]);
        assert!(matches!(
            read_paf_oriented(&p, &other, &other, 0),
            Err(PafError::LengthMismatch { .. })
        ));
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn short_alignments_are_dropped_on_read_but_names_still_checked() {
        let d = std::env::temp_dir().join(format!("dottir_paf_minlen_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        let p = d.join("a.paf");
        std::fs::write(
            &p,
            "q\t1000\t0\t900\t+\tt\t1000\t0\t880\t850\t900\t60\n\
             q\t1000\t0\t50\t+\tt\t1000\t500\t560\t50\t60\t60\n",
        )
        .unwrap();
        let (q, t) = (recs(&[("q", 1000)]), recs(&[("t", 1000)]));
        assert_eq!(read_paf_file(&p, &q, &t, 0).unwrap().len(), 2);
        // The second is 60 long on the target side.
        assert_eq!(read_paf_file(&p, &q, &t, 60).unwrap().len(), 2);
        assert_eq!(read_paf_file(&p, &q, &t, 61).unwrap().len(), 1);
        std::fs::write(&p, "q\t1000\t0\t50\t+\tzz\t1000\t0\t50\t50\t50\t60\n").unwrap();
        assert!(read_paf_file(&p, &q, &t, 1000).is_err(), "unknown name");
        let _ = std::fs::remove_dir_all(&d);
    }

    #[test]
    fn reads_columns_and_tags_and_round_trips() {
        let q = recs(&[("q1", 1000), ("q2", 500)]);
        let t = recs(&[("t1", 2000)]);
        let a = read_paf(PAF.as_bytes(), &q, &t).unwrap();
        assert_eq!(a.len(), 2);
        assert_eq!((a[0].q_rec, a[0].q_start, a[0].q_end), (0, 10, 500));
        assert!(a[0].primary && !a[0].reverse && a[0].divergence.is_none());
        assert!(!a[1].primary && a[1].reverse);
        assert!((a[1].identity() - 0.9877).abs() < 1e-4);
        let mut out = Vec::new();
        write_paf(&mut out, &a, &q, &t).unwrap();
        let back = read_paf(&out[..], &q, &t).unwrap();
        assert_eq!(back, a);
    }

    #[test]
    fn errors_name_the_problem() {
        let q = recs(&[("q1", 1000)]);
        let t = recs(&[("t1", 2000)]);
        let e = read_paf(PAF.as_bytes(), &q, &t).unwrap_err();
        assert!(
            matches!(e, PafError::UnknownRecords { count: 1, .. }),
            "{e}"
        );
        let wrong = recs(&[("q1", 999), ("q2", 500)]);
        let e = read_paf(PAF.as_bytes(), &wrong, &t).unwrap_err();
        assert!(matches!(e, PafError::LengthMismatch { line: 2, .. }), "{e}");
        let e = read_paf("q1\t1000\t10\n".as_bytes(), &q, &t).unwrap_err();
        assert!(matches!(e, PafError::Malformed { line: 1, .. }), "{e}");
        let oob = "q1\t1000\t10\t1500\t+\tt1\t2000\t0\t10\t5\t10\t0\n";
        assert!(read_paf(oob.as_bytes(), &q, &t).is_err());
    }
}
