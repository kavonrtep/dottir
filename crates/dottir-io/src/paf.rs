//! PAF (minimap2 pairwise mapping format) reader and writer for genome
//! mode (ADR 0008). Plain or gzipped input. Only the 12 mandatory columns
//! and the `tp:A` / `de:f` tags are interpreted; other tags are ignored.
//!
//! There is no noodles PAF crate; the format is 12 tab-separated columns.

use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

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
pub fn read_paf_file(
    path: &Path,
    query: &[GenomeRecord],
    target: &[GenomeRecord],
) -> Result<Vec<Alignment>, PafError> {
    read_paf(open_maybe_gz(path)?, query, target)
}

/// Read a PAF for a plot with `vertical` records on the query side and
/// `horizontal` ones on the target side. A PAF made the other way round
/// (its query names unknown, but valid when swapped) is read swapped; the
/// flag says so.
pub fn read_paf_oriented(
    path: &Path,
    vertical: &[GenomeRecord],
    horizontal: &[GenomeRecord],
) -> Result<(Vec<Alignment>, bool), PafError> {
    match read_paf_file(path, vertical, horizontal) {
        Ok(a) => Ok((a, false)),
        Err(first @ PafError::UnknownRecords { .. }) => {
            match read_paf_file(path, horizontal, vertical) {
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
    mut reader: R,
    query: &[GenomeRecord],
    target: &[GenomeRecord],
) -> Result<Vec<Alignment>, PafError> {
    let mut qn = Names::new("query", query);
    let mut tn = Names::new("target", target);
    let mut out = Vec::new();
    let mut buf = String::new();
    let mut line = 0u64;
    loop {
        buf.clear();
        if reader.read_line(&mut buf)? == 0 {
            break;
        }
        line += 1;
        let l = buf.trim_end_matches(['\n', '\r']);
        if l.is_empty() || l.starts_with('#') {
            continue;
        }
        let bad = |msg: &str| PafError::Malformed {
            line,
            msg: msg.to_string(),
        };
        let cols: Vec<&str> = l.split('\t').collect();
        if cols.len() < 12 {
            return Err(bad(&format!(
                "{} columns, expected at least 12",
                cols.len()
            )));
        }
        let num = |i: usize| -> Result<u64, PafError> {
            cols[i]
                .parse()
                .map_err(|_| bad(&format!("column {} is not a number: {:?}", i + 1, cols[i])))
        };
        let (q_len, q_start, q_end) = (num(1)?, num(2)?, num(3)?);
        let reverse = match cols[4] {
            "+" => false,
            "-" => true,
            s => return Err(bad(&format!("strand must be + or -, got {s:?}"))),
        };
        let (t_len, t_start, t_end) = (num(6)?, num(7)?, num(8)?);
        let (matches, block_len) = (num(9)?, num(10)?);
        let mapq = num(11)?.min(255) as u8;
        if q_start > q_end || q_end > q_len || t_start > t_end || t_end > t_len {
            return Err(bad("coordinates out of range"));
        }
        let mut primary = true;
        let mut divergence = None;
        for tag in &cols[12..] {
            if let Some(v) = tag.strip_prefix("tp:A:") {
                primary = matches!(v, "P" | "I");
            } else if let Some(v) = tag.strip_prefix("de:f:") {
                divergence = Some(v.parse().map_err(|_| bad(&format!("bad de:f tag {v:?}")))?);
            }
        }
        let q = qn.resolve(line, cols[0], q_len)?;
        let t = tn.resolve(line, cols[5], t_len)?;
        let (Some(q_rec), Some(t_rec)) = (q, t) else {
            continue;
        };
        out.push(Alignment {
            q_rec,
            q_start,
            q_end,
            t_rec,
            t_start,
            t_end,
            reverse,
            matches,
            block_len,
            mapq,
            primary,
            divergence,
        });
    }
    qn.finish()?;
    tn.finish()?;
    Ok(out)
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
