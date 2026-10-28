//! `--export-alignment-at`: the residues around a dotplot position, as a
//! FASTA pair, Stockholm block or plain text (spec §4.4.5, plan D5).

use anyhow::{bail, Context, Result};
use clap::{Args, ValueEnum};
use dottir_io::alignment::{slice_pair, SliceFormat};
use dottir_io::Sequence;

#[derive(Args, Debug, Clone)]
pub struct AlignExportArgs {
    /// Write the residues around query position Q and subject position S
    /// to stdout. Each is 1-based along the whole input (`1500`) or within
    /// a record (`chr2:1500`). With no `-o`, nothing is plotted.
    #[arg(long, value_name = "Q,S")]
    pub export_alignment_at: Option<String>,

    /// Residues on each side of the position.
    #[arg(
        long,
        value_name = "N",
        default_value_t = 30,
        requires = "export_alignment_at"
    )]
    pub alignment_window: usize,

    /// Output format of `--export-alignment-at`.
    #[arg(long, value_enum, default_value_t = AlignFormat::Fasta, requires = "export_alignment_at")]
    pub alignment_format: AlignFormat,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum AlignFormat {
    Fasta,
    Stockholm,
    Text,
}

impl AlignFormat {
    fn slice_format(self) -> SliceFormat {
        match self {
            AlignFormat::Fasta => SliceFormat::FastaPair,
            AlignFormat::Stockholm => SliceFormat::Stockholm,
            AlignFormat::Text => SliceFormat::Plain,
        }
    }
}

/// Record index and 0-based position within it for a 1-based `spec`.
fn resolve(seq: &Sequence, spec: &str, what: &str) -> Result<(usize, usize)> {
    let spec = spec.trim();
    // Record names may contain ':' (drill-down regions), so split last.
    let (rec, pos) = match spec.rsplit_once(':') {
        Some((name, pos)) => {
            let rec = seq
                .records
                .iter()
                .position(|r| r.id == name)
                .with_context(|| format!("{what}: no record named {name:?}"))?;
            (Some(rec), pos)
        }
        None => (None, spec),
    };
    let pos: usize = pos
        .replace([',', '_'], "")
        .parse()
        .with_context(|| format!("{what}: bad position {pos:?}"))?;
    if pos == 0 {
        bail!("{what}: positions are 1-based");
    }
    match rec {
        Some(i) => {
            let len = seq.records[i].len();
            if pos > len {
                bail!(
                    "{what}: {pos} is past the end of {} ({len} bp)",
                    seq.records[i].id
                );
            }
            Ok((i, pos - 1))
        }
        None => {
            let (r, local) = seq
                .record_at(pos - 1)
                .with_context(|| format!("{what}: {pos} is past the end ({} bp)", seq.len()))?;
            let i = seq
                .records
                .iter()
                .position(|x| std::ptr::eq(x, r))
                .expect("own record");
            Ok((i, local))
        }
    }
}

/// The exported text. The window stays inside the records holding the
/// two positions.
pub fn export(
    args: &AlignExportArgs,
    at: &str,
    query: &Sequence,
    subject: &Sequence,
) -> Result<String> {
    let (q, s) = at
        .split_once(',')
        .with_context(|| format!("--export-alignment-at wants Q,S, got {at:?}"))?;
    let (qi, qpos) = resolve(query, q, "query")?;
    let (si, spos) = resolve(subject, s, "subject")?;
    let (qr, sr) = (&query.records[qi], &subject.records[si]);
    let sliced = slice_pair(
        &query.bytes()[qr.range.clone()],
        &subject.bytes()[sr.range.clone()],
        qpos,
        spos,
        args.alignment_window,
        args.alignment_format.slice_format(),
        &qr.id,
        &sr.id,
    );
    Ok(sliced.text)
}

#[cfg(test)]
mod tests {
    use super::*;
    use dottir_io::fasta::parse_fasta;

    fn seq(text: &str) -> Sequence {
        Sequence::from_records(parse_fasta(text).unwrap(), None)
    }

    fn args(window: usize, format: AlignFormat) -> AlignExportArgs {
        AlignExportArgs {
            export_alignment_at: None,
            alignment_window: window,
            alignment_format: format,
        }
    }

    #[test]
    fn positions_by_record_or_global_stay_in_their_record() {
        let q = seq(">a\nAAAAAAAAAA\n>b\nCCCCGGGGTT\n");
        let s = seq(">x\nCCCCGGGGTT\n");
        let a = args(2, AlignFormat::Text);
        let by_name = export(&a, "b:5,x:5", &q, &s).unwrap();
        let global = export(&a, "15,5", &q, &s).unwrap();
        assert_eq!(by_name, global);
        assert!(by_name.contains("CCGGG"), "{by_name}");
        // At a record start the window does not reach into record a.
        let edge = export(&a, "b:1,x:1", &q, &s).unwrap();
        assert!(edge.contains("CCC") && !edge.contains('A'), "{edge}");
        let fasta = export(&args(2, AlignFormat::Fasta), "b:5,x:5", &q, &s).unwrap();
        assert!(fasta.starts_with(">b:3..7"), "{fasta}");
    }

    #[test]
    fn bad_positions_are_errors() {
        let q = seq(">a\nACGT\n");
        let a = args(2, AlignFormat::Text);
        assert!(export(&a, "0,1", &q, &q).is_err());
        assert!(export(&a, "5,1", &q, &q).is_err());
        assert!(export(&a, "zz:1,1", &q, &q).is_err());
        assert!(export(&a, "1", &q, &q).is_err());
    }
}
