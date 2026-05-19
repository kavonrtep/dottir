//! `dottir find-peaks` subcommand.
//!
//! Reads a `dottir periodogram` TSV (either the main periodogram
//! output or its `--fft` companion — auto-detected by the column
//! header) and writes a classified-peaks TSV: one row per detected
//! fundamental / harmonic / sub-repeat, per FASTA record.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};

use dottir_core::{
    find_peaks_in_periodogram, find_peaks_in_spectrum, Peak, PeakKind, PeaksConfig, Spectrum,
    SubrepeatConfig,
};

#[derive(Parser, Debug)]
#[command(version)]
pub struct FindPeaksArgs {
    /// Input TSV — either a `dottir periodogram` output or its
    /// `--fft` companion. The format is auto-detected from the
    /// column header.
    #[arg(value_name = "INPUT")]
    pub input: PathBuf,

    /// Output TSV. Default: stdout.
    #[arg(short = 'o', long, value_name = "PATH")]
    pub output: Option<PathBuf>,

    // ── Periodogram-only knobs ───────────────────────────────────
    /// Periodogram only: column to rank peaks by. Default `signal_mean`
    /// — robust against the analytical z-score bias that uniformly
    /// inflates noise across all offsets. Pass `z_score` only when
    /// the source used `--z-score empirical`.
    #[arg(long, value_enum, default_value_t = RankByArg::SignalMean)]
    pub rank_by: RankByArg,

    /// Minimum score for a peak candidate. Defaults:
    /// 10 for `--rank-by signal_mean` (the default), 5 for
    /// `--rank-by z_score`.
    #[arg(long, value_name = "FLOAT")]
    pub min_score: Option<f64>,

    /// Drop peaks at `period > boundary_fraction × max_period` —
    /// kernel edge artifacts saturate at `k ≈ N/2`. Default 0.9;
    /// pass 1.0 to disable.
    #[arg(long, default_value_t = 0.9, value_name = "0.0-1.0")]
    pub boundary_fraction: f64,

    // ── Shared classification knobs ──────────────────────────────
    /// Relative tolerance for integer-ratio harmonic detection.
    /// Default 0.02 (2%).
    #[arg(long, default_value_t = 0.02, value_name = "FLOAT")]
    pub harmonic_tolerance: f64,

    /// Largest integer harmonic to consider. Default 30. Long
    /// tandem arrays show clean harmonics out to n>20.
    #[arg(long, default_value_t = 30, value_name = "INT")]
    pub max_harmonic_n: u32,

    /// Keep only fundamentals with at least this many detected
    /// harmonics. Default 1 — a "fundamental" with zero harmonics
    /// is almost always noise.
    #[arg(long, default_value_t = 1, value_name = "INT")]
    pub min_harmonics: u32,

    /// Limit output to the top-N highest-scoring peaks per record.
    /// 0 = no limit (default).
    #[arg(long, default_value_t = 0, value_name = "INT")]
    pub top_n: usize,

    /// Include harmonics in the output. By default only
    /// fundamentals (and sub-repeats, if `--subrepeats`) are
    /// emitted — harmonics clutter the view for normal use.
    #[arg(long)]
    pub show_harmonics: bool,

    // ── Sub-repeat scan (periodogram only) ───────────────────────
    /// Enable sub-repeat detection: for each fundamental, scan
    /// integer divisors at fundamental/n with ±tolerance bp and
    /// surface as `subrepeat` any peak found there.
    #[arg(long)]
    pub subrepeats: bool,

    /// Minimum score for a sub-repeat (typically lower than
    /// `--min-score`). Default 5.
    #[arg(long, default_value_t = 5.0, value_name = "FLOAT")]
    pub subrepeat_min_score: f64,

    /// Largest integer divisor to scan for sub-repeats. Default 6.
    #[arg(long, default_value_t = 6, value_name = "INT")]
    pub max_divisor: u32,

    /// `±` bp tolerance when matching sub-repeat candidates around
    /// `fundamental / n`. Default 2.
    #[arg(long, default_value_t = 2, value_name = "INT")]
    pub period_tolerance: u32,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum RankByArg {
    SignalMean,
    ZScore,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum InputFormat {
    Periodogram,
    Fft,
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

pub fn run(args: FindPeaksArgs) -> Result<()> {
    let (fmt, header, rows) =
        load_tsv(&args.input).with_context(|| format!("reading {}", args.input.display()))?;

    let min_score = args.min_score.unwrap_or(match args.rank_by {
        RankByArg::SignalMean => 10.0,
        RankByArg::ZScore => 5.0,
    });
    let cfg = PeaksConfig {
        min_score,
        harmonic_tolerance: args.harmonic_tolerance,
        max_harmonic_n: args.max_harmonic_n,
        min_harmonics: args.min_harmonics,
        boundary_fraction: args.boundary_fraction,
        subrepeats: if args.subrepeats && fmt == InputFormat::Periodogram {
            Some(SubrepeatConfig {
                min_score: args.subrepeat_min_score,
                max_divisor: args.max_divisor,
                period_tolerance: args.period_tolerance,
            })
        } else {
            None
        },
    };

    let mut out: Box<dyn Write> = match args.output.as_ref() {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).with_context(|| format!("creating {}", p.display()))?,
        )),
        None => Box::new(BufWriter::new(std::io::stdout().lock())),
    };

    write_header(&mut out, &args, fmt, &cfg)?;
    writeln!(
        out,
        "record_id\tperiod_bp\tscore\tkind\tparent\tharmonic_n\tdivisor_n\tn_harmonics\tharmonics"
    )?;

    let records = group_by_record(rows, &header)?;
    for (rid, record_rows) in records {
        let peaks = match fmt {
            InputFormat::Periodogram => {
                process_periodogram_record(&record_rows, &header, args.rank_by, &cfg)?
            }
            InputFormat::Fft => process_fft_record(&record_rows, &header, &cfg)?,
        };
        // Default output: fundamentals + sub-repeats; harmonics
        // hidden unless --show-harmonics. Sub-repeats are emitted
        // right after their parent fundamental for visual scanning.
        let visible: Vec<&Peak> = if args.show_harmonics {
            peaks.iter().collect()
        } else {
            peaks
                .iter()
                .filter(|p| !matches!(p.kind, PeakKind::Harmonic))
                .collect()
        };
        // Re-order: each fundamental immediately followed by its
        // sub-repeats. Harmonics (if shown) preserve their natural
        // score-descending order, appended after the last fundamental.
        let mut ordered: Vec<&Peak> = Vec::with_capacity(visible.len());
        let fundamentals: Vec<&Peak> = visible
            .iter()
            .copied()
            .filter(|p| p.kind == PeakKind::Fundamental)
            .collect();
        for fund in &fundamentals {
            ordered.push(*fund);
            for sub in visible
                .iter()
                .copied()
                .filter(|p| p.kind == PeakKind::Subrepeat && p.parent_period == Some(fund.period))
            {
                ordered.push(sub);
            }
        }
        if args.show_harmonics {
            for harm in visible
                .iter()
                .copied()
                .filter(|p| p.kind == PeakKind::Harmonic)
            {
                ordered.push(harm);
            }
        }
        let limited: Vec<&Peak> = if args.top_n > 0 {
            ordered.into_iter().take(args.top_n).collect()
        } else {
            ordered
        };
        for p in &limited {
            emit_peak(&mut out, &rid, p)?;
        }
    }

    out.flush().context("flushing output")?;
    if let Some(p) = args.output.as_ref() {
        tracing::info!("wrote {}", p.display());
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// TSV loading + format detection
// ---------------------------------------------------------------------------

/// One TSV data row as `(column_name -> value)`.
type Row = std::collections::HashMap<String, String>;

/// Load a TSV, skipping `#` comment lines. Returns
/// `(format, column_header, rows)`.
fn load_tsv(path: &std::path::Path) -> Result<(InputFormat, Vec<String>, Vec<Row>)> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let mut header: Vec<String> = Vec::new();
    let mut rows: Vec<Row> = Vec::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        if header.is_empty() {
            header = line.split('\t').map(|s| s.to_string()).collect();
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != header.len() {
            continue; // tolerate ragged rows
        }
        let mut row: Row = std::collections::HashMap::with_capacity(header.len());
        for (k, v) in header.iter().zip(parts.iter()) {
            row.insert(k.clone(), v.to_string());
        }
        rows.push(row);
    }
    let fmt = detect_format(&header)?;
    Ok((fmt, header, rows))
}

fn detect_format(header: &[String]) -> Result<InputFormat> {
    let cols: std::collections::HashSet<&str> = header.iter().map(|s| s.as_str()).collect();
    let periodogram_required = ["record_id", "k", "signal_mean", "z_score"];
    let fft_required = [
        "record_id",
        "bin",
        "period_residues",
        "amplitude",
        "peak_rank",
    ];
    if periodogram_required.iter().all(|c| cols.contains(c)) {
        Ok(InputFormat::Periodogram)
    } else if fft_required.iter().all(|c| cols.contains(c)) {
        Ok(InputFormat::Fft)
    } else {
        anyhow::bail!(
            "unrecognised TSV columns: {:?}\n\
             expected periodogram ({:?}) or FFT ({:?})",
            header,
            periodogram_required,
            fft_required,
        )
    }
}

/// Group rows by `record_id`, preserving input order both within
/// and across records.
fn group_by_record(rows: Vec<Row>, _header: &[String]) -> Result<Vec<(String, Vec<Row>)>> {
    let mut order: Vec<String> = Vec::new();
    let mut groups: std::collections::HashMap<String, Vec<Row>> = std::collections::HashMap::new();
    for row in rows {
        let rid = row
            .get("record_id")
            .cloned()
            .ok_or_else(|| anyhow::anyhow!("row missing record_id"))?;
        if !groups.contains_key(&rid) {
            order.push(rid.clone());
        }
        groups.entry(rid).or_default().push(row);
    }
    Ok(order
        .into_iter()
        .map(|rid| {
            let rows = groups.remove(&rid).unwrap();
            (rid, rows)
        })
        .collect())
}

// ---------------------------------------------------------------------------
// Per-record processing
// ---------------------------------------------------------------------------

fn process_periodogram_record(
    rows: &[Row],
    _header: &[String],
    rank_by: RankByArg,
    cfg: &PeaksConfig,
) -> Result<Vec<Peak>> {
    let col = match rank_by {
        RankByArg::SignalMean => "signal_mean",
        RankByArg::ZScore => "z_score",
    };
    // Build the per-offset score array. The TSV is already sorted
    // by k within a record (it's how dottir periodogram writes it),
    // and the offsets are dense — `k = min_offset + i`.
    if rows.is_empty() {
        return Ok(Vec::new());
    }
    let mut ks: Vec<u32> = Vec::with_capacity(rows.len());
    let mut scores: Vec<f64> = Vec::with_capacity(rows.len());
    for r in rows {
        let k: u32 = r
            .get("k")
            .ok_or_else(|| anyhow::anyhow!("row missing k"))?
            .parse()
            .context("parsing k")?;
        let s: f64 = r
            .get(col)
            .ok_or_else(|| anyhow::anyhow!("row missing {col}"))?
            .parse()
            .context(format!("parsing {col}"))?;
        ks.push(k);
        scores.push(s);
    }
    let min_offset = ks[0];
    // Sanity check that offsets are dense.
    for (i, &k) in ks.iter().enumerate() {
        if k != min_offset + i as u32 {
            anyhow::bail!(
                "non-dense offsets in periodogram TSV at row {i}: \
                 expected k={}, got k={k}",
                min_offset + i as u32
            );
        }
    }
    Ok(find_peaks_in_periodogram(&scores, min_offset, cfg)?)
}

fn process_fft_record(rows: &[Row], _header: &[String], cfg: &PeaksConfig) -> Result<Vec<Peak>> {
    if rows.is_empty() {
        return Ok(Vec::new());
    }
    // Build a Spectrum struct from the FFT TSV rows. The TSV is
    // already sorted by bin (dottir periodogram's --fft writer
    // iterates 0..bin_count).
    let n = rows.len();
    let mut amplitude: Vec<f64> = Vec::with_capacity(n);
    let mut peak_ranks: Vec<Option<u32>> = Vec::with_capacity(n);
    let mut last_bin: i64 = -1;
    for r in rows {
        let bin: u32 = r
            .get("bin")
            .ok_or_else(|| anyhow::anyhow!("row missing bin"))?
            .parse()
            .context("parsing bin")?;
        if bin as i64 != last_bin + 1 {
            anyhow::bail!(
                "non-dense bins in FFT TSV: expected {}, got {}",
                last_bin + 1,
                bin
            );
        }
        last_bin = bin as i64;
        let amp: f64 = r
            .get("amplitude")
            .ok_or_else(|| anyhow::anyhow!("row missing amplitude"))?
            .parse()
            .context("parsing amplitude")?;
        amplitude.push(amp);
        let rank_str = r.get("peak_rank").cloned().unwrap_or_default();
        let rank = if rank_str.trim().is_empty() {
            None
        } else {
            Some(rank_str.parse::<u32>().context("parsing peak_rank")?)
        };
        peak_ranks.push(rank);
    }
    // padded_len = 2 × (bin_count - 1) where bin_count = n. The
    // Spectrum::period helper computes `padded_len / bin`, so we
    // need the right value. For a real-FFT output of length L
    // padded to M, the spectrum has M/2 + 1 bins → bin_count = M/2 + 1
    // → M = 2 × (bin_count - 1).
    let padded_len = 2 * (n.saturating_sub(1)).max(1);
    let spectrum = Spectrum {
        input_len: padded_len, // unused for find_peaks_in_spectrum
        padded_len,
        amplitude,
        peak_ranks,
    };
    Ok(find_peaks_in_spectrum(&spectrum, cfg)?)
}

// ---------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------

fn write_header(
    out: &mut dyn Write,
    args: &FindPeaksArgs,
    fmt: InputFormat,
    cfg: &PeaksConfig,
) -> std::io::Result<()> {
    writeln!(out, "# dottir find-peaks v{}", env!("CARGO_PKG_VERSION"))?;
    writeln!(out, "# input: {}", args.input.display())?;
    writeln!(
        out,
        "# input_format: {}",
        match fmt {
            InputFormat::Periodogram => "periodogram",
            InputFormat::Fft => "fft",
        }
    )?;
    if fmt == InputFormat::Periodogram {
        writeln!(
            out,
            "# rank_by: {}",
            match args.rank_by {
                RankByArg::SignalMean => "signal_mean",
                RankByArg::ZScore => "z_score",
            }
        )?;
        writeln!(out, "# boundary_fraction: {}", cfg.boundary_fraction)?;
    }
    writeln!(out, "# min_score: {}", cfg.min_score)?;
    writeln!(out, "# harmonic_tolerance: {}", cfg.harmonic_tolerance)?;
    writeln!(out, "# max_harmonic_n: {}", cfg.max_harmonic_n)?;
    writeln!(out, "# min_harmonics: {}", cfg.min_harmonics)?;
    if let Some(sub) = cfg.subrepeats {
        writeln!(out, "# subrepeats: true")?;
        writeln!(out, "# subrepeat_min_score: {}", sub.min_score)?;
        writeln!(out, "# max_divisor: {}", sub.max_divisor)?;
        writeln!(out, "# period_tolerance: {}", sub.period_tolerance)?;
    }
    if args.top_n > 0 {
        writeln!(out, "# top_n: {}", args.top_n)?;
    }
    writeln!(out, "# show_harmonics: {}", args.show_harmonics)?;
    Ok(())
}

fn emit_peak(out: &mut dyn Write, rid: &str, p: &Peak) -> std::io::Result<()> {
    let period_s = if p.period.fract().abs() < 1e-9 {
        format!("{}", p.period as i64)
    } else {
        format!("{:.2}", p.period)
    };
    let kind_s = match p.kind {
        PeakKind::Fundamental => "fundamental",
        PeakKind::Harmonic => "harmonic",
        PeakKind::Subrepeat => "subrepeat",
    };
    let parent_s = p
        .parent_period
        .map(|v| {
            if v.fract().abs() < 1e-9 {
                format!("{}", v as i64)
            } else {
                format!("{:.2}", v)
            }
        })
        .unwrap_or_default();
    let harm_n_s = p.harmonic_n.map(|n| n.to_string()).unwrap_or_default();
    let div_n_s = p.divisor_n.map(|n| n.to_string()).unwrap_or_default();
    let n_harm_s = if p.kind == PeakKind::Fundamental {
        p.n_harmonics.to_string()
    } else {
        String::new()
    };
    let harm_list_s = if p.kind == PeakKind::Fundamental {
        p.harmonics
            .iter()
            .map(|n| n.to_string())
            .collect::<Vec<_>>()
            .join(",")
    } else {
        String::new()
    };
    writeln!(
        out,
        "{rid}\t{period_s}\t{:.4}\t{kind_s}\t{parent_s}\t{harm_n_s}\t{div_n_s}\t{n_harm_s}\t{harm_list_s}",
        p.score
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn make_tsv(name: &str, content: &str) -> std::path::PathBuf {
        // Per-test name avoids parallel-test races on a shared path.
        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "dottir_find_peaks_test_{}_{}.tsv",
            std::process::id(),
            name
        ));
        let mut f = File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn detects_periodogram_format() {
        let path = make_tsv(
            "periodogram",
            "# header\nrecord_id\tk\traw_sum\tsignal_sum\tsignal_mean\tz_score\n\
             r1\t3\t0\t0\t0.0\t0.0\n",
        );
        let (fmt, _, _) = load_tsv(&path).unwrap();
        assert_eq!(fmt, InputFormat::Periodogram);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn detects_fft_format() {
        let path = make_tsv(
            "fft",
            "# header\nrecord_id\tbin\tfrequency\tperiod_residues\tamplitude\tpeak_rank\n\
             r1\t0\t0.0\tinf\t0.0\t\n",
        );
        let (fmt, _, _) = load_tsv(&path).unwrap();
        assert_eq!(fmt, InputFormat::Fft);
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn rejects_unknown_format() {
        let path = make_tsv("unknown", "record_id\tfoo\tbar\nr1\t1\t2\n");
        let result = load_tsv(&path);
        assert!(result.is_err());
        std::fs::remove_file(&path).ok();
    }
}
