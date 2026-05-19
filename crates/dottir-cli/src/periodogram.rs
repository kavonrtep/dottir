//! `dottir periodogram` subcommand.
//!
//! Loads a DNA FASTA, runs a self-comparison periodogram per record,
//! writes a single TSV with the per-offset signal and z-score.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Parser, ValueEnum};

use dottir_core::{
    analytical_null, analytical_z_scores, compute_periodogram, compute_periodogram_parallel,
    empirical_null_stats, AnalyticalNull, Periodogram, PeriodogramConfig, ScoreMatrix, Sensitivity,
};
use dottir_io::fasta;

#[derive(Parser, Debug)]
#[command(version)]
pub struct PeriodogramArgs {
    /// DNA FASTA. Multi-record inputs produce one periodogram block
    /// per record in the same output TSV.
    #[arg(value_name = "INPUT")]
    pub input: PathBuf,

    /// Output TSV. Default: stdout.
    #[arg(short = 'o', long, value_name = "PATH")]
    pub output: Option<PathBuf>,

    /// DNA match score for the built-in identity-style matrix.
    /// Combine with `--mismatch-score` to build a custom matrix.
    #[arg(long, default_value_t = 5, value_name = "INT")]
    pub match_score: i32,

    /// DNA mismatch score. Defaults to dotter's `-4`.
    #[arg(long, default_value_t = -4, value_name = "INT")]
    pub mismatch_score: i32,

    /// Window size W in residues. Default: per-record Karlin/Altschul
    /// estimate, clamped to `[3, 50]`.
    #[arg(short = 'W', long, value_name = "INT")]
    pub window: Option<u32>,

    /// Pixel factor for the `min(255, score * pixel_fac / W)` step.
    /// `0` = per-record Karlin auto (`0.2 * 256 / E[M]`).
    #[arg(long, default_value_t = 0, value_name = "INT")]
    pub pixel_fac: u32,

    /// Greyramp noise floor (0..=255). Pixel values ≤ this contribute
    /// 0 to the periodogram. Default matches GUI greyramp (40).
    #[arg(long, default_value_t = 40, value_name = "0-255")]
    pub greyramp_white: u8,

    /// Greyramp saturation (0..=255). Pixel values ≥ this contribute
    /// 255. Default matches GUI greyramp (100).
    #[arg(long, default_value_t = 100, value_name = "0-255")]
    pub greyramp_black: u8,

    /// Smallest offset k to report. Default 3 (skips main diagonal and
    /// trivial k=1, k=2 homopolymer spikes).
    #[arg(long, default_value_t = 3, value_name = "INT")]
    pub min_offset: u32,

    /// Largest offset k to report per record. Default: floor(N / 2).
    #[arg(long, value_name = "INT")]
    pub max_offset: Option<u32>,

    /// z-score policy.
    ///
    /// * `auto` (default): analytical when sensitivity is identity
    ///   (`--greyramp-white 0 --greyramp-black 255`), empirical
    ///   otherwise.
    /// * `analytical`: closed-form, cheap, approximate. Documented
    ///   in `dottir_core::periodogram::analytical_z_scores`.
    /// * `empirical`: shuffle-based, slower, principled.
    /// * `off`: emit `nan` in the z_score column.
    #[arg(long, value_enum, default_value_t = ZScoreModeArg::Auto)]
    pub z_score: ZScoreModeArg,

    /// Number of shuffles for empirical z-score.
    #[arg(long, default_value_t = 200, value_name = "INT")]
    pub z_shuffles: u32,

    /// Seed for empirical-mode shuffling. Same seed → identical output.
    #[arg(long, default_value_t = 0, value_name = "INT")]
    pub seed: u64,

    /// Memory cap per record. Bytes; suffix-less integer.
    /// Default 1 GiB. The streaming algorithm is O(N) memory, so this
    /// only guards against pathologically large records.
    #[arg(long, default_value_t = 1 << 30, value_name = "BYTES")]
    pub memory_limit: u64,

    /// Rayon worker threads. `1` forces single-threaded; `0` uses the
    /// rayon default (one per CPU). Set this before any rayon call.
    #[arg(long, default_value_t = 0, value_name = "INT")]
    pub threads: usize,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
pub enum ZScoreModeArg {
    Auto,
    Analytical,
    Empirical,
    Off,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum ResolvedZMode {
    Analytical,
    Empirical,
    Off,
}

pub fn run(args: PeriodogramArgs) -> Result<()> {
    // Configure rayon pool BEFORE any compute. `build_global` is a
    // process-once call; it errors if called twice, so swallow that
    // (subsequent invocations within tests etc. just keep the
    // existing pool).
    if args.threads > 0 {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global();
    }

    let matrix = ScoreMatrix::custom_dna(args.match_score, args.mismatch_score);

    let sensitivity = Sensitivity {
        white: args.greyramp_white,
        black: args.greyramp_black,
    };
    let resolved_z = resolve_z_mode(args.z_score, sensitivity);

    let cfg = PeriodogramConfig {
        matrix: matrix.clone(),
        window_size: args.window,
        pixel_fac: args.pixel_fac,
        sensitivity,
        min_offset: args.min_offset,
        max_offset: args.max_offset,
        memory_limit_bytes: args.memory_limit,
    };

    tracing::info!("reading {}", args.input.display());
    let loaded = fasta::load_fasta_file(&args.input)
        .with_context(|| format!("reading {}", args.input.display()))?;
    if loaded.records.is_empty() {
        anyhow::bail!("input FASTA contains no records");
    }

    let mut out: Box<dyn Write> = match args.output.as_ref() {
        Some(p) => Box::new(BufWriter::new(
            File::create(p).with_context(|| format!("creating {}", p.display()))?,
        )),
        None => Box::new(BufWriter::new(std::io::stdout().lock())),
    };

    write_global_header(&mut out, &args.input, &args, sensitivity, resolved_z)?;
    writeln!(
        out,
        "record_id\tk\traw_sum\tsignal_sum\tsignal_mean\tz_score"
    )?;

    let force_serial = args.threads == 1;
    for record in &loaded.records {
        process_record(
            &mut out,
            &record.id,
            &record.sequence,
            &cfg,
            resolved_z,
            args.z_shuffles,
            args.seed,
            force_serial,
        )
        .with_context(|| format!("processing record {}", record.id))?;
    }

    out.flush().context("flushing periodogram output")?;
    if let Some(p) = args.output.as_ref() {
        tracing::info!("wrote {}", p.display());
    }
    Ok(())
}

fn resolve_z_mode(arg: ZScoreModeArg, sensitivity: Sensitivity) -> ResolvedZMode {
    match arg {
        ZScoreModeArg::Off => ResolvedZMode::Off,
        ZScoreModeArg::Analytical => ResolvedZMode::Analytical,
        ZScoreModeArg::Empirical => ResolvedZMode::Empirical,
        ZScoreModeArg::Auto => {
            if sensitivity.is_identity() {
                ResolvedZMode::Analytical
            } else {
                ResolvedZMode::Empirical
            }
        }
    }
}

fn write_global_header<W: Write>(
    out: &mut W,
    input: &std::path::Path,
    args: &PeriodogramArgs,
    sensitivity: Sensitivity,
    z_mode: ResolvedZMode,
) -> Result<()> {
    writeln!(out, "# dottir periodogram v{}", env!("CARGO_PKG_VERSION"))?;
    writeln!(out, "# input: {}", input.display())?;
    writeln!(
        out,
        "# matrix: dna match={} mismatch={}",
        args.match_score, args.mismatch_score
    )?;
    match args.window {
        Some(w) => writeln!(out, "# window: {w} (explicit)")?,
        None => writeln!(out, "# window: auto-Karlin per record")?,
    }
    match args.pixel_fac {
        0 => writeln!(out, "# pixel_fac: auto-Karlin per record")?,
        n => writeln!(out, "# pixel_fac: {n} (explicit)")?,
    }
    writeln!(
        out,
        "# sensitivity: white={} black={}",
        sensitivity.white, sensitivity.black
    )?;
    writeln!(out, "# min_offset: {}", args.min_offset)?;
    if let Some(m) = args.max_offset {
        writeln!(out, "# max_offset: {m} (explicit)")?;
    } else {
        writeln!(out, "# max_offset: floor(N/2) per record")?;
    }
    let z_label = match z_mode {
        ResolvedZMode::Analytical => "analytical",
        ResolvedZMode::Empirical => "empirical",
        ResolvedZMode::Off => "off",
    };
    writeln!(
        out,
        "# z_score: {z_label} (requested {:?}, shuffles={}, seed={})",
        args.z_score, args.z_shuffles, args.seed
    )?;
    writeln!(out, "# memory_limit: {} bytes", args.memory_limit)?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn process_record<W: Write>(
    out: &mut W,
    record_id: &str,
    seq: &[u8],
    cfg: &PeriodogramConfig,
    z_mode: ResolvedZMode,
    z_shuffles: u32,
    seed: u64,
    force_serial: bool,
) -> Result<()> {
    let periodogram = if force_serial {
        compute_periodogram(seq, cfg)
    } else {
        compute_periodogram_parallel(seq, cfg)
    }
    .context("compute_periodogram failed")?;

    let n_buckets = periodogram.signal_sum.len();
    let z = match z_mode {
        ResolvedZMode::Off => vec![f64::NAN; n_buckets],
        ResolvedZMode::Analytical => {
            let null = analytical_null(&cfg.matrix, &periodogram.residue_freqs);
            analytical_z_scores(&periodogram, &null)
        }
        ResolvedZMode::Empirical => {
            let stats = empirical_null_stats(seq, cfg, z_shuffles, seed)
                .context("empirical_null_stats failed")?;
            periodogram
                .signal_sum
                .iter()
                .zip(stats.iter())
                .map(|(&sum, &(mean, std))| {
                    if std <= f64::EPSILON {
                        0.0
                    } else {
                        (sum as f64 - mean) / std
                    }
                })
                .collect()
        }
    };

    // Per-record header line.
    let analytical = analytical_null(&cfg.matrix, &periodogram.residue_freqs);
    write_record_header(out, record_id, &periodogram, analytical, z_mode, z_shuffles)?;

    // Data rows.
    let rows = periodogram
        .raw_sum
        .iter()
        .zip(periodogram.signal_sum.iter())
        .zip(periodogram.n_pairs.iter())
        .zip(z.iter())
        .enumerate();
    for (i, (((raw, sig), np), z_val)) in rows {
        let k = periodogram.min_offset + i as u32;
        let mean = if *np == 0 {
            0.0
        } else {
            *sig as f64 / *np as f64
        };
        writeln!(out, "{record_id}\t{k}\t{raw}\t{sig}\t{mean:.6}\t{z_val:.4}")?;
    }
    Ok(())
}

fn write_record_header<W: Write>(
    out: &mut W,
    record_id: &str,
    p: &Periodogram,
    null: AnalyticalNull,
    z_mode: ResolvedZMode,
    z_shuffles: u32,
) -> Result<()> {
    let z_label = match z_mode {
        ResolvedZMode::Analytical => "analytical",
        ResolvedZMode::Empirical => "empirical",
        ResolvedZMode::Off => "off",
    };
    write!(
        out,
        "# record_id={record_id} length={} window={} pixel_fac={} \
         freq_A={:.3} freq_C={:.3} freq_G={:.3} freq_T={:.3} \
         null_mean={:.4} null_var={:.4} z={z_label}",
        p.seq_len,
        p.window_size,
        p.pixel_fac,
        p.residue_freqs[0],
        p.residue_freqs[1],
        p.residue_freqs[2],
        p.residue_freqs[3],
        null.mean_per_pair,
        null.var_per_pair,
    )?;
    if matches!(z_mode, ResolvedZMode::Empirical) {
        write!(out, " z_shuffles={z_shuffles}")?;
    }
    writeln!(out)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn auto_resolves_to_analytical_for_identity_sensitivity() {
        assert_eq!(
            resolve_z_mode(ZScoreModeArg::Auto, Sensitivity::identity()),
            ResolvedZMode::Analytical
        );
    }

    #[test]
    fn auto_resolves_to_empirical_when_greyramp_clips() {
        assert_eq!(
            resolve_z_mode(ZScoreModeArg::Auto, Sensitivity::gui_default()),
            ResolvedZMode::Empirical
        );
    }

    #[test]
    fn explicit_overrides_auto_resolution() {
        // Even with identity sensitivity, --z-score empirical sticks.
        assert_eq!(
            resolve_z_mode(ZScoreModeArg::Empirical, Sensitivity::identity()),
            ResolvedZMode::Empirical
        );
        // And vice versa.
        assert_eq!(
            resolve_z_mode(ZScoreModeArg::Analytical, Sensitivity::gui_default()),
            ResolvedZMode::Analytical
        );
    }
}
