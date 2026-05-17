//! `dottir` CLI binary — batch mode (spec §4.4.8).
//!
//! Computes a dotplot from a (query, subject) FASTA pair and writes a
//! greyscale PNG + a `.params.toml` sidecar. Mirrors C dotter CLI option
//! names where reasonable.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{Parser, Subcommand, ValueEnum};

use dottir_core::{
    compute_dotplot, BlastMode, PlotConfig, ScoreMatrix, Strand, Triangle,
    PIXELMAP_FORMAT_VERSION,
};
use dottir_io::{
    fasta,
    params::{
        sha256_bytes, ParamsSidecar, DottirInfo, HostInfo, InputInfo,
        KarlinInfo, PlotParamsInfo,
    },
    png_export,
};

#[derive(Parser, Debug)]
#[command(
    name = "dottir",
    version,
    about = "Modern Rust reimplementation of Dotter (Sonnhammer & Durbin 1995)"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Compute a dotplot and write a PNG + params sidecar.
    Batch(BatchArgs),
}

#[derive(clap::Args, Debug)]
struct BatchArgs {
    /// Query FASTA (horizontal axis).
    #[arg(value_name = "QUERY")]
    query: PathBuf,
    /// Subject FASTA (vertical axis).
    #[arg(value_name = "SUBJECT")]
    subject: PathBuf,
    /// Output PNG path. The params sidecar is written to
    /// `<output>.params.toml` alongside it.
    #[arg(short = 'o', long, value_name = "PATH")]
    output: PathBuf,

    /// BLAST mode.
    #[arg(long, value_enum, default_value_t = ModeArg::Blastn)]
    mode: ModeArg,

    /// Built-in score matrix name. Default: BLOSUM62 for protein modes,
    /// DNA+5/-4 for BLASTN.
    #[arg(long)]
    matrix: Option<String>,

    /// Window size W. Default: Karlin/Altschul estimate.
    #[arg(short = 'W', long, value_name = "INT")]
    window: Option<u32>,

    /// Zoom factor (pixels per matrix block).
    #[arg(short = 'z', long, default_value_t = 1)]
    zoom: u32,

    /// Pixel factor (multiplier in `min(255, score * pixel_fac / W)`).
    #[arg(long, default_value_t = 50)]
    pixel_fac: u32,

    /// Strand selection. BLASTP ignores this.
    #[arg(long, value_enum, default_value_t = StrandArg::Both)]
    strand: StrandArg,

    /// Compute only the forward (watson) strand. Equivalent to
    /// `--strand forward`. Original Dotter flag: `--watson-only`.
    #[arg(long, default_value_t = false, conflicts_with_all = ["crick_only"])]
    watson_only: bool,

    /// Compute only the reverse (crick) strand. Equivalent to
    /// `--strand reverse`. Original Dotter flag: `--crick-only`.
    #[arg(long, default_value_t = false)]
    crick_only: bool,

    /// Pre-process: reverse-complement the query before computation.
    /// Original Dotter flag: `-r`.
    #[arg(short = 'r', long, default_value_t = false)]
    reverse_query: bool,

    /// Pre-process: reverse-complement the subject before computation.
    /// Original Dotter flag: `-v`.
    #[arg(short = 'v', long, default_value_t = false)]
    reverse_subject: bool,

    /// Compute as a self-comparison (query == subject).
    #[arg(long, default_value_t = false)]
    self_comparison: bool,

    /// Self-comparison triangle.
    #[arg(long, value_enum, default_value_t = TriangleArg::Both)]
    triangle: TriangleArg,

    /// Disable the self-comparison mirror post-process.
    #[arg(long, default_value_t = false)]
    disable_mirror: bool,

    /// Memory limit for the pixelmap, in bytes.
    #[arg(long, default_value_t = 512 * 1024 * 1024)]
    memory_limit_bytes: u64,

    /// If set, automatically pick `zoom` so the largest output dimension
    /// is at most this. Avoids surprise OOMs on large inputs
    /// (spec §4.4.8).
    #[arg(long, value_name = "PIXELS")]
    auto_zoom: Option<u32>,

    /// Skip writing the `.params.toml` sidecar.
    #[arg(long, default_value_t = false)]
    no_sidecar: bool,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum ModeArg {
    Blastn,
    Blastp,
    Blastx,
}

impl ModeArg {
    fn to_core(self) -> BlastMode {
        match self {
            ModeArg::Blastn => BlastMode::Blastn,
            ModeArg::Blastp => BlastMode::Blastp,
            ModeArg::Blastx => BlastMode::Blastx,
        }
    }
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum StrandArg {
    Forward,
    Reverse,
    Both,
}

impl StrandArg {
    fn to_core(self) -> Strand {
        match self {
            StrandArg::Forward => Strand::Forward,
            StrandArg::Reverse => Strand::Reverse,
            StrandArg::Both => Strand::Both,
        }
    }
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum TriangleArg {
    Both,
    Upper,
    Lower,
}

impl TriangleArg {
    fn to_core(self) -> Triangle {
        match self {
            TriangleArg::Both => Triangle::Both,
            TriangleArg::Upper => Triangle::Upper,
            TriangleArg::Lower => Triangle::Lower,
        }
    }
}

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();
    let cli = Cli::parse();
    match cli.command {
        Command::Batch(args) => run_batch(args),
    }
}

fn run_batch(args: BatchArgs) -> Result<()> {
    let mode = args.mode.to_core();
    let matrix = pick_matrix(args.matrix.as_deref(), mode)?;

    tracing::info!("reading query  {}", args.query.display());
    let query_loaded = fasta::load_fasta_file(&args.query)
        .with_context(|| format!("reading query {}", args.query.display()))?;
    let query = dottir_io::Sequence::from_records(
        query_loaded.records.clone(),
        Some(args.query.clone()),
    );

    tracing::info!("reading subject {}", args.subject.display());
    let subject_loaded = fasta::load_fasta_file(&args.subject)
        .with_context(|| format!("reading subject {}", args.subject.display()))?;
    let subject = dottir_io::Sequence::from_records(
        subject_loaded.records.clone(),
        Some(args.subject.clone()),
    );

    // auto-zoom: pick zoom so max(qlen, slen) / zoom <= auto_zoom.
    let zoom = match args.auto_zoom {
        Some(target) => {
            let max_dim = query.len().max(subject.len()) as u32;
            ((max_dim + target - 1) / target).max(1)
        }
        None => args.zoom,
    };
    if zoom != args.zoom {
        tracing::info!("auto_zoom chose zoom={zoom} (--zoom={} ignored)", args.zoom);
    }

    // --watson-only / --crick-only override --strand. Mutually
    // exclusive (clap's `conflicts_with` already enforced that).
    let strand = if args.watson_only {
        Strand::Forward
    } else if args.crick_only {
        Strand::Reverse
    } else {
        args.strand.to_core()
    };
    let mut cfg = PlotConfig {
        mode,
        matrix,
        window_size: args.window,
        zoom,
        pixel_fac: args.pixel_fac,
        strand,
        self_comparison: args.self_comparison,
        triangle: args.triangle.to_core(),
        disable_mirror: args.disable_mirror,
        memory_limit_bytes: args.memory_limit_bytes,
        separate_strand_channels: false,
        reverse_query: args.reverse_query,
        reverse_subject: args.reverse_subject,
    };
    // BLASTP only supports Forward.
    if mode == BlastMode::Blastp {
        cfg.strand = Strand::Forward;
    }

    tracing::info!(
        "computing {}×{} dotplot, mode={:?}, W={:?}, zoom={}",
        query.len(),
        subject.len(),
        mode,
        args.window,
        zoom
    );
    let plot = compute_dotplot(query.bytes(), subject.bytes(), &cfg)
        .context("compute_dotplot failed")?;
    tracing::info!(
        "wrote pixelmap {}x{} ({} bytes)",
        plot.width,
        plot.height,
        plot.pixels.len()
    );

    // Provenance text chunks for the PNG (spec §4.4.4).
    let dottir_version = env!("CARGO_PKG_VERSION");
    let resolved = format!(
        "mode={:?};matrix={};W={};zoom={};pixel_fac={};strand={:?}",
        cfg.mode, cfg.matrix.name, plot.params.window_size, cfg.zoom,
        cfg.pixel_fac, cfg.strand
    );
    let text_chunks = [
        ("dottir-version", dottir_version),
        ("dottir-pixelmap-format-version", &PIXELMAP_FORMAT_VERSION.to_string()),
        ("dottir-parameters", &resolved),
    ];
    let text_chunk_refs: Vec<(&str, &str)> = text_chunks
        .iter()
        .map(|(a, b)| (*a, *b as &str))
        .collect();

    png_export::write_grayscale_png(
        &args.output,
        plot.width,
        plot.height,
        &plot.pixels,
        &text_chunk_refs,
    )?;
    tracing::info!("wrote {}", args.output.display());

    if !args.no_sidecar {
        let sidecar_path = sidecar_path(&args.output);
        let sidecar = ParamsSidecar {
            dottir: DottirInfo {
                version: dottir_version.to_string(),
                git_sha: option_env!("GIT_SHA").map(|s| s.to_string()),
                pixelmap_format_version: PIXELMAP_FORMAT_VERSION,
            },
            query: input_info(
                &args.query,
                &query_loaded.bytes,
                &query.records,
                query.bytes(),
            )?,
            subject: input_info(
                &args.subject,
                &subject_loaded.bytes,
                &subject.records,
                subject.bytes(),
            )?,
            plot: PlotParamsInfo {
                mode: format!("{:?}", cfg.mode),
                matrix: cfg.matrix.name.clone(),
                strand: format!("{:?}", cfg.strand),
                window_size: plot.params.window_size,
                zoom: cfg.zoom,
                pixel_fac: cfg.pixel_fac,
                self_comparison: cfg.self_comparison,
                karlin: plot.params.karlin.map(|k| KarlinInfo {
                    lambda: k.lambda,
                    k: k.k,
                    h: k.h,
                    expected_residue_score: k.expected_residue_score,
                    expected_msp_score: k.expected_msp_score,
                    predicted_msp_length: k.predicted_msp_length,
                }),
                width: plot.width,
                height: plot.height,
            },
            host: HostInfo {
                hostname: dottir_io::params::hostname(),
                os: std::env::consts::OS.to_string(),
                timestamp_utc: now_iso_utc(),
            },
        };
        std::fs::write(&sidecar_path, sidecar.to_toml()?)?;
        tracing::info!("wrote sidecar {}", sidecar_path.display());
    }

    Ok(())
}

fn pick_matrix(name: Option<&str>, mode: BlastMode) -> Result<ScoreMatrix> {
    match (name, mode) {
        (None, BlastMode::Blastn) => Ok(ScoreMatrix::dna_identity()),
        (None, _) => Ok(ScoreMatrix::blosum62()),
        (Some(n), BlastMode::Blastn) => {
            if n.eq_ignore_ascii_case("DNA") || n.eq_ignore_ascii_case("DNA+5/-4") {
                Ok(ScoreMatrix::dna_identity())
            } else {
                anyhow::bail!("matrix {n:?} not valid for BLASTN — use DNA+5/-4")
            }
        }
        (Some(n), _) => ScoreMatrix::by_name(n)
            .ok_or_else(|| anyhow::anyhow!("unknown matrix {n:?}")),
    }
}

fn input_info(
    path: &Path,
    bytes: &[u8],
    recs: &[dottir_io::RecordSpan],
    sequence: &[u8],
) -> Result<InputInfo> {
    Ok(InputInfo {
        path: path.display().to_string(),
        sha256: sha256_bytes(bytes),
        size_bytes: bytes.len() as u64,
        n_records: recs.len(),
        total_residues: sequence.len(),
    })
}

fn sidecar_path(output: &Path) -> PathBuf {
    let mut s = output.as_os_str().to_owned();
    s.push(".params.toml");
    PathBuf::from(s)
}

fn now_iso_utc() -> String {
    use std::time::{SystemTime, UNIX_EPOCH};
    let secs = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0);
    // Cheap and correct enough for a sidecar (a few-second drift is fine).
    // Avoids pulling in chrono.
    iso_from_unix(secs as i64)
}

fn iso_from_unix(secs: i64) -> String {
    // Days from epoch to Y-M-D using a public-domain Howard Hinnant
    // algorithm (date.h). Civil time = UTC; no timezone handling.
    let z = secs.div_euclid(86_400);
    let secs_of_day = secs.rem_euclid(86_400);
    let (h, m, s) = (
        (secs_of_day / 3600) as u32,
        ((secs_of_day / 60) % 60) as u32,
        (secs_of_day % 60) as u32,
    );
    let z = z + 719_468; // shift epoch to 0000-03-01
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m_civil = (if mp < 10 { mp + 3 } else { mp - 9 }) as u32;
    let y_civil = y + if m_civil <= 2 { 1 } else { 0 };
    format!("{y_civil:04}-{m_civil:02}-{d:02}T{h:02}:{m:02}:{s:02}Z")
}
