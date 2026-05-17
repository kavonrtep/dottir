//! `dottir-gui` — interactive frontend (egui/eframe).
//!
//! Phase 5 MVP per IMPLEMENTATION_PLAN.md §5:
//!
//! * Top-level eframe app skeleton.
//! * Dotplot canvas as a textured quad with pan + scroll-zoom.
//! * Greyramp panel — black/white sliders + swap/reset, applied as a
//!   256-byte LUT to the displayed image. The underlying pixelmap is
//!   not recomputed (spec §4.2.1).
//! * Crosshair: click sets it; arrow keys nudge by 1; Shift = ×10,
//!   Ctrl = ×100.
//! * Status bar with synchronised (q, s) coordinates and the pixelmap
//!   value under the crosshair.
//! * File menu: Open Query / Open Subject via the native file picker
//!   (rfd); recompute when both are loaded.
//! * Settings panel: matrix, window-size override, zoom, pixel-fac,
//!   strand, self-comparison toggle.
//! * Command-line pre-load: `dottir-gui QUERY SUBJECT [options]`
//!   mirrors the original Dotter invocation. When both sequences are
//!   supplied the dotplot is computed on first frame.

#![forbid(unsafe_code)]

mod app;
mod session;

use std::path::PathBuf;

use anyhow::Result;
use app::{DottirApp, StartupConfig};
use clap::Parser;
use dottir_core::{BlastMode, Strand};

/// Mirror the original `dotter [options] <horizontal_sequence>
/// <vertical_sequence>` CLI where it makes sense for the GUI.
///
/// The two positional arguments are FASTA paths for the query
/// (horizontal axis) and subject (vertical axis). All flags are
/// optional; anything not given falls back to the GUI's defaults and
/// can be changed from the Settings window.
#[derive(Parser, Debug)]
#[command(
    name = "dottir-gui",
    version,
    about = "Interactive dotplot viewer (Rust reimplementation of Dotter)",
    long_about = "Interactive dotplot viewer.\n\
                  \n\
                  Run with no arguments to open an empty window and use \
                  File → Open to load FASTA files, or pass two FASTA paths \
                  to pre-load them and compute the dotplot on startup."
)]
struct Cli {
    /// Query FASTA path — horizontal axis. Optional.
    #[arg(value_name = "QUERY")]
    query: Option<PathBuf>,

    /// Subject FASTA path — vertical axis. Optional.
    #[arg(value_name = "SUBJECT")]
    subject: Option<PathBuf>,

    /// Sliding window size. Omit to let Karlin/Altschul pick.
    /// Original Dotter equivalent: `-W <int>`.
    #[arg(short = 'W', long, value_name = "INT")]
    window: Option<u32>,

    /// Computation zoom factor (pixels per matrix block).
    /// Original Dotter equivalent: `-z <int>`.
    #[arg(short = 'z', long, default_value_t = 1)]
    zoom: u32,

    /// Pixel factor (multiplier in `min(255, score * pixel_fac / W)`).
    /// Original Dotter equivalent: `-p <int>`.
    #[arg(short = 'p', long, default_value_t = 50)]
    pixel_fac: u32,

    /// BLAST mode.
    #[arg(long, value_enum, default_value_t = ModeArg::Blastn)]
    mode: ModeArg,

    /// Built-in score matrix name. BLASTN defaults to DNA+5/-4;
    /// protein modes default to BLOSUM62.
    #[arg(long, value_name = "NAME")]
    matrix: Option<String>,

    /// Strand selection (BLASTN only).
    #[arg(long, value_enum, default_value_t = StrandArg::Both)]
    strand: StrandArg,

    /// Compute as a self-comparison (query == subject).
    #[arg(long, default_value_t = false)]
    self_comparison: bool,

    /// Memory cap for the pixelmap, in MiB.
    /// Original Dotter equivalent: `-m <float>`.
    #[arg(short = 'm', long, value_name = "MiB", default_value_t = 512)]
    memory_mib: u32,
}

#[derive(Copy, Clone, Debug, clap::ValueEnum)]
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

#[derive(Copy, Clone, Debug, clap::ValueEnum)]
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

impl Cli {
    fn into_startup(self) -> StartupConfig {
        let mode = self.mode.to_core();
        StartupConfig {
            query: self.query,
            subject: self.subject,
            mode,
            matrix_name: self.matrix.unwrap_or_else(|| match mode {
                BlastMode::Blastn => "DNA+5/-4".into(),
                _ => "BLOSUM62".into(),
            }),
            window_size: self.window,
            zoom: self.zoom.max(1),
            pixel_fac: self.pixel_fac.max(1),
            strand: self.strand.to_core(),
            self_comparison: self.self_comparison,
            memory_limit_bytes: (self.memory_mib as u64) * 1024 * 1024,
        }
    }
}

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();

    let cli = Cli::parse();
    let startup = cli.into_startup();

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("dottir")
            .with_inner_size([1200.0, 800.0])
            .with_min_inner_size([640.0, 480.0]),
        ..Default::default()
    };

    eframe::run_native(
        "dottir",
        options,
        Box::new(move |cc| Ok(Box::new(DottirApp::new(cc, startup)))),
    )
    .map_err(|e| anyhow::anyhow!("eframe failed: {e}"))
}
