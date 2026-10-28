//! `dottir-gui` — interactive frontend (egui/eframe).
//!
//! Phase 5 MVP per docs/IMPLEMENTATION_PLAN.md §5:
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

mod annotation_overlay;
mod app;
mod compute_worker;
mod genome_view;
mod record_sort;
mod session;

use std::path::PathBuf;

use anyhow::Result;
use app::{DottirApp, StartupConfig};
use clap::Parser;
use dottir_core::{BlastMode, Strand};
use genome_view::GenomeInputs;

/// Mirror the original `dotter [options] <horizontal_sequence>
/// [<vertical_sequence>]` CLI where it makes sense for the GUI.
///
/// One positional path triggers a self-comparison; two paths load a
/// query (horizontal) / subject (vertical) pair. Each path is a FASTA
/// file, or a GFF3 file carrying its own sequences after a `##FASTA`
/// directive — the latter needs no separate `--gff-*` flag. All flags
/// are optional; anything not given falls back to the GUI's defaults
/// and can be changed from the Settings window.
#[derive(Parser, Debug)]
#[command(
    name = "dottir-gui",
    version,
    about = "Interactive dotplot viewer (Rust reimplementation of Dotter)",
    long_about = "Interactive dotplot viewer.\n\
                  \n\
                  Run with no arguments to open an empty window and use \
                  File → Open to load sequences, or pass one or two paths \
                  to pre-load them and compute the dotplot on startup.\n\
                  \n\
                  A path may be a FASTA file, or a GFF3 file with an \
                  embedded ##FASTA section — the latter supplies both the \
                  sequence and its annotation overlay from one file."
)]
struct Cli {
    /// Query sequence — horizontal axis. FASTA, or a GFF3 with an
    /// embedded `##FASTA` section (its features are overlaid on this
    /// axis). Optional.
    #[arg(value_name = "QUERY")]
    query: Option<PathBuf>,

    /// Subject sequence — vertical axis. FASTA or GFF3-with-FASTA.
    /// Omit (with QUERY given) for a self-comparison.
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

    /// Memory cap for the pixelmap, in MiB.
    /// Original Dotter equivalent: `-m <float>`.
    #[arg(short = 'm', long, value_name = "MiB", default_value_t = 512)]
    memory_mib: u32,

    /// Annotation file (GFF3 or BED; format auto-detected) for the
    /// query (horizontal) axis. `--bed-query` is an accepted alias.
    /// Overrides features embedded in a GFF3 query input.
    #[arg(long = "gff-query", visible_alias = "bed-query", value_name = "PATH")]
    annot_query: Option<PathBuf>,

    /// Annotation file (GFF3 or BED) for the subject (vertical) axis.
    /// `--bed-subject` is an accepted alias. Overrides features
    /// embedded in a GFF3 subject input.
    #[arg(
        long = "gff-subject",
        visible_alias = "bed-subject",
        value_name = "PATH"
    )]
    annot_subject: Option<PathBuf>,

    /// Annotation file (GFF3 or BED) applied to both axes — the natural
    /// choice for a self-comparison. `--bed` is an accepted alias.
    /// Per-axis flags override this.
    #[arg(long = "gff", visible_alias = "bed", value_name = "PATH")]
    annot_self: Option<PathBuf>,

    /// Open QUERY and SUBJECT as a genome comparison (ADR 0008): QUERY
    /// is the reference (horizontal), SUBJECT the assembly drawn against
    /// it (omit for a self-comparison). Alignments come from an aligner
    /// or from `--paf`.
    #[arg(long)]
    genome: bool,

    /// Genome comparison from this PAF (implies `--genome`).
    #[arg(long, value_name = "FILE")]
    paf: Option<PathBuf>,

    /// Genome mode aligner: `auto` (mm2plus, else minimap2), a name on
    /// PATH, or a path.
    #[arg(long, value_name = "NAME|PATH", default_value = "auto")]
    aligner: String,

    /// Genome mode minimap2 preset (`-x`).
    #[arg(long, value_name = "PRESET", default_value = "asm20")]
    preset: String,

    /// Extra aligner arguments, space-separated.
    #[arg(long, value_name = "ARGS", allow_hyphen_values = true)]
    aligner_args: Option<String>,

    /// Genome mode: keep all copies of repeats (`-P`).
    #[arg(long, alias = "self-all-copies")]
    all_copies: bool,

    /// Reorder and reverse-complement the subject's (`vertical`) or the
    /// query's (`horizontal`) records along the other sequence after
    /// loading (the D-GENIES contig sort; DNA, two inputs). Also in
    /// *View → Sort records*.
    #[arg(long, value_enum, value_name = "AXIS")]
    sort_records: Option<SortArg>,

    /// Draw the alignments of a PAF file over the plot (HSPs). Records
    /// not plotted are skipped; either PAF orientation works.
    #[arg(long, value_name = "PAF")]
    hsp_paf: Option<PathBuf>,

    /// How HSPs are drawn.
    #[arg(long, value_name = "MODE", value_parser = ["line", "score", "greyscale", "off"])]
    hsp_mode: Option<String>,
}

#[derive(Copy, Clone, Debug, clap::ValueEnum)]
enum SortArg {
    Vertical,
    Horizontal,
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
        if self.genome || self.paf.is_some() {
            if let Some(reference) = self.query.clone() {
                return StartupConfig {
                    genome: Some(GenomeInputs {
                        reference,
                        query: self.subject,
                        paf: self.paf,
                        aligner: self.aligner,
                        preset: self.preset,
                        aligner_args: self.aligner_args.unwrap_or_default(),
                        all_copies: self.all_copies,
                    }),
                    ..StartupConfig::default()
                };
            }
        }
        let mode = self.mode.to_core();
        // One positional → self-comparison (subject mirrors query).
        // Two positionals → pairwise as given.
        let (subject, self_comparison) = match (self.query.as_ref(), self.subject) {
            (Some(q), None) => (Some(q.clone()), true),
            (_, sub) => (sub, false),
        };
        StartupConfig {
            query: self.query,
            subject,
            mode,
            matrix_name: self.matrix.unwrap_or_else(|| match mode {
                BlastMode::Blastn => "DNA+5/-4".into(),
                _ => "BLOSUM62".into(),
            }),
            window_size: self.window,
            zoom: self.zoom.max(1),
            pixel_fac: self.pixel_fac.max(1),
            strand: self.strand.to_core(),
            self_comparison,
            memory_limit_bytes: (self.memory_mib as u64) * 1024 * 1024,
            // `--gff`/`--bed` applies to both axes; per-axis flags win.
            annot_query: self.annot_query.or_else(|| self.annot_self.clone()),
            annot_subject: self.annot_subject.or(self.annot_self),
            genome: None,
            hsp_paf: self.hsp_paf,
            hsp_mode: self
                .hsp_mode
                .as_deref()
                .and_then(dottir_io::hsp::HspMode::parse),
            sort_records: self.sort_records.map(|a| match a {
                SortArg::Vertical => dottir_core::genome::SortAxis::Vertical,
                SortArg::Horizontal => dottir_core::genome::SortAxis::Horizontal,
            }),
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
