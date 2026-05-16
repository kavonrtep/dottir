//! `dottir-gui` — interactive frontend.
//!
//! **Status (Phase 5):** the egui/eframe runtime is deferred until the
//! workspace MSRV bumps past Rust 1.75. The egui ≥ 0.27 dependency tree
//! pulls `winit` → `toml_edit-0.25` → `indexmap-2.14`, all of which
//! require edition2024 (Rust 1.85+). See `docs/adr/0003-gui-msrv.md`
//! for the decision.
//!
//! Until then this binary acts as a "compute-and-write-PNG" headless
//! frontend so users can still drive `dottir-core` against a FASTA pair
//! without the egui runtime. It is functionally a thinner version of
//! `dottir-cli batch`.

use std::path::PathBuf;

use anyhow::{Context, Result};

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();

    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("dottir-gui {} (GUI deferred — see docs/adr/0003-gui-msrv.md)", env!("CARGO_PKG_VERSION"));
        eprintln!();
        eprintln!("Headless fallback usage:");
        eprintln!("  dottir-gui <QUERY.fa> <SUBJECT.fa> <OUT.png>");
        eprintln!();
        eprintln!("For full CLI options use the dottir-cli binary:");
        eprintln!("  dottir batch ...");
        std::process::exit(2);
    }
    let query = PathBuf::from(&args[1]);
    let subject = PathBuf::from(&args[2]);
    let output = PathBuf::from(&args[3]);

    use dottir_core::{compute_dotplot, PlotConfig, ScoreMatrix};
    use dottir_io::{fasta, png_export};

    let q_recs = fasta::read_fasta_file(&query)
        .with_context(|| format!("reading query {}", query.display()))?;
    let s_recs = fasta::read_fasta_file(&subject)
        .with_context(|| format!("reading subject {}", subject.display()))?;
    let q_seq = fasta::concatenate(&q_recs);
    let s_seq = fasta::concatenate(&s_recs);
    eprintln!(
        "loaded {} ({} residues) and {} ({} residues)",
        query.display(), q_seq.len(),
        subject.display(), s_seq.len()
    );

    let cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    let plot = compute_dotplot(&q_seq, &s_seq, &cfg)?;
    eprintln!(
        "computed {}×{} pixelmap (W={}, λ={:?})",
        plot.width,
        plot.height,
        plot.params.window_size,
        plot.params.karlin.map(|k| k.lambda)
    );
    png_export::write_grayscale_png(
        &output,
        plot.width,
        plot.height,
        &plot.pixels,
        &[("dottir-gui", env!("CARGO_PKG_VERSION"))],
    )?;
    eprintln!("wrote {}", output.display());
    Ok(())
}
