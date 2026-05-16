//! `dottir` CLI binary. Phase 4 fills this in; for now it is a stub so the
//! workspace builds.

use anyhow::Result;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(
    name = "dottir",
    version,
    about = "Modern Rust reimplementation of Dotter (Sonnhammer & Durbin 1995)"
)]
struct Cli {
    /// Subcommand. Only `batch` is planned for v1.0.
    #[command(subcommand)]
    command: Option<Command>,
}

#[derive(clap::Subcommand, Debug)]
enum Command {
    /// Headless batch mode: compute a dotplot and write image/sidecar.
    Batch,
}

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();
    let cli = Cli::parse();
    match cli.command {
        None => {
            eprintln!("dottir {} — CLI is not yet implemented (Phase 4)", env!("CARGO_PKG_VERSION"));
            eprintln!("See IMPLEMENTATION_PLAN.md");
        }
        Some(Command::Batch) => {
            eprintln!("dottir batch: not yet implemented (Phase 4)");
        }
    }
    Ok(())
}
