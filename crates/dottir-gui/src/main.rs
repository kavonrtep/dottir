//! `dottir-gui` binary. Phase 5 fills this in with egui/eframe.

use anyhow::Result;

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();
    eprintln!("dottir-gui {} — GUI is not yet implemented (Phase 5)", env!("CARGO_PKG_VERSION"));
    eprintln!("See IMPLEMENTATION_PLAN.md");
    Ok(())
}
