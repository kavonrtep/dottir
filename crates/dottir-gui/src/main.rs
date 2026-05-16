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
//!
//! Out of scope for Phase 5 (tracked separately): GFF3 / PAF tracks
//! (ADR 0004 follow-up), spawn-sub-dotter, alignment-view dock,
//! breakline rendering for multi-record FASTA.

#![forbid(unsafe_code)]

mod app;

use anyhow::Result;
use app::DottirApp;

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();

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
        Box::new(|cc| Ok(Box::new(DottirApp::new(cc)))),
    )
    .map_err(|e| anyhow::anyhow!("eframe failed: {e}"))
}
