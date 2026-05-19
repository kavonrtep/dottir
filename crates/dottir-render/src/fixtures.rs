//! Synthetic `DotPlot` fixtures with known geometric structure.
//!
//! Each fixture is built by feeding a deterministic sequence into
//! `dottir_core::compute_dotplot`, so the rendering tests exercise
//! the same code path that the GUI does. The "known geometry" is
//! what the metrics in [`crate::metrics`] check against.

use dottir_core::{compute_dotplot, DotPlot, PlotConfig, ScoreMatrix, Strand};

/// Self-comparison of `unit` repeated `n_units` times.
///
/// Expected geometric structure in the raster:
///
/// * Parallel diagonals at offsets `y - x = k * unit.len()` for
///   integer `k` in `-(W-1)..(W-1)` where `W = unit.len() * n_units`.
/// * Each diagonal is at most one pixel wide (window equals unit
///   length, so the recurrence only saturates exactly on-period).
///
/// Any rendering policy that preserves these two properties at all
/// physical scales is "faithful". The current GUI policy fails on
/// the second property at non-integer scales — that's the bug we
/// are measuring.
pub fn tandem_repeat(unit: &[u8], n_units: usize) -> DotPlot {
    assert!(!unit.is_empty());
    assert!(n_units >= 2);
    let seq: Vec<u8> = unit
        .iter()
        .cycle()
        .take(unit.len() * n_units)
        .copied()
        .collect();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.zoom = 1;
    cfg.window_size = Some(unit.len() as u32);
    cfg.self_comparison = true;
    compute_dotplot(&seq, &seq, &cfg).expect("tandem_repeat fixture: compute_dotplot failed")
}
