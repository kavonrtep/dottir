//! `dottir-core` — algorithmic core of dottir (a Rust reimplementation of
//! Dotter, Sonnhammer & Durbin 1995).
//!
//! This crate is I/O-free and GUI-free by design: it can be embedded in
//! notebooks, other Rust tools, or (eventually) `pyo3` bindings. All file
//! handling lives in the `dottir-io` crate; all rendering in `dottir-gui`.
//!
//! ## Scope by phase
//!
//! - **Phase 0** (this version): Karlin/Altschul statistics and the built-in
//!   score matrices. See [`karlin`] and [`matrix`].
//! - **Phase 1+**: the sliding-window dotplot kernel. The public entry point
//!   [`compute_dotplot`] is reserved but currently returns
//!   [`DottirError::NotImplemented`]; the supporting modules
//!   ([`score_vec`], [`sliding`], [`pixel`], [`antidiag`]) are skeletons.
//!
//! ## Determinism
//!
//! Per spec §4.1.11 and CLAUDE.md, identical inputs MUST produce
//! byte-identical pixelmaps across runs and thread counts. The core
//! deliberately avoids `HashMap` iteration in hot paths.

// Default-deny unsafe; locally-allowed only where the std lib's
// equivalent safe API is still nightly. See `pixel::PixelMap::view_mut`
// and `pixel::PixelMap::into_vec` — both transmute between
// layout-identical types (`u8` ↔ `AtomicU8`) and document why they
// are sound.
#![deny(unsafe_code)]
#![warn(missing_debug_implementations)]

pub mod alphabet;
pub mod antidiag;
pub mod error;
pub mod find_peaks;
pub mod karlin;
pub mod matrix;
pub mod periodogram;
pub mod pixel;
pub mod plot;
pub mod ridges;
pub mod score_vec;
pub mod sliding;
pub mod spectrum;
pub mod translation;

pub use error::DottirError;
pub use find_peaks::{
    auto_threshold_mad, find_peaks_in_periodogram, find_peaks_in_spectrum, HarmonicDirection, Peak,
    PeakKind, PeaksConfig, SubrepeatConfig,
};
pub use karlin::{karlin_window_size, KarlinConfig, KarlinResult};
pub use matrix::{BlastMode, ScoreMatrix};
pub use periodogram::{
    analytical_null, analytical_z_scores, compute_periodogram, compute_periodogram_parallel,
    empirical_null_stats, AnalyticalNull, Periodogram, PeriodogramConfig, Sensitivity,
};
pub use plot::{
    compute_dotplot, pick_auto_zoom, reverse_complement, snap_zoom_to_period_divisor, DotPlot,
    PlotConfig, PlotParams, Strand, Triangle,
};
pub use ridges::{extract_ridges, extract_ridges_from_pixels, Ridge, RidgeDirection, RidgeParams};
pub use spectrum::{spectrum_of_signal, DetrendMode, Spectrum, SpectrumConfig, WindowFn};

/// Bumped whenever the algorithmic contract changes such that previously
/// pinned golden pixelmaps must be regenerated. See CLAUDE.md.
pub const PIXELMAP_FORMAT_VERSION: u32 = 0;
