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

#![forbid(unsafe_code)]
#![warn(missing_debug_implementations)]

pub mod alphabet;
pub mod antidiag;
pub mod error;
pub mod karlin;
pub mod matrix;
pub mod pixel;
pub mod plot;
pub mod score_vec;
pub mod sliding;

pub use error::DottirError;
pub use karlin::{karlin_window_size, KarlinConfig, KarlinResult};
pub use matrix::{BlastMode, ScoreMatrix};
pub use plot::{compute_dotplot, DotPlot, PlotConfig, PlotParams, Strand};

/// Bumped whenever the algorithmic contract changes such that previously
/// pinned golden pixelmaps must be regenerated. See CLAUDE.md.
pub const PIXELMAP_FORMAT_VERSION: u32 = 0;
