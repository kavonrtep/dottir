//! Genome mode (ADR 0008): whole-assembly dotplots drawn from alignments
//! produced by an external aligner (PAF), at D-GENIES scale.
//!
//! This module is I/O-free. `dottir-io` reads PAF into [`Alignment`]s and
//! record lists from FASTA indexes; everything here works on those.
//!
//! * [`alignment`] — the alignment model, identity, chunk-cut joining;
//! * [`layout`] — record order and orientation → axis coordinates.

pub mod alignment;
pub mod layout;

pub use alignment::{Alignment, GenomeRecord};
pub use layout::{AxisLayout, AxisSegment};
