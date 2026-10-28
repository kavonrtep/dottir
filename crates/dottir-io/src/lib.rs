//! I/O for dottir: FASTA, GFF3, PAF, matrix files, `.dot`, image exports.
//!
//! `dottir-core` is intentionally I/O-free; everything that touches the
//! filesystem, the network, or an output image lives here.
//!
//! ## Phase 4 surface
//!
//! * [`fasta`] — minimal FASTA reader (plain + gzipped).
//! * [`input`] — positional sequence input: FASTA, or GFF3 carrying
//!   its own sequences after `##FASTA`.
//! * [`png_export`] — greyscale 8-bit PNG with `tEXt` provenance chunks.
//! * [`params`] — TOML sidecar describing the inputs and parameters.
//! * [`paf`], [`faidx`], [`cache`] — genome mode (ADR 0008): PAF reader,
//!   FASTA index and region reader, user cache directory.
//!
//! ## Future phases
//!
//! GFF3 (`noodles-gff`), `.dot` binary format,
//! and SVG/PDF export are not yet wired — they land in Phases 4-extra
//! and 6.

#![forbid(unsafe_code)]

pub mod alignment;
pub mod annotation;
pub mod cache;
pub mod dot_format;
pub mod faidx;
pub mod fasta;
pub mod genome_export;
pub mod input;
pub mod large;
pub mod paf;
pub mod params;
pub mod png_export;
pub mod sequence;
pub mod svg_export;
pub mod text_overlay;

pub use annotation::{AnnotError, AnnotSet, AnnotSource, Feature, Strand, TYPE_KEY};
pub use input::{load_sequence_input, InputError, InputKind, SeqInput};
pub use sequence::{detect_alphabet, DetectedAlphabet, RecordSpan, Sequence};
