//! I/O for dottir: FASTA, GFF3, PAF, matrix files, `.dot`, image exports.
//!
//! `dottir-core` is intentionally I/O-free; everything that touches the
//! filesystem, the network, or an output image lives here.
//!
//! ## Phase 4 surface
//!
//! * [`fasta`] — minimal FASTA reader (plain + gzipped).
//! * [`png_export`] — greyscale 8-bit PNG with `tEXt` provenance chunks.
//! * [`params`] — TOML sidecar describing the inputs and parameters.
//!
//! ## Future phases
//!
//! GFF3 (`noodles-gff`), PAF (`noodles-paf`), `.dot` binary format,
//! and SVG/PDF export are not yet wired — they land in Phases 4-extra
//! and 6.

#![forbid(unsafe_code)]

pub mod fasta;
pub mod params;
pub mod png_export;
