//! I/O for dottir: FASTA, GFF3, PAF, matrix files, `.dot`, image exports.
//!
//! `dottir-core` is intentionally I/O-free; everything that touches the
//! filesystem, the network, or an output image lives here. Phases 1+ will
//! add the actual readers/writers.

#![forbid(unsafe_code)]
