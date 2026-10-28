//! Genome mode (ADR 0008): whole-assembly dotplots drawn from alignments
//! produced by an external aligner (PAF), at D-GENIES scale.
//!
//! This module is I/O-free. `dottir-io` reads PAF into [`Alignment`]s and
//! record lists from FASTA indexes; everything here works on those.
//!
//! * [`alignment`] — the alignment model, identity, chunk-cut joining;
//! * [`filter`] — length / identity / primary filters, D-GENIES noise;
//! * [`order`] — query record sorting and orientation (D-GENIES);
//! * [`layout`] — record order and orientation → axis coordinates;
//! * [`raster`] — order-independent rasterisation.
//!
//! Filters, sorting and the small-record mix blocks follow D-GENIES
//! (Cabanettes & Klopp 2018, `dgenies/lib/paf.py`, GPL-3.0).

/// Version of the genome-mode raster. Bump on any change to the pixels
/// produced for the same alignments and parameters.
pub const GENOME_FORMAT_VERSION: u32 = 1;

pub mod alignment;
pub mod filter;
pub mod layout;
pub mod order;
pub mod raster;

pub use alignment::{Alignment, GenomeRecord};
pub use filter::FilterConfig;
pub use layout::{AxisLayout, AxisSegment};
pub use order::{sort_query, QueryOrder};
pub use raster::{rasterize, GenomeRaster, RasterConfig};
