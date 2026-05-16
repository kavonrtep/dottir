//! Top-level dotplot computation.
//!
//! Phase 0 ships the public signatures; Phase 1 wires the BLASTN forward
//! kernel through them; Phase 2 adds reverse-strand, BLASTP/BLASTX, and
//! self-comparison.

use crate::error::DottirError;
use crate::karlin::KarlinResult;
use crate::matrix::{BlastMode, ScoreMatrix};

/// Which strand(s) of the *subject* sequence to compute against the query.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    Both,
}

/// Choice of triangle in self-comparison mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Triangle {
    Upper,
    Lower,
}

/// Everything needed to compute one dotplot. See spec §6.2.
#[derive(Debug, Clone)]
pub struct PlotConfig {
    pub mode: BlastMode,
    pub matrix: ScoreMatrix,
    /// `None` = use Karlin/Altschul estimate.
    pub window_size: Option<u32>,
    /// `zoom × zoom` block size per output pixel. Spec §4.1.5.
    pub zoom: u32,
    /// Multiplier in the final scale step `min(255, score * pixel_fac / W)`.
    pub pixel_fac: u32,
    pub strand: Strand,
    pub self_comparison: bool,
    pub triangle: Triangle,
    pub disable_mirror: bool,
    /// Refuse to allocate a pixelmap larger than this many bytes.
    /// Default 0.5 GiB; spec §4.5.2.
    pub memory_limit_bytes: u64,
    /// If true, keep forward and reverse dot channels separate in the
    /// returned [`DotPlot`] (spec §4.4.3 — inverted-repeat highlighting).
    pub separate_strand_channels: bool,
}

impl PlotConfig {
    /// Sensible BLASTN defaults: forward+reverse, zoom 1, pixel_fac matching
    /// C dotter (50), 0.5 GiB cap.
    pub fn default_blastn(matrix: ScoreMatrix) -> Self {
        Self {
            mode: BlastMode::Blastn,
            matrix,
            window_size: None,
            zoom: 1,
            pixel_fac: 50,
            strand: Strand::Both,
            self_comparison: false,
            triangle: Triangle::Upper,
            disable_mirror: false,
            memory_limit_bytes: 512 * 1024 * 1024,
            separate_strand_channels: false,
        }
    }

    /// Sensible BLASTP defaults.
    pub fn default_blastp(matrix: ScoreMatrix) -> Self {
        Self {
            mode: BlastMode::Blastp,
            matrix,
            window_size: None,
            zoom: 1,
            pixel_fac: 50,
            strand: Strand::Forward,
            self_comparison: false,
            triangle: Triangle::Upper,
            disable_mirror: false,
            memory_limit_bytes: 512 * 1024 * 1024,
            separate_strand_channels: false,
        }
    }
}

/// Resolved parameters for a completed dotplot. Mirrors `KarlinResult` plus
/// the inputs that the user might not have specified.
#[derive(Debug, Clone, Copy)]
pub struct PlotParams {
    pub window_size: u32,
    pub zoom: u32,
    pub pixel_fac: u32,
    pub karlin: Option<KarlinResult>,
}

/// A computed dotplot. `pixels.len() == width * height`; `forward_pixels`
/// and `reverse_pixels` are populated only when
/// [`PlotConfig::separate_strand_channels`] is set.
#[derive(Debug, Clone)]
pub struct DotPlot {
    pub width: u32,
    pub height: u32,
    pub pixels: Vec<u8>,
    pub forward_pixels: Option<Vec<u8>>,
    pub reverse_pixels: Option<Vec<u8>>,
    pub params: PlotParams,
}

/// Compute a dotplot. Stub for Phase 0; populated by Phase 1+.
pub fn compute_dotplot(
    _query: &[u8],
    _subject: &[u8],
    _config: &PlotConfig,
) -> Result<DotPlot, DottirError> {
    Err(DottirError::NotImplemented(
        "compute_dotplot — Phase 1 will land the BLASTN forward kernel",
    ))
}
