//! Top-level dotplot computation.
//!
//! Phase 0 ships the public signatures; Phase 1 wires the BLASTN forward
//! kernel through them; Phase 2 adds reverse-strand, BLASTP/BLASTX, and
//! self-comparison.

use crate::alphabet::encode;
use crate::error::DottirError;
use crate::karlin::{karlin_window_size, KarlinResult};
use crate::matrix::{BlastMode, ScoreMatrix};
use crate::pixel::{image_dimension, PixelMap};
use crate::score_vec::ScoreVec;
use crate::sliding::{sliding_window_pass, Direction};

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

/// Compute a dotplot for a `(query, subject)` pair.
///
/// **Phase 1 scope**: BLASTN with `strand = Forward`, no self-comparison.
/// Other modes return [`DottirError::NotImplemented`] and will be filled in
/// by Phase 2 (reverse strand, BLASTP, BLASTX, self-comparison).
///
/// # Example
///
/// ```
/// use dottir_core::{compute_dotplot, ScoreMatrix, PlotConfig, Strand};
///
/// let q = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".repeat(2);
/// let s = q.clone();
/// let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
/// cfg.strand = Strand::Forward;
/// cfg.window_size = Some(8);
/// cfg.zoom = 1;
/// let plot = compute_dotplot(&q, &s, &cfg).unwrap();
/// // The main diagonal should be heavily lit on a self-similar input.
/// let diag_pixel = plot.pixels[plot.width as usize / 2 * plot.width as usize + plot.width as usize / 2];
/// assert!(diag_pixel > 0);
/// ```
pub fn compute_dotplot(
    query: &[u8],
    subject: &[u8],
    config: &PlotConfig,
) -> Result<DotPlot, DottirError> {
    if query.is_empty() || subject.is_empty() {
        return Err(DottirError::EmptySequence);
    }
    if config.zoom == 0 {
        return Err(DottirError::InvalidConfig("zoom must be >= 1".into()));
    }
    if config.pixel_fac == 0 {
        return Err(DottirError::InvalidConfig("pixel_fac must be >= 1".into()));
    }
    if config.self_comparison {
        return Err(DottirError::NotImplemented(
            "self_comparison — Phase 2",
        ));
    }
    if config.mode != BlastMode::Blastn {
        return Err(DottirError::NotImplemented(
            "BLASTP / BLASTX — Phase 2",
        ));
    }
    if matches!(config.strand, Strand::Reverse | Strand::Both) {
        return Err(DottirError::NotImplemented(
            "reverse / both strand — Phase 2",
        ));
    }

    // Karlin window estimate, or honour the override.
    let (window, karlin_result) = match config.window_size {
        Some(w) => (w, None),
        None => {
            let r = karlin_window_size(&config.matrix, query, subject, config.mode)?;
            (r.window_size, Some(r))
        }
    };
    if window < 1 {
        return Err(DottirError::InvalidConfig(
            "window size must be >= 1".into(),
        ));
    }

    // Encode sequences using the alphabet for the mode.
    let alpha = config.mode.alphabet();
    let q_encoded = encode(query, alpha);
    let s_encoded = encode(subject, alpha);

    // Output dimensions.
    let width = image_dimension(q_encoded.len(), config.zoom);
    let height = image_dimension(s_encoded.len(), config.zoom);
    let mut pixelmap =
        PixelMap::new_checked(width, height, config.memory_limit_bytes)?;

    // Precomputed score vector keyed on the query.
    let score_vec = ScoreVec::build(&config.matrix, &q_encoded);

    // Phase 1: single forward pass.
    sliding_window_pass(
        &score_vec,
        &s_encoded,
        window,
        config.zoom,
        config.pixel_fac,
        Direction::Forward,
        &mut pixelmap,
    );

    Ok(DotPlot {
        width,
        height,
        pixels: pixelmap.into_vec(),
        forward_pixels: None,
        reverse_pixels: None,
        params: PlotParams {
            window_size: window,
            zoom: config.zoom,
            pixel_fac: config.pixel_fac,
            karlin: karlin_result,
        },
    })
}
