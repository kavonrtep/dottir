//! Top-level dotplot computation.
//!
//! Phase 0 ships the public signatures; Phase 1 wires the BLASTN forward
//! kernel through them; Phase 2 adds reverse-strand, BLASTP/BLASTX, and
//! self-comparison.

use crate::alphabet::{complement_encoded, encode, reverse_complement_dna};
use crate::error::DottirError;
use crate::karlin::{karlin_window_size, KarlinResult};
use crate::matrix::{BlastMode, ScoreMatrix};
use crate::pixel::{image_dimension, merge_max_into, PixelMap};
use crate::score_vec::ScoreVec;
use crate::sliding::{sliding_window_pass_chunked, Direction};

/// Which strand(s) of the *subject* sequence to compute against the query.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    Both,
}

/// Which triangles of a self-comparison plot end up populated.
///
/// The kernel only fills the lower triangle when `self_comparison` is set
/// (matching C dotter's `qmax = sIdx + 1` cap). After the kernel,
/// `mirror_self_comparison` post-processes per this enum:
///
/// * [`Triangle::Both`]: mirror lower into upper — the default, gives a
///   symmetric plot. Matches C `DOTTER_TRIANGLE_BOTH`.
/// * [`Triangle::Upper`]: copy lower into upper, then zero the lower.
///   Matches C `DOTTER_TRIANGLE_UPPER`.
/// * [`Triangle::Lower`]: leave the lower as computed, leave the upper
///   blank. Matches C `DOTTER_TRIANGLE_LOWER` ("done by default").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Triangle {
    Both,
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
    /// Pre-process: reverse-complement the query before computation.
    /// Equivalent to the original Dotter's `-r` flag (spec §4.1.10).
    /// Only meaningful for BLASTN; ignored for BLASTP.
    pub reverse_query: bool,
    /// Pre-process: reverse-complement the subject before computation.
    /// Equivalent to the original Dotter's `-v` flag.
    pub reverse_subject: bool,
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
            triangle: Triangle::Both,
            disable_mirror: false,
            memory_limit_bytes: 512 * 1024 * 1024,
            separate_strand_channels: false,
            reverse_query: false,
            reverse_subject: false,
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
            triangle: Triangle::Both,
            disable_mirror: false,
            memory_limit_bytes: 512 * 1024 * 1024,
            separate_strand_channels: false,
            reverse_query: false,
            reverse_subject: false,
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
/// Supports BLASTN (forward, reverse, both, with optional separate
/// strand channels and self-comparison mirroring) and BLASTP. BLASTX
/// (three-frame translation) is Phase 2-extra and still returns
/// [`DottirError::NotImplemented`].
///
/// # Example — BLASTN self-comparison, both strands
///
/// ```
/// use dottir_core::{compute_dotplot, ScoreMatrix, PlotConfig, Strand};
///
/// let q = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".repeat(2);
/// let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
/// cfg.strand = Strand::Both;
/// cfg.window_size = Some(8);
/// cfg.zoom = 1;
/// cfg.self_comparison = true;
/// let plot = compute_dotplot(&q, &q, &cfg).unwrap();
/// // Main diagonal is symmetric.
/// let n = plot.width as usize;
/// for i in 8..n {
///     assert_eq!(plot.pixels[i * n + i], plot.pixels[i * n + i]);
/// }
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
    if config.mode == BlastMode::Blastx {
        return Err(DottirError::NotImplemented(
            "BLASTX three-frame translation — Phase 2-extra",
        ));
    }
    if config.self_comparison && query.len() != subject.len() {
        return Err(DottirError::InvalidConfig(
            "self_comparison requires query.len() == subject.len()".into(),
        ));
    }
    // BLASTP only does a single forward pass; reverse-strand options
    // are meaningless on proteins.
    if config.mode == BlastMode::Blastp
        && matches!(config.strand, Strand::Reverse | Strand::Both)
    {
        return Err(DottirError::InvalidConfig(
            "BLASTP only supports Strand::Forward (proteins have no reverse strand)".into(),
        ));
    }

    // Karlin window estimate or override.
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

    let alpha = config.mode.alphabet();
    // Spec §4.1.10: `-r` / `-v` reverse-complement the corresponding
    // axis sequence before computation. Only meaningful for BLASTN
    // (protein has no reverse strand). We do this on raw ASCII bytes
    // via reverse_complement_dna, then encode.
    let query_buf: std::borrow::Cow<[u8]> = if config.reverse_query
        && config.mode == BlastMode::Blastn
    {
        std::borrow::Cow::Owned(reverse_complement_dna(query))
    } else {
        std::borrow::Cow::Borrowed(query)
    };
    let subject_buf: std::borrow::Cow<[u8]> = if config.reverse_subject
        && config.mode == BlastMode::Blastn
    {
        std::borrow::Cow::Owned(reverse_complement_dna(subject))
    } else {
        std::borrow::Cow::Borrowed(subject)
    };
    let q_encoded = encode(&query_buf, alpha);
    let s_encoded_forward = encode(&subject_buf, alpha);

    let width = image_dimension(q_encoded.len(), config.zoom);
    let height = image_dimension(s_encoded_forward.len(), config.zoom);

    let score_vec = ScoreVec::build(&config.matrix, &q_encoded);

    // Per spec §4.4.3, separate-strand channels keep forward and reverse
    // hits distinct (used for inverted-repeat highlighting in the GUI).
    // Otherwise both passes max-merge into a single pixelmap.
    //
    // After Phase A1 these are *atomic* PixelMaps — at most one per
    // strand for the entire pass, shared across all rayon workers.
    let forward_map = PixelMap::new_checked(width, height, config.memory_limit_bytes)?;
    let mut reverse_map: Option<PixelMap> = None;

    // Decide which passes to run. For BLASTN, Strand::Both = forward +
    // reverse; Strand::Forward = forward only; Strand::Reverse = reverse
    // only. For BLASTP, always Forward.
    let do_forward = matches!(config.strand, Strand::Forward | Strand::Both)
        || config.mode == BlastMode::Blastp;
    let do_reverse = config.mode == BlastMode::Blastn
        && matches!(config.strand, Strand::Reverse | Strand::Both);

    if do_forward {
        run_pass(
            &score_vec,
            &s_encoded_forward,
            window,
            config.zoom,
            config.pixel_fac,
            Direction::Forward,
            config.self_comparison,
            &forward_map,
        );
    }

    if do_reverse {
        // Match C dotter's reverse pass: build the score vector against
        // the *complement* of the query (C uses the ntob_compl[] table
        // for the same effect), then scan the subject backwards.
        let q_complement = complement_encoded(&q_encoded);
        let score_vec_rev = ScoreVec::build(&config.matrix, &q_complement);
        let target_map: &PixelMap = if config.separate_strand_channels {
            reverse_map = Some(PixelMap::new_checked(
                width,
                height,
                config.memory_limit_bytes,
            )?);
            reverse_map.as_ref().unwrap()
        } else {
            &forward_map
        };
        run_pass(
            &score_vec_rev,
            &s_encoded_forward,
            window,
            config.zoom,
            config.pixel_fac,
            Direction::Reverse,
            config.self_comparison,
            target_map,
        );
    }

    // Compute is done — drop into raw bytes for post-processing. All
    // workers have completed (rayon::for_each is synchronous), so the
    // PixelMaps' atomic cells are now exclusive and consume cleanly.
    let mut forward_pixels = forward_map.into_vec();
    let mut reverse_pixels = reverse_map.map(|m| m.into_vec());

    // Self-comparison post-processing (spec §4.1.8).
    if config.self_comparison {
        mirror_self_comparison(
            &mut forward_pixels,
            width as usize,
            height as usize,
            config.triangle,
            config.disable_mirror,
        );
        if let Some(ref mut r) = reverse_pixels {
            mirror_self_comparison(
                r,
                width as usize,
                height as usize,
                config.triangle,
                config.disable_mirror,
            );
        }
    }

    let (combined, fwd_split, rev_split) =
        if config.separate_strand_channels && reverse_pixels.is_some() {
            // Combined channel: max-merge of forward + reverse for
            // downstream code that wants the unified view.
            let mut combined = forward_pixels.clone();
            if let Some(ref r) = reverse_pixels {
                merge_max_into(&mut combined, r);
            }
            (combined, Some(forward_pixels), reverse_pixels)
        } else {
            (forward_pixels, None, None)
        };

    Ok(DotPlot {
        width,
        height,
        pixels: combined,
        forward_pixels: fwd_split,
        reverse_pixels: rev_split,
        params: PlotParams {
            window_size: window,
            zoom: config.zoom,
            pixel_fac: config.pixel_fac,
            karlin: karlin_result,
        },
    })
}

/// Reverse-complement a DNA sequence in place (helper for callers who want
/// to pre-flip a sequence before passing to [`compute_dotplot`]; equivalent
/// to the original Dotter's `-r` / `-v` CLI flags per spec §4.1.10).
pub fn reverse_complement(dna: &[u8]) -> Vec<u8> {
    reverse_complement_dna(dna)
}

/// Run a single sliding-window pass into `out`, chunking on the subject
/// axis with rayon when the `rayon` feature is enabled.
///
/// Phase A1 model: there is exactly **one** [`PixelMap`] allocated for
/// the whole pass (regardless of thread count). Workers share `&out`
/// and emit pixels through [`PixelMap::max_merge`], which uses an
/// atomic `compare_exchange_weak` loop. Integer `max` is associative
/// and commutative, so the result is byte-identical to a serial pass
/// — verified by `tests/parallel_determinism.rs` at thread counts
/// 1, 2, 4, 8 (spec §4.1.11).
///
/// Memory budget: `O(W × H) + O(n_threads × qlen)` for the shared
/// pixelmap plus each worker's local ping-pong sum buffers
/// (`sum1`, `sum2` are i32). The previous design allocated
/// `n_chunks + 1` full pixelmaps; that hard-violated
/// `memory_limit_bytes` on multi-threaded runs and is fixed here.
#[allow(clippy::too_many_arguments)]
fn run_pass(
    score_vec: &ScoreVec,
    subject: &[u8],
    window: u32,
    zoom: u32,
    pixel_fac: u32,
    direction: Direction,
    self_comp: bool,
    out: &PixelMap,
) {
    #[cfg(feature = "rayon")]
    {
        use rayon::prelude::*;

        let slen = subject.len();
        let n_threads = rayon::current_num_threads().max(1);
        // Skip parallelisation if the input is small. The chunk overhead
        // (rayon scheduling + warm-up of W extra subject positions per
        // chunk) dominates and serial is faster.
        let min_for_parallel = (window as usize).saturating_mul(64).max(2048);
        if n_threads <= 1 || slen < min_for_parallel {
            sliding_window_pass_chunked(
                score_vec, subject, window, zoom, pixel_fac, direction,
                self_comp, 0..slen, out,
            );
            return;
        }

        // Aim for ~4× as many chunks as threads, so straggling workers
        // don't dominate runtime. Self-comparison's lower-triangle cap
        // means later subject positions are heavier than earlier ones,
        // and finer chunks balance the work better.
        let target_chunks = (n_threads * 4).max(2);
        let chunk_size = slen.div_ceil(target_chunks).max(window as usize);

        let chunks: Vec<std::ops::Range<usize>> = (0..slen)
            .step_by(chunk_size)
            .map(|lo| lo..(lo + chunk_size).min(slen))
            .collect();

        chunks.par_iter().for_each(|range| {
            sliding_window_pass_chunked(
                score_vec,
                subject,
                window,
                zoom,
                pixel_fac,
                direction,
                self_comp,
                range.clone(),
                out,
            );
        });
    }

    #[cfg(not(feature = "rayon"))]
    {
        sliding_window_pass_chunked(
            score_vec,
            subject,
            window,
            zoom,
            pixel_fac,
            direction,
            self_comp,
            0..subject.len(),
            out,
        );
    }
}

/// Post-process a self-comparison pixelmap according to [`Triangle`].
/// The kernel only filled the *lower* triangle (q < s, i.e. row > col)
/// thanks to the `qmax = s + 1` cap, so:
///
/// * [`Triangle::Both`]: copy lower into upper — full symmetric plot.
/// * [`Triangle::Upper`]: copy lower into upper, then zero the lower.
/// * [`Triangle::Lower`]: do nothing (already filled).
///
/// `disable_mirror = true` short-circuits to a no-op regardless of mode,
/// honouring spec §4.1.8 / the original `--disable-mirror` flag.
fn mirror_self_comparison(
    data: &mut [u8],
    width: usize,
    height: usize,
    triangle: Triangle,
    disable_mirror: bool,
) {
    if disable_mirror {
        return;
    }
    debug_assert_eq!(data.len(), width * height);
    let stride = width;
    let dim = stride.min(height);
    for s in 0..dim {
        for q in 0..s {
            // (row=s, col=q): row > col → LOWER triangle (filled by kernel).
            // (row=q, col=s): row < col → UPPER triangle (empty).
            let lower_idx = s * stride + q;
            let upper_idx = q * stride + s;
            match triangle {
                Triangle::Both => {
                    data[upper_idx] = data[lower_idx];
                }
                Triangle::Upper => {
                    data[upper_idx] = data[lower_idx];
                    data[lower_idx] = 0;
                }
                Triangle::Lower => {
                    // No-op: the kernel already filled the lower triangle.
                }
            }
        }
    }
}
