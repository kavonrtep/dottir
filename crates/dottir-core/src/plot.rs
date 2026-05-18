//! Top-level dotplot computation.
//!
//! Phase 0 ships the public signatures; Phase 1 wires the BLASTN forward
//! kernel through them; Phase 2 adds reverse-strand, BLASTP/BLASTX, and
//! self-comparison.

use crate::alphabet::{complement_encoded, encode, reverse_complement_dna};
use crate::error::DottirError;
use crate::karlin::{karlin_window_size, KarlinResult};
use crate::matrix::{BlastMode, ScoreMatrix};
use crate::pixel::{image_dimension, PixelMap, PixelView};
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
    ///
    /// `0` means **auto** — derive from Karlin's expected residue score
    /// in MSP via dotter's formula `round(0.2 * 256 / exp_res_score)`
    /// (see `dotplot.c:854`). This positions a Karlin-expected match
    /// residue at ~1/5 of the displayable range and pushes background
    /// noise toward white. The resolved value is reported back via
    /// [`PlotParams::pixel_fac`].
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
    /// Sensible BLASTN defaults: forward+reverse, zoom 1, auto pixel_fac
    /// (matches dotter — derived from Karlin's expected residue score),
    /// 0.5 GiB cap.
    pub fn default_blastn(matrix: ScoreMatrix) -> Self {
        Self {
            mode: BlastMode::Blastn,
            matrix,
            window_size: None,
            zoom: 1,
            pixel_fac: 0, // auto
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
            pixel_fac: 0, // auto
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

/// Pick a `zoom` value such that the largest output pixelmap dimension
/// (`max(qlen, slen).div_ceil(zoom)`) is at most `target_max_dim`.
///
/// Returns at least `1`. `target_max_dim` is clamped to `>= 1` to avoid
/// a divide-by-zero edge case. Used by the CLI's `--auto-zoom` flag and
/// the GUI's load-time auto-zoom (spec §4.4.8 — "avoid surprise OOMs on
/// large inputs").
///
/// # Examples
///
/// ```
/// use dottir_core::pick_auto_zoom;
/// // Tiny inputs already fit — zoom stays at 1.
/// assert_eq!(pick_auto_zoom(100, 100, 4096), 1);
/// // 27 000 residues squashed into 4096 px each axis → zoom 7.
/// assert_eq!(pick_auto_zoom(27_000, 27_000, 4096), 7);
/// // Asymmetric: pick zoom for the larger axis.
/// assert_eq!(pick_auto_zoom(1_000, 100_000, 4096), 25);
/// ```
pub fn pick_auto_zoom(qlen: usize, slen: usize, target_max_dim: u32) -> u32 {
    let target = target_max_dim.max(1) as u64;
    let max_dim = qlen.max(slen) as u64;
    if max_dim == 0 {
        return 1;
    }
    (max_dim.div_ceil(target)).max(1).min(u32::MAX as u64) as u32
}

/// Snap an automatically chosen zoom to a divisor of the supplied repeat
/// periods when a nearby divisor exists.
///
/// This is a visualization aid for tandem-repeat inputs: if the repeat
/// period is not divisible by the compute zoom, identical monomer tiles can
/// land on different pixel phases. The function is intentionally conservative:
/// it prefers divisors at or above `base_zoom` so auto-fit does not silently
/// allocate a larger pixelmap. If no coarser nearby divisor exists, it may pick
/// a finer divisor inside the tolerance window.
pub fn snap_zoom_to_period_divisor(base_zoom: u32, periods: &[usize], tolerance: f64) -> u32 {
    let base = base_zoom.max(1);
    if periods.is_empty() || !tolerance.is_finite() || tolerance < 1.0 {
        return base;
    }

    let mut common = 0usize;
    for &period in periods.iter().filter(|&&p| p > 1) {
        common = if common == 0 {
            period
        } else {
            gcd_usize(common, period)
        };
    }
    if common == 0 || common % base as usize == 0 {
        return base;
    }

    let min_zoom = ((base as f64) / tolerance).ceil().max(1.0) as u32;
    let max_zoom = ((base as f64) * tolerance).floor().max(base as f64) as u32;

    let mut coarser: Option<u32> = None;
    let mut finer: Option<u32> = None;
    let limit = (common as f64).sqrt() as usize;
    for d in 1..=limit {
        if common % d != 0 {
            continue;
        }
        for cand in [d, common / d] {
            if cand == 0 || cand > u32::MAX as usize {
                continue;
            }
            let z = cand as u32;
            if z < min_zoom || z > max_zoom {
                continue;
            }
            if z >= base {
                if coarser.is_none_or(|best| z < best) {
                    coarser = Some(z);
                }
            } else if finer.is_none_or(|best| z > best) {
                finer = Some(z);
            }
        }
    }

    coarser.or(finer).unwrap_or(base)
}

fn gcd_usize(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
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

/// A computed dotplot. `pixels.len() == width * height`.
///
/// **Channel semantics:**
/// - When [`Self::reverse_pixels`] is `None` (the common case), `pixels` is
///   the only channel — the merged "combined" view if both strands ran,
///   or the single strand's output if only one ran.
/// - When `reverse_pixels` is `Some` (set only by callers that pass
///   [`PlotConfig::separate_strand_channels`] = `true`), `pixels` holds
///   the **forward** channel and `reverse_pixels` holds the reverse
///   channel. Use [`Self::combined`] to get a merged view on demand —
///   the merge is *not* precomputed, saving ~33 % of the peak memory
///   when separate channels are requested.
#[derive(Debug, Clone)]
pub struct DotPlot {
    pub width: u32,
    pub height: u32,
    pub pixels: Vec<u8>,
    pub reverse_pixels: Option<Vec<u8>>,
    pub params: PlotParams,
}

impl DotPlot {
    /// Borrow the primary channel. Always equivalent to `&self.pixels`.
    /// Useful at call sites that want to signal "I treat this as the
    /// combined view" — when separate channels are off, `pixels` *is*
    /// the merged view; when on, callers should usually go through
    /// [`Self::combined`] instead.
    #[inline]
    pub fn primary(&self) -> &[u8] {
        &self.pixels
    }

    /// Return a merged forward+reverse view.
    ///
    /// - When `reverse_pixels.is_none()`: returns `pixels` borrowed
    ///   (no allocation) — `pixels` already *is* the merged view.
    /// - When `reverse_pixels.is_some()`: allocates and returns a fresh
    ///   buffer that is the element-wise max of `pixels` (forward) and
    ///   `reverse_pixels`.
    ///
    /// This replaces the eagerly-stored `forward + reverse + combined`
    /// triple from earlier versions of the type: the combined channel
    /// is built lazily here only when a caller asks.
    pub fn combined(&self) -> std::borrow::Cow<'_, [u8]> {
        match &self.reverse_pixels {
            None => std::borrow::Cow::Borrowed(&self.pixels),
            Some(rev) => {
                let mut out = self.pixels.clone();
                crate::pixel::merge_max_into(&mut out, rev);
                std::borrow::Cow::Owned(out)
            }
        }
    }
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
    // `config.pixel_fac == 0` means "auto-derive from Karlin's
    // expected_residue_score" — handled below once Karlin has run.
    if config.mode == BlastMode::Blastx {
        return compute_blastx(query, subject, config);
    }
    if config.self_comparison && query.len() != subject.len() {
        return Err(DottirError::InvalidConfig(
            "self_comparison requires query.len() == subject.len()".into(),
        ));
    }
    // BLASTP only does a single forward pass; reverse-strand options
    // are meaningless on proteins.
    if config.mode == BlastMode::Blastp && matches!(config.strand, Strand::Reverse | Strand::Both) {
        return Err(DottirError::InvalidConfig(
            "BLASTP only supports Strand::Forward (proteins have no reverse strand)".into(),
        ));
    }

    // Karlin window estimate (and/or `expected_residue_score` for the
    // auto pixel-factor). We need Karlin whenever EITHER the window
    // size or pixel_fac is on auto. Matches dotter `dotplot.c:1104-1106`
    // which calls `winsizeFromlambdak` even when the user fixed W in
    // order to get `expResScore`.
    let need_karlin = config.window_size.is_none() || config.pixel_fac == 0;
    let karlin_result = if need_karlin {
        Some(karlin_window_size(
            &config.matrix,
            query,
            subject,
            config.mode,
        )?)
    } else {
        None
    };
    let window = match config.window_size {
        Some(w) => w,
        None => {
            karlin_result
                .as_ref()
                .expect("karlin_result is Some when window_size is None")
                .window_size
        }
    };
    if window < 1 {
        return Err(DottirError::InvalidConfig(
            "window size must be >= 1".into(),
        ));
    }
    // Resolve `pixel_fac == 0` to dotter's auto formula
    // (dotplot.c:854): `0.2 * NUM_COLORS / expected_residue_score`. The
    // intent is that a Karlin-expected MSP residue scores at 1/5 of the
    // displayable range, leaving headroom for exceptional matches and
    // pushing background noise toward white.
    let pixel_fac = if config.pixel_fac == 0 {
        let exp_res = karlin_result
            .as_ref()
            .expect("karlin_result is Some when pixel_fac == 0")
            .expected_residue_score;
        if !(exp_res.is_finite() && exp_res > 0.0) {
            return Err(DottirError::InvalidConfig(format!(
                "auto pixel_fac requires positive Karlin expected_residue_score, got {exp_res:.3}"
            )));
        }
        ((0.2 * 256.0 / exp_res).round() as u32).max(1)
    } else {
        config.pixel_fac
    };

    let alpha = config.mode.alphabet();
    // Spec §4.1.10: `-r` / `-v` reverse-complement the corresponding
    // axis sequence before computation. Only meaningful for BLASTN
    // (protein has no reverse strand). We do this on raw ASCII bytes
    // via reverse_complement_dna, then encode.
    let query_buf: std::borrow::Cow<[u8]> =
        if config.reverse_query && config.mode == BlastMode::Blastn {
            std::borrow::Cow::Owned(reverse_complement_dna(query))
        } else {
            std::borrow::Cow::Borrowed(query)
        };
    let subject_buf: std::borrow::Cow<[u8]> =
        if config.reverse_subject && config.mode == BlastMode::Blastn {
            std::borrow::Cow::Owned(reverse_complement_dna(subject))
        } else {
            std::borrow::Cow::Borrowed(subject)
        };
    let q_encoded = encode(&query_buf, alpha);
    let s_encoded_forward = encode(&subject_buf, alpha);

    let width = image_dimension(q_encoded.len(), config.zoom);
    let height = image_dimension(s_encoded_forward.len(), config.zoom);

    // Decide which passes to run. For BLASTN, Strand::Both = forward +
    // reverse; Strand::Forward = forward only; Strand::Reverse = reverse
    // only. For BLASTP, always Forward.
    let do_forward =
        matches!(config.strand, Strand::Forward | Strand::Both) || config.mode == BlastMode::Blastp;
    let do_reverse =
        config.mode == BlastMode::Blastn && matches!(config.strand, Strand::Reverse | Strand::Both);

    // Honest up-front memory budget. The peak holding cost is one
    // pixelmap per *retained* channel — 2 when the caller asked for
    // separate forward/reverse channels and both strands run, 1
    // otherwise. The kernel's ping-pong sum buffers are i32 × qlen ×
    // n_threads ≈ negligible compared to width×height. Checking once
    // here is the source of truth; `PixelMap::new_checked`'s own
    // check becomes a secondary guard for direct callers.
    let per_channel = (width as u64) * (height as u64);
    let channels: u32 = if do_reverse && config.separate_strand_channels {
        2
    } else {
        1
    };
    let requested = per_channel.saturating_mul(channels as u64);
    if requested > config.memory_limit_bytes {
        return Err(DottirError::OutOfMemory {
            requested,
            per_channel,
            channels,
            limit: config.memory_limit_bytes,
        });
    }

    let score_vec = ScoreVec::build(&config.matrix, &q_encoded);

    // Per spec §4.4.3, separate-strand channels keep forward and reverse
    // hits distinct (used for inverted-repeat highlighting in the GUI).
    // Otherwise both passes max-merge into a single pixelmap.
    //
    // After Phase A1 these are owner-side `PixelMap`s — at most one per
    // strand for the entire pass; each pass borrows a transient
    // [`PixelView`] from its target map (atomic write lens shared
    // across all rayon workers).
    let mut forward_map = PixelMap::new_checked(width, height, config.memory_limit_bytes)?;
    let mut reverse_map: Option<PixelMap> = None;

    if do_forward {
        let view = forward_map.view_mut();
        run_pass(
            &score_vec,
            &s_encoded_forward,
            window,
            config.zoom,
            pixel_fac,
            Direction::Forward,
            config.self_comparison,
            &view,
        );
    }

    if do_reverse {
        // Match C dotter's reverse pass: build the score vector against
        // the *complement* of the query (C uses the ntob_compl[] table
        // for the same effect), then scan the subject backwards.
        let q_complement = complement_encoded(&q_encoded);
        let score_vec_rev = ScoreVec::build(&config.matrix, &q_complement);
        if config.separate_strand_channels {
            let mut rm = PixelMap::new_checked(width, height, config.memory_limit_bytes)?;
            {
                let view = rm.view_mut();
                run_pass(
                    &score_vec_rev,
                    &s_encoded_forward,
                    window,
                    config.zoom,
                    pixel_fac,
                    Direction::Reverse,
                    config.self_comparison,
                    &view,
                );
            }
            reverse_map = Some(rm);
        } else {
            let view = forward_map.view_mut();
            run_pass(
                &score_vec_rev,
                &s_encoded_forward,
                window,
                config.zoom,
                pixel_fac,
                Direction::Reverse,
                config.self_comparison,
                &view,
            );
        }
    }

    // Compute is done — drop into raw bytes for post-processing. The
    // views' lifetimes ended with each `run_pass` call, so the
    // PixelMaps are now exclusively owned and `into_vec` is zero-copy.
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

    // Channel-storage decision:
    // - separate=true with both strands: keep them split (forward in
    //   `pixels`, reverse in `reverse_pixels`). The combined view is
    //   built lazily via `DotPlot::combined()` if a caller asks.
    // - otherwise: a single channel, with `reverse_pixels = None`.
    //   When `do_reverse` ran without separate channels, it merged
    //   into the same map as forward at compute time (see the run_pass
    //   target_map decision above), so `forward_pixels` already holds
    //   the combined bytes.
    let (pixels, reverse_pixels) = if config.separate_strand_channels && reverse_pixels.is_some() {
        (forward_pixels, reverse_pixels)
    } else {
        (forward_pixels, None)
    };

    Ok(DotPlot {
        width,
        height,
        pixels,
        reverse_pixels,
        params: PlotParams {
            window_size: window,
            zoom: config.zoom,
            pixel_fac,
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
    out: &PixelView<'_>,
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
                score_vec,
                subject,
                window,
                zoom,
                pixel_fac,
                direction,
                self_comp,
                0..slen,
                out,
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

/// BLASTX: translate the DNA query in three forward reading frames,
/// run the protein kernel against each, and max-merge into a single
/// pixelmap. Subject is already in protein space. Self-comparison and
/// reverse-strand options don't apply to BLASTX (the C dotter
/// short-circuits the strand machinery for this mode).
fn compute_blastx(
    query: &[u8],
    subject: &[u8],
    config: &PlotConfig,
) -> Result<DotPlot, DottirError> {
    if config.matrix.kind != crate::alphabet::AlphabetKind::Protein {
        return Err(DottirError::InvalidMatrix(
            "BLASTX requires a protein score matrix (BLOSUM/PAM)".into(),
        ));
    }
    // Karlin estimate uses the translated frame-0 sequence vs the
    // protein subject — matches the C dotter's BLASTX choice of
    // pepQSeqLen for Karlin (it calls winsizeFromlambdak on frame 0).
    let frame0 = crate::translation::translate_frame(query, 0);
    if frame0.is_empty() {
        return Err(DottirError::InvalidConfig(
            "BLASTX query is too short to translate (<3 bases)".into(),
        ));
    }
    // Same auto-Karlin logic as compute_dotplot: run Karlin when
    // either window_size OR pixel_fac is on auto.
    let need_karlin = config.window_size.is_none() || config.pixel_fac == 0;
    let karlin_result = if need_karlin {
        Some(karlin_window_size(
            &config.matrix,
            &frame0,
            subject,
            BlastMode::Blastp,
        )?)
    } else {
        None
    };
    let window = match config.window_size {
        Some(w) => w,
        None => {
            karlin_result
                .as_ref()
                .expect("karlin_result is Some when window_size is None")
                .window_size
        }
    };
    if window < 1 {
        return Err(DottirError::InvalidConfig(
            "window size must be >= 1".into(),
        ));
    }
    let pixel_fac = if config.pixel_fac == 0 {
        let exp_res = karlin_result
            .as_ref()
            .expect("karlin_result is Some when pixel_fac == 0")
            .expected_residue_score;
        if !(exp_res.is_finite() && exp_res > 0.0) {
            return Err(DottirError::InvalidConfig(format!(
                "auto pixel_fac requires positive Karlin expected_residue_score, got {exp_res:.3}"
            )));
        }
        ((0.2 * 256.0 / exp_res).round() as u32).max(1)
    } else {
        config.pixel_fac
    };

    // Each frame has length ⌈(qlen − f) / 3⌉; the C dotter sizes the
    // image to pepQSeqLen = qlen / 3. Match that. Frames whose
    // translated length is shorter just emit nothing past their tip.
    let pepqlen = query.len() / 3;
    let alpha = crate::alphabet::AlphabetKind::Protein;
    let s_encoded = crate::alphabet::encode(subject, alpha);
    let width = image_dimension(pepqlen, config.zoom);
    let height = image_dimension(subject.len(), config.zoom);
    let mut pixmap = PixelMap::new_checked(width, height, config.memory_limit_bytes)?;

    for frame_offset in 0..3 {
        let translated = crate::translation::translate_frame(query, frame_offset);
        if translated.is_empty() {
            continue;
        }
        let translated_padded = if translated.len() < pepqlen {
            // Pad with SENTINEL so the score vector has consistent
            // width across frames. The padding columns score 0 against
            // any subject residue (the unknown-residue zero row in
            // ScoreVec handles it).
            let mut v = translated;
            v.resize(pepqlen, crate::alphabet::SENTINEL);
            v
        } else {
            translated
        };
        let q_encoded = crate::alphabet::encode(&translated_padded, alpha);
        let score_vec = ScoreVec::build(&config.matrix, &q_encoded);
        let view = pixmap.view_mut();
        run_pass(
            &score_vec,
            &s_encoded,
            window,
            config.zoom,
            pixel_fac,
            Direction::Forward,
            false, // self-comparison not supported for BLASTX
            &view,
        );
    }

    Ok(DotPlot {
        width,
        height,
        pixels: pixmap.into_vec(),
        reverse_pixels: None,
        params: PlotParams {
            window_size: window,
            zoom: config.zoom,
            pixel_fac,
            karlin: karlin_result,
        },
    })
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
