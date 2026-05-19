//! Windowed, greyramp-sensitized self-periodogram for DNA sequences.
//!
//! Conceptually: take the dotter self-comparison dotplot at zoom=1,
//! apply a greyramp-style sensitivity ramp to every pixel, and sum
//! the result along each anti-diagonal at offset `k = s - q`. The
//! output is a function `P(k)` whose peaks correspond exactly to the
//! visible diagonals in the dotplot view — same window kernel, same
//! noise floor, same saturation.
//!
//! Implementation streams the same sum recurrence used by the
//! dotplot kernel ([`crate::sliding`]) but writes into a 1-D
//! offset-indexed accumulator instead of an N × N pixmap. Memory is
//! O(N) regardless of input size, so the periodogram scales to
//! arbitrary records without needing a memory budget for the
//! pixmap itself.
//!
//! See `docs/dottir_specification.md` (future periodogram section)
//! and the design discussion that produced this module.
//!
//! # Quick example
//!
//! ```
//! use dottir_core::matrix::ScoreMatrix;
//! use dottir_core::periodogram::{
//!     compute_periodogram, PeriodogramConfig, Sensitivity,
//! };
//!
//! // Pure tandem repeat of period 7.
//! let unit = b"ACGTACA";
//! let seq: Vec<u8> = unit.iter().cycle().take(unit.len() * 30).copied().collect();
//!
//! let cfg = PeriodogramConfig {
//!     matrix: ScoreMatrix::dna_identity(),
//!     window_size: Some(7),
//!     pixel_fac: 50,
//!     sensitivity: Sensitivity::identity(),
//!     min_offset: 3,
//!     max_offset: Some(50),
//!     memory_limit_bytes: 1 << 30,
//! };
//! let p = compute_periodogram(&seq, &cfg).unwrap();
//! // Peak at k = unit.len() (and its multiples).
//! let argmax = (0..p.signal_sum.len())
//!     .max_by_key(|&i| p.signal_sum[i]).unwrap();
//! assert_eq!(cfg.min_offset + argmax as u32, unit.len() as u32);
//! ```

use crate::alphabet::{encode, AlphabetKind, SENTINEL};
use crate::error::DottirError;
use crate::karlin::{karlin_window_size, KarlinResult};
use crate::matrix::{BlastMode, ScoreMatrix};
use crate::score_vec::ScoreVec;

// ---------------------------------------------------------------------------
// Sensitivity ramp
// ---------------------------------------------------------------------------

/// Greyramp-style sensitivity ramp for periodogram contributions.
///
/// Given a per-pixel value `p` in `0..=255` (after the kernel's
/// `min(255, score * pixel_fac / W)` scaling), the contribution to
/// the periodogram is:
///
/// ```text
///     0                                     if p <= white
///     255                                   if p >= black
///     ((p - white) * 255) / (black - white) otherwise
/// ```
///
/// So `white` is the noise floor (pixels at or below it contribute
/// nothing) and `black` is saturation (pixels at or above it
/// contribute fully). This mirrors how the GUI greyramp shapes
/// display intensity — the same `(white, black)` calibration the
/// user tunes on the dotplot view transfers directly to the
/// periodogram.
///
/// [`Self::identity`] passes raw pixels through unchanged
/// (`white = 0, black = 255`) — useful when you want the unfiltered
/// windowed sum.
#[derive(Clone, Copy, Debug)]
pub struct Sensitivity {
    pub white: u8,
    pub black: u8,
}

impl Sensitivity {
    /// `white = 0, black = 255`: every pixel value `p` contributes
    /// exactly `p`. Used by the fast analytical z-score path
    /// (closed-form formulas require no clipping).
    pub const fn identity() -> Self {
        Self {
            white: 0,
            black: 255,
        }
    }

    /// GUI default — `white = 40, black = 100`. Matches
    /// `dottir-gui::Greyramp::default()`.
    pub const fn gui_default() -> Self {
        Self {
            white: 40,
            black: 100,
        }
    }

    /// True iff this ramp passes raw pixels through unchanged.
    #[inline]
    pub fn is_identity(&self) -> bool {
        self.white == 0 && self.black == 255
    }

    /// Map a raw pixel value to its periodogram contribution.
    #[inline]
    pub fn signal(&self, pixel: u8) -> u8 {
        if pixel <= self.white {
            return 0;
        }
        if pixel >= self.black {
            return 255;
        }
        let num = (pixel as u32 - self.white as u32) * 255;
        let denom = (self.black as u32 - self.white as u32).max(1);
        (num / denom) as u8
    }
}

// ---------------------------------------------------------------------------
// Config + result
// ---------------------------------------------------------------------------

/// Configuration for [`compute_periodogram`].
#[derive(Clone, Debug)]
pub struct PeriodogramConfig {
    /// DNA score matrix (4×4 over ACGT or larger).
    pub matrix: ScoreMatrix,
    /// Sliding window size in residues. `None` = auto-derive from
    /// Karlin/Altschul statistics (clamped to `[3, 50]`).
    pub window_size: Option<u32>,
    /// Multiplier in the per-pixel scale step
    /// `min(255, score * pixel_fac / W)`. `0` = auto-derive from
    /// Karlin's `expected_residue_score` (matches dotter:
    /// `0.2 * 256 / E[M]`).
    pub pixel_fac: u32,
    /// Sensitivity ramp applied to each pixel before summation.
    pub sensitivity: Sensitivity,
    /// Smallest offset `k` to report. Typically `3` to skip the
    /// trivial main-diagonal and homopolymer spikes at `k = 1, 2`.
    pub min_offset: u32,
    /// Largest offset to report. `None` = `floor(N / 2)`.
    pub max_offset: Option<u32>,
    /// Hard cap on the total bytes the computation may allocate.
    /// Streaming algorithm needs only O(N) memory; this guards
    /// against pathologically large inputs.
    pub memory_limit_bytes: u64,
}

impl PeriodogramConfig {
    /// Sensible BLASTN-style defaults: DNA identity matrix, auto
    /// window + pixel_fac, sensitivity matching the GUI greyramp
    /// default, min offset 3, 1 GiB memory cap.
    pub fn default_blastn() -> Self {
        Self {
            matrix: ScoreMatrix::dna_identity(),
            window_size: None,
            pixel_fac: 0,
            sensitivity: Sensitivity::gui_default(),
            min_offset: 3,
            max_offset: None,
            memory_limit_bytes: 1 << 30,
        }
    }
}

/// Result of [`compute_periodogram`]. Per-offset vectors are indexed
/// such that bucket `i` corresponds to offset `min_offset + i`.
#[derive(Clone, Debug)]
pub struct Periodogram {
    pub min_offset: u32,
    /// Sum of `min(255, score * pixel_fac / W)` along diagonal k.
    /// Sensitivity-independent — useful for diagnostics and the
    /// analytical z-score path.
    pub raw_sum: Vec<u64>,
    /// Sum of sensitivity-mapped pixel values along diagonal k.
    /// Primary signal: peaks here track what the dotplot view
    /// displays after the user's greyramp settings.
    pub signal_sum: Vec<u64>,
    /// Number of (q, s) pairs *visited* per offset. For a
    /// length-N sequence with window W:
    /// `n_pairs[k - min_offset] = max(0, N - W - k)`. Used as the
    /// denominator for `signal_mean = signal_sum / n_pairs`.
    pub n_pairs: Vec<u32>,
    /// Window size actually used (resolved from Karlin if auto).
    pub window_size: u32,
    /// Pixel factor actually used (resolved from Karlin if auto).
    pub pixel_fac: u32,
    /// Per-residue frequencies over ACGT, normalised to sum to 1.
    /// Bytes outside ACGT (N, soft-masked, gaps, etc.) are dropped
    /// before normalisation. Used by [`analytical_null`].
    pub residue_freqs: [f64; 4],
    /// Length of the input sequence (residues). Convenience for
    /// callers building TSV headers.
    pub seq_len: u32,
}

// ---------------------------------------------------------------------------
// compute_periodogram
// ---------------------------------------------------------------------------

/// Compute the windowed, sensitivity-shaped self-periodogram of a
/// DNA sequence.
///
/// Single-threaded; for the parallel driver see
/// [`compute_periodogram_parallel`]. Returns `Err` for empty input,
/// invalid window size, or a `(N, window, max_offset)` combination
/// that would overflow [`PeriodogramConfig::memory_limit_bytes`].
pub fn compute_periodogram(
    seq: &[u8],
    config: &PeriodogramConfig,
) -> Result<Periodogram, DottirError> {
    let r = resolve(seq, config)?;

    let mut raw_sum = vec![0u64; r.n_buckets];
    let mut signal_sum = vec![0u64; r.n_buckets];

    kernel_pass(
        &r.score_vec,
        &r.seq_encoded,
        r.window,
        r.pixel_fac,
        config.sensitivity,
        config.min_offset,
        r.max_offset,
        0..seq.len(),
        &mut raw_sum,
        &mut signal_sum,
    );

    Ok(finalise(config, seq.len(), &r, raw_sum, signal_sum))
}

/// Rayon-parallel variant. Splits the subject (and-query) axis into
/// chunks with W-position warm-up overlap, runs the kernel on each
/// chunk into a local accumulator, then sums the locals. Output is
/// byte-identical to [`compute_periodogram`] (sums are
/// associative + commutative).
///
/// Falls back to serial when `n_threads <= 1` or the input is too
/// small for the chunking overhead to be worth it.
#[cfg(feature = "rayon")]
pub fn compute_periodogram_parallel(
    seq: &[u8],
    config: &PeriodogramConfig,
) -> Result<Periodogram, DottirError> {
    use rayon::prelude::*;

    let r = resolve(seq, config)?;
    let n = seq.len();
    let n_threads = rayon::current_num_threads().max(1);
    // Heuristic mirroring `plot::run_pass`: skip parallel for inputs
    // where chunk overhead dominates.
    let min_for_parallel = (r.window as usize).saturating_mul(64).max(2048);
    if n_threads <= 1 || n < min_for_parallel {
        let mut raw_sum = vec![0u64; r.n_buckets];
        let mut signal_sum = vec![0u64; r.n_buckets];
        kernel_pass(
            &r.score_vec,
            &r.seq_encoded,
            r.window,
            r.pixel_fac,
            config.sensitivity,
            config.min_offset,
            r.max_offset,
            0..n,
            &mut raw_sum,
            &mut signal_sum,
        );
        return Ok(finalise(config, n, &r, raw_sum, signal_sum));
    }

    let target_chunks = (n_threads * 4).max(2);
    let chunk_size = n.div_ceil(target_chunks).max(r.window as usize);
    let chunks: Vec<std::ops::Range<usize>> = (0..n)
        .step_by(chunk_size)
        .map(|lo| lo..(lo + chunk_size).min(n))
        .collect();

    let (raw_sum, signal_sum) = chunks
        .into_par_iter()
        .map(|range| {
            let mut local_raw = vec![0u64; r.n_buckets];
            let mut local_sig = vec![0u64; r.n_buckets];
            kernel_pass(
                &r.score_vec,
                &r.seq_encoded,
                r.window,
                r.pixel_fac,
                config.sensitivity,
                config.min_offset,
                r.max_offset,
                range,
                &mut local_raw,
                &mut local_sig,
            );
            (local_raw, local_sig)
        })
        .reduce(
            || (vec![0u64; r.n_buckets], vec![0u64; r.n_buckets]),
            |(mut a_raw, mut a_sig), (b_raw, b_sig)| {
                for i in 0..a_raw.len() {
                    a_raw[i] += b_raw[i];
                    a_sig[i] += b_sig[i];
                }
                (a_raw, a_sig)
            },
        );

    Ok(finalise(config, n, &r, raw_sum, signal_sum))
}

#[cfg(not(feature = "rayon"))]
pub fn compute_periodogram_parallel(
    seq: &[u8],
    config: &PeriodogramConfig,
) -> Result<Periodogram, DottirError> {
    compute_periodogram(seq, config)
}

// ---------------------------------------------------------------------------
// Internal: shared resolution + final assembly
// ---------------------------------------------------------------------------

struct Resolved {
    score_vec: ScoreVec,
    seq_encoded: Vec<u8>,
    residue_freqs: [f64; 4],
    window: u32,
    pixel_fac: u32,
    n_buckets: usize,
    max_offset: u32,
}

fn resolve(seq: &[u8], config: &PeriodogramConfig) -> Result<Resolved, DottirError> {
    if seq.is_empty() {
        return Err(DottirError::EmptySequence);
    }
    if config.matrix.kind != AlphabetKind::Dna {
        return Err(DottirError::InvalidConfig(format!(
            "periodogram requires a DNA matrix, got {:?}",
            config.matrix.kind
        )));
    }
    if config.min_offset == 0 {
        return Err(DottirError::InvalidConfig(
            "min_offset must be >= 1 (k=0 is the trivial main diagonal)".into(),
        ));
    }

    // Auto window + pixel_fac via Karlin (self-comparison: query == subject).
    let need_karlin = config.window_size.is_none() || config.pixel_fac == 0;
    let karlin: Option<KarlinResult> = if need_karlin {
        Some(karlin_window_size(
            &config.matrix,
            seq,
            seq,
            BlastMode::Blastn,
        )?)
    } else {
        None
    };
    let window = match config.window_size {
        Some(w) => w,
        None => karlin.as_ref().unwrap().window_size,
    };
    if window < 1 {
        return Err(DottirError::InvalidConfig(
            "window size must be >= 1".into(),
        ));
    }
    let pixel_fac = if config.pixel_fac == 0 {
        let e = karlin.as_ref().unwrap().expected_residue_score;
        if !(e.is_finite() && e > 0.0) {
            return Err(DottirError::InvalidConfig(format!(
                "auto pixel_fac requires positive expected_residue_score, got {e}"
            )));
        }
        ((0.2 * 256.0 / e).round() as u32).max(1)
    } else {
        config.pixel_fac
    };

    let n = seq.len() as u32;
    let max_offset = config.max_offset.unwrap_or(n / 2);
    if max_offset < config.min_offset {
        return Err(DottirError::InvalidConfig(format!(
            "max_offset {} < min_offset {}",
            max_offset, config.min_offset
        )));
    }
    let n_buckets = (max_offset - config.min_offset + 1) as usize;

    // Memory check: two row buffers (i32 × N) + two output vecs (u64 × n_buckets).
    let needed = (n as u64).saturating_mul(4 * 2) + (n_buckets as u64).saturating_mul(8 * 2);
    if needed > config.memory_limit_bytes {
        return Err(DottirError::OutOfMemory {
            requested: needed,
            per_channel: needed,
            channels: 1,
            limit: config.memory_limit_bytes,
        });
    }

    let seq_encoded = encode(seq, AlphabetKind::Dna);
    let residue_freqs = residue_freqs_dna(&seq_encoded);
    let score_vec = ScoreVec::build(&config.matrix, &seq_encoded);
    Ok(Resolved {
        score_vec,
        seq_encoded,
        residue_freqs,
        window,
        pixel_fac,
        n_buckets,
        max_offset,
    })
}

fn finalise(
    config: &PeriodogramConfig,
    seq_len: usize,
    r: &Resolved,
    raw_sum: Vec<u64>,
    signal_sum: Vec<u64>,
) -> Periodogram {
    let n = seq_len as i64;
    let w = r.window as i64;
    let n_pairs: Vec<u32> = (config.min_offset..=r.max_offset)
        .map(|k| (n - w - k as i64).max(0) as u32)
        .collect();
    Periodogram {
        min_offset: config.min_offset,
        raw_sum,
        signal_sum,
        n_pairs,
        window_size: r.window,
        pixel_fac: r.pixel_fac,
        residue_freqs: r.residue_freqs,
        seq_len: seq_len as u32,
    }
}

/// Compute ACGT frequencies (normalised to sum to 1) from an encoded
/// sequence. Bytes outside `0..=3` are skipped (N, ambiguity codes,
/// sentinel bytes). Returns `[0, 0, 0, 0]` if no scorable residues
/// were seen.
fn residue_freqs_dna(encoded: &[u8]) -> [f64; 4] {
    let mut counts = [0u64; 4];
    let mut total = 0u64;
    for &b in encoded {
        if b == SENTINEL {
            continue;
        }
        let i = b as usize;
        if i < 4 {
            counts[i] += 1;
            total += 1;
        }
    }
    if total == 0 {
        return [0.0; 4];
    }
    let inv = 1.0 / total as f64;
    [
        counts[0] as f64 * inv,
        counts[1] as f64 * inv,
        counts[2] as f64 * inv,
        counts[3] as f64 * inv,
    ]
}

// ---------------------------------------------------------------------------
// Streaming kernel
// ---------------------------------------------------------------------------

/// Streaming variant of the dotter self-comparison kernel. Walks
/// the sum recurrence over the (q, s) pairs but, instead of
/// max-merging into a pixmap, accumulates each pixel's scaled value
/// into `raw_sum[k - min_offset]` and its sensitivity-mapped value
/// into `signal_sum[k - min_offset]`, where `k = s - q`.
///
/// `s_emit_range` specifies the subject positions whose pixels
/// should be emitted (matches the dotplot kernel's chunked emit
/// region for rayon driver). Warm-up positions before the range
/// are walked silently to populate the ping-pong buffers.
#[allow(clippy::too_many_arguments)]
fn kernel_pass(
    score_vec: &ScoreVec,
    seq_encoded: &[u8],
    window: u32,
    pixel_fac: u32,
    sensitivity: Sensitivity,
    min_offset: u32,
    max_offset: u32,
    s_emit_range: std::ops::Range<usize>,
    raw_sum: &mut [u64],
    signal_sum: &mut [u64],
) {
    let qlen = score_vec.qlen;
    let n = seq_encoded.len();
    let w = window as usize;
    if w == 0 || qlen == 0 || n == 0 || s_emit_range.is_empty() {
        return;
    }

    let mut sum1 = vec![0i32; qlen];
    let mut sum2 = vec![0i32; qlen];

    // Extend the iteration start by `w` positions so the chunk's
    // ping-pong buffers are fully populated by the time we enter
    // the emit range. Matches `sliding_window_pass_chunked`.
    let s_iter_start = s_emit_range.start.saturating_sub(w);
    let s_iter_end = s_emit_range.end;

    let mut iter_idx: usize = 0;
    for s in s_iter_start..s_iter_end {
        let s_row = score_vec.subject_row(seq_encoded[s]) as usize;
        let add_row = score_vec.row(s_row);
        let del_row: &[i32] = if iter_idx < w {
            score_vec.row(score_vec.unknown_row() as usize)
        } else {
            let prev_s = s - w;
            let prev_row = score_vec.subject_row(seq_encoded[prev_s]) as usize;
            score_vec.row(prev_row)
        };

        let (oldsum, newsum): (&[i32], &mut [i32]) = if iter_idx & 1 == 1 {
            (sum2.as_slice(), sum1.as_mut_slice())
        } else {
            (sum1.as_slice(), sum2.as_mut_slice())
        };

        // Initialise newsum[0] from the first column of addrow.
        newsum[0] = add_row[0];
        // Warm-up q region: no delrow subtraction.
        let q_warmup_end = w.min(qlen);
        for q in 1..q_warmup_end {
            newsum[q] = oldsum[q - 1] + add_row[q];
        }

        if w < qlen {
            let s_valid = iter_idx >= w;
            // Self-comparison cap matches the dotplot kernel: only
            // visit q in [W, s+1) so the lower triangle is processed.
            // Symmetric (k > 0) is sufficient — upper triangle would
            // double-count every offset.
            let q_end = qlen.min(s + 1);

            for q in w..q_end {
                let new_val = oldsum[q - 1] + add_row[q] - del_row[q - w];
                newsum[q] = new_val;

                if !s_valid || new_val <= 0 {
                    continue;
                }
                if !s_emit_range.contains(&s) {
                    continue;
                }
                let k_signed = s as i64 - q as i64;
                if k_signed < min_offset as i64 || k_signed > max_offset as i64 {
                    continue;
                }
                let k = k_signed as u32;
                let scaled_signed = new_val as i64 * pixel_fac as i64 / window as i64;
                let scaled = scaled_signed.clamp(0, 255) as u8;
                let idx = (k - min_offset) as usize;
                raw_sum[idx] += scaled as u64;
                signal_sum[idx] += sensitivity.signal(scaled) as u64;
            }
        }

        iter_idx += 1;
    }
}

// ---------------------------------------------------------------------------
// z-score: analytical null
// ---------------------------------------------------------------------------

/// Closed-form per-pair mean and variance under a shuffle null
/// based on residue composition. Cheap (O(4²)) and approximate
/// (treats overlapping windows and clipping as independent).
/// Sufficient for picking peaks; underestimates the true variance
/// for non-identity sensitivity, so the auto z-score policy
/// switches to [`empirical_null_stats`] when the sensitivity ramp
/// clips.
#[derive(Clone, Copy, Debug)]
pub struct AnalyticalNull {
    /// `E[M(X, Y)]` where X, Y are independent residues drawn from
    /// the per-record frequency distribution.
    pub mean_per_pair: f64,
    /// `Var[M(X, Y)]`.
    pub var_per_pair: f64,
}

/// Compute the per-pair mean / variance of the score matrix under
/// the supplied residue-frequency null. Only the 4×4 ACGT block
/// is consulted (DNA-only).
pub fn analytical_null(matrix: &ScoreMatrix, freqs: &[f64; 4]) -> AnalyticalNull {
    let mut mean = 0.0;
    let mut sec = 0.0;
    for i in 0..4 {
        for j in 0..4 {
            let p = freqs[i] * freqs[j];
            let s = matrix.get(i, j) as f64;
            mean += p * s;
            sec += p * s * s;
        }
    }
    AnalyticalNull {
        mean_per_pair: mean,
        var_per_pair: (sec - mean * mean).max(0.0),
    }
}

/// Analytical z-score per offset for the `raw_sum` channel.
///
/// Approximation:
/// * `E[P_raw(k)]   ≈ n_pairs[k] * pixel_fac * mean_per_pair`
/// * `Var[P_raw(k)] ≈ pixel_fac² * n_pairs[k] * var_per_pair`
///
/// The pixel-fac dependence comes from the per-pixel scaling
/// `scaled = round(raw * pixel_fac / W)`; the `W²` from windowed
/// summation cancels against the `1/W²` from the scaling.
///
/// **Caveats:** the kernel's `new_val > 0` gate clips negative
/// scores, so for matrices with `mean_per_pair < 0` (DNA identity
/// is `-1.75` under uniform composition) the analytical mean
/// overestimates and the true distribution is skewed. Use
/// [`empirical_null_stats`] when accuracy matters.
pub fn analytical_z_scores(p: &Periodogram, null: &AnalyticalNull) -> Vec<f64> {
    let pf = p.pixel_fac as f64;
    let mu_per_pair = pf * null.mean_per_pair;
    let var_per_pair = pf * pf * null.var_per_pair;
    p.raw_sum
        .iter()
        .zip(p.n_pairs.iter())
        .map(|(&sum, &np)| {
            if np == 0 {
                return 0.0;
            }
            let np = np as f64;
            let mean = np * mu_per_pair;
            let sd = (np * var_per_pair).sqrt();
            if sd <= f64::EPSILON {
                0.0
            } else {
                (sum as f64 - mean) / sd
            }
        })
        .collect()
}

// ---------------------------------------------------------------------------
// z-score: empirical null via sequence shuffles
// ---------------------------------------------------------------------------

/// Per-offset `(mean, std)` of the periodogram's `signal_sum`
/// under `shuffles` Fisher-Yates shuffles of the input. Compute
/// `z(k) = (signal_sum[k] - mean[k]) / std[k]` (caller-side, to
/// avoid handing back a third vector).
///
/// Uses an internal `xorshift64` PRNG seeded from `seed`; passing
/// the same seed twice produces byte-identical output, which makes
/// regression-testing the empirical path tractable.
///
/// Cost: `shuffles` × the cost of [`compute_periodogram`]. The
/// shuffle itself is O(N) per iteration; the periodogram dominates.
pub fn empirical_null_stats(
    seq: &[u8],
    config: &PeriodogramConfig,
    shuffles: u32,
    seed: u64,
) -> Result<Vec<(f64, f64)>, DottirError> {
    if shuffles == 0 {
        return Err(DottirError::InvalidConfig(
            "empirical null requires shuffles >= 1".into(),
        ));
    }
    // One pass to size the output vec + ensure config validity.
    let probe = compute_periodogram(seq, config)?;
    let n_buckets = probe.signal_sum.len();

    let mut sum = vec![0f64; n_buckets];
    let mut sum_sq = vec![0f64; n_buckets];

    let mut rng = Xorshift64::new(seed);
    // Re-use one buffer across shuffles.
    let mut shuffled = seq.to_vec();
    for _ in 0..shuffles {
        rng.shuffle(&mut shuffled);
        let p = compute_periodogram(&shuffled, config)?;
        for i in 0..n_buckets {
            let x = p.signal_sum[i] as f64;
            sum[i] += x;
            sum_sq[i] += x * x;
        }
    }
    let n = shuffles as f64;
    Ok((0..n_buckets)
        .map(|i| {
            let mean = sum[i] / n;
            let var = (sum_sq[i] / n - mean * mean).max(0.0);
            (mean, var.sqrt())
        })
        .collect())
}

/// 64-bit xorshift PRNG. Sufficient for shuffling a DNA sequence
/// (no cryptographic claims, no concern over modulo bias at this
/// scale). Seed `0` is reset to a non-zero default since xorshift
/// gets stuck at 0.
#[derive(Clone, Copy, Debug)]
struct Xorshift64(u64);

impl Xorshift64 {
    fn new(seed: u64) -> Self {
        // Xorshift64 is stuck at the all-zeros state; pick a
        // canonical non-zero default when the user passes 0.
        Self(if seed == 0 {
            0x9E37_79B9_7F4A_7C15
        } else {
            seed
        })
    }

    #[inline]
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// Fisher-Yates shuffle of `slice`. Modulo bias at this scale
    /// (slice lengths up to 2³² and a 64-bit PRNG) is < 1 part in
    /// 2³²; negligible for periodogram stats.
    fn shuffle<T>(&mut self, slice: &mut [T]) {
        for i in (1..slice.len()).rev() {
            let j = (self.next_u64() as usize) % (i + 1);
            slice.swap(i, j);
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a length-`n` tandem repeat of `unit`. Truncates if the
    /// repeat doesn't divide evenly.
    fn tandem(unit: &[u8], n: usize) -> Vec<u8> {
        unit.iter().cycle().take(n).copied().collect()
    }

    fn cfg(window: u32, max_offset: u32) -> PeriodogramConfig {
        PeriodogramConfig {
            matrix: ScoreMatrix::dna_identity(),
            window_size: Some(window),
            // Fixed pixel_fac so tests are independent of Karlin (which
            // would refuse to run on a perfectly periodic sequence
            // because `lambda` is ill-defined).
            pixel_fac: 50,
            sensitivity: Sensitivity::identity(),
            min_offset: 3,
            max_offset: Some(max_offset),
            memory_limit_bytes: 1 << 30,
        }
    }

    #[test]
    fn sensitivity_identity_passes_through() {
        let s = Sensitivity::identity();
        for v in [0u8, 1, 40, 100, 200, 255] {
            assert_eq!(s.signal(v), v, "identity should be a no-op at v={v}");
        }
    }

    #[test]
    fn sensitivity_clips_and_ramps() {
        let s = Sensitivity::gui_default(); // white=40, black=100
        let denom = 60u32; // black - white
        assert_eq!(s.signal(0), 0);
        assert_eq!(s.signal(40), 0);
        assert_eq!(s.signal(41), (255 / denom) as u8);
        assert_eq!(s.signal(70), (30 * 255 / denom) as u8); // ~127
        assert_eq!(s.signal(100), 255);
        assert_eq!(s.signal(255), 255);
    }

    #[test]
    fn empty_input_errors() {
        let err = compute_periodogram(&[], &cfg(7, 20)).unwrap_err();
        assert!(matches!(err, DottirError::EmptySequence));
    }

    #[test]
    fn peak_lands_at_repeat_period() {
        // Period-7 tandem repeat, 30 copies (210 residues).
        let seq = tandem(b"ACGTACA", 210);
        let p = compute_periodogram(&seq, &cfg(7, 80)).unwrap();
        let argmax = p
            .signal_sum
            .iter()
            .enumerate()
            .max_by_key(|(_, &v)| v)
            .unwrap()
            .0;
        let k = p.min_offset + argmax as u32;
        assert_eq!(k, 7, "expected peak at the unit length 7, got k={k}");
    }

    #[test]
    fn peak_height_grows_with_period_multiples() {
        // For a perfect period-P repeat, P, 2P, 3P, ... should also
        // peak (anti-diagonals at multiples of the period). The
        // dominant peak is still P, but multiples are clearly above
        // background.
        let p_unit = 5;
        let seq = tandem(b"ACGTA", p_unit as usize * 60);
        let p = compute_periodogram(&seq, &cfg(p_unit, 30)).unwrap();
        let off = |k: u32| -> u64 { p.signal_sum[(k - p.min_offset) as usize] };
        // Peak at P.
        assert!(off(5) > off(3));
        assert!(off(5) > off(4));
        assert!(off(5) > off(6));
        // 10 and 15 are above 11 and 16 (non-multiples).
        assert!(off(10) > off(11));
        assert!(off(15) > off(16));
    }

    #[test]
    fn random_sequence_has_no_dominant_peak() {
        // A sequence shuffled deterministically should produce a
        // periodogram dominated by noise, no single peak more than
        // 3× the median signal.
        let mut seq = tandem(b"ACGT", 4000);
        // Fisher-Yates with the local Xorshift to keep this hermetic.
        let mut rng = Xorshift64::new(42);
        rng.shuffle(&mut seq);
        let p = compute_periodogram(&seq, &cfg(8, 100)).unwrap();
        let mut signals: Vec<u64> = p.signal_sum.clone();
        signals.sort_unstable();
        let median = signals[signals.len() / 2];
        let max = *signals.last().unwrap();
        assert!(
            (max as f64) < (median as f64) * 3.0 + 1.0,
            "shuffled sequence should not have peaks > 3× median; max={max} median={median}"
        );
    }

    #[test]
    fn parallel_matches_serial() {
        let seq = tandem(b"ACGTAGT", 3000);
        let cfg = cfg(7, 200);
        let serial = compute_periodogram(&seq, &cfg).unwrap();
        let parallel = compute_periodogram_parallel(&seq, &cfg).unwrap();
        assert_eq!(serial.raw_sum, parallel.raw_sum);
        assert_eq!(serial.signal_sum, parallel.signal_sum);
        assert_eq!(serial.n_pairs, parallel.n_pairs);
    }

    #[test]
    fn identity_sensitivity_is_passthrough() {
        let seq = tandem(b"ACGTACA", 210);
        let p = compute_periodogram(&seq, &cfg(7, 50)).unwrap();
        // With Sensitivity::identity() (white=0, black=255) every
        // pixel value `p` contributes exactly `p`, so signal_sum
        // matches raw_sum element-wise.
        assert_eq!(p.raw_sum, p.signal_sum);
    }

    #[test]
    fn sensitivity_preserves_peak_and_kills_subthreshold() {
        // GUI default sensitivity (white=40, black=100) stretches the
        // 40..100 range across the full 0..255 output. That filters
        // sub-noise contributions to 0 and saturates strong matches —
        // exactly the dotplot-display denoising we want. The peak at
        // the true period should still be the argmax.
        let seq = tandem(b"ACGTACA", 210);
        let mut c = cfg(7, 50);
        c.sensitivity = Sensitivity::gui_default();
        let p = compute_periodogram(&seq, &c).unwrap();
        let argmax = p
            .signal_sum
            .iter()
            .enumerate()
            .max_by_key(|(_, &v)| v)
            .unwrap()
            .0;
        assert_eq!(p.min_offset + argmax as u32, 7);
    }

    #[test]
    fn analytical_null_matches_hand_calc() {
        // Uniform composition + DNA identity matrix:
        //   E[M]   = 5 * Σ f² - 4 * Σ_{a≠b} f_a f_b
        //          = 5 * 4 * (1/16) - 4 * 12 * (1/16)
        //          = 5/4 - 3 = -1.75
        //   Var[M] = E[M²] - E[M]²
        //   E[M²]  = 25 * 4 * 1/16 + 16 * 12 * 1/16 = 25/4 + 12 = 18.25
        //   Var    = 18.25 - 1.75² = 15.1875
        let freqs = [0.25, 0.25, 0.25, 0.25];
        let n = analytical_null(&ScoreMatrix::dna_identity(), &freqs);
        assert!((n.mean_per_pair - -1.75).abs() < 1e-9);
        assert!((n.var_per_pair - 15.1875).abs() < 1e-9);
    }

    #[test]
    fn empirical_null_seeded_is_reproducible() {
        let seq = tandem(b"ACGTACA", 200);
        let cfg = cfg(7, 30);
        let a = empirical_null_stats(&seq, &cfg, 10, 42).unwrap();
        let b = empirical_null_stats(&seq, &cfg, 10, 42).unwrap();
        assert_eq!(a.len(), b.len());
        for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
            assert!(
                (x.0 - y.0).abs() < 1e-9 && (x.1 - y.1).abs() < 1e-9,
                "mismatch at i={i}: {x:?} vs {y:?}"
            );
        }
    }

    #[test]
    fn n_pairs_formula() {
        let seq = tandem(b"ACGTACA", 100);
        let p = compute_periodogram(&seq, &cfg(7, 40)).unwrap();
        // n_pairs[k - 3] = max(0, N - W - k) = max(0, 100 - 7 - k)
        for (i, &np) in p.n_pairs.iter().enumerate() {
            let k = p.min_offset + i as u32;
            let expected = (100i64 - 7 - k as i64).max(0) as u32;
            assert_eq!(np, expected, "k={k}");
        }
    }

    #[test]
    fn memory_limit_rejects_oversize() {
        let seq = tandem(b"ACGT", 100_000);
        let mut cfg = cfg(8, 50_000);
        cfg.memory_limit_bytes = 1024;
        let err = compute_periodogram(&seq, &cfg).unwrap_err();
        assert!(matches!(err, DottirError::OutOfMemory { .. }));
    }
}
