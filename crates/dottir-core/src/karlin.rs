//! Karlin/Altschul statistics: λ, K, H and the derived window size.
//!
//! This module is a faithful port of `dotterApp/dotterKarlin.c` (Sonnhammer
//! & Durbin 1995, derived from Stephen Altschul's BLAST). The numerical
//! procedure is preserved bit-for-bit where possible:
//!
//! * Bisection on λ with 25 fixed iterations.
//! * `MAXIT = 20`, `SUMLIMIT = 0.01` for the K series.
//! * `K` "fudge" to 0.1 when the geometric correction ratio is too close
//!   to 1.
//!
//! Loop structures are translated as-is, with pointer arithmetic replaced by
//! plain slice indexing. Per CLAUDE.md and spec §4.1.11, future changes
//! here MUST be guarded by golden tests; non-obvious deviations require an
//! ADR.
//!
//! References:
//!
//! * Karlin, S. & Altschul, S. F. (1990) "Methods for Assessing the
//!   Statistical Significance of Molecular Sequence Features by Using
//!   General Scoring Schemes." PNAS 87, 2264–2268.
//! * `dotterApp/dotterKarlin.c:144` (`karlin`), `:343` (`winsizeFromlambdak`).
//!
//! ## Lints
//!
//! `needless_range_loop` is allowed for the doubly-nested matrix walks —
//! the C source uses the same `for (i = 0; i < n; ++i)` shape and a
//! pair of `matrix.get(i, j)` calls reads more naturally than the
//! enumerated-iterator alternative when both indices are used.

#![allow(clippy::needless_range_loop)]

use crate::alphabet::{encode_dna, encode_protein, SENTINEL};
use crate::error::DottirError;
use crate::matrix::{BlastMode, ScoreMatrix};

/// Fixed iteration cap on the K series (C: `MAXIT`).
const MAXIT: usize = 20;
/// Convergence threshold on the partial sum (C: `SUMLIMIT`).
const SUMLIMIT: f64 = 0.01;
/// Fixed iteration count for the λ bisection (C: hard-coded 25).
const LAMBDA_ITERATIONS: usize = 25;
/// Nominal dot-matrix size used to compute expected MSP score
/// (C: local `n = 100` in `winsizeFromlambdak`).
const NOMINAL_MATRIX_DIM: u32 = 100;

/// Tunables for the window-size estimate.
#[derive(Debug, Clone, Copy)]
pub struct KarlinConfig {
    /// Lower clamp on the resulting window size (default 3). Spec §4.1.1.
    pub min_window: u32,
    /// Upper clamp on the resulting window size (default 50). Spec §4.1.1.
    pub max_window: u32,
    /// Nominal MSP matrix side (default 100). The expected MSP score is
    /// `(log(n*n) + log(K)) / λ`. Exposed for tests / experimentation.
    pub nominal_matrix_dim: u32,
}

impl Default for KarlinConfig {
    fn default() -> Self {
        Self {
            min_window: 3,
            max_window: 50,
            nominal_matrix_dim: NOMINAL_MATRIX_DIM,
        }
    }
}

/// Output of [`karlin_window_size`]. All four numbers are exposed so the GUI
/// can show them and so golden tests can pin them.
#[derive(Debug, Clone, Copy)]
pub struct KarlinResult {
    pub lambda: f64,
    pub k: f64,
    pub h: f64,
    /// Predicted residue score in an MSP, i.e. `Σ q_i q_j s_ij exp(λ s_ij)`.
    pub expected_residue_score: f64,
    /// Expected MSP score for a `n × n` random matrix, where `n` is
    /// [`KarlinConfig::nominal_matrix_dim`].
    pub expected_msp_score: f64,
    /// `expected_msp_score / expected_residue_score` rounded; pre-clamp.
    pub predicted_msp_length: u32,
    /// The clamped window size returned to the caller.
    pub window_size: u32,
}

/// Compute the Karlin/Altschul window size for a query/subject pair.
///
/// `mode` selects the alphabet and karlin bucket size. Sequences are encoded
/// with the same ASCII→index tables as C dotter; any byte not in the
/// alphabet is dropped from the residue count (spec §4.1.1).
///
/// # Example — DNA (BLASTN)
///
/// ```
/// use dottir_core::karlin::karlin_window_size;
/// use dottir_core::matrix::{BlastMode, ScoreMatrix};
///
/// let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
/// let subject = query;
/// let matrix = ScoreMatrix::dna_identity();
/// let r = karlin_window_size(&matrix, query, subject, BlastMode::Blastn).unwrap();
/// // Window is always in [3, 50] by default (spec §4.1.1).
/// assert!((3..=50).contains(&r.window_size));
/// // λ, K, H are positive for well-formed inputs.
/// assert!(r.lambda > 0.0 && r.k > 0.0 && r.h > 0.0);
/// ```
///
/// # Example — protein (BLASTP)
///
/// ```
/// use dottir_core::karlin::karlin_window_size;
/// use dottir_core::matrix::{BlastMode, ScoreMatrix};
///
/// let q = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVK";
/// let s = b"MAATKRIIRQRYTIKHYVTRLREHIDHEEQVRKDLDEHKHRADRMLEELAGAILAAEH";
/// let r = karlin_window_size(
///     &ScoreMatrix::blosum62(),
///     q, s,
///     BlastMode::Blastp,
/// ).unwrap();
/// assert!(r.lambda > 0.1 && r.lambda < 1.0);
/// ```
pub fn karlin_window_size(
    matrix: &ScoreMatrix,
    query: &[u8],
    subject: &[u8],
    mode: BlastMode,
) -> Result<KarlinResult, DottirError> {
    karlin_window_size_with(matrix, query, subject, mode, KarlinConfig::default())
}

/// Same as [`karlin_window_size`] but with overridable [`KarlinConfig`].
pub fn karlin_window_size_with(
    matrix: &ScoreMatrix,
    query: &[u8],
    subject: &[u8],
    mode: BlastMode,
    cfg: KarlinConfig,
) -> Result<KarlinResult, DottirError> {
    if matrix.kind != mode.alphabet() {
        return Err(DottirError::InvalidMatrix(format!(
            "matrix alphabet {:?} doesn't match BLAST mode {:?}",
            matrix.kind, mode
        )));
    }
    let abetsize = mode.karlin_size();
    if matrix.size() < abetsize {
        return Err(DottirError::InvalidMatrix(format!(
            "matrix size {} < karlin size {}",
            matrix.size(),
            abetsize
        )));
    }
    if cfg.min_window == 0 || cfg.min_window > cfg.max_window {
        return Err(DottirError::InvalidConfig(format!(
            "min_window={} max_window={} invalid",
            cfg.min_window, cfg.max_window
        )));
    }

    let (fq_q, _qlen) = residue_frequencies(query, mode, abetsize)?;
    let (fq_s, _slen) = residue_frequencies(subject, mode, abetsize)?;

    // Lowest / highest score across the abetsize-by-abetsize submatrix.
    let mut lows: i64 = 0;
    let mut highs: i64 = 0;
    for i in 0..abetsize {
        for j in 0..abetsize {
            let s = matrix.get(i, j) as i64;
            if s < lows {
                lows = s;
            }
            if s > highs {
                highs = s;
            }
        }
    }

    let range = (highs - lows) as usize;
    let mut prob = vec![0.0_f64; range + 1];
    for i in 0..abetsize {
        for j in 0..abetsize {
            let s = matrix.get(i, j) as i64;
            prob[(s - lows) as usize] += fq_q[i] * fq_s[j];
        }
    }

    let (lambda, k, h) = karlin(lows, highs, &mut prob)?;

    // Expected per-residue score in an MSP: Σ q_i q_j s_ij exp(λ s_ij).
    let mut exp_res_score = 0.0_f64;
    let mut sum = 0.0_f64;
    for i in 0..abetsize {
        for j in 0..abetsize {
            let s = matrix.get(i, j) as f64;
            let qij = fq_q[i] * fq_s[j] * (lambda * s).exp();
            sum += qij;
            exp_res_score += qij * s;
        }
    }
    let _ = sum; // C only warns when |sum-1| > 1e-4; we don't gate on it.

    let n = cfg.nominal_matrix_dim as f64;
    let exp_msp_score = ((n * n).ln() + k.ln()) / lambda;
    if exp_res_score <= 0.0 {
        return Err(DottirError::KarlinFailure(format!(
            "expected residue score in MSP is non-positive ({exp_res_score:.5})"
        )));
    }
    let predicted = (exp_msp_score / exp_res_score + 0.5).floor() as i64;
    let predicted_u = predicted.max(0) as u32;
    let window = predicted_u.clamp(cfg.min_window, cfg.max_window);

    Ok(KarlinResult {
        lambda,
        k,
        h,
        expected_residue_score: exp_res_score,
        expected_msp_score: exp_msp_score,
        predicted_msp_length: predicted_u,
        window_size: window,
    })
}

/// Count scorable residues in a sequence and normalise to frequencies.
///
/// Mirrors `winsizeFromlambdak`'s residue-counting loop: only bytes whose
/// alphabet index is `< abetsize` contribute (so 'N' is excluded for DNA;
/// for protein the C code's quirk of bucketing unknowns into the '*' slot
/// is preserved by encoding `*` to index 23 in [`crate::alphabet`]).
fn residue_frequencies(
    seq: &[u8],
    mode: BlastMode,
    abetsize: usize,
) -> Result<(Vec<f64>, usize), DottirError> {
    if seq.is_empty() {
        return Err(DottirError::EmptySequence);
    }
    let mut counts = vec![0_u64; abetsize];
    let mut total = 0_usize;
    for &b in seq {
        let idx = match mode {
            BlastMode::Blastn => encode_dna(b),
            BlastMode::Blastp | BlastMode::Blastx => encode_protein(b),
        };
        if idx != SENTINEL && (idx as usize) < abetsize {
            counts[idx as usize] += 1;
            total += 1;
        }
    }
    if total == 0 {
        return Err(DottirError::NoScorableResidues);
    }
    let denom = total as f64;
    Ok((
        counts.into_iter().map(|c| c as f64 / denom).collect(),
        total,
    ))
}

/// The Karlin/Altschul λ/K/H solver.
///
/// Faithful port of `dotterKarlin.c:144`. The input `pr` is a probability
/// vector of length `range + 1` indexed by `(score − low)`. On return the
/// vector has been normalised in place.
fn karlin(low: i64, high: i64, pr: &mut [f64]) -> Result<(f64, f64, f64), DottirError> {
    if low >= 0 {
        return Err(DottirError::KarlinFailure(
            "no negative score in substitution matrix".into(),
        ));
    }
    let range = (high - low) as usize;
    debug_assert_eq!(pr.len(), range + 1);

    // Find the rightmost positive-score probability. If everything from -low+1
    // upwards is zero, there's no positive score and Karlin can't proceed.
    {
        let neg_low = (-low) as usize;
        let mut i = range;
        while i > neg_low && pr[i] == 0.0 {
            i -= 1;
        }
        if i <= neg_low {
            return Err(DottirError::KarlinFailure(
                "positive score impossible for this matrix / composition".into(),
            ));
        }
    }

    // Sum + non-negativity check.
    let mut sum_raw = 0.0_f64;
    for &v in pr.iter() {
        if v < 0.0 {
            return Err(DottirError::KarlinFailure(
                "negative score probability".into(),
            ));
        }
        sum_raw += v;
    }

    // Normalised copy `p[i] = pr[i] / sum_raw`, accumulated `Sum = Σ i*p[i] + low`.
    let mut p = vec![0.0_f64; range + 1];
    let mut s_expected = low as f64;
    for i in 0..=range {
        p[i] = pr[i] / sum_raw;
        s_expected += (i as f64) * p[i];
    }
    if s_expected >= 0.0 {
        return Err(DottirError::KarlinFailure(format!(
            "non-negative expected score: {s_expected:.3}"
        )));
    }

    // -------------------------------------------------------------- λ ---
    // Find an `up` bracket such that Σ p_i exp(up·i) > 1.
    let mut up = 0.5_f64;
    loop {
        up *= 2.0;
        let mut s = 0.0_f64;
        for i in 0..=range {
            // C: sum += *ptr1++ * exp(up*i), where i runs from `low` to
            // `high`. We reproduce by using `low + i` as the exponent index.
            s += p[i] * ((up * (low as f64 + i as f64)).exp());
        }
        if s >= 1.0 {
            break;
        }
    }

    // 25 bisection iterations seeking λ such that Σ p_i exp(λ·i) = 1.
    let mut lambda = 0.0_f64;
    for _ in 0..LAMBDA_ITERATIONS {
        let mid = (lambda + up) / 2.0;
        let mut s = 0.0_f64;
        for i in 0..=range {
            s += p[i] * ((mid * (low as f64 + i as f64)).exp());
        }
        if s > 1.0 {
            up = mid;
        } else {
            lambda = mid;
        }
    }
    let beta = lambda.exp();

    // -------------------------------------------------------------- H ---
    // Relative entropy H = λ · Σ p_i · i · exp(λ·i).
    let mut av = 0.0_f64;
    for i in 0..=range {
        let score = low as f64 + i as f64;
        av += p[i] * score * (lambda * score).exp();
    }
    let h = lambda * av;

    // -------------------------------------------------------------- K ---
    // Short-circuit cases from the C source.
    if low == -1 || high == 1 {
        let mut k = if high == 1 {
            av
        } else {
            s_expected * s_expected / av
        };
        k *= 1.0 - 1.0 / beta;
        return Ok((lambda, k, h));
    }

    // The K series. C uses a `MAXIT*(range+1)` scratch buffer indexed
    // through `lo`, `hi` running offsets. We use the same convention.
    let mut big_sum = 0.0_f64;
    let mut lo: i64 = 0;
    let mut hi: i64 = 0;
    let mut big_p = vec![0.0_f64; MAXIT * (range + 1)];
    big_p[0] = 1.0;
    let mut s = 1.0_f64;
    let mut oldsum = 1.0_f64;
    let mut oldsum2 = 1.0_f64;
    let mut j = 0_usize;
    while j < MAXIT && s > SUMLIMIT {
        let mut first: i64 = range as i64;
        let mut last: i64 = range as i64;
        // Walk the running window: lo += low, hi += high (low is negative).
        lo += low;
        hi += high;
        let p_top = hi - lo; // index of the last element written this iteration

        // Inner DP: for each output position from p_top down to 0,
        //   big_p[k] = Σ big_p[k - first..=k - first - (last-first)] * p[first..=last]
        // The C code uses two descending pointers; we keep the same structure.
        let mut k = p_top;
        while k >= 0 {
            let mut s_inner = 0.0_f64;
            for i in first..=last {
                let bp_idx = k - i;
                if bp_idx < 0 {
                    break;
                }
                s_inner += big_p[bp_idx as usize] * p[i as usize];
            }
            big_p[k as usize] = s_inner;
            if first != 0 {
                first -= 1;
            }
            if k <= range as i64 {
                last -= 1;
            }
            k -= 1;
        }

        // Aggregate the freshly-written big_p window into the new partial
        // sum. The C code uses two pointer-advancing loops which read
        //   * for score i in [lo, -1):  weight = beta^i (= exp(λ·i))
        //   * for score i in [0, hi]:   weight = 1
        // and reads big_p at offset `i - lo`. We translate to direct
        // indexing.
        let mut s_new = 0.0_f64;
        // Initialise weight = beta^(lo-1) and multiply by beta *before* use,
        // so that at the first iteration weight = beta^lo. The seemingly
        // redundant `powi(beta, lo-1) · beta` (rather than `powi(beta, lo)`)
        // is exactly what the C reference does, and any reformulation
        // diverges from C at the ULP level.
        let mut weight = powi(beta, (lo - 1) as i32);
        let mut i = lo;
        while i != 0 {
            weight *= beta;
            s_new += big_p[(i - lo) as usize] * weight;
            i += 1;
        }
        while i <= hi {
            s_new += big_p[(i - lo) as usize];
            i += 1;
        }

        // C's `oldsum2 = oldsum` happens at end of body, before the
        // for-update (`oldsum = sum, Sum += sum /= ++j`). So at the end of
        // the body, `oldsum` still holds the *previous* iteration's
        // pre-divide sum, and `s_new` is the *current* iteration's
        // pre-divide sum. We snapshot oldsum into oldsum2, then move
        // s_new into oldsum, then divide and accumulate.
        oldsum2 = oldsum;
        oldsum = s_new;
        j += 1;
        s = s_new / j as f64;
        big_sum += s;
    }

    // Geometric-progression correction. If the ratio is very close to 1 the
    // series hasn't converged usefully; the C code fudges K to 0.1.
    let ratio = oldsum / oldsum2;
    let k = if ratio >= 1.0 - SUMLIMIT * 0.001 {
        0.1
    } else {
        // Extend the series with geometric correction.
        let mut s_corr = s;
        let mut sum = big_sum;
        let mut oldsum_corr = oldsum;
        let mut j_corr = j;
        while s_corr > SUMLIMIT * 0.01 {
            oldsum_corr *= ratio;
            j_corr += 1;
            s_corr = oldsum_corr / j_corr as f64;
            sum += s_corr;
        }

        // Determine the gcd `g` of all score differences that have non-zero
        // probability (the lattice spacing). Then
        //   K = g · exp(-2·Sum) / (av · etop(λ · g))
        // where `etop(x) = -expm1(-x)`.
        let mut first_nonzero_idx: i64 = -1;
        for i in 0..=range {
            if p[i] != 0.0 {
                first_nonzero_idx = low + i as i64;
                break;
            }
        }
        let mut g = -first_nonzero_idx;
        if g <= 0 {
            return Err(DottirError::KarlinFailure(
                "could not determine score lattice spacing".into(),
            ));
        }
        // Walk the rest looking for the gcd.
        let mut cur = first_nonzero_idx;
        let high_i = high;
        while cur < high_i && g > 1 {
            cur += 1;
            if p[(cur - low) as usize] != 0.0 {
                g = gcd(g, cur);
            }
        }
        let g_f = g as f64;
        // C: `K = (g · exp(-2·Sum)) / (av · etop(λ·g))` where
        // `etop(x) = -fct_expm1(-x)`. The exact order of operations
        // (and our [`fct_expm1`] using the C polynomial rather than
        // libm's correctly-rounded `exp_m1`) is what gives ULP-identical
        // results against the C reference.
        let etop = -fct_expm1(-lambda * g_f);
        g_f * (-2.0 * sum).exp() / (av * etop)
    };

    Ok((lambda, k, h))
}

/// `exp(x) - 1` evaluated the same way as `dotterKarlin.c::fct_expm1`:
/// 12-term Horner-folded Taylor series for `|x| <= 0.33`, otherwise
/// `exp(x) - 1`. We use this rather than [`f64::exp_m1`] because the
/// libm implementation is correctly rounded but a different bit pattern
/// from the C reference — and Karlin K is sensitive to that last ULP.
//
// `rustfmt::skip`: the 13-deep Horner nest below triggers exponential
// time in rustfmt 1.9.0 (the function never returns). The expression
// is a verbatim port of the C reference and we preserve its shape
// anyway, so skipping the formatter here is both a workaround and the
// right call on its own merits.
#[rustfmt::skip]
fn fct_expm1(x: f64) -> f64 {
    let absx = if x < 0.0 { -x } else { x };
    if absx > 0.33 {
        return x.exp() - 1.0;
    }
    if absx < 1.0e-16 {
        return x;
    }
    x * (1.0
        + x * (0.5
            + x * (1.0 / 6.0
                + x * (1.0 / 24.0
                    + x * (1.0 / 120.0
                        + x * (1.0 / 720.0
                            + x * (1.0 / 5040.0
                                + x * (1.0 / 40320.0
                                    + x * (1.0 / 362880.0
                                        + x * (1.0 / 3628800.0
                                            + x * (1.0 / 39916800.0
                                                + x * (1.0 / 479001600.0
                                                    + x / 6227020800.0))))))))))))
}

#[inline]
fn gcd(a: i64, b: i64) -> i64 {
    let mut a = a.unsigned_abs();
    let mut b = b.unsigned_abs();
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a as i64
}

/// `x^n` via the C dotter's `fct_powi` (integer exponentiation by squaring).
/// Reproduced exactly so the K-series uses the same intermediate values.
fn powi(mut x: f64, n: i32) -> f64 {
    let mut y = 1.0_f64;
    let mut i = n.unsigned_abs() as i64;
    while i > 0 {
        if i & 1 == 1 {
            y *= x;
        }
        x *= x;
        i /= 2;
    }
    if n >= 0 {
        y
    } else {
        1.0 / y
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_identity_blastn_window() {
        // A roughly balanced random DNA pair; result should be within the
        // [min,max] clamp and reproducible across runs.
        let q = b"ACGTACGTACGTACGTAACGTNNNNCCCCAAAAGGTGTGTAGCTAGCT".repeat(40);
        let s = b"GTGTACGAGCATCGTCTACTGAGCTACGTATCGATCGTAGCTACGATG".repeat(40);
        let m = ScoreMatrix::dna_identity();
        let res = karlin_window_size(&m, &q, &s, BlastMode::Blastn).unwrap();
        assert!(
            res.lambda > 0.0,
            "lambda must be positive, got {}",
            res.lambda
        );
        assert!(res.k > 0.0, "k must be positive, got {}", res.k);
        assert!(res.h > 0.0, "h must be positive, got {}", res.h);
        assert!(
            (3..=50).contains(&res.window_size),
            "window {} out of [3,50]",
            res.window_size
        );
    }

    #[test]
    fn protein_blosum62_window() {
        // Use uniform protein composition; matches the Altschul default
        // assumption fairly well and exercises the K-series.
        let q = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKAL\
                  RSALEFNAHVDEMVRLRREVGNQLEELQNRLREYIQRDHRGHEALQQYRVKQVHLDQEEIA";
        let s = b"MAATKRIIRQRYTIKHYVTRLREHIDHEEQVRKDLDEHKHRADRMLEELAGAILAAEHRLRD\
                  AREAFEQLLDKLEEHLRYAEELQEKFAKLERELAEHRLEEIEGRLAQAEEEFVEQHRRLENEL";
        let m = ScoreMatrix::blosum62();
        let res = karlin_window_size(&m, q, s, BlastMode::Blastp).unwrap();
        assert!(
            res.lambda > 0.1 && res.lambda < 1.0,
            "lambda={}",
            res.lambda
        );
        assert!(res.k > 0.0 && res.k < 1.0, "k={}", res.k);
        assert!(
            (3..=50).contains(&res.window_size),
            "window {} out of [3,50]",
            res.window_size
        );
    }

    #[test]
    fn empty_sequence_errors() {
        let m = ScoreMatrix::dna_identity();
        let err = karlin_window_size(&m, b"", b"ACGT", BlastMode::Blastn).unwrap_err();
        assert!(matches!(err, DottirError::EmptySequence));
    }

    #[test]
    fn all_unscorable_errors() {
        // NNNN has zero scorable residues under the DNA Karlin alphabet
        // (size 4, 'N' bucket excluded).
        let m = ScoreMatrix::dna_identity();
        let err = karlin_window_size(&m, b"NNNN", b"ACGT", BlastMode::Blastn).unwrap_err();
        assert!(matches!(err, DottirError::NoScorableResidues));
    }

    #[test]
    fn matrix_mode_mismatch_errors() {
        // Pass a DNA matrix in blastp mode.
        let m = ScoreMatrix::dna_identity();
        let err =
            karlin_window_size(&m, b"ACDEFGHIK", b"ACDEFGHIK", BlastMode::Blastp).unwrap_err();
        assert!(matches!(err, DottirError::InvalidMatrix(_)));
    }

    #[test]
    fn gcd_helper() {
        assert_eq!(gcd(12, 8), 4);
        assert_eq!(gcd(-12, 8), 4);
        assert_eq!(gcd(7, 0), 7);
        assert_eq!(gcd(0, 7), 7);
        assert_eq!(gcd(1, 100), 1);
    }

    #[test]
    fn powi_helper_matches_pow() {
        assert!((powi(2.0, 10) - 1024.0).abs() < 1e-12);
        assert!((powi(0.5, -3) - 8.0).abs() < 1e-12);
        assert!((powi(1.5, 0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn determinism_same_inputs_same_outputs() {
        // Spec §4.1.11: byte-identical results across runs. We can't span
        // processes here, but we *can* assert bit-equality on repeated
        // invocations within a process.
        let q = b"ACGTACGTACGTACGTACGT".repeat(40);
        let s = b"GTGTACGAGCATCGTCTACT".repeat(40);
        let m = ScoreMatrix::dna_identity();
        let r1 = karlin_window_size(&m, &q, &s, BlastMode::Blastn).unwrap();
        let r2 = karlin_window_size(&m, &q, &s, BlastMode::Blastn).unwrap();
        assert_eq!(r1.lambda.to_bits(), r2.lambda.to_bits());
        assert_eq!(r1.k.to_bits(), r2.k.to_bits());
        assert_eq!(r1.h.to_bits(), r2.h.to_bits());
        assert_eq!(r1.window_size, r2.window_size);
    }
}
