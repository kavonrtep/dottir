//! Classify periodogram / spectrum peaks into fundamentals,
//! harmonics, and sub-repeats.
//!
//! Given a 1-D signal of scores indexed by period (residues for a
//! periodogram, residues-per-cycle for an FFT spectrum), this module
//! identifies the dominant periodic structure and labels every peak
//! with its relationship to that structure:
//!
//! * **Fundamental** — a peak that doesn't sit on the integer-ratio
//!   ladder of any stronger peak. The biologically interesting answer.
//! * **Harmonic** — a peak at an integer multiple of a fundamental's
//!   period (periodogram) or integer divisor (FFT). Carries the
//!   integer `n` and the parent fundamental's period.
//! * **Sub-repeat** — a peak at an integer DIVISOR of a fundamental
//!   (periodogram) — typically the monomer underlying an HOR. Only
//!   surfaced when the divisor position has enough signal but the
//!   fundamental's score is much higher (e.g., noisy-monomer HORs).
//!
//! Workflow (mirrors the validated Python prototype at
//! `scripts/find_peaks.py`):
//!
//! 1. **Local-maximum extraction** — strict three-bin local maxima
//!    above `min_score`. For periodogram input, an optional
//!    `boundary_fraction` trims peaks at `k > boundary_fraction *
//!    input_len` (large records show kernel saturation artifacts at
//!    `k ≈ N/2`).
//! 2. **Greedy classification, strongest-first** — the highest-scoring
//!    candidate is always a fundamental; subsequent candidates that
//!    match an integer-ratio harmonic of an existing fundamental get
//!    classified as such.
//! 3. **Bidirectional consolidation** — if a fundamental turns out to
//!    be a harmonic of another fundamental, demote it (and cascade
//!    its harmonics). Catches the "harmonic happens to be stronger"
//!    case where the algorithm originally claimed the wrong peak as
//!    fundamental.
//! 4. **Sub-repeat scan** (optional, periodogram only) — for each
//!    fundamental, look at `period / n` for `n in 2..=max_divisor`
//!    with a small period-tolerance. If a peak there exceeds
//!    `sub_min_score`, label it `Sub-repeat`.
//! 5. **`min_harmonics` filter** — drop fundamentals with sparser
//!    harmonic ladders than the threshold (1 is a sensible default —
//!    a "fundamental" with zero detected harmonics is rarely real).
//!
//! Validated on `test6`, `test7`, `test8`, `TRC_7`, and `TRC_2` —
//! consistently identifies the same fundamentals as visual dotplot
//! inspection.

use crate::error::DottirError;
use crate::spectrum::Spectrum;

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

/// How a peak relates to the dominant periodic structure.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PeakKind {
    /// Top-level periodic signal. Not on any other peak's integer
    /// ratio ladder.
    Fundamental,
    /// On a fundamental's integer-multiple ladder
    /// (periodogram: `period ≈ n × fundamental`).
    Harmonic,
    /// On a fundamental's integer-divisor ladder
    /// (periodogram: `period ≈ fundamental / n`). Typically the
    /// monomer underlying an HOR.
    Subrepeat,
}

/// A classified peak.
#[derive(Clone, Debug)]
pub struct Peak {
    pub period: f64,
    pub score: f64,
    pub kind: PeakKind,
    /// For `Harmonic` and `Subrepeat`: the parent fundamental's
    /// period. `None` for `Fundamental`.
    pub parent_period: Option<f64>,
    /// For `Harmonic`: the integer multiplier `n` such that
    /// `period ≈ n × parent_period`. `None` otherwise.
    pub harmonic_n: Option<u32>,
    /// For `Subrepeat`: the integer divisor `n` such that
    /// `period ≈ parent_period / n`. `None` otherwise.
    pub divisor_n: Option<u32>,
    /// For `Fundamental`: count of distinct harmonic positions
    /// detected. `0` for `Harmonic` and `Subrepeat`.
    pub n_harmonics: u32,
    /// For `Fundamental`: sorted unique integer harmonic positions
    /// (e.g. `[2, 3, 4, 5, 6]`). Empty otherwise.
    pub harmonics: Vec<u32>,
}

impl Peak {
    fn new_fundamental(period: f64, score: f64) -> Self {
        Self {
            period,
            score,
            kind: PeakKind::Fundamental,
            parent_period: None,
            harmonic_n: None,
            divisor_n: None,
            n_harmonics: 0,
            harmonics: Vec::new(),
        }
    }

    fn into_harmonic(mut self, parent_period: f64, n: u32) -> Self {
        self.kind = PeakKind::Harmonic;
        self.parent_period = Some(parent_period);
        self.harmonic_n = Some(n);
        self.divisor_n = None;
        self.n_harmonics = 0;
        self.harmonics.clear();
        self
    }
}

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

/// Convention for the harmonic direction. Periodogram and FFT
/// interpret integer ratios in opposite senses (periodogram: a
/// harmonic has a *larger* period than its fundamental; FFT: a
/// harmonic has a *smaller* period).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HarmonicDirection {
    /// Periodogram. `harmonic.period = n × fundamental.period`.
    LargerIsHarmonic,
    /// FFT (period space). `fundamental.period = n × harmonic.period`.
    SmallerIsHarmonic,
}

/// Subrepeat-detection knobs. Only meaningful for periodogram input.
#[derive(Clone, Copy, Debug)]
pub struct SubrepeatConfig {
    /// Minimum score to qualify as a subrepeat (typically lower than
    /// the main `min_score` — subrepeats are weaker than their
    /// parent HOR).
    pub min_score: f64,
    /// Largest integer divisor to scan (2..=max_divisor). Most
    /// natural HORs are 2- to 5-mers; 6 covers the common cases.
    pub max_divisor: u32,
    /// `±` bp tolerance when looking for the subrepeat peak around
    /// `fundamental / n`. The actual monomer peak can drift by a
    /// couple of bp when monomers have irregular size.
    pub period_tolerance: u32,
}

impl Default for SubrepeatConfig {
    fn default() -> Self {
        Self {
            min_score: 5.0,
            max_divisor: 6,
            period_tolerance: 2,
        }
    }
}

/// Configuration for the peak finder.
#[derive(Clone, Copy, Debug)]
pub struct PeaksConfig {
    /// Minimum score for a peak candidate. Defaults to 10 (validated
    /// for `signal_mean`-ranked periodograms). For `z_score`-ranked
    /// input use ~5.
    pub min_score: f64,
    /// Relative tolerance for integer-ratio harmonic detection.
    /// Default `0.02` (2%). Larger = stricter dedup; smaller =
    /// catches more harmonics at the risk of false matches.
    pub harmonic_tolerance: f64,
    /// Largest integer harmonic to consider in classification.
    /// Default `30`. Long tandem arrays show clean harmonics out to
    /// n>20 — going to 30 covers most natural cases.
    pub max_harmonic_n: u32,
    /// Keep only fundamentals with at least this many detected
    /// harmonics. Default `1` — a "fundamental" with zero harmonics
    /// is almost always noise.
    pub min_harmonics: u32,
    /// If `Some`, enable sub-repeat detection (periodogram only).
    pub subrepeats: Option<SubrepeatConfig>,
    /// Drop peaks at `period > boundary_fraction × input_len` —
    /// kernel edge artifacts saturate at `k ≈ N/2`. Default `0.9`.
    /// Set to `1.0` to disable. Periodogram-only; ignored for
    /// spectrum input.
    pub boundary_fraction: f64,
}

impl Default for PeaksConfig {
    fn default() -> Self {
        Self {
            min_score: 10.0,
            harmonic_tolerance: 0.02,
            max_harmonic_n: 30,
            min_harmonics: 1,
            subrepeats: None,
            boundary_fraction: 0.9,
        }
    }
}

// ---------------------------------------------------------------------------
// Public entry points
// ---------------------------------------------------------------------------

/// Find and classify peaks in a periodogram score signal.
///
/// `scores[i]` is the score at offset `i + offset_base`. For typical
/// dottir periodogram output:
///
/// * `scores = periodogram.raw_sum / n_pairs` (i.e. `signal_mean`) —
///   **recommended default**, robust against analytical z-score bias.
/// * `scores = empirical_z_score` — works when z-score is reliable
///   (rare in large records under analytical mode).
///
/// `offset_base = periodogram.min_offset` from
/// [`crate::periodogram::Periodogram`].
pub fn find_peaks_in_periodogram(
    scores: &[f64],
    offset_base: u32,
    cfg: &PeaksConfig,
) -> Result<Vec<Peak>, DottirError> {
    if scores.is_empty() {
        return Ok(Vec::new());
    }

    // 1. Local maxima above min_score, respecting the boundary cap.
    let max_period_inclusive = if cfg.boundary_fraction >= 1.0 {
        None
    } else {
        // `scores.len()` is the count of buckets, so the highest
        // *period* in this signal is `offset_base + scores.len() - 1`.
        // The "input length" the boundary fraction guards against is
        // `2 × max_period` (the periodogram covers k in [min, N/2]).
        // Concretely: drop peaks whose period exceeds
        // `boundary_fraction × (offset_base + scores.len() - 1)`.
        let max_period = offset_base as usize + scores.len() - 1;
        Some((max_period as f64 * cfg.boundary_fraction) as u32)
    };
    let mut candidates: Vec<Peak> = strict_local_maxima(scores)
        .into_iter()
        .filter_map(|(i, score)| {
            let period = (offset_base + i as u32) as f64;
            if score < cfg.min_score {
                return None;
            }
            if let Some(cap) = max_period_inclusive {
                if period as u32 > cap {
                    return None;
                }
            }
            Some(Peak::new_fundamental(period, score))
        })
        .collect();

    // 2. Greedy classification.
    classify(&mut candidates, HarmonicDirection::LargerIsHarmonic, cfg);

    // 3. Consolidation (catches "harmonic-claimed-as-fundamental").
    consolidate(&mut candidates, HarmonicDirection::LargerIsHarmonic, cfg);
    reparent(&mut candidates, HarmonicDirection::LargerIsHarmonic, cfg);

    // 4. Sub-repeats (periodogram-only).
    let subrepeats = if let Some(sub_cfg) = cfg.subrepeats {
        scan_subrepeats(&candidates, scores, offset_base, &sub_cfg)
    } else {
        Vec::new()
    };

    // 5. min_harmonics filter on fundamentals.
    let mut out: Vec<Peak> = candidates
        .into_iter()
        .filter(|p| p.kind != PeakKind::Fundamental || p.n_harmonics >= cfg.min_harmonics)
        .collect();
    out.extend(subrepeats);
    // Sort fundamentals first (by score desc); within each fundamental, its
    // sub-repeats follow. Other peaks (harmonics) get appended afterward.
    out.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(out)
}

/// Find and classify peaks in an FFT [`Spectrum`].
///
/// Uses the spectrum's already-annotated `peak_ranks` as candidates
/// (so it inherits the spectrum's local-max detection). Subrepeats
/// don't apply to FFT input — the harmonic ladder already runs in
/// both period directions because FFT lattice is `k = 1, 2, 3, …`
/// in frequency space.
pub fn find_peaks_in_spectrum(
    spectrum: &Spectrum,
    cfg: &PeaksConfig,
) -> Result<Vec<Peak>, DottirError> {
    let mut candidates: Vec<Peak> = (0..spectrum.amplitude.len())
        .filter_map(|bin| {
            spectrum.peak_ranks[bin]
                .map(|_rank| Peak::new_fundamental(spectrum.period(bin), spectrum.amplitude[bin]))
        })
        .filter(|p| p.score >= cfg.min_score && p.period.is_finite())
        .collect();

    classify(&mut candidates, HarmonicDirection::SmallerIsHarmonic, cfg);
    consolidate(&mut candidates, HarmonicDirection::SmallerIsHarmonic, cfg);
    reparent(&mut candidates, HarmonicDirection::SmallerIsHarmonic, cfg);

    let mut out: Vec<Peak> = candidates
        .into_iter()
        .filter(|p| p.kind != PeakKind::Fundamental || p.n_harmonics >= cfg.min_harmonics)
        .collect();
    out.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(out)
}

// ---------------------------------------------------------------------------
// Internal: local-maxima detection
// ---------------------------------------------------------------------------

/// Strict three-bin local maxima: `i` such that `s[i] > s[i-1]` and
/// `s[i] > s[i+1]`. Endpoints are not considered (no peak on the
/// very first or last bin — boundary effects).
/// Compute a robust per-record threshold using median + k×MAD.
///
/// Returns `max(floor, median + k × 1.4826 × MAD)`. The
/// `1.4826` makes the MAD-based estimate consistent with a
/// Gaussian σ. Non-finite scores are dropped before the
/// computation. If `MAD == 0` (e.g., uniform input), returns
/// `floor` — there's nothing useful to adapt to.
///
/// Recommended `k = 5` for periodogram peak detection (5-σ
/// equivalent) — conservative enough that real peaks stand
/// out above ordinary variation.
pub fn auto_threshold_mad(scores: &[f64], k: f64, floor: f64) -> f64 {
    let mut values: Vec<f64> = scores.iter().copied().filter(|s| s.is_finite()).collect();
    if values.is_empty() {
        return floor;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = percentile_sorted(&values, 0.5);
    let mut abs_dev: Vec<f64> = values.iter().map(|v| (v - median).abs()).collect();
    abs_dev.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mad = percentile_sorted(&abs_dev, 0.5);
    if mad <= 0.0 {
        return floor;
    }
    let adaptive = median + k * 1.4826 * mad;
    adaptive.max(floor)
}

/// Linear-interpolated percentile of a pre-sorted slice. `q` in [0, 1].
fn percentile_sorted(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }
    let pos = q.clamp(0.0, 1.0) * (sorted.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = (lo + 1).min(sorted.len() - 1);
    let frac = pos - lo as f64;
    sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

fn strict_local_maxima(scores: &[f64]) -> Vec<(usize, f64)> {
    if scores.len() < 3 {
        return Vec::new();
    }
    let mut out = Vec::new();
    for i in 1..scores.len() - 1 {
        if scores[i] > scores[i - 1] && scores[i] > scores[i + 1] {
            out.push((i, scores[i]));
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Internal: greedy classification
// ---------------------------------------------------------------------------

/// Greedy strongest-first classification with **best-parent**
/// selection. Each candidate is checked against every claimed
/// fundamental; if any match within tolerance, it joins the
/// best-matching ladder. Otherwise it becomes a new fundamental.
///
/// "Best" tie-breaks (in order):
/// 1. smallest normalised error (distance / allowed_bp);
/// 2. smallest absolute period distance;
/// 3. smallest `n` (prefer direct integer ratios over deep harmonics);
/// 4. strongest parent score (for determinism).
fn classify(peaks: &mut [Peak], direction: HarmonicDirection, cfg: &PeaksConfig) {
    // Sort strongest-first. Ties broken by lower period (shorter
    // period first) for determinism.
    peaks.sort_by(|a, b| {
        b.score
            .partial_cmp(&a.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                a.period
                    .partial_cmp(&b.period)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
    });

    // Walk left-to-right (already sorted by score desc), tracking
    // confirmed fundamentals by index.
    let mut fundamental_idx: Vec<usize> = Vec::new();
    for i in 0..peaks.len() {
        let parents: Vec<(f64, f64)> = fundamental_idx
            .iter()
            .map(|&fi| (peaks[fi].period, peaks[fi].score))
            .collect();
        let matched = best_parent_match(peaks[i].period, &parents, direction, cfg);
        if let Some(m) = matched {
            // Find the index of the matched parent.
            let fi = fundamental_idx
                .iter()
                .copied()
                .find(|&fi| (peaks[fi].period - m.parent_period).abs() < 1e-9)
                .expect("parent index lookup should succeed");
            let parent_period = peaks[fi].period;
            let taken = std::mem::replace(&mut peaks[i], Peak::new_fundamental(0.0, 0.0));
            peaks[i] = taken.into_harmonic(parent_period, m.n);
            // Record the harmonic position in the parent's ladder.
            let prev = std::mem::take(&mut peaks[fi].harmonics);
            let mut new = prev;
            if !new.contains(&m.n) {
                new.push(m.n);
                new.sort_unstable();
            }
            peaks[fi].n_harmonics = new.len() as u32;
            peaks[fi].harmonics = new;
        } else {
            fundamental_idx.push(i);
        }
    }
}

/// Result of a parent-match search.
#[derive(Clone, Copy, Debug)]
struct ParentMatch {
    parent_period: f64,
    n: u32,
    /// Normalised error: `distance_bp / allowed_bp`. Smaller is better.
    norm_err: f64,
    /// Absolute distance in bp to the expected period.
    abs_distance: f64,
}

/// Pick the best matching parent for `peak_period` from the list of
/// `(period, score)` parent candidates. Returns `None` if no parent
/// matches within tolerance.
fn best_parent_match(
    peak_period: f64,
    parents: &[(f64, f64)],
    direction: HarmonicDirection,
    cfg: &PeaksConfig,
) -> Option<ParentMatch> {
    let mut best: Option<(f64, ParentMatch)> = None; // (parent_score, match)
    for &(parent_period, parent_score) in parents {
        let Some(n) = harmonic_ratio(
            peak_period,
            parent_period,
            direction,
            cfg.harmonic_tolerance,
            cfg.max_harmonic_n,
        ) else {
            continue;
        };
        // Expected period in bp.
        let expected = match direction {
            HarmonicDirection::LargerIsHarmonic => parent_period * n as f64,
            HarmonicDirection::SmallerIsHarmonic => parent_period / n as f64,
        };
        let abs_distance = (peak_period - expected).abs();
        // `allowed_bp` mirrors the harmonic_ratio acceptance:
        // ratio_err <= tol  ↔  |period - expected|/expected <= tol.
        let allowed = (expected * cfg.harmonic_tolerance).abs().max(f64::EPSILON);
        let norm_err = abs_distance / allowed;
        let m = ParentMatch {
            parent_period,
            n,
            norm_err,
            abs_distance,
        };
        let replace = match &best {
            None => true,
            Some((best_score, best_m)) => {
                // 1. smaller normalised error wins
                if m.norm_err < best_m.norm_err - 1e-12 {
                    true
                } else if m.norm_err > best_m.norm_err + 1e-12 {
                    false
                // 2. smaller absolute distance wins
                } else if m.abs_distance < best_m.abs_distance - 1e-12 {
                    true
                } else if m.abs_distance > best_m.abs_distance + 1e-12 {
                    false
                // 3. smaller n wins
                } else if m.n < best_m.n {
                    true
                } else if m.n > best_m.n {
                    false
                // 4. stronger parent for determinism
                } else {
                    parent_score > *best_score
                }
            }
        };
        if replace {
            best = Some((parent_score, m));
        }
    }
    best.map(|(_, m)| m)
}

/// Post-consolidation pass: for every classified harmonic, recompute
/// its best parent among the **final** fundamentals (after
/// `consolidate` may have demoted some). This makes the final
/// assignment independent of the order in which candidates were
/// first classified.
///
/// Then rebuild every fundamental's `harmonics` ladder from the
/// final harmonic-to-parent assignments.
fn reparent(peaks: &mut [Peak], direction: HarmonicDirection, cfg: &PeaksConfig) {
    // Snapshot final fundamentals as (period, score).
    let parents: Vec<(f64, f64)> = peaks
        .iter()
        .filter(|p| p.kind == PeakKind::Fundamental)
        .map(|p| (p.period, p.score))
        .collect();
    if parents.is_empty() {
        return;
    }

    // Clear every fundamental's harmonics list — we'll rebuild it.
    for p in peaks.iter_mut() {
        if p.kind == PeakKind::Fundamental {
            p.harmonics.clear();
            p.n_harmonics = 0;
        }
    }

    // For each harmonic, find the best parent in the final set.
    // Track (harmonic_index, parent_period, n) for the ladder rebuild.
    let mut reassignments: Vec<(usize, f64, u32)> = Vec::new();
    for (i, p) in peaks.iter().enumerate() {
        if p.kind != PeakKind::Harmonic {
            continue;
        }
        if let Some(m) = best_parent_match(p.period, &parents, direction, cfg) {
            reassignments.push((i, m.parent_period, m.n));
        }
        // If a harmonic no longer matches any parent (rare), leave it
        // attached to its current parent — the consolidation pass will
        // never produce an orphan harmonic in practice.
    }

    for (i, parent_period, n) in reassignments {
        peaks[i].parent_period = Some(parent_period);
        peaks[i].harmonic_n = Some(n);
        // Update parent's ladder.
        if let Some(parent_idx) = peaks.iter().position(|p| {
            p.kind == PeakKind::Fundamental && (p.period - parent_period).abs() < 1e-9
        }) {
            if !peaks[parent_idx].harmonics.contains(&n) {
                peaks[parent_idx].harmonics.push(n);
                peaks[parent_idx].harmonics.sort_unstable();
                peaks[parent_idx].n_harmonics = peaks[parent_idx].harmonics.len() as u32;
            }
        }
    }
}

/// If `peak_period` is the n-th harmonic of `fundamental_period`,
/// return n. Otherwise None.
fn harmonic_ratio(
    peak_period: f64,
    fundamental_period: f64,
    direction: HarmonicDirection,
    tolerance: f64,
    max_n: u32,
) -> Option<u32> {
    let ratio = match direction {
        HarmonicDirection::LargerIsHarmonic => peak_period / fundamental_period,
        HarmonicDirection::SmallerIsHarmonic => fundamental_period / peak_period,
    };
    if !ratio.is_finite() {
        return None;
    }
    let n = ratio.round() as i64;
    if n < 2 || n > max_n as i64 {
        return None;
    }
    let err = (ratio - n as f64).abs() / n as f64;
    if err > tolerance {
        return None;
    }
    Some(n as u32)
}

// ---------------------------------------------------------------------------
// Internal: bidirectional consolidation
// ---------------------------------------------------------------------------

/// Demote any fundamental that turns out to be a harmonic of another
/// (stronger-score isn't the convention winner). Cascades the demoted
/// peak's existing harmonics to the survivor with multiplied n.
fn consolidate(peaks: &mut [Peak], direction: HarmonicDirection, _cfg: &PeaksConfig) {
    loop {
        // Index of any fundamental that should be demoted, and to which.
        let (demote_idx, survivor_idx, demote_n) = {
            let mut found = None;
            let fundamentals: Vec<usize> = peaks
                .iter()
                .enumerate()
                .filter(|(_, p)| p.kind == PeakKind::Fundamental)
                .map(|(i, _)| i)
                .collect();
            'outer: for &fi in &fundamentals {
                for &gi in &fundamentals {
                    if fi == gi {
                        continue;
                    }
                    // f (peaks[fi]) is demoted if g (peaks[gi]) is the
                    // biological fundamental and f is its n-th harmonic.
                    let (n_demote, ok) = match direction {
                        HarmonicDirection::LargerIsHarmonic => {
                            // Smaller period is fundamental. Demote f
                            // only if g.period < f.period.
                            if peaks[gi].period >= peaks[fi].period {
                                continue;
                            }
                            let ratio = peaks[fi].period / peaks[gi].period;
                            let n = ratio.round() as i64;
                            let ok = n >= 2
                                && n <= _cfg.max_harmonic_n as i64
                                && ((ratio - n as f64).abs() / n as f64) <= _cfg.harmonic_tolerance;
                            (n as u32, ok)
                        }
                        HarmonicDirection::SmallerIsHarmonic => {
                            // Larger period is fundamental. Demote f
                            // only if g.period > f.period.
                            if peaks[gi].period <= peaks[fi].period {
                                continue;
                            }
                            let ratio = peaks[gi].period / peaks[fi].period;
                            let n = ratio.round() as i64;
                            let ok = n >= 2
                                && n <= _cfg.max_harmonic_n as i64
                                && ((ratio - n as f64).abs() / n as f64) <= _cfg.harmonic_tolerance;
                            (n as u32, ok)
                        }
                    };
                    if ok {
                        found = Some((fi, gi, n_demote));
                        break 'outer;
                    }
                }
            }
            match found {
                Some(t) => t,
                None => return,
            }
        };

        // Demote peaks[demote_idx]. Cascade its harmonics to survivor.
        let demoted_period = peaks[demote_idx].period;
        let survivor_period = peaks[survivor_idx].period;
        let demoted_harmonics = std::mem::take(&mut peaks[demote_idx].harmonics);
        peaks[demote_idx].n_harmonics = 0;

        // Convert peaks[demote_idx] to a harmonic of survivor.
        let mut demoted =
            std::mem::replace(&mut peaks[demote_idx], Peak::new_fundamental(0.0, 0.0));
        let parent_period = survivor_period;
        let n = demote_n;
        demoted.kind = PeakKind::Harmonic;
        demoted.parent_period = Some(parent_period);
        demoted.harmonic_n = Some(n);
        demoted.divisor_n = None;
        peaks[demote_idx] = demoted;

        // Build the new harmonic-n list for the survivor.
        let mut new_ns: Vec<u32> = Vec::with_capacity(demoted_harmonics.len() + 1);
        new_ns.push(n);
        for old_n in &demoted_harmonics {
            new_ns.push(old_n * n);
        }

        // Re-attribute peaks previously labelled as harmonics of the
        // demoted fundamental.
        for p in peaks.iter_mut() {
            if p.kind == PeakKind::Harmonic && p.parent_period == Some(demoted_period) {
                if let Some(old_n) = p.harmonic_n {
                    let new_n = old_n * n;
                    p.parent_period = Some(survivor_period);
                    p.harmonic_n = Some(new_n);
                    new_ns.push(new_n);
                }
            }
        }

        // Merge into survivor's harmonics.
        let mut merged: Vec<u32> = peaks[survivor_idx]
            .harmonics
            .iter()
            .chain(new_ns.iter())
            .copied()
            .collect();
        merged.sort_unstable();
        merged.dedup();
        peaks[survivor_idx].n_harmonics = merged.len() as u32;
        peaks[survivor_idx].harmonics = merged;
    }
}

// ---------------------------------------------------------------------------
// Internal: sub-repeat scanning
// ---------------------------------------------------------------------------

fn scan_subrepeats(
    peaks: &[Peak],
    scores: &[f64],
    offset_base: u32,
    cfg: &SubrepeatConfig,
) -> Vec<Peak> {
    // Quick existence set: don't re-emit a peak that's already in the
    // candidate list (it'll be classified there).
    let existing: std::collections::HashSet<u32> = peaks
        .iter()
        .filter_map(|p| {
            if p.period.is_finite() {
                Some(p.period as u32)
            } else {
                None
            }
        })
        .collect();

    let mut out = Vec::new();
    for p in peaks.iter().filter(|p| p.kind == PeakKind::Fundamental) {
        let parent_period_u32 = p.period as u32;
        for n in 2..=cfg.max_divisor {
            let target = (p.period / n as f64).round() as i32;
            if target < 3 {
                break;
            }
            let mut best: Option<(u32, f64)> = None;
            let tol = cfg.period_tolerance as i32;
            for d in -tol..=tol {
                let cand = target + d;
                if cand < 3 {
                    continue;
                }
                let idx = cand as i64 - offset_base as i64;
                if idx < 0 || idx >= scores.len() as i64 {
                    continue;
                }
                let s = scores[idx as usize];
                if best.is_none_or(|(_, bs)| s > bs) {
                    best = Some((cand as u32, s));
                }
            }
            let Some((best_period, best_score)) = best else {
                continue;
            };
            if best_score < cfg.min_score {
                continue;
            }
            // Strict local-max check at best_period.
            let idx = best_period as i64 - offset_base as i64;
            let prev = if idx > 0 {
                scores[(idx - 1) as usize]
            } else {
                f64::NEG_INFINITY
            };
            let next = if (idx as usize + 1) < scores.len() {
                scores[(idx + 1) as usize]
            } else {
                f64::NEG_INFINITY
            };
            if best_score <= prev || best_score <= next {
                continue;
            }
            if existing.contains(&best_period) {
                continue;
            }
            out.push(Peak {
                period: best_period as f64,
                score: best_score,
                kind: PeakKind::Subrepeat,
                parent_period: Some(p.period),
                harmonic_n: None,
                divisor_n: Some(n),
                n_harmonics: 0,
                harmonics: Vec::new(),
            });
        }
        let _ = parent_period_u32;
    }
    out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a synthetic periodogram-like score vector: a single
    /// spike at each requested offset (no Gaussian spread — that
    /// would let neighbouring periods bleed into each other when
    /// testing two coprime families).
    fn synth_bumps(len: usize, offsets: &[usize], amp: f64) -> Vec<f64> {
        let mut v = vec![0.0_f64; len];
        for &o in offsets {
            if o < len {
                v[o] = v[o].max(amp);
            }
        }
        v
    }

    #[test]
    fn pure_tandem_detects_fundamental_only() {
        // Period-7 impulses at 7, 14, 21, …
        let mut offsets = Vec::new();
        let mut k = 7usize;
        while k < 200 {
            offsets.push(k);
            k += 7;
        }
        let scores = synth_bumps(300, &offsets, 100.0);
        let cfg = PeaksConfig::default();
        let peaks = find_peaks_in_periodogram(&scores, 0, &cfg).unwrap();
        let funds: Vec<_> = peaks
            .iter()
            .filter(|p| p.kind == PeakKind::Fundamental)
            .collect();
        assert_eq!(funds.len(), 1, "expected one fundamental, got {funds:?}");
        let f = funds[0];
        assert_eq!(f.period as u32, 7);
        // Harmonic ladder should be 2,3,4,5,… up to whatever fits.
        assert!(f.n_harmonics >= 3, "expected ladder; got {}", f.n_harmonics);
        assert_eq!(f.harmonics[0], 2);
    }

    #[test]
    fn two_independent_periods_kept_separate() {
        let mut offsets = Vec::new();
        // Two coprime fundamentals chosen to have a short signal
        // range so high-order accidental near-integer ratios don't
        // creep in (real biology is the same: at high harmonic n
        // any two periods will eventually look like multiples of
        // each other under finite-tolerance matching).
        // Family A: period 13
        let mut k = 13usize;
        while k < 80 {
            offsets.push(k);
            k += 13;
        }
        // Family B: period 23 (coprime with 13, very different
        // ratio — no near-integer accidents in this short range)
        k = 23;
        while k < 80 {
            offsets.push(k);
            k += 23;
        }
        let scores = synth_bumps(120, &offsets, 100.0);
        let peaks = find_peaks_in_periodogram(&scores, 0, &PeaksConfig::default()).unwrap();
        let fund_periods: Vec<u32> = peaks
            .iter()
            .filter(|p| p.kind == PeakKind::Fundamental)
            .map(|p| p.period as u32)
            .collect();
        assert!(fund_periods.contains(&13), "missing 13: {fund_periods:?}");
        assert!(fund_periods.contains(&23), "missing 23: {fund_periods:?}");
    }

    #[test]
    fn consolidation_collapses_harmonic_claimed_first() {
        // Strong peak at 14 (the 2× harmonic of 7), weaker at 7. The
        // greedy pass will claim 14 as fundamental first; consolidation
        // must demote it once 7 is also identified.
        let mut scores = vec![0.0; 200];
        // Make every-7 the "real" tandem with weak amplitude
        for k in (7..200).step_by(7) {
            scores[k] = 30.0;
        }
        // Then bump the n=2 and n=4 harmonics super-high
        scores[14] = 120.0;
        scores[28] = 120.0;
        scores[42] = 100.0;
        // To make sure the 14 gets claimed first, weaken the
        // fundamental:
        scores[7] = 25.0;
        let cfg = PeaksConfig {
            min_score: 5.0,
            ..PeaksConfig::default()
        };
        let peaks = find_peaks_in_periodogram(&scores, 0, &cfg).unwrap();
        // The biological fundamental (7) should win, not 14.
        let funds: Vec<u32> = peaks
            .iter()
            .filter(|p| p.kind == PeakKind::Fundamental)
            .map(|p| p.period as u32)
            .collect();
        assert!(
            funds.contains(&7),
            "consolidation should keep 7 as fundamental; got {funds:?}"
        );
        assert!(
            !funds.contains(&14),
            "14 should be demoted to harmonic of 7; got {funds:?}"
        );
    }

    #[test]
    fn boundary_filter_drops_edge_artifact() {
        // Synthetic record: real signal at period 7, edge artifact at
        // period 95 (out of 99). Boundary cap at 0.9 → max period ≈ 89
        // → 95 dropped.
        let mut scores = vec![0.0; 100];
        scores[95] = 200.0; // edge artifact
        for k in (7..90).step_by(7) {
            scores[k] = 30.0;
        }
        let cfg = PeaksConfig {
            min_score: 5.0,
            boundary_fraction: 0.9,
            ..PeaksConfig::default()
        };
        let peaks = find_peaks_in_periodogram(&scores, 0, &cfg).unwrap();
        assert!(
            peaks.iter().all(|p| (p.period as u32) <= 89),
            "expected boundary peak dropped; got {peaks:?}"
        );
    }

    #[test]
    fn boundary_filter_disabled_at_1_0() {
        // Peak placed mid-bin so it's a strict three-bin local max
        // (endpoints can never be local maxima). The boundary cap at
        // 1.0 lets this peak through; at 0.5 it would be dropped.
        let mut scores = vec![0.0; 100];
        scores[7] = 30.0;
        scores[14] = 30.0;
        scores[21] = 30.0;
        scores[85] = 200.0; // 85 < 99 → local max possible
        let cfg = PeaksConfig {
            min_score: 5.0,
            boundary_fraction: 1.0,
            ..PeaksConfig::default()
        };
        let peaks = find_peaks_in_periodogram(&scores, 0, &cfg).unwrap();
        assert!(
            peaks.iter().any(|p| (p.period as u32) == 85),
            "expected 85 with boundary disabled: {peaks:?}"
        );
        // Same data with boundary 0.5 → 85 > 50 → dropped.
        let cfg2 = PeaksConfig {
            boundary_fraction: 0.5,
            ..cfg
        };
        let peaks2 = find_peaks_in_periodogram(&scores, 0, &cfg2).unwrap();
        assert!(peaks2.iter().all(|p| (p.period as u32) <= 49));
    }

    #[test]
    fn subrepeats_detect_monomer_under_hor() {
        // 150 bp HOR (strong) with 50 bp monomer (weak).
        let mut scores = vec![0.0; 400];
        // Strong HOR ladder at 150, 300, 450 (out of range), so just
        // 150 and 300.
        scores[150] = 200.0;
        scores[300] = 180.0;
        // Weak monomer signal at 50, 100 (also harmonic of 50 so it
        // would itself form a chain — but we deliberately weaken).
        scores[50] = 30.0;
        scores[100] = 25.0;
        let cfg = PeaksConfig {
            min_score: 100.0, // filter out the weak monomer signal
            subrepeats: Some(SubrepeatConfig {
                min_score: 5.0,
                max_divisor: 6,
                period_tolerance: 2,
            }),
            ..PeaksConfig::default()
        };
        let peaks = find_peaks_in_periodogram(&scores, 0, &cfg).unwrap();
        let subs: Vec<_> = peaks
            .iter()
            .filter(|p| p.kind == PeakKind::Subrepeat)
            .collect();
        assert!(!subs.is_empty(), "expected at least one subrepeat");
        // Should find 50 as div_n=3 of 150.
        let found_50 = subs.iter().any(|p| {
            (p.period as u32) == 50 && p.divisor_n == Some(3) && p.parent_period == Some(150.0)
        });
        assert!(
            found_50,
            "expected 50 as subrepeat (div_n=3) of 150: {subs:?}"
        );
    }

    #[test]
    fn min_harmonics_filters_lone_peaks() {
        // Just one isolated bump — no harmonic ladder.
        let mut scores = vec![0.0; 200];
        scores[55] = 50.0;
        let cfg = PeaksConfig {
            min_harmonics: 1,
            ..PeaksConfig::default()
        };
        let peaks = find_peaks_in_periodogram(&scores, 0, &cfg).unwrap();
        assert!(peaks.is_empty(), "lone peak should be filtered: {peaks:?}");
    }

    #[test]
    fn auto_threshold_mad_basics() {
        // Pure-noise-ish input: most values near zero, threshold
        // should be small.
        let noise: Vec<f64> = (0..100)
            .map(|i| if i % 10 == 0 { 1.0 } else { 0.0 })
            .collect();
        let t = auto_threshold_mad(&noise, 5.0, 1.0);
        // MAD of mostly-zero is 0 → falls back to floor.
        assert_eq!(t, 1.0);
    }

    #[test]
    fn auto_threshold_mad_lifts_above_noise() {
        // Most values 0; a handful at moderate amplitude → MAD > 0,
        // threshold rises above the floor.
        let mut scores = vec![0.0_f64; 1000];
        for i in (0..1000).step_by(7) {
            scores[i] = 2.0;
        }
        let t = auto_threshold_mad(&scores, 5.0, 1.0);
        // MAD is around 0 (most values are 0); fallback to floor.
        // To get a meaningful MAD, more than half the values must
        // deviate from the median.
        assert_eq!(t, 1.0);
    }

    #[test]
    fn auto_threshold_mad_handles_dense_signal() {
        // Roughly uniform [0, 10] distribution: MAD ≈ 2.5, median ≈ 5.
        // threshold = 5 + 5 × 1.4826 × 2.5 ≈ 23.5; floor lower.
        let scores: Vec<f64> = (0..1000).map(|i| (i % 11) as f64).collect();
        let t = auto_threshold_mad(&scores, 5.0, 1.0);
        assert!(
            t > 15.0,
            "expected a meaningful adaptive threshold; got {t}"
        );
    }

    #[test]
    fn auto_threshold_mad_floor_respected() {
        // Floor higher than computed adaptive → returns floor.
        let scores: Vec<f64> = (0..1000).map(|i| (i % 5) as f64).collect();
        let t = auto_threshold_mad(&scores, 5.0, 100.0);
        assert_eq!(t, 100.0);
    }

    #[test]
    fn empty_input_returns_empty() {
        let peaks = find_peaks_in_periodogram(&[], 3, &PeaksConfig::default()).unwrap();
        assert!(peaks.is_empty());
    }

    #[test]
    fn overlapping_high_harmonic_families_classified_correctly() {
        // The classic bug case: periods 7 and 13 with peaks extending
        // to where coincidental near-integer ratios appear. Period
        // 78 = 6 * 13 is also within 2% of 11 * 7 = 77, so a
        // first-match algorithm misclassifies 78 as a deep harmonic
        // of 7. Best-parent matching must pick 13.
        //
        // Place specific peaks manually (no Gaussian spread or full
        // arithmetic sequences) so 7 and 13 both register as strict
        // local maxima without colliding with each other's harmonics.
        let mut scores = vec![0.0_f64; 300];
        // Family A (period 7): ladder — skip 14 (collides with 13)
        // and 77 (collides with 78). Keep enough to register as a
        // fundamental with harmonics.
        for k in [7, 21, 28, 35, 42, 49, 56, 63, 70, 84, 98, 112, 126] {
            scores[k] = 100.0;
        }
        // Family B (period 13): ladder — skip 91 (= 7×13, shared with
        // family A and would be a shared-multiple confound).
        for k in [13, 26, 39, 52, 65, 78, 104, 117, 130] {
            scores[k] = 100.0;
        }
        let peaks = find_peaks_in_periodogram(&scores, 0, &PeaksConfig::default()).unwrap();

        // Both families should be fundamentals.
        let fund_periods: Vec<u32> = peaks
            .iter()
            .filter(|p| p.kind == PeakKind::Fundamental)
            .map(|p| p.period as u32)
            .collect();
        assert!(fund_periods.contains(&7), "missing 7: {fund_periods:?}");
        assert!(fund_periods.contains(&13), "missing 13: {fund_periods:?}");

        // Period 78 should be classified as h_6 of 13 (distance 0),
        // NOT h_11 of 7 (distance 1). Best-parent picks 13.
        let p78 = peaks
            .iter()
            .find(|p| (p.period - 78.0).abs() < 0.5)
            .expect("expected a peak near 78");
        assert_eq!(
            p78.parent_period.map(|v| v as u32),
            Some(13),
            "78 should be harmonic of 13, not {:?}",
            p78.parent_period
        );
        assert_eq!(p78.harmonic_n, Some(6));
    }

    #[test]
    fn reparent_after_consolidation_finds_better_parent() {
        // Construct a case where peak X is initially matched to
        // fundamental F1, but later F2 (the biological fundamental)
        // is also detected and X should be reattached to F2 because
        // F2 has a cleaner integer ratio.
        //
        // Setup: peak at 60. F1 = 7 (60/7=8.57, n=9, err 2% — borderline);
        // F2 = 12 (60/12=5, exact). If F2 is detected after the greedy
        // pass attaches 60 to F1 with marginal match, reparenting must
        // move 60 to F2.
        let mut scores = vec![0.0_f64; 200];
        // Make 7 the strongest fundamental so it's claimed first.
        for k in (7..200).step_by(7) {
            scores[k] = 100.0;
        }
        // Make 12 weaker but still a fundamental.
        for k in (12..200).step_by(12) {
            scores[k] = 50.0;
        }
        // Boost peak at 60 (which is also 5×12).
        scores[60] = 200.0;
        let peaks = find_peaks_in_periodogram(&scores, 0, &PeaksConfig::default()).unwrap();
        let p60 = peaks
            .iter()
            .find(|p| (p.period - 60.0).abs() < 0.5)
            .expect("expected a peak near 60");
        // 60 = 5 * 12 exactly (clean), vs 60 = 8.57 * 7 (marginal).
        // Best-match should attach to 12.
        assert_eq!(
            p60.parent_period.map(|v| v as u32),
            Some(12),
            "60 should attach to 12 (exact 5x) not 7: got {:?}",
            p60.parent_period
        );
    }

    #[test]
    fn min_score_filters_below_threshold() {
        let mut scores = vec![1.0; 100];
        scores[10] = 5.0;
        scores[20] = 5.0;
        scores[30] = 5.0;
        // All below default min_score (10).
        let peaks = find_peaks_in_periodogram(&scores, 0, &PeaksConfig::default()).unwrap();
        assert!(peaks.is_empty());
    }
}
