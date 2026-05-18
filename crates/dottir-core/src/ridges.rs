//! Vector ridge extraction — visualisation aid for the GUI.
//!
//! Walks the pixmap along both diagonal directions and finds
//! coherent runs of "lit" cells (≥ `threshold`). Each run becomes a
//! [`Ridge`] with start, end, direction, and peak strength. The GUI
//! draws these as anti-aliased line segments over the raster, which
//! masks the per-window intensity oscillation that imperfect-
//! homology diagonals show under a stretched greyramp without
//! altering the underlying pixmap (spec §4.1 stays exact).
//!
//! Ridge extraction is a pure function over the same pixmap the
//! kernel produced. No spec deviation; the raster is still the
//! scientific truth — the overlay is presentational.
//!
//! ## Algorithm
//!
//! For each of the two diagonal directions:
//! 1. For every diagonal (parameterised by `d = q − s` for forward,
//!    `d = q + s` for reverse) walk one cell at a time.
//! 2. Track a single in-progress ridge per diagonal. A lit cell
//!    starts it (or extends the current). A non-lit cell increments
//!    a gap counter; once `gap_count > max_gap` the ridge ends and
//!    is emitted iff its total lit length ≥ `min_length`.
//! 3. Diagonal end closes any active ridge.
//!
//! Complexity: O(W × H). One pass per direction. The kernel that
//! produced the pixmap is far more expensive, so this extraction
//! is essentially free on top.

use crate::DotPlot;

/// Parameters that tune ridge sensitivity. `threshold` is the pixmap
/// value above which a cell counts as "lit"; in the GUI we tie it
/// to `greyramp.white` so the overlay follows the user's noise floor.
/// `min_length` filters out short coincidental runs; `max_gap`
/// bridges brief drops in homology (a single mismatch cluster can
/// kill 2-3 consecutive windows even on a 95 % conserved diagonal).
#[derive(Clone, Copy, Debug)]
pub struct RidgeParams {
    pub threshold: u8,
    pub min_length: u32,
    pub max_gap: u32,
}

impl Default for RidgeParams {
    fn default() -> Self {
        Self {
            threshold: 40,
            min_length: 8,
            max_gap: 2,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RidgeDirection {
    Forward,
    Reverse,
}

/// One detected ridge — a coherent diagonal run in the pixmap.
/// Coordinates are in pixmap-pixel space (`0..plot.width`,
/// `0..plot.height`). The GUI translates to screen coords by
/// subtracting the current `view_offset`.
#[derive(Clone, Copy, Debug)]
pub struct Ridge {
    pub start: (u32, u32),
    pub end: (u32, u32),
    pub direction: RidgeDirection,
    /// Peak pixmap value seen along the ridge — useful if a caller
    /// wants to fade weak ridges or color-code by strength.
    pub strength: u8,
}

/// Convenience wrapper around [`extract_ridges_from_pixels`] that
/// reads `width`/`height` from a [`DotPlot`] and operates on its
/// `pixels` field.
pub fn extract_ridges(plot: &DotPlot, params: &RidgeParams) -> Vec<Ridge> {
    extract_ridges_from_pixels(&plot.pixels, plot.width, plot.height, params)
}

/// Find coherent diagonal runs in a row-major pixmap.
///
/// Returns ridges in *both* directions: forward (slope +1 in
/// (q, s) space) and reverse (slope −1). The GUI can colour the
/// two channels differently when both strands are visible.
pub fn extract_ridges_from_pixels(
    pixels: &[u8],
    width: u32,
    height: u32,
    params: &RidgeParams,
) -> Vec<Ridge> {
    let mut ridges = Vec::new();
    let w = width as usize;
    let h = height as usize;
    if pixels.len() != w * h || w == 0 || h == 0 {
        return ridges;
    }
    let thr = params.threshold;
    let min_len = params.min_length as usize;
    let max_gap = params.max_gap as usize;

    // Forward diagonals: q − s = d. The diagonal at d ≥ 0 starts at
    // (d, 0); at d < 0 it starts at (0, −d). Each step is (+1, +1)
    // until we leave the grid.
    let h_i = h as isize;
    let w_i = w as isize;
    for d in -(h_i - 1)..=(w_i - 1) {
        let (q0, s0) = if d >= 0 {
            (d as usize, 0)
        } else {
            (0, (-d) as usize)
        };
        let len = (w - q0).min(h - s0);
        let coords = (0..len).map(|k| (q0 + k, s0 + k));
        scan_diagonal(
            pixels,
            w,
            coords,
            RidgeDirection::Forward,
            thr,
            min_len,
            max_gap,
            &mut ridges,
        );
    }

    // Reverse diagonals: q + s = d. The diagonal at d < h starts at
    // (0, d); at d ≥ h it starts at (d − (h − 1), h − 1). Each step
    // is (+1, −1) until we leave the grid.
    for d in 0..(w + h - 1) {
        let (q0, s0) = if d < h {
            (0usize, d)
        } else {
            (d - (h - 1), h - 1)
        };
        let len = (w - q0).min(s0 + 1);
        let coords = (0..len).map(|k| (q0 + k, s0 - k));
        scan_diagonal(
            pixels,
            w,
            coords,
            RidgeDirection::Reverse,
            thr,
            min_len,
            max_gap,
            &mut ridges,
        );
    }

    ridges
}

#[allow(clippy::too_many_arguments)]
fn scan_diagonal(
    pixels: &[u8],
    w: usize,
    coords: impl Iterator<Item = (usize, usize)>,
    dir: RidgeDirection,
    thr: u8,
    min_len: usize,
    max_gap: usize,
    out: &mut Vec<Ridge>,
) {
    struct Active {
        start_q: u32,
        start_s: u32,
        end_q: u32,
        end_s: u32,
        max_val: u8,
        lit_count: usize,
        gap_count: usize,
    }
    let mut active: Option<Active> = None;

    let flush = |a: Active, out: &mut Vec<Ridge>| {
        if a.lit_count >= min_len {
            out.push(Ridge {
                start: (a.start_q, a.start_s),
                end: (a.end_q, a.end_s),
                direction: dir,
                strength: a.max_val,
            });
        }
    };

    for (q, s) in coords {
        let v = pixels[s * w + q];
        if v >= thr {
            match active.as_mut() {
                None => {
                    active = Some(Active {
                        start_q: q as u32,
                        start_s: s as u32,
                        end_q: q as u32,
                        end_s: s as u32,
                        max_val: v,
                        lit_count: 1,
                        gap_count: 0,
                    });
                }
                Some(r) => {
                    r.end_q = q as u32;
                    r.end_s = s as u32;
                    if v > r.max_val {
                        r.max_val = v;
                    }
                    r.lit_count += 1;
                    r.gap_count = 0;
                }
            }
        } else if let Some(r) = active.as_mut() {
            r.gap_count += 1;
            if r.gap_count > max_gap {
                if let Some(a) = active.take() {
                    flush(a, out);
                }
            }
        }
    }
    if let Some(a) = active {
        flush(a, out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a pixmap by drawing horizontal/diagonal/single-cell
    /// strokes from a "pattern" expression — keeps the tests
    /// readable. Cells default to 0; lit cells take a fixed value.
    fn pixmap(w: u32, h: u32) -> Vec<u8> {
        vec![0u8; (w * h) as usize]
    }
    fn set(p: &mut [u8], w: u32, q: u32, s: u32, v: u8) {
        p[(s * w + q) as usize] = v;
    }

    #[test]
    fn single_forward_diagonal_emits_one_ridge() {
        let w = 20;
        let h = 20;
        let mut p = pixmap(w, h);
        // Lit diagonal at d = 0, cells (3, 3) .. (15, 15).
        for k in 3..=15 {
            set(&mut p, w, k, k, 200);
        }
        let r = extract_ridges_from_pixels(
            &p,
            w,
            h,
            &RidgeParams {
                threshold: 40,
                min_length: 5,
                max_gap: 0,
            },
        );
        assert_eq!(r.len(), 1, "expected exactly one ridge, got {r:?}");
        let ridge = r[0];
        assert_eq!(ridge.direction, RidgeDirection::Forward);
        assert_eq!(ridge.start, (3, 3));
        assert_eq!(ridge.end, (15, 15));
        assert_eq!(ridge.strength, 200);
    }

    #[test]
    fn parallel_forward_diagonals_emit_two_ridges() {
        let w = 20;
        let h = 20;
        let mut p = pixmap(w, h);
        for k in 0..10 {
            set(&mut p, w, k, k, 200); // d = 0
            set(&mut p, w, 3 + k, k, 200); // d = 3
        }
        let r = extract_ridges_from_pixels(
            &p,
            w,
            h,
            &RidgeParams {
                threshold: 40,
                min_length: 5,
                max_gap: 0,
            },
        );
        assert_eq!(r.len(), 2, "expected two ridges, got {r:?}");
        // Both forward; we don't promise order, so collect starts.
        let starts: Vec<(u32, u32)> = r.iter().map(|x| x.start).collect();
        assert!(starts.contains(&(0, 0)));
        assert!(starts.contains(&(3, 0)));
    }

    #[test]
    fn short_run_below_min_length_is_dropped() {
        let w = 20;
        let h = 20;
        let mut p = pixmap(w, h);
        for k in 0..4 {
            set(&mut p, w, k, k, 200);
        }
        let r = extract_ridges_from_pixels(
            &p,
            w,
            h,
            &RidgeParams {
                threshold: 40,
                min_length: 5,
                max_gap: 0,
            },
        );
        assert!(r.is_empty(), "short run should be filtered, got {r:?}");
    }

    #[test]
    fn small_gap_bridged_when_max_gap_is_one() {
        let w = 20;
        let h = 20;
        let mut p = pixmap(w, h);
        // Long run with a single 1-cell gap at k = 5.
        for k in 0..10 {
            if k != 5 {
                set(&mut p, w, k, k, 200);
            }
        }
        let r = extract_ridges_from_pixels(
            &p,
            w,
            h,
            &RidgeParams {
                threshold: 40,
                min_length: 8,
                max_gap: 1,
            },
        );
        assert_eq!(r.len(), 1, "1-cell gap should bridge, got {r:?}");
        assert_eq!(r[0].start, (0, 0));
        assert_eq!(r[0].end, (9, 9));
    }

    #[test]
    fn large_gap_splits_into_two_ridges() {
        let w = 20;
        let h = 20;
        let mut p = pixmap(w, h);
        // Two runs separated by a 3-cell gap; max_gap = 2 → splits.
        for k in 0..6 {
            set(&mut p, w, k, k, 200);
        }
        for k in 9..15 {
            set(&mut p, w, k, k, 200);
        }
        let r = extract_ridges_from_pixels(
            &p,
            w,
            h,
            &RidgeParams {
                threshold: 40,
                min_length: 5,
                max_gap: 2,
            },
        );
        assert_eq!(r.len(), 2, "wide gap should split into 2, got {r:?}");
    }

    #[test]
    fn reverse_diagonal_detected_with_reverse_direction() {
        let w = 20;
        let h = 20;
        let mut p = pixmap(w, h);
        // Reverse diagonal q + s = 18 — cells (3, 15) .. (15, 3).
        for k in 0..=12 {
            set(&mut p, w, 3 + k, 15 - k, 200);
        }
        let r = extract_ridges_from_pixels(
            &p,
            w,
            h,
            &RidgeParams {
                threshold: 40,
                min_length: 5,
                max_gap: 0,
            },
        );
        assert!(
            r.iter().any(|x| x.direction == RidgeDirection::Reverse),
            "expected a reverse ridge, got {r:?}"
        );
    }

    #[test]
    fn below_threshold_emits_nothing() {
        let w = 10;
        let h = 10;
        let mut p = pixmap(w, h);
        for k in 0..8 {
            set(&mut p, w, k, k, 30); // below threshold of 40
        }
        let r = extract_ridges_from_pixels(
            &p,
            w,
            h,
            &RidgeParams {
                threshold: 40,
                min_length: 5,
                max_gap: 0,
            },
        );
        assert!(r.is_empty());
    }

    #[test]
    fn empty_or_mismatched_input_returns_empty() {
        let r = extract_ridges_from_pixels(&[], 0, 0, &RidgeParams::default());
        assert!(r.is_empty());
        let r = extract_ridges_from_pixels(&[0u8; 5], 10, 10, &RidgeParams::default());
        assert!(r.is_empty(), "pixels.len() != w*h should return empty");
    }
}
