//! Sliding-window sum recurrence (spec §4.1.4) — the dotplot inner loop.
//!
//! Faithful port of `doCalculateImage` at
//! `third_party/seqtools/dotterApp/dotplot.c:1308`. Single-threaded for
//! Phase 1; Phase 3 layers rayon on top by chunking the subject axis.
//!
//! ## What this computes
//!
//! For a fixed subject position `s` and query position `q`, define
//! ```text
//!     score(q, s) = Σ_{k = 0..W} matrix[ subject[s-W+1+k], query[q-W+1+k] ]
//! ```
//! the sum of substitution-matrix scores along the diagonal landing at
//! `(q, s)`, over a window of length `W`. The recurrence is
//! ```text
//!     score(q, s) = score(q-1, s-1) + add(s, q) - del(s, q)
//! ```
//! where `add(s, q) = matrix[subject[s], query[q]]` and
//! `del(s, q) = matrix[subject[s-W], query[q-W]]`. In the warm-up region
//! (q < W or s < W), `del` is replaced by zero.
//!
//! Two ping-pong row buffers (`sum1`, `sum2`) carry `score(_, s-1)` into
//! the next iteration without per-iteration allocation.

use crate::antidiag::keep_dot;
use crate::pixel::PixelMap;
use crate::score_vec::ScoreVec;

/// Direction along the subject axis. Forward iterates `s = 0..slen`,
/// reverse iterates `s = slen-1..=0`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Direction {
    Forward,
    Reverse,
}

impl Direction {
    /// `+1` for forward, `-1` for reverse — mirroring C's `incrementVal`.
    #[inline]
    pub fn step(self) -> i32 {
        match self {
            Direction::Forward => 1,
            Direction::Reverse => -1,
        }
    }
}

/// One full pass of the sliding-window sum recurrence over a (query,
/// subject) pair, max-merging into `out`.
///
/// * `score_vec` — precomputed `(n+1) × qlen` row-major scores. Row `n`
///   is the synthetic zero row used for unknown subject residues.
/// * `subject_encoded` — slen-length subject sequence as alphabet indices
///   (use [`ScoreVec::subject_row`] to map each byte to a valid row).
/// * `window` — sliding window size `W`. Must be `>= 1`.
/// * `zoom` — `zoom × zoom` block per output pixel.
/// * `pixel_fac` — multiplier in `min(255, score * pixel_fac / W)`.
/// * `direction` — forward or reverse along the subject axis.
/// * `out` — pixelmap to max-merge into.
///
/// # Determinism
///
/// Single-threaded; iteration order is `s = 0..slen` (or its reverse),
/// `q = 0..pepqlen`. With max-merge being associative+commutative, parallel
/// chunking by subject in Phase 3 will preserve byte-identical output.
pub fn sliding_window_pass(
    score_vec: &ScoreVec,
    subject_encoded: &[u8],
    window: u32,
    zoom: u32,
    pixel_fac: u32,
    direction: Direction,
    out: &mut PixelMap,
) {
    let qlen = score_vec.qlen;
    let slen = subject_encoded.len();
    let w = window as usize;
    let win2 = w / 2;
    let zoom_us = zoom as usize;

    if w == 0 || qlen == 0 || slen == 0 {
        return;
    }

    // Ping-pong buffers for the previous and current row of partial sums.
    let mut sum1 = vec![0_i32; qlen];
    let mut sum2 = vec![0_i32; qlen];

    // Iteration bounds per direction:
    //   forward: s = 0, 1, ..., slen - 1     (s_in_valid_range when s >= W)
    //   reverse: s = slen-1, ..., 0          (s_in_valid_range when s <  slen - W)
    let (s_start, s_end_excl, step): (i64, i64, i64) = match direction {
        Direction::Forward => (0, slen as i64, 1),
        Direction::Reverse => (slen as i64 - 1, -1, -1),
    };

    // Iteration counter for the ping-pong parity. C uses `sIdx & 1`, which
    // assumes sIdx is non-negative; for the reverse pass it still works
    // because slen-1 has the same parity as 0..slen-1 in the same order.
    // We use an explicit counter for clarity and determinism.
    let mut iter_idx: usize = 0;

    let mut s_signed = s_start;
    while s_signed != s_end_excl {
        let s = s_signed as usize;
        let s_row = score_vec.subject_row(subject_encoded[s]) as usize;
        let add_row = score_vec.row(s_row);

        // Row that's leaving the window. For forward: s-W (rows `W` steps
        // ago). For reverse: s+W (rows `W` steps ago in iteration order,
        // which is `W` ahead in array order).
        let del_row: &[i32] = match direction {
            Direction::Forward => {
                if s >= w {
                    let prev_s = s - w;
                    let prev_row = score_vec.subject_row(subject_encoded[prev_s]) as usize;
                    score_vec.row(prev_row)
                } else {
                    score_vec.row(score_vec.unknown_row() as usize)
                }
            }
            Direction::Reverse => {
                if s + w < slen {
                    let prev_s = s + w;
                    let prev_row = score_vec.subject_row(subject_encoded[prev_s]) as usize;
                    score_vec.row(prev_row)
                } else {
                    score_vec.row(score_vec.unknown_row() as usize)
                }
            }
        };

        // Pick ping/pong rows. `oldsum` holds the previous iteration's
        // partial sums; `newsum` will be overwritten.
        let (oldsum, newsum): (&[i32], &mut [i32]) = if iter_idx & 1 == 1 {
            (sum2.as_slice(), sum1.as_mut_slice())
        } else {
            (sum1.as_slice(), sum2.as_mut_slice())
        };

        // Initialise newsum[0] from the first column of the addrow.
        newsum[0] = add_row[0];

        // Warm-up region: q in [1, min(W, qlen)). No delrow subtraction.
        let q_warmup_end = w.min(qlen);
        for q in 1..q_warmup_end {
            newsum[q] = oldsum[q - 1] + add_row[q];
        }

        // Steady state: q in [W, qlen). Full recurrence; this is where
        // pixels are emitted.
        if w < qlen {
            // For the steady state we need to know whether s is in the
            // "valid range" (full window built up):
            //   forward: valid when s >= W
            //   reverse: valid when s < slen - W
            let s_valid = match direction {
                Direction::Forward => s >= w,
                Direction::Reverse => s + w < slen,
            };

            for q in w..qlen {
                let new_val = oldsum[q - 1] + add_row[q] - del_row[q - w];
                newsum[q] = new_val;

                if !s_valid || new_val <= 0 {
                    continue;
                }
                // Pixel emission. Coordinates mirror C dotplot.c:1391–1424.
                emit_pixel(
                    new_val, q, s, window, win2, zoom_us, pixel_fac, direction, out,
                );
            }
        }

        iter_idx += 1;
        s_signed += step;
    }
}

/// Emit one pixel candidate at sub-pixel `(q_idx, s_idx)`. Applies the
/// anti-diagonal suppression rule (spec §4.1.6) and the
/// `min(255, score * pixel_fac / W)` scale step, max-merging into `out`.
#[inline]
#[allow(clippy::too_many_arguments)]
fn emit_pixel(
    score: i32,
    q_idx: usize,
    s_idx: usize,
    window: u32,
    win2: usize,
    zoom: usize,
    pixel_fac: u32,
    direction: Direction,
    out: &mut PixelMap,
) {
    // C uses `incrementVal*win2` for the subject offset. Forward: subtract
    // win2; reverse: add win2 (since incrementVal = -1).
    let dotposq = (q_idx - win2) / zoom;
    let dotposs_centered = match direction {
        Direction::Forward => s_idx - win2,
        Direction::Reverse => s_idx + win2,
    };
    let dotposs = dotposs_centered / zoom;

    if dotposq >= out.width() || dotposs >= out.height() {
        return;
    }

    let q_local = (q_idx - win2) - dotposq * zoom;
    let s_local = dotposs_centered - dotposs * zoom;

    let reverse = matches!(direction, Direction::Reverse);
    if !keep_dot(zoom as u32, q_local as u32, s_local as u32, reverse) {
        return;
    }

    // The C does: `val = newsum * pixelFac / slidingWinSize`, then
    // `(val > 255 ? 255 : (unsigned char)val)`. Negative `val` cannot
    // reach here because of the `*newsum > 0` gate; we still clamp for
    // safety against signed overflow in pathological inputs.
    let scaled_signed = score as i64 * pixel_fac as i64 / window as i64;
    let scaled = scaled_signed.clamp(0, 255) as u8;
    out.max_merge(dotposq, dotposs, scaled);
}
