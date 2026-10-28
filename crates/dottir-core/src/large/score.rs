//! Exact sliding-window scoring along anchored diagonals.
//!
//! Each cell value is the same `W`-residue window sum the exact engine
//! computes (`sliding.rs`): forward `Σ_k M(q-k, s-k)`, reverse
//! `Σ_k Mc(q-k, s+k)` with `Mc` scoring the complemented query. Emission
//! uses the exact engine's validity bounds, pixel placement, anti-diagonal
//! rule and scaling, so every pixel written here is `<=` the exact engine's
//! pixel and equal to the exact value of each scored cell.

use std::ops::Range;
use std::sync::atomic::{AtomicU64, Ordering};

use super::depth::{ClassThresholds, RepeatClass};
use super::hits::Run;
use crate::antidiag::keep_dot;
use crate::pixel::PixelView;

/// Residue-pair score table over encoded codes, with every code `>= n`
/// (N, sentinel) scoring 0 — the same convention as `ScoreVec`.
#[derive(Debug, Clone)]
pub(crate) struct PairTable {
    n: usize,
    fwd: Vec<i32>,
    rev: Vec<i32>,
}

impl PairTable {
    pub(crate) fn new(matrix: &crate::matrix::ScoreMatrix) -> Self {
        let n = matrix.size();
        let m = n + 1;
        let mut fwd = vec![0; m * m];
        let mut rev = vec![0; m * m];
        for s in 0..n {
            for q in 0..n {
                fwd[s * m + q] = matrix.get(s, q);
                let qc = crate::alphabet::complement_dna_encoded(q as u8) as usize;
                rev[s * m + q] = if qc < n { matrix.get(s, qc) } else { 0 };
            }
        }
        Self { n, fwd, rev }
    }

    #[inline]
    fn idx(&self, s: u8, q: u8) -> usize {
        (s as usize).min(self.n) * (self.n + 1) + (q as usize).min(self.n)
    }
    #[inline]
    fn fwd(&self, s: u8, q: u8) -> i32 {
        self.fwd[self.idx(s, q)]
    }
    #[inline]
    fn rev(&self, s: u8, q: u8) -> i32 {
        self.rev[self.idx(s, q)]
    }
}

/// Output geometry and the per-(strand, class) pixel views.
pub(crate) struct Canvas<'a> {
    pub q0: usize,
    pub s0: usize,
    pub zoom: usize,
    pub width: usize,
    pub height: usize,
    /// `views[class]` for the strand being scored.
    pub fwd: &'a [PixelView<'a>],
    pub rev: &'a [PixelView<'a>],
    /// Coverage counters `[strand * COUNT + class][pixel]`, when requested.
    pub coverage: Option<&'a [Vec<AtomicU64>]>,
}

/// Collects the cells of one diagonal pixel by pixel. Consecutive cells
/// of a diagonal mostly fall in the same pixel, and within one pixel the
/// anti-diagonal rule gives the same answer for every cell of the diagonal
/// (forward: `s_local - q_local` is constant; reverse: `s_local + q_local`
/// is). So per pixel crossing the accumulator keeps a running max score and
/// a cell count per class, and touches the shared pixel arrays only when
/// the diagonal leaves the pixel. The result is identical to emitting
/// every cell, because the scaling is monotone.
struct PixelAcc<'c, 'a> {
    ctx: &'c ScoreCtx<'c>,
    canvas: &'c Canvas<'a>,
    reverse: bool,
    /// Base ranges of the current pixel; empty until the first cell.
    q_px: Range<usize>,
    s_px: Range<usize>,
    /// Canvas pixel index of the current pixel, `None` when off canvas.
    pixel: Option<(usize, usize)>,
    keep: bool,
    best: [i32; RepeatClass::COUNT],
    count: [u64; RepeatClass::COUNT],
}

impl<'c, 'a> PixelAcc<'c, 'a> {
    fn new(ctx: &'c ScoreCtx<'c>, canvas: &'c Canvas<'a>, reverse: bool) -> Self {
        Self {
            ctx,
            canvas,
            reverse,
            q_px: 0..0,
            s_px: 0..0,
            pixel: None,
            keep: false,
            best: [0; RepeatClass::COUNT],
            count: [0; RepeatClass::COUNT],
        }
    }

    /// A cell with a positive window `score` at sequence position (qi, si).
    #[inline]
    fn cell(&mut self, score: i32, qi: usize, si: usize) {
        let win2 = self.ctx.window / 2;
        let qc = qi - win2;
        let sc = if self.reverse { si + win2 } else { si - win2 };
        if !(self.q_px.contains(&qc) && self.s_px.contains(&sc)) {
            self.flush();
            self.enter(qc, sc);
        }
        if self.pixel.is_none() {
            return;
        }
        let c = self
            .ctx
            .thresholds
            .classify(self.ctx.q_depth[qc].max(self.ctx.s_depth[sc]))
            .index();
        self.best[c] = self.best[c].max(score);
        self.count[c] += 1;
    }

    /// Start the pixel containing (qc, sc) — the mirror of
    /// `sliding::emit_pixel`, offset by the canvas origin.
    fn enter(&mut self, qc: usize, sc: usize) {
        let canvas = self.canvas;
        let zoom = canvas.zoom;
        let (x_abs, y_abs) = (qc / zoom, sc / zoom);
        self.q_px = x_abs * zoom..(x_abs + 1) * zoom;
        self.s_px = y_abs * zoom..(y_abs + 1) * zoom;
        self.pixel = None;
        // Canvas origins are multiples of the zoom.
        if qc < canvas.q0 || sc < canvas.s0 {
            return;
        }
        let x = x_abs - canvas.q0 / zoom;
        let y = y_abs - canvas.s0 / zoom;
        if x >= canvas.width || y >= canvas.height {
            return;
        }
        let q_local = (qc - x_abs * zoom) as u32;
        let s_local = (sc - y_abs * zoom) as u32;
        self.keep = keep_dot(zoom as u32, q_local, s_local, self.reverse);
        self.pixel = Some((x, y));
    }

    /// Write the current pixel. Coverage counts every positive cell,
    /// including those the anti-diagonal rule drops from the score
    /// channels: it measures similarity, not line drawing.
    fn flush(&mut self) {
        let Some((x, y)) = self.pixel else {
            return;
        };
        let canvas = self.canvas;
        let views = if self.reverse { canvas.rev } else { canvas.fwd };
        for c in 0..RepeatClass::COUNT {
            if self.count[c] == 0 {
                continue;
            }
            if let Some(cov) = canvas.coverage {
                cov[self.reverse as usize * RepeatClass::COUNT + c][y * canvas.width + x]
                    .fetch_add(self.count[c], Ordering::Relaxed);
            }
            if self.keep {
                let scaled = (self.best[c] as i64 * self.ctx.pixel_fac as i64
                    / self.ctx.window as i64)
                    .clamp(0, 255) as u8;
                views[c].max_merge(x, y, scaled);
            }
            self.best[c] = 0;
            self.count[c] = 0;
        }
    }
}

pub(crate) struct ScoreCtx<'a> {
    pub table: &'a PairTable,
    pub q_enc: &'a [u8],
    pub s_enc: &'a [u8],
    pub q_depth: &'a [u16],
    pub s_depth: &'a [u16],
    pub thresholds: ClassThresholds,
    pub window: usize,
    pub pixel_fac: u32,
    pub span: usize,
    pub extend: usize,
    pub lower_only: bool,
}

/// Score the cells covered by `runs` whose query position lies in `own`;
/// returns the number of cells scored. A run covers every window that
/// overlaps one of its anchor seeds, plus `extend` cells beyond. Runs on the
/// same diagonal are merged first, so each cell is scored once; ownership by
/// query position makes the result independent of how the query is cut into
/// chunks, provided `runs` holds every anchor within reach of `own`.
pub(crate) fn score_owned(
    ctx: &ScoreCtx<'_>,
    canvas: &Canvas<'_>,
    runs: &[Run],
    own: &Range<usize>,
) -> u64 {
    let w = ctx.window;
    let qlen = ctx.q_enc.len();
    let slen = ctx.s_enc.len();
    if w == 0 || w >= qlen {
        return 0;
    }
    let win2 = w / 2;
    // Cells whose pixel column can land on the canvas.
    let q_lo = own.start.max(canvas.q0 + win2).max(w);
    let q_hi = (own.end.min(qlen) - 1).min(canvas.q0 + canvas.width * canvas.zoom + win2);
    if q_lo > q_hi {
        return 0;
    }
    let mut cells = 0u64;
    let mut score = |reverse: bool, diag: u64, lo: usize, hi: usize| {
        let (lo, hi) = (lo.max(q_lo), hi.min(q_hi));
        if lo > hi {
            return;
        }
        cells += if reverse {
            score_reverse(ctx, canvas, diag as i64, lo, hi, slen)
        } else {
            score_forward(ctx, canvas, diag as i64 - qlen as i64, lo, hi, slen)
        };
    };
    let mut cur: Option<(bool, u64, usize, usize)> = None;
    for run in runs {
        let lo = (run.first_p as usize).saturating_sub(ctx.extend);
        let hi = run.last_p as usize + ctx.span + w - 2 + ctx.extend;
        match &mut cur {
            Some((r, d, _, c_hi)) if *r == run.reverse && *d == run.diag && lo <= *c_hi + 1 => {
                *c_hi = (*c_hi).max(hi);
            }
            _ => {
                if let Some((r, d, l, h)) = cur.take() {
                    score(r, d, l, h);
                }
                cur = Some((run.reverse, run.diag, lo, hi));
            }
        }
    }
    if let Some((r, d, l, h)) = cur {
        score(r, d, l, h);
    }
    cells
}

fn score_forward(
    ctx: &ScoreCtx<'_>,
    canvas: &Canvas<'_>,
    d: i64,
    mut lo: usize,
    mut hi: usize,
    slen: usize,
) -> u64 {
    let w = ctx.window;
    let win2 = w / 2;
    if ctx.lower_only && d < 0 {
        return 0;
    }
    // s = q + d must satisfy w <= s <= slen - 1, and its pixel row must be
    // on the canvas.
    let s_lo = (w as i64).max((canvas.s0 + win2) as i64);
    let s_hi = (slen as i64 - 1).min((canvas.s0 + canvas.height * canvas.zoom + win2) as i64);
    lo = lo.max((s_lo - d).max(0) as usize);
    let hi_s = s_hi - d;
    if hi_s < 0 {
        return 0;
    }
    hi = hi.min(hi_s as usize);
    if lo > hi {
        return 0;
    }
    let t = ctx.table;
    let (q, s) = (ctx.q_enc, ctx.s_enc);
    let at = |qi: usize| (qi as i64 + d) as usize;

    let mut sum: i32 = (0..w).map(|k| t.fwd(s[at(lo) - k], q[lo - k])).sum();
    let mut cells = 0u64;
    let mut acc = PixelAcc::new(ctx, canvas, false);
    let mut qi = lo;
    loop {
        if sum > 0 {
            acc.cell(sum, qi, at(qi));
        }
        cells += 1;
        if qi == hi {
            acc.flush();
            break;
        }
        qi += 1;
        let si = at(qi);
        sum += t.fwd(s[si], q[qi]) - t.fwd(s[si - w], q[qi - w]);
    }
    cells
}

fn score_reverse(
    ctx: &ScoreCtx<'_>,
    canvas: &Canvas<'_>,
    c: i64,
    mut lo: usize,
    mut hi: usize,
    slen: usize,
) -> u64 {
    let w = ctx.window;
    let win2 = w / 2;
    // s = c - q must satisfy 0 <= s <= slen - 1 - w, and the pixel row
    // (s + win2) / zoom must be on the canvas.
    let s_max = (slen as i64 - 1 - w as i64)
        .min((canvas.s0 + canvas.height * canvas.zoom) as i64 - win2 as i64 - 1);
    let s_min = (canvas.s0 as i64 - win2 as i64).max(0);
    if s_max < s_min {
        return 0;
    }
    lo = lo.max((c - s_max).max(0) as usize);
    let hi_c = c - s_min;
    if hi_c < 0 {
        return 0;
    }
    hi = hi.min(hi_c as usize);
    if ctx.lower_only {
        // q <= s = c - q
        if c < 0 {
            return 0;
        }
        hi = hi.min((c / 2) as usize);
    }
    if lo > hi {
        return 0;
    }
    let t = ctx.table;
    let (q, s) = (ctx.q_enc, ctx.s_enc);
    let at = |qi: usize| (c - qi as i64) as usize;

    let mut sum: i32 = (0..w).map(|k| t.rev(s[at(lo) + k], q[lo - k])).sum();
    let mut cells = 0u64;
    let mut acc = PixelAcc::new(ctx, canvas, true);
    let mut qi = lo;
    loop {
        if sum > 0 {
            acc.cell(sum, qi, at(qi));
        }
        cells += 1;
        if qi == hi {
            acc.flush();
            break;
        }
        qi += 1;
        let si = at(qi);
        sum += t.rev(s[si], q[qi]) - t.rev(s[si + w], q[qi - w]);
    }
    cells
}
