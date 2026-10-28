//! Exact sliding-window scoring along anchored diagonals.
//!
//! Each cell value is the same `W`-residue window sum the exact engine
//! computes (`sliding.rs`): forward `Σ_k M(q-k, s-k)`, reverse
//! `Σ_k Mc(q-k, s+k)` with `Mc` scoring the complemented query. Emission
//! uses the exact engine's validity bounds, pixel placement, anti-diagonal
//! rule and scaling, so every pixel written here is `<=` the exact engine's
//! pixel and equal to the exact value of each scored cell.

use std::sync::atomic::{AtomicU32, Ordering};

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
    pub coverage: Option<&'a [Vec<AtomicU32>]>,
}

/// Batches coverage increments: consecutive cells of a diagonal usually
/// land in the same pixel, so the atomic add happens once per pixel
/// crossing rather than once per cell.
struct CovAcc<'a> {
    buf: Option<&'a [Vec<AtomicU32>]>,
    channel: usize,
    pixel: usize,
    n: u32,
}

impl<'a> CovAcc<'a> {
    fn new(buf: Option<&'a [Vec<AtomicU32>]>) -> Self {
        Self {
            buf,
            channel: 0,
            pixel: 0,
            n: 0,
        }
    }

    #[inline]
    fn add(&mut self, channel: usize, pixel: usize) {
        if self.buf.is_none() {
            return;
        }
        if self.n > 0 && (channel, pixel) == (self.channel, self.pixel) {
            self.n += 1;
            return;
        }
        self.flush();
        (self.channel, self.pixel, self.n) = (channel, pixel, 1);
    }

    fn flush(&mut self) {
        if let (Some(buf), true) = (self.buf, self.n > 0) {
            buf[self.channel][self.pixel].fetch_add(self.n, Ordering::Relaxed);
            self.n = 0;
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

/// Score one run; returns the number of cells scored. The cell interval
/// covers every window overlapping an anchor seed, `extend` beyond, and
/// through to any fill partner in a neighbouring chunk.
pub(crate) fn score_run(ctx: &ScoreCtx<'_>, canvas: &Canvas<'_>, run: &Run) -> u64 {
    let w = ctx.window;
    let qlen = ctx.q_enc.len();
    let slen = ctx.s_enc.len();
    if w == 0 || w >= qlen {
        return 0;
    }
    let win2 = w / 2;
    let mut lo = (run.first_p as usize).saturating_sub(ctx.extend);
    let mut hi = run.last_p as usize + ctx.span + w - 2 + ctx.extend;
    if let Some(l) = run.lo_fill {
        lo = lo.min(l as usize);
    }
    if let Some(h) = run.hi_fill {
        hi = hi.max(h as usize);
    }
    // Cells whose pixel column can land on the canvas.
    lo = lo.max(canvas.q0 + win2).max(w);
    hi = hi
        .min(qlen - 1)
        .min(canvas.q0 + canvas.width * canvas.zoom + win2);

    if run.reverse {
        score_reverse(ctx, canvas, run.diag as i64, lo, hi, slen)
    } else {
        score_forward(ctx, canvas, run.diag as i64 - qlen as i64, lo, hi, slen)
    }
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
    let mut cov = CovAcc::new(canvas.coverage);
    let mut qi = lo;
    loop {
        let si = at(qi);
        if sum > 0 {
            emit(ctx, canvas, &mut cov, sum, qi, si, false);
        }
        cells += 1;
        if qi == hi {
            cov.flush();
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
    let mut cov = CovAcc::new(canvas.coverage);
    let mut qi = lo;
    loop {
        let si = at(qi);
        if sum > 0 {
            emit(ctx, canvas, &mut cov, sum, qi, si, true);
        }
        cells += 1;
        if qi == hi {
            cov.flush();
            break;
        }
        qi += 1;
        let si = at(qi);
        sum += t.rev(s[si], q[qi]) - t.rev(s[si + w], q[qi - w]);
    }
    cells
}

/// Mirror of `sliding::emit_pixel`, offset by the canvas origin and routed
/// to the class channel of the cell. Coverage counts every positive cell
/// of the pixel, including those the anti-diagonal rule drops from the
/// score channels: it measures similarity, not line drawing.
#[inline]
fn emit(
    ctx: &ScoreCtx<'_>,
    canvas: &Canvas<'_>,
    cov: &mut CovAcc<'_>,
    score: i32,
    qi: usize,
    si: usize,
    reverse: bool,
) {
    let win2 = ctx.window / 2;
    let zoom = canvas.zoom;
    let qc = qi - win2;
    let sc = if reverse { si + win2 } else { si - win2 };
    if qc < canvas.q0 || sc < canvas.s0 {
        return;
    }
    let x_abs = qc / zoom;
    let y_abs = sc / zoom;
    let x = x_abs - canvas.q0 / zoom;
    let y = y_abs - canvas.s0 / zoom;
    if x >= canvas.width || y >= canvas.height {
        return;
    }
    let class: RepeatClass = ctx
        .thresholds
        .classify(ctx.q_depth[qc].max(ctx.s_depth[sc]));
    cov.add(
        reverse as usize * RepeatClass::COUNT + class.index(),
        y * canvas.width + x,
    );
    let q_local = (qc - x_abs * zoom) as u32;
    let s_local = (sc - y_abs * zoom) as u32;
    if !keep_dot(zoom as u32, q_local, s_local, reverse) {
        return;
    }
    let scaled = (score as i64 * ctx.pixel_fac as i64 / ctx.window as i64).clamp(0, 255) as u8;
    let views = if reverse { canvas.rev } else { canvas.fwd };
    views[class.index()].max_merge(x, y, scaled);
}
