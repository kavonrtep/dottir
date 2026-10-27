//! Query streaming, the two-hit diagonal filter, and merging of surviving
//! anchors into diagonal runs.
//!
//! Hits are packed as `(diag << 32) | p` where `p` is the query seed start
//! and `diag` is `s - p + qlen` (forward) or `p + s + span - 1` (reverse —
//! the anti-diagonal `q + s` of the matched cells). Sorting the packed
//! values groups hits by diagonal, ordered along it, which is all the
//! two-hit test and the merge need.

use std::ops::Range;

use super::depth::{ClassThresholds, RepeatClass};
use super::index::SeedIndex;
use super::seed::for_each_window;

/// A run of two-hit anchors on one (anti-)diagonal.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Run {
    pub reverse: bool,
    pub diag: u64,
    /// Query start of the first anchor seed.
    pub first_p: u32,
    /// Query start of the last anchor seed.
    pub last_p: u32,
    pub n_anchors: u32,
    pub q_rec: u32,
    pub s_rec: u32,
    pub class: RepeatClass,
    /// Query start of a two-hit partner just outside the chunk core that
    /// would have merged with the first (`lo_fill`) / last (`hi_fill`)
    /// anchor. Scoring through to it makes the scored cell set independent
    /// of where chunk boundaries fall.
    pub lo_fill: Option<u32>,
    pub hi_fill: Option<u32>,
}

/// Everything a chunk worker reads.
pub(crate) struct HitCtx<'a> {
    pub index: &'a SeedIndex,
    pub q_enc: &'a [u8],
    pub q_breaks: &'a [usize],
    pub s_breaks: &'a [usize],
    /// Per-base depth for run classes; empty slices classify as unique.
    pub q_depth: &'a [u16],
    pub s_depth: &'a [u16],
    pub thresholds: ClassThresholds,
    /// Subject seed starts considered.
    pub s_range: Range<usize>,
    pub max_occ: u32,
    pub two_hit_dist: u32,
    pub forward: bool,
    pub reverse: bool,
    /// Self-comparison lower triangle: drop hits that can only produce
    /// cells with `q > s`.
    pub lower_only: bool,
}

#[derive(Debug, Default, Clone, Copy)]
pub(crate) struct HitStats {
    pub raw_hits: u64,
    pub capped_positions: u64,
    pub anchors: u64,
}

pub(crate) struct ChunkHits {
    pub runs: Vec<Run>,
    /// Query seed starts in the chunk core whose lookup was capped.
    pub capped: Vec<u32>,
    pub stats: HitStats,
}

/// Seeds starting in `stream` are looked up; only anchors whose start lies
/// in `core` are kept (the margin supplies two-hit partners).
pub(crate) fn chunk_hits(ctx: &HitCtx<'_>, core: Range<usize>, stream: Range<usize>) -> ChunkHits {
    let seed = ctx.index.seed();
    let span = seed.span() as u64;
    let qlen = ctx.q_enc.len() as u64;
    let mut fwd_hits: Vec<u64> = Vec::new();
    let mut rev_hits: Vec<u64> = Vec::new();
    let mut capped = Vec::new();
    let mut stats = HitStats::default();

    for_each_window(seed, ctx.q_enc, ctx.q_breaks, stream, |p, fwd, rc| {
        let p64 = p as u64;
        let mut was_capped = false;
        if ctx.forward {
            let key = seed.key(fwd);
            if ctx.index.count(key) > ctx.max_occ {
                was_capped = true;
            } else {
                for &s in ctx.index.bucket_in(key, &ctx.s_range) {
                    let s = s as u64;
                    if ctx.lower_only && s < p64 {
                        continue;
                    }
                    fwd_hits.push(((s + qlen - p64) << 32) | p64);
                }
            }
        }
        if ctx.reverse {
            let key = seed.key(rc);
            if ctx.index.count(key) > ctx.max_occ {
                was_capped = true;
            } else {
                for &s in ctx.index.bucket_in(key, &ctx.s_range) {
                    let s = s as u64;
                    // Cells of this hit have q + s = p + s + span - 1; any
                    // with q <= s need p <= s + span.
                    if ctx.lower_only && p64 > s + span {
                        continue;
                    }
                    rev_hits.push(((p64 + s + span - 1) << 32) | p64);
                }
            }
        }
        if was_capped && core.contains(&p) {
            capped.push(p as u32);
        }
    });
    stats.capped_positions = capped.len() as u64;
    stats.raw_hits = (fwd_hits.len() + rev_hits.len()) as u64;

    let mut runs = Vec::new();
    for (reverse, mut hits) in [(false, fwd_hits), (true, rev_hits)] {
        if hits.is_empty() {
            continue;
        }
        hits.sort_unstable();
        stats.anchors += two_hit_runs(ctx, &hits, reverse, &core, &mut runs);
    }
    ChunkHits {
        runs,
        capped,
        stats,
    }
}

#[inline]
fn diag_of(h: u64) -> u64 {
    h >> 32
}
#[inline]
fn pos_of(h: u64) -> u32 {
    h as u32
}

#[inline]
fn record_of(breaks: &[usize], pos: usize) -> u32 {
    breaks.partition_point(|&b| b <= pos) as u32
}

/// Apply the two-hit test to sorted `hits` and append merged runs; returns
/// the number of anchors kept.
fn two_hit_runs(
    ctx: &HitCtx<'_>,
    hits: &[u64],
    reverse: bool,
    core: &Range<usize>,
    out: &mut Vec<Run>,
) -> u64 {
    let dist = ctx.two_hit_dist;
    let span = ctx.index.seed().span() as u64;
    let half = span / 2;
    let qlen = ctx.q_enc.len() as u64;
    let near = |a: u64, b: u64| diag_of(a) == diag_of(b) && pos_of(b) - pos_of(a) <= dist;
    let subject_of = |h: u64| -> usize {
        let (d, p) = (diag_of(h), pos_of(h) as u64);
        (if reverse {
            d + 1 - span - p
        } else {
            d + p - qlen
        }) as usize
    };
    let records = |h: u64| {
        (
            record_of(ctx.q_breaks, pos_of(h) as usize),
            record_of(ctx.s_breaks, subject_of(h)),
        )
    };
    let in_core = |h: u64| core.contains(&(pos_of(h) as usize));
    // A margin hit adjacent to a run end, in the same record pair, is an
    // anchor that the neighbouring chunk merges with this run.
    let fill = |h: u64, other: Option<u64>| -> Option<u32> {
        let o = other?;
        let adjacent = if pos_of(o) < pos_of(h) {
            near(o, h)
        } else {
            near(h, o)
        };
        (adjacent && !in_core(o) && records(o) == records(h)).then_some(pos_of(o))
    };

    let mut kept = 0u64;
    let mut cur: Option<(Run, usize)> = None; // run + index of its last hit
    let mut votes = [0u32; RepeatClass::COUNT];

    let mut flush = |cur: &mut Option<(Run, usize)>, votes: &mut [u32; RepeatClass::COUNT]| {
        if let Some((mut r, last)) = cur.take() {
            // Majority class; ties go to the more repetitive class.
            let best = (0..RepeatClass::COUNT)
                .rev()
                .max_by_key(|&i| votes[i])
                .unwrap_or(0);
            r.class = RepeatClass::from_index(best).unwrap_or(RepeatClass::Unique);
            r.hi_fill = fill(hits[last], hits.get(last + 1).copied());
            out.push(r);
        }
        *votes = [0; RepeatClass::COUNT];
    };

    for i in 0..hits.len() {
        let h = hits[i];
        if !in_core(h) {
            continue;
        }
        let passes =
            (i > 0 && near(hits[i - 1], h)) || (i + 1 < hits.len() && near(h, hits[i + 1]));
        if !passes {
            continue;
        }
        kept += 1;
        let d = diag_of(h);
        let p = pos_of(h);
        let (q_rec, s_rec) = records(h);
        let qc = (p as u64 + half) as usize;
        let sc = (subject_of(h) as u64 + half) as usize;
        let depth = |d: &[u16], i: usize| d.get(i).copied().unwrap_or(0);
        let class = ctx
            .thresholds
            .classify(depth(ctx.q_depth, qc).max(depth(ctx.s_depth, sc)));

        let extend = cur.as_ref().is_some_and(|(r, _)| {
            r.diag == d && p - r.last_p <= dist && r.q_rec == q_rec && r.s_rec == s_rec
        });
        if extend {
            let (r, last) = cur.as_mut().expect("checked above");
            r.last_p = p;
            r.n_anchors += 1;
            *last = i;
        } else {
            flush(&mut cur, &mut votes);
            let lo_fill = if i > 0 {
                fill(h, Some(hits[i - 1]))
            } else {
                None
            };
            cur = Some((
                Run {
                    reverse,
                    diag: d,
                    first_p: p,
                    last_p: p,
                    n_anchors: 1,
                    q_rec,
                    s_rec,
                    class: RepeatClass::Unique,
                    lo_fill,
                    hi_fill: None,
                },
                i,
            ));
        }
        votes[class.index()] += 1;
    }
    flush(&mut cur, &mut votes);
    kept
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::{encode, reverse_complement_dna, AlphabetKind};
    use crate::large::seed::SpacedSeed;

    pub(crate) fn lcg_dna(n: usize, mut state: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(state >> 62) as usize]
            })
            .collect()
    }

    fn run_all(q: &[u8], s: &[u8], dist: u32) -> ChunkHits {
        let seed = SpacedSeed::default();
        let qe = encode(q, AlphabetKind::Dna);
        let se = encode(s, AlphabetKind::Dna);
        let index = SeedIndex::build(&seed, &se, &[]).unwrap();
        let qd = vec![0u16; q.len()];
        let sd = vec![0u16; s.len()];
        let ctx = HitCtx {
            index: &index,
            q_enc: &qe,
            q_breaks: &[],
            s_breaks: &[],
            q_depth: &qd,
            s_depth: &sd,
            thresholds: ClassThresholds::default(),
            s_range: 0..s.len(),
            max_occ: 1000,
            two_hit_dist: dist,
            forward: true,
            reverse: true,
            lower_only: false,
        };
        chunk_hits(&ctx, 0..q.len(), 0..q.len())
    }

    #[test]
    fn random_hits_are_filtered() {
        let q = lcg_dna(200_000, 11);
        let s = lcg_dna(200_000, 12);
        let out = run_all(&q, &s, 100);
        // ~2 × 4e10 / 1.7e7 ≈ 4800 raw hits expected by chance.
        assert!(out.stats.raw_hits > 1000, "{:?}", out.stats);
        assert!(
            out.stats.anchors * 50 < out.stats.raw_hits,
            "two-hit filter kept too much: {:?}",
            out.stats
        );
    }

    #[test]
    fn diverged_copy_forms_runs_on_both_strands() {
        let unit = lcg_dna(2000, 5);
        // ~20 % substitutions, deterministic.
        let mut state = 99u64;
        let mutated: Vec<u8> = unit
            .iter()
            .map(|&b| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                if (state >> 59) < 6 {
                    b"ACGT"[((state >> 40) & 3) as usize]
                } else {
                    b
                }
            })
            .collect();
        let mut q = lcg_dna(3000, 1);
        q.splice(500..500, unit.iter().copied());
        let mut s = lcg_dna(3000, 2);
        s.splice(1000..1000, mutated.iter().copied());
        s.splice(4000..4000, reverse_complement_dna(&mutated));
        let out = run_all(&q, &s, 100);
        let fwd: Vec<_> = out.runs.iter().filter(|r| !r.reverse).collect();
        let rev: Vec<_> = out.runs.iter().filter(|r| r.reverse).collect();
        assert!(!fwd.is_empty() && !rev.is_empty());
        // forward copy lies on s - q = 500
        assert!(fwd.iter().all(|r| r.diag == 500 + q.len() as u64));
        let covered: u32 = fwd.iter().map(|r| r.last_p - r.first_p + 19).sum();
        assert!(covered > 1000, "forward coverage {covered}");
        // reverse copy: q in 500..2500 pairs with s in 4000..6000 reversed
        // ⇒ q + s = 500 + 5999 = 6499.
        assert!(rev.iter().all(|r| r.diag == 6499));
    }
}
