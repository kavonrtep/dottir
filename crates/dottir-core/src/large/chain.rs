//! Basic co-linear chaining of diagonal runs into synteny segments.
//!
//! A display aid, not an aligner: chains carry coordinates and a score,
//! never a base-level alignment. The DP is the minimap2 (Li 2018) shape,
//! simplified: `f(i) = w_i + max(0, max_j f(j) - β(j, i))` over a bounded
//! look-back, with `β` penalising diagonal drift and gap length.

use super::depth::RepeatClass;
use super::Segment;

#[derive(Debug, Clone, PartialEq)]
pub struct ChainConfig {
    /// Segments more repetitive than this are not chained.
    pub max_class: RepeatClass,
    /// Largest query or subject gap bridged between two segments.
    pub max_gap: u32,
    /// Predecessors examined per segment.
    pub lookback: usize,
    /// Allowed overlap between consecutive segments, per axis.
    pub overlap_tol: u32,
    /// Penalty per base of `|dq - ds|`.
    pub drift_penalty: f64,
    /// Penalty per base of `min(dq, ds)`.
    pub gap_penalty: f64,
    /// Chains scoring below this are dropped.
    pub min_score: f64,
    /// Chains with fewer two-hit anchors are dropped. With seed sampling
    /// the two-hit distance grows, and a lone pair of chance hits can then
    /// span enough bases to pass `min_score`.
    pub min_anchors: u32,
    /// A chain is dropped when at least this share of both its query and
    /// its subject span lies inside a better chain of the same record pair
    /// and strand (minimap2's `--mask-level`). Side chains along small
    /// indels of a synteny block would otherwise be reported separately.
    pub mask_level: f64,
}

impl Default for ChainConfig {
    fn default() -> Self {
        Self {
            max_class: RepeatClass::Low,
            max_gap: 10_000,
            lookback: 50,
            overlap_tol: 50,
            drift_penalty: 0.1,
            gap_penalty: 0.01,
            min_score: 300.0,
            min_anchors: 3,
            mask_level: 0.5,
        }
    }
}

/// A chain of segments on one strand of one record pair. Coordinates are
/// 0-based, half-open, in the concatenated sequences.
#[derive(Debug, Clone, PartialEq)]
pub struct Chain {
    pub reverse: bool,
    pub q_rec: u32,
    pub s_rec: u32,
    pub q_start: u32,
    pub q_end: u32,
    pub s_start: u32,
    pub s_end: u32,
    pub score: f64,
    pub n_segments: u32,
    pub n_anchors: u32,
}

/// Chain `segments`. Output order: record pair, strand, then query start.
pub fn chain_segments(segments: &[Segment], cfg: &ChainConfig) -> Vec<Chain> {
    let mut idx: Vec<usize> = (0..segments.len())
        .filter(|&i| segments[i].class <= cfg.max_class)
        .collect();
    let key = |s: &Segment| {
        (
            s.q_rec, s.s_rec, s.reverse, s.q_start, s.s_start, s.q_end, s.s_end,
        )
    };
    idx.sort_by_key(|&i| key(&segments[i]));

    let mut chains = Vec::new();
    let mut a = 0;
    while a < idx.len() {
        let g = &segments[idx[a]];
        let mut b = a + 1;
        while b < idx.len() && {
            let h = &segments[idx[b]];
            (h.q_rec, h.s_rec, h.reverse) == (g.q_rec, g.s_rec, g.reverse)
        } {
            b += 1;
        }
        let group: Vec<&Segment> = idx[a..b].iter().map(|&i| &segments[i]).collect();
        chain_group(&group, cfg, &mut chains);
        a = b;
    }
    chains.sort_by_key(|c| (c.q_rec, c.s_rec, c.reverse, c.q_start, c.s_start));
    chains
}

/// Gap `(dq, ds)` from `j` to `i` along the strand direction; negative
/// values are overlaps.
#[inline]
fn gaps(j: &Segment, i: &Segment) -> (i64, i64) {
    let dq = i.q_start as i64 - j.q_end as i64;
    let ds = if i.reverse {
        j.s_start as i64 - i.s_end as i64
    } else {
        i.s_start as i64 - j.s_end as i64
    };
    (dq, ds)
}

fn chain_group(group: &[&Segment], cfg: &ChainConfig, out: &mut Vec<Chain>) {
    let n = group.len();
    let weight = |s: &Segment| (s.q_end - s.q_start) as f64;
    let mut f = vec![0f64; n];
    let mut pred: Vec<Option<usize>> = vec![None; n];
    let tol = cfg.overlap_tol as i64;
    let max_gap = cfg.max_gap as i64;

    for i in 0..n {
        let si = group[i];
        let wi = weight(si);
        f[i] = wi;
        for j in (i.saturating_sub(cfg.lookback)..i).rev() {
            let sj = group[j];
            let (dq, ds) = gaps(sj, si);
            if dq < -tol || ds < -tol || dq > max_gap || ds > max_gap {
                continue;
            }
            // Strictly co-linear progress on both axes.
            let progresses = si.q_end > sj.q_end
                && if si.reverse {
                    si.s_start < sj.s_start
                } else {
                    si.s_end > sj.s_end
                };
            if !progresses {
                continue;
            }
            let overlap = (-dq.min(ds)).max(0) as f64;
            let drift = (dq - ds).abs() as f64;
            let gap = dq.min(ds).max(0) as f64;
            let cand = f[j] + wi - overlap - cfg.drift_penalty * drift - cfg.gap_penalty * gap;
            if cand > f[i] {
                f[i] = cand;
                pred[i] = Some(j);
            }
        }
    }

    // Greedy backtrack from the best end points; each segment is used once.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&x, &y| f[y].total_cmp(&f[x]).then(x.cmp(&y)));
    let mut used = vec![false; n];
    let mut kept: Vec<(u32, u32, u32, u32)> = Vec::new();
    for &end in &order {
        if used[end] {
            continue;
        }
        let mut members = vec![];
        let mut cur = Some(end);
        let mut base = 0.0;
        while let Some(c) = cur {
            if used[c] {
                base = f[c];
                break;
            }
            used[c] = true;
            members.push(c);
            cur = pred[c];
        }
        let score = f[end] - base;
        if score < cfg.min_score {
            continue;
        }
        let segs: Vec<&Segment> = members.iter().map(|&m| group[m]).collect();
        let n_anchors: u32 = segs.iter().map(|s| s.n_anchors).sum();
        if n_anchors < cfg.min_anchors {
            continue;
        }
        let first = segs[0];
        let q_start = segs.iter().map(|s| s.q_start).min().unwrap_or(0);
        let q_end = segs.iter().map(|s| s.q_end).max().unwrap_or(0);
        let s_start = segs.iter().map(|s| s.s_start).min().unwrap_or(0);
        let s_end = segs.iter().map(|s| s.s_end).max().unwrap_or(0);
        // Chains come out best first, so `kept` holds only better ones.
        let masked = kept.iter().any(|&(q0, q1, s0, s1)| {
            covered(q_start, q_end, q0, q1) >= cfg.mask_level
                && covered(s_start, s_end, s0, s1) >= cfg.mask_level
        });
        if masked {
            continue;
        }
        kept.push((q_start, q_end, s_start, s_end));
        out.push(Chain {
            reverse: first.reverse,
            q_rec: first.q_rec,
            s_rec: first.s_rec,
            q_start,
            q_end,
            s_start,
            s_end,
            score,
            n_segments: segs.len() as u32,
            n_anchors,
        });
    }
}

/// Share of `[a0, a1)` inside `[b0, b1)`.
fn covered(a0: u32, a1: u32, b0: u32, b1: u32) -> f64 {
    let inter = a1.min(b1).saturating_sub(a0.max(b0));
    inter as f64 / (a1 - a0).max(1) as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seg(reverse: bool, q: u32, s: u32, len: u32) -> Segment {
        let (s_start, s_end) = (s, s + len);
        Segment {
            reverse,
            q_rec: 0,
            s_rec: 0,
            q_start: q,
            q_end: q + len,
            s_start,
            s_end,
            n_anchors: len / 20,
            class: RepeatClass::Unique,
        }
    }

    #[test]
    fn chains_colinear_forward_and_skips_outlier() {
        let segs = vec![
            seg(false, 0, 1000, 400),
            seg(false, 500, 1510, 400),
            seg(false, 1000, 2000, 400),
            seg(false, 600, 9000, 100), // off-diagonal, too short alone
        ];
        let chains = chain_segments(&segs, &ChainConfig::default());
        assert_eq!(chains.len(), 1);
        let c = &chains[0];
        assert_eq!(
            (c.q_start, c.q_end, c.s_start, c.s_end),
            (0, 1400, 1000, 2400)
        );
        assert_eq!(c.n_segments, 3);
    }

    #[test]
    fn chains_reverse_strand_with_decreasing_subject() {
        let segs = vec![
            seg(true, 0, 5000, 400),
            seg(true, 450, 4500, 400),
            seg(true, 900, 4000, 400),
        ];
        let chains = chain_segments(&segs, &ChainConfig::default());
        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].n_segments, 3);
        assert!(chains[0].reverse);
    }

    #[test]
    fn repeats_and_record_pairs_are_kept_apart() {
        let mut a = seg(false, 0, 1000, 400);
        let mut b = seg(false, 500, 1500, 400);
        b.s_rec = 1;
        let mut c = seg(false, 1000, 2000, 400);
        c.class = RepeatClass::High;
        a.n_anchors = 10;
        let chains = chain_segments(&[a, b, c], &ChainConfig::default());
        assert_eq!(chains.len(), 2);
        assert!(chains.iter().all(|c| c.n_segments == 1));
    }
}

#[cfg(test)]
mod mask_tests {
    use super::*;

    fn seg(q: u32, s: u32, len: u32) -> Segment {
        Segment {
            reverse: false,
            q_rec: 0,
            s_rec: 0,
            q_start: q,
            q_end: q + len,
            s_start: s,
            s_end: s + len,
            n_anchors: len / 20,
            class: RepeatClass::Unique,
        }
    }

    #[test]
    fn side_chains_inside_a_block_are_masked() {
        // A 3-segment block, plus a 700 bp run on a diagonal shifted by a
        // small indel that overlaps the second segment too much to join it.
        let segs = vec![
            seg(0, 0, 2000),
            seg(2050, 2060, 2000),
            seg(4100, 4115, 2000),
            seg(2400, 2390, 700),
        ];
        let chains = chain_segments(&segs, &ChainConfig::default());
        assert_eq!(chains.len(), 1, "{chains:?}");
        assert_eq!(chains[0].n_segments, 3);
        // Without masking the side run is its own chain.
        let cfg = ChainConfig {
            mask_level: 1.01,
            ..ChainConfig::default()
        };
        assert_eq!(chain_segments(&segs, &cfg).len(), 2);
    }

    #[test]
    fn two_anchor_chains_are_dropped() {
        let mut lone = seg(0, 0, 400);
        lone.n_anchors = 2;
        assert!(chain_segments(&[lone.clone()], &ChainConfig::default()).is_empty());
        lone.n_anchors = 3;
        assert_eq!(chain_segments(&[lone], &ChainConfig::default()).len(), 1);
    }

    #[test]
    fn duplicated_blocks_elsewhere_are_kept() {
        // Same query block matching two subject places (a segdup): only
        // the query overlaps, so both chains survive.
        let segs = vec![seg(0, 0, 3000), seg(0, 50_000, 3000)];
        assert_eq!(chain_segments(&segs, &ChainConfig::default()).len(), 2);
    }
}
