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
        let first = segs[0];
        out.push(Chain {
            reverse: first.reverse,
            q_rec: first.q_rec,
            s_rec: first.s_rec,
            q_start: segs.iter().map(|s| s.q_start).min().unwrap_or(0),
            q_end: segs.iter().map(|s| s.q_end).max().unwrap_or(0),
            s_start: segs.iter().map(|s| s.s_start).min().unwrap_or(0),
            s_end: segs.iter().map(|s| s.s_end).max().unwrap_or(0),
            score,
            n_segments: segs.len() as u32,
            n_anchors: segs.iter().map(|s| s.n_anchors).sum(),
        });
    }
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
