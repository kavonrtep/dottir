//! Alignments between two assemblies, as read from PAF.

/// A sequence record of an assembly: name and length only. Genome mode
/// never holds the residues.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomeRecord {
    pub name: String,
    pub len: u64,
}

/// One alignment. Coordinates are record-local, 0-based, end-exclusive,
/// on the forward strand of each record (PAF convention). On the reverse
/// strand, query `q_start` pairs with target `t_end`.
#[derive(Debug, Clone, PartialEq)]
pub struct Alignment {
    /// Index into the query record list.
    pub q_rec: u32,
    pub q_start: u64,
    pub q_end: u64,
    /// Index into the target (reference) record list.
    pub t_rec: u32,
    pub t_start: u64,
    pub t_end: u64,
    pub reverse: bool,
    /// PAF column 10: residue matches.
    pub matches: u64,
    /// PAF column 11: alignment block length, gaps included.
    pub block_len: u64,
    pub mapq: u8,
    /// `tp:A:P` or no `tp` tag.
    pub primary: bool,
    /// `de:f`: gap-compressed per-base divergence, when the aligner
    /// computed base-level alignments (`-c`).
    pub divergence: Option<f32>,
}

impl Alignment {
    /// Identity in `[0, 1]`: `1 − de` when present, else matches / block
    /// length (spec §4.7.4).
    pub fn identity(&self) -> f64 {
        match self.divergence {
            Some(de) => (1.0 - f64::from(de)).clamp(0.0, 1.0),
            None if self.block_len == 0 => 0.0,
            None => (self.matches as f64 / self.block_len as f64).clamp(0.0, 1.0),
        }
    }

    /// Query span length.
    pub fn q_len(&self) -> u64 {
        self.q_end - self.q_start
    }

    /// Target span length.
    pub fn t_len(&self) -> u64 {
        self.t_end - self.t_start
    }

    /// Length used by length filters: the longer of the two spans.
    pub fn len(&self) -> u64 {
        self.q_len().max(self.t_len())
    }

    /// True for a zero-length alignment.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Target coordinate paired with `q_start` (`t_end` on the reverse
    /// strand).
    pub fn t_at_q_start(&self) -> u64 {
        if self.reverse {
            self.t_end
        } else {
            self.t_start
        }
    }

    /// The same alignment with query and target exchanged.
    pub fn swapped(&self) -> Self {
        Self {
            q_rec: self.t_rec,
            q_start: self.t_start,
            q_end: self.t_end,
            t_rec: self.q_rec,
            t_start: self.q_start,
            t_end: self.q_end,
            ..self.clone()
        }
    }

    /// Target coordinate paired with `q_end`.
    pub fn t_at_q_end(&self) -> u64 {
        if self.reverse {
            self.t_start
        } else {
            self.t_end
        }
    }
}

/// Join alignments split where the query was cut into chunks for the
/// aligner (ADR 0008 §3). Two alignments are joined when they lie on the
/// same record pair and strand, one ends within `tol` of a cut and the
/// other starts within `tol` after it, and the target gap between them is
/// within `tol` in the strand's direction. `cuts` are query-record-local
/// positions per query record (sorted). Output is sorted.
pub fn join_at_cuts(mut alns: Vec<Alignment>, cuts: &[Vec<u64>], tol: u64) -> Vec<Alignment> {
    sort_alignments(&mut alns);
    let near_cut = |rec: u32, pos: u64| -> Option<u64> {
        let c = cuts.get(rec as usize)?;
        let i = c.partition_point(|&x| x + tol < pos);
        c.get(i).copied().filter(|&x| x.abs_diff(pos) <= tol)
    };
    // Candidates ending at a cut, keyed by (q_rec, t_rec, reverse, cut).
    let mut open: std::collections::BTreeMap<(u32, u32, bool, u64), Vec<usize>> =
        std::collections::BTreeMap::new();
    for (i, a) in alns.iter().enumerate() {
        if let Some(c) = near_cut(a.q_rec, a.q_end) {
            open.entry((a.q_rec, a.t_rec, a.reverse, c))
                .or_default()
                .push(i);
        }
    }
    let mut absorbed = vec![false; alns.len()];
    let mut merged_into: Vec<Option<usize>> = vec![None; alns.len()];
    for j in 0..alns.len() {
        let b = &alns[j];
        let Some(c) = near_cut(b.q_rec, b.q_start) else {
            continue;
        };
        let Some(cands) = open.get(&(b.q_rec, b.t_rec, b.reverse, c)) else {
            continue;
        };
        // The best partner: co-linear with the smallest target gap.
        let best = cands
            .iter()
            .copied()
            .filter(|&i| i != j && merged_into[i].is_none())
            .filter(|&i| {
                let a = &alns[i];
                let gap = if b.reverse {
                    a.t_start as i128 - b.t_end as i128
                } else {
                    b.t_start as i128 - a.t_end as i128
                };
                a.q_end <= b.q_start + tol && gap >= -(tol as i128) && gap <= tol as i128
            })
            .min_by_key(|&i| {
                let a = &alns[i];
                if b.reverse {
                    a.t_start.abs_diff(b.t_end)
                } else {
                    b.t_start.abs_diff(a.t_end)
                }
            });
        if let Some(i) = best {
            merged_into[i] = Some(j);
            absorbed[j] = true;
        }
    }
    // Follow each chain head (an alignment nobody merged into) to its end.
    let mut out = Vec::with_capacity(alns.len());
    for i in 0..alns.len() {
        if absorbed[i] {
            continue;
        }
        let mut acc = alns[i].clone();
        let mut k = i;
        while let Some(j) = merged_into[k] {
            let b = &alns[j];
            acc.q_end = acc.q_end.max(b.q_end);
            acc.t_start = acc.t_start.min(b.t_start);
            acc.t_end = acc.t_end.max(b.t_end);
            acc.matches += b.matches;
            acc.block_len += b.block_len;
            acc.mapq = acc.mapq.min(b.mapq);
            acc.primary &= b.primary;
            acc.divergence = match (acc.divergence, b.divergence) {
                (Some(x), Some(y)) => Some(x.max(y)),
                _ => None,
            };
            k = j;
        }
        out.push(acc);
    }
    sort_alignments(&mut out);
    out
}

/// Canonical order: by query record, query start, target record, target
/// start, strand, then the remaining fields. Rendering never depends on
/// it, but exports and tests do.
pub fn sort_alignments(alns: &mut [Alignment]) {
    alns.sort_by(|a, b| {
        (
            a.q_rec, a.q_start, a.t_rec, a.t_start, a.reverse, a.q_end, a.t_end,
        )
            .cmp(&(
                b.q_rec, b.q_start, b.t_rec, b.t_start, b.reverse, b.q_end, b.t_end,
            ))
            .then(a.matches.cmp(&b.matches))
            .then(a.block_len.cmp(&b.block_len))
    });
}

/// Overlap fraction (of the shorter query span) above which a lower-
/// scoring alignment is secondary to a primary, as minimap2's default
/// `--mask-level`.
pub const MASK_LEVEL: f64 = 0.5;

/// Mark primaries the way minimap2 does without `-P`, for PAFs that have
/// none (`-P` reports every chain as `tp:A:S`): in descending score
/// (matches, then block length, then position), an alignment is primary
/// unless at least `mask_level` of the shorter query span overlaps an
/// earlier primary of the same query record.
pub fn assign_primary(alns: &mut [Alignment], mask_level: f64) {
    use std::collections::BTreeMap;
    let mut by_score: Vec<usize> = (0..alns.len()).collect();
    by_score.sort_by(|&i, &j| {
        let (a, b) = (&alns[i], &alns[j]);
        (b.matches, b.block_len)
            .cmp(&(a.matches, a.block_len))
            .then(
                (
                    a.q_rec, a.q_start, a.q_end, a.t_rec, a.t_start, a.t_end, a.reverse,
                )
                    .cmp(&(
                        b.q_rec, b.q_start, b.q_end, b.t_rec, b.t_start, b.t_end, b.reverse,
                    )),
            )
    });
    // Per query record: primaries keyed by (start, index) → end, and the
    // longest primary, which bounds how far left an overlap can start.
    type Primaries = (BTreeMap<(u64, usize), u64>, u64);
    let mut primaries: BTreeMap<u32, Primaries> = BTreeMap::new();
    for i in by_score {
        let a = &alns[i];
        let (set, longest) = primaries.entry(a.q_rec).or_default();
        let len = a.q_len();
        let secondary = set
            .range(..(a.q_end, 0))
            .rev()
            .take_while(|&(&(start, _), _)| start + *longest > a.q_start)
            .any(|(&(start, _), &end)| {
                let overlap = end.min(a.q_end).saturating_sub(start.max(a.q_start));
                overlap > 0 && overlap as f64 >= mask_level * len.min(end - start) as f64
            });
        if !secondary {
            set.insert((a.q_start, i), a.q_end);
            *longest = (*longest).max(len);
        }
        alns[i].primary = !secondary;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primaries_follow_minimap2_mask_level() {
        // A long best hit, a repeat copy inside it, a hit half inside it,
        // one barely overlapping, and a hit on another stretch.
        let mut alns = vec![
            aln((1000, 2000), (0, 100), false),
            aln((0, 10_000), (0, 10_000), false),
            aln((9_500, 10_500), (50_000, 51_000), true),
            aln((9_900, 12_000), (60_000, 62_100), false),
            aln((20_000, 21_000), (70_000, 71_000), false),
        ];
        for a in &mut alns {
            a.primary = false;
        }
        assign_primary(&mut alns, MASK_LEVEL);
        let p: Vec<bool> = alns.iter().map(|a| a.primary).collect();
        assert_eq!(p, [false, true, false, true, true]);
        // Input order does not matter.
        let mut rev: Vec<Alignment> = alns.iter().rev().cloned().collect();
        assign_primary(&mut rev, MASK_LEVEL);
        assert_eq!(rev.iter().rev().map(|a| a.primary).collect::<Vec<_>>(), p);
    }

    fn aln(q: (u64, u64), t: (u64, u64), reverse: bool) -> Alignment {
        Alignment {
            q_rec: 0,
            q_start: q.0,
            q_end: q.1,
            t_rec: 0,
            t_start: t.0,
            t_end: t.1,
            reverse,
            matches: (q.1 - q.0) * 9 / 10,
            block_len: q.1 - q.0,
            mapq: 60,
            primary: true,
            divergence: None,
        }
    }

    #[test]
    fn identity_prefers_divergence() {
        let mut a = aln((0, 100), (0, 100), false);
        assert!((a.identity() - 0.9).abs() < 1e-9);
        a.divergence = Some(0.02);
        assert!((a.identity() - 0.98).abs() < 1e-6);
    }

    #[test]
    fn forward_pieces_join_across_a_cut() {
        let cuts = vec![vec![1000]];
        let alns = vec![
            aln((1000, 1800), (5998, 6800), false),
            aln((200, 997), (5200, 5997), false),
        ];
        let out = join_at_cuts(alns, &cuts, 20);
        assert_eq!(out.len(), 1);
        assert_eq!((out[0].q_start, out[0].q_end), (200, 1800));
        assert_eq!((out[0].t_start, out[0].t_end), (5200, 6800));
        assert_eq!(out[0].block_len, 1597);
    }

    #[test]
    fn reverse_pieces_join_across_a_cut() {
        // Query 200..1000 pairs with target 6000..5200 (descending).
        let cuts = vec![vec![1000]];
        let alns = vec![
            aln((200, 1000), (5200, 6000), true),
            aln((1000, 1800), (4400, 5200), true),
        ];
        let out = join_at_cuts(alns, &cuts, 20);
        assert_eq!(out.len(), 1);
        assert_eq!((out[0].t_start, out[0].t_end), (4400, 6000));
    }

    #[test]
    fn pieces_not_colinear_or_not_at_a_cut_stay_apart() {
        let cuts = vec![vec![1000]];
        // Target jumps 50 kb: a rearrangement at the cut, not a split.
        let jump = vec![
            aln((200, 1000), (5200, 6000), false),
            aln((1000, 1800), (56000, 56800), false),
        ];
        assert_eq!(join_at_cuts(jump, &cuts, 20).len(), 2);
        // Adjacent but far from any cut.
        let far = vec![
            aln((2000, 3000), (0, 1000), false),
            aln((3000, 4000), (1000, 2000), false),
        ];
        assert_eq!(join_at_cuts(far, &cuts, 20).len(), 2);
        // Opposite strands.
        let strands = vec![
            aln((200, 1000), (5200, 6000), false),
            aln((1000, 1800), (6000, 6800), true),
        ];
        assert_eq!(join_at_cuts(strands, &cuts, 20).len(), 2);
    }

    #[test]
    fn a_three_piece_alignment_over_two_cuts_joins_fully() {
        let cuts = vec![vec![1000, 2000]];
        let alns = vec![
            aln((2000, 2500), (2000, 2500), false),
            aln((500, 1000), (500, 1000), false),
            aln((1000, 2000), (1000, 2000), false),
        ];
        let out = join_at_cuts(alns, &cuts, 20);
        assert_eq!(out.len(), 1);
        assert_eq!((out[0].q_start, out[0].q_end), (500, 2500));
    }
}
