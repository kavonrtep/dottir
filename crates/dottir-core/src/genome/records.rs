//! Record sorting for ordinary dotplots: the D-GENIES order of genome
//! mode ([`sort_query`]), computed from large-engine chains instead of
//! aligner output.

use std::ops::Range;

use super::alignment::{Alignment, GenomeRecord};
use super::layout::AxisLayout;
use super::order::{sort_query, QueryOrder};
use crate::large::Chain;

/// The dotplot axis whose records are reordered.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SortAxis {
    /// Subject records, along the query.
    Vertical,
    /// Query records, along the subject.
    Horizontal,
}

/// A record of a concatenated dotplot input: name and buffer range.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InputRecord {
    pub name: String,
    pub range: Range<usize>,
}

/// Order and orient the records of `axis` along the other axis from
/// `chains` (query = horizontal, subject = vertical, as in the large
/// engine). Chains stand in for alignments: only their extents and
/// strands are used, as D-GENIES uses only alignment lengths and
/// positions to sort.
pub fn sort_by_chains(
    chains: &[Chain],
    horizontal: &[InputRecord],
    vertical: &[InputRecord],
    axis: SortAxis,
) -> QueryOrder {
    let recs = |r: &[InputRecord]| -> Vec<GenomeRecord> {
        r.iter()
            .map(|r| GenomeRecord {
                name: r.name.clone(),
                len: r.range.len() as u64,
            })
            .collect()
    };
    let (h, v) = (recs(horizontal), recs(vertical));
    let alns: Vec<Alignment> = chains
        .iter()
        .filter_map(|c| {
            let hr = horizontal.get(c.q_rec as usize)?;
            let vr = vertical.get(c.s_rec as usize)?;
            let (ts, te) = (
                c.q_start as usize - hr.range.start,
                c.q_end as usize - hr.range.start,
            );
            let (qs, qe) = (
                c.s_start as usize - vr.range.start,
                c.s_end as usize - vr.range.start,
            );
            let (ql, tl) = ((qe - qs) as u64, (te - ts) as u64);
            let a = Alignment {
                q_rec: c.s_rec,
                q_start: qs as u64,
                q_end: qe as u64,
                t_rec: c.q_rec,
                t_start: ts as u64,
                t_end: te as u64,
                reverse: c.reverse,
                matches: ql.min(tl),
                block_len: ql.max(tl),
                mapq: 255,
                primary: true,
                divergence: None,
            };
            Some(match axis {
                SortAxis::Vertical => a,
                SortAxis::Horizontal => a.swapped(),
            })
        })
        .collect();
    let refs: Vec<&Alignment> = alns.iter().collect();
    let (query, target) = match axis {
        SortAxis::Vertical => (&v, &h),
        SortAxis::Horizontal => (&h, &v),
    };
    let layout = AxisLayout::new(target, None, &[], 0.0);
    sort_query(&refs, query, target, &layout)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(name: &str, range: Range<usize>) -> InputRecord {
        InputRecord {
            name: name.into(),
            range,
        }
    }

    fn chain(q_rec: u32, s_rec: u32, q: (u32, u32), s: (u32, u32), reverse: bool) -> Chain {
        Chain {
            reverse,
            q_rec,
            s_rec,
            q_start: q.0,
            q_end: q.1,
            s_start: s.0,
            s_end: s.1,
            score: 1000.0,
            n_segments: 1,
            n_anchors: 10,
        }
    }

    #[test]
    fn subject_contigs_follow_the_query_and_vice_versa() {
        // Query: one 10 kb chromosome. Subject: c0 (3 kb) = chr 6–9 kb
        // reversed, c1 (4 kb) = chr 0–4 kb, c2 unaligned.
        let h = [rec("chr", 0..10_000)];
        let v = [
            rec("c0", 0..3000),
            rec("c1", 3000..7000),
            rec("c2", 7000..8000),
        ];
        let chains = [
            chain(0, 0, (6000, 9000), (0, 3000), true),
            chain(0, 1, (0, 4000), (3000, 7000), false),
        ];
        let o = sort_by_chains(&chains, &h, &v, SortAxis::Vertical);
        assert_eq!(o.order, [1, 0, 2]);
        assert_eq!(o.flips, [true, false, false]);
        assert!(o.placements[2].is_none());

        // The other way: chromosomes placed along one long contig.
        let h = [rec("a", 0..4000), rec("b", 4000..7000)];
        let v = [rec("ctg", 0..10_000)];
        let chains = [
            chain(0, 0, (0, 4000), (6000, 10_000), false),
            chain(1, 0, (4000, 7000), (0, 3000), true),
        ];
        let o = sort_by_chains(&chains, &h, &v, SortAxis::Horizontal);
        assert_eq!(o.order, [1, 0]);
        assert_eq!(o.flips, [false, true]);
    }
}
