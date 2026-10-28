//! Query record ordering and orientation along the reference (spec
//! §4.7.6), following D-GENIES `Paf.sort` and
//! `Paf.is_contig_well_oriented`.

use super::alignment::{Alignment, GenomeRecord};
use super::filter::line_len;
use super::layout::AxisLayout;

/// Where a query record was placed.
#[derive(Debug, Clone, PartialEq)]
pub struct Placement {
    /// Target record with the highest gravity.
    pub t_rec: u32,
    /// Gravity-weighted mean target midpoint, in target axis coordinates.
    pub position: f64,
    /// Draw the record reverse-complemented.
    pub reverse: bool,
    /// Extent of the alignments to `t_rec`, record-local (forward
    /// coordinates of the query record).
    pub q_span: (u64, u64),
    pub t_span: (u64, u64),
}

/// The result of sorting: a record order for [`AxisLayout::new`], flips,
/// and per-record placements (`None` = no alignment).
#[derive(Debug, Clone, PartialEq)]
pub struct QueryOrder {
    pub order: Vec<u32>,
    pub flips: Vec<bool>,
    pub placements: Vec<Option<Placement>>,
}

struct Line {
    t_rec: u32,
    /// Target midpoint (axis coordinates).
    mid_t: f64,
    /// Query midpoint (record-local).
    mid_q: f64,
    len: f64,
    reverse: bool,
    q: (u64, u64),
    t: (u64, u64),
}

/// Sort query records along the target. Each alignment weighs
/// `(1 + line length)²` ("gravity"); a query record goes to the target
/// record with the largest total gravity, at the gravity-weighted mean
/// of its alignments' target midpoints there. Records without
/// alignments go last. Ties are broken by record name.
///
/// Orientation, from the best block's long alignments (longer than 10 %
/// of the longest and at least 1 % of the shorter record): with two or
/// more, the query midpoints should rise along the target (mean step
/// direction > −0.1); with one, it should be on the forward strand.
pub fn sort_query(
    alns: &[&Alignment],
    query: &[GenomeRecord],
    target: &[GenomeRecord],
    target_layout: &AxisLayout,
) -> QueryOrder {
    let mut by_query: Vec<Vec<Line>> = (0..query.len()).map(|_| Vec::new()).collect();
    for a in alns {
        let (x0, x1) = (
            target_layout.to_axis(a.t_rec, a.t_start),
            target_layout.to_axis(a.t_rec, a.t_end),
        );
        by_query[a.q_rec as usize].push(Line {
            t_rec: a.t_rec,
            mid_t: (x0.min(x1) as f64 + x0.max(x1) as f64) / 2.0,
            mid_q: (a.q_start as f64 + a.q_end as f64) / 2.0,
            len: line_len(a),
            reverse: a.reverse,
            q: (a.q_start, a.q_end),
            t: (a.t_start, a.t_end),
        });
    }
    let placements: Vec<Option<Placement>> = by_query
        .iter()
        .enumerate()
        .map(|(q, lines)| place(lines, query[q].len, target, target_layout))
        .collect();
    let mut order: Vec<u32> = (0..query.len() as u32).collect();
    order.sort_by(|&a, &b| {
        let pos = |r: u32| {
            placements[r as usize]
                .as_ref()
                .map_or(f64::INFINITY, |p| p.position)
        };
        pos(a)
            .total_cmp(&pos(b))
            .then_with(|| query[a as usize].name.cmp(&query[b as usize].name))
    });
    let flips = placements
        .iter()
        .map(|p| p.as_ref().is_some_and(|p| p.reverse))
        .collect();
    QueryOrder {
        order,
        flips,
        placements,
    }
}

fn place(
    lines: &[Line],
    q_len: u64,
    target: &[GenomeRecord],
    target_layout: &AxisLayout,
) -> Option<Placement> {
    if lines.is_empty() {
        return None;
    }
    let gravity = |l: &Line| (1.0 + l.len).powi(2);
    let mut per_target: std::collections::BTreeMap<u32, f64> = Default::default();
    for l in lines {
        *per_target.entry(l.t_rec).or_default() += gravity(l);
    }
    // Highest gravity; ties go to the target drawn first.
    let (&t_rec, _) = per_target.iter().max_by(|(ta, ga), (tb, gb)| {
        ga.total_cmp(gb).then_with(|| {
            target_layout
                .segment(**tb)
                .start
                .cmp(&target_layout.segment(**ta).start)
        })
    })?;
    let mut block: Vec<&Line> = lines.iter().filter(|l| l.t_rec == t_rec).collect();
    let (sum, weight) = block.iter().fold((0.0, 0.0), |(s, w), l| {
        (s + l.mid_t * gravity(l), w + gravity(l))
    });
    let position = sum / weight;
    let span = |f: fn(&Line) -> (u64, u64)| {
        block.iter().fold((u64::MAX, 0), |(lo, hi), l| {
            (lo.min(f(l).0), hi.max(f(l).1))
        })
    };
    let (q_span, t_span) = (span(|l| l.q), span(|l| l.t));

    block.sort_by(|a, b| b.len.total_cmp(&a.len).then(a.mid_t.total_cmp(&b.mid_t)));
    let max_len = block[0].len;
    let min_rec = q_len.min(target[t_rec as usize].len) as f64;
    let n = block
        .iter()
        .take_while(|l| l.len > 0.1 * max_len && l.len >= 0.01 * min_rec)
        .count();
    let mut selected = block[..n].to_vec();
    let well_oriented = match selected.len() {
        0 => true,
        1 => !selected[0].reverse,
        _ => {
            selected.sort_by(|a, b| a.mid_t.total_cmp(&b.mid_t));
            let steps: Vec<f64> = selected
                .windows(2)
                .map(|w| if w[1].mid_q > w[0].mid_q { 1.0 } else { -1.0 })
                .collect();
            steps.iter().sum::<f64>() / steps.len() as f64 > -0.1
        }
    };
    Some(Placement {
        t_rec,
        position,
        reverse: !well_oriented,
        q_span,
        t_span,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(name: &str, len: u64) -> GenomeRecord {
        GenomeRecord {
            name: name.into(),
            len,
        }
    }

    fn aln(q_rec: u32, q: (u64, u64), t_rec: u32, t: (u64, u64), reverse: bool) -> Alignment {
        Alignment {
            q_rec,
            q_start: q.0,
            q_end: q.1,
            t_rec,
            t_start: t.0,
            t_end: t.1,
            reverse,
            matches: q.1 - q.0,
            block_len: q.1 - q.0,
            mapq: 60,
            primary: true,
            divergence: None,
        }
    }

    #[test]
    fn contigs_are_placed_along_the_target_and_reoriented() {
        let target = vec![rec("chr1", 10_000), rec("chr2", 10_000)];
        let tl = AxisLayout::new(&target, None, &[], 0.0);
        // c0 → chr2 start; c1 → chr1 end, reversed (two long lines whose
        // query midpoints fall along the target); c2 → chr1 start, with
        // a weak hit on chr2; c3 unaligned.
        let query = vec![
            rec("c0", 2000),
            rec("c1", 3000),
            rec("c2", 2000),
            rec("c3", 500),
        ];
        let alns = [
            aln(0, (0, 2000), 1, (0, 2000), false),
            aln(1, (0, 1500), 0, (8500, 10_000), true),
            aln(1, (1500, 3000), 0, (7000, 8500), true),
            aln(2, (0, 2000), 0, (0, 2000), false),
            aln(2, (0, 300), 1, (5000, 5300), false),
        ];
        let refs: Vec<&Alignment> = alns.iter().collect();
        let o = sort_query(&refs, &query, &target, &tl);
        assert_eq!(o.order, vec![2, 1, 0, 3]);
        assert_eq!(o.flips, vec![false, true, false, false]);
        let p = o.placements[2].as_ref().unwrap();
        assert_eq!((p.t_rec, p.q_span, p.t_span), (0, (0, 2000), (0, 2000)));
        assert_eq!(o.placements[1].as_ref().unwrap().t_span, (7000, 10_000));
        assert!(o.placements[3].is_none());
    }

    #[test]
    fn a_single_reverse_line_flips_and_ties_sort_by_name() {
        let target = vec![rec("chr1", 10_000)];
        let tl = AxisLayout::new(&target, None, &[], 0.0);
        let query = vec![rec("b", 1000), rec("a", 1000)];
        let alns = [
            aln(0, (0, 1000), 0, (0, 1000), true),
            aln(1, (0, 1000), 0, (0, 1000), false),
        ];
        let refs: Vec<&Alignment> = alns.iter().collect();
        let o = sort_query(&refs, &query, &target, &tl);
        assert_eq!(o.order, vec![1, 0]);
        assert_eq!(o.flips, vec![true, false]);
    }
}
