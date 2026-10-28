//! Filters, sorting and axis layouts for one genome plot, shared by the
//! batch CLI and the GUI.

use super::alignment::{Alignment, GenomeRecord};
use super::filter::{keep_mask, FilterConfig};
use super::layout::{AxisLayout, DEFAULT_SMALL_FRACTION};
use super::order::{sort_query, Placement};
use super::records::SortAxis;

/// What to draw and where.
#[derive(Debug, Clone, PartialEq)]
pub struct Arranged {
    /// Per alignment: passes the filters.
    pub keep: Vec<bool>,
    /// Noise length threshold, when the noise filter applied.
    pub noise: Option<f64>,
    /// Horizontal axis (reference, alignment targets): file order unless
    /// sorted along the query.
    pub x: AxisLayout,
    /// Vertical axis (alignment queries).
    pub y: AxisLayout,
    /// Per vertical record: its best reference record and extents.
    pub placements: Vec<Option<Placement>>,
}

impl Arranged {
    /// The alignments that pass the filters.
    pub fn drawn<'a>(&'a self, alns: &'a [Alignment]) -> impl Iterator<Item = &'a Alignment> {
        alns.iter()
            .zip(&self.keep)
            .filter(|(_, &k)| k)
            .map(|(a, _)| a)
    }
}

/// Filter `alns` and lay out both axes. `sort` orders and orients one
/// axis's records along the other as D-GENIES does, from the
/// noise-filtered alignments (spec §4.7.6): [`SortAxis::Vertical`] the
/// query records along the reference, [`SortAxis::Horizontal`] the
/// reference records along the query. `flips` (per vertical record, may
/// be empty) reverses records on top of that.
pub fn arrange(
    alns: &[Alignment],
    x_records: &[GenomeRecord],
    y_records: &[GenomeRecord],
    filter: &FilterConfig,
    sort: Option<SortAxis>,
    flips: &[bool],
) -> Arranged {
    let (keep, noise) = keep_mask(alns, filter);
    let flipped = |base: &[bool]| -> Vec<bool> {
        (0..y_records.len())
            .map(|r| base.get(r).copied().unwrap_or(false) ^ flips.get(r).copied().unwrap_or(false))
            .collect()
    };
    let sort_keep = || {
        keep_mask(
            alns,
            &FilterConfig {
                remove_noise: true,
                ..filter.clone()
            },
        )
        .0
    };
    let file_order_x = || AxisLayout::new(x_records, None, &[], DEFAULT_SMALL_FRACTION);
    let (x, y, placements) = match sort {
        Some(SortAxis::Vertical) => {
            let x = file_order_x();
            let sort_keep = sort_keep();
            let used: Vec<&Alignment> = select(alns, &sort_keep).collect();
            let o = sort_query(&used, y_records, x_records, &x);
            let y = AxisLayout::new(
                y_records,
                Some(&o.order),
                &flipped(&o.flips),
                DEFAULT_SMALL_FRACTION,
            );
            (x, y, o.placements)
        }
        Some(SortAxis::Horizontal) => {
            let y = AxisLayout::new(y_records, None, &flipped(&[]), DEFAULT_SMALL_FRACTION);
            let sort_keep = sort_keep();
            let swapped: Vec<Alignment> =
                select(alns, &sort_keep).map(Alignment::swapped).collect();
            let used: Vec<&Alignment> = swapped.iter().collect();
            let o = sort_query(&used, x_records, y_records, &y);
            let x = AxisLayout::new(x_records, Some(&o.order), &o.flips, DEFAULT_SMALL_FRACTION);
            let used: Vec<&Alignment> = select(alns, &keep).collect();
            let placements = sort_query(&used, y_records, x_records, &x).placements;
            (x, y, placements)
        }
        None => {
            let x = file_order_x();
            let used: Vec<&Alignment> = select(alns, &keep).collect();
            let placements = sort_query(&used, y_records, x_records, &x).placements;
            let y = AxisLayout::new(y_records, None, &flipped(&[]), DEFAULT_SMALL_FRACTION);
            (x, y, placements)
        }
    };
    Arranged {
        keep,
        noise,
        x,
        y,
        placements,
    }
}

fn select<'a>(alns: &'a [Alignment], keep: &'a [bool]) -> impl Iterator<Item = &'a Alignment> {
    alns.iter().zip(keep).filter(|(_, &k)| k).map(|(a, _)| a)
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

    fn aln(q_rec: u32, t_rec: u32, t: (u64, u64), reverse: bool) -> Alignment {
        Alignment {
            q_rec,
            q_start: 0,
            q_end: t.1 - t.0,
            t_rec,
            t_start: t.0,
            t_end: t.1,
            reverse,
            matches: t.1 - t.0,
            block_len: t.1 - t.0,
            mapq: 60,
            primary: true,
            divergence: None,
        }
    }

    #[test]
    fn sorting_and_user_flips_combine() {
        let x = [rec("chr1", 100_000)];
        let y = [rec("a", 20_000), rec("b", 30_000)];
        // b sits before a on chr1; a is reversed.
        let alns = [
            aln(0, 0, (60_000, 80_000), true),
            aln(1, 0, (0, 30_000), false),
        ];
        let cfg = FilterConfig::default();
        let s = arrange(&alns, &x, &y, &cfg, Some(SortAxis::Vertical), &[]);
        let order: Vec<u32> = s.y.segments.iter().map(|g| g.rec).collect();
        assert_eq!(order, [1, 0]);
        assert!(s.y.segment(0).reverse && !s.y.segment(1).reverse);
        // A user flip of a undoes the sort's flip; of b adds one.
        let s = arrange(&alns, &x, &y, &cfg, Some(SortAxis::Vertical), &[true, true]);
        assert!(!s.y.segment(0).reverse && s.y.segment(1).reverse);
        let u = arrange(&alns, &x, &y, &cfg, None, &[false, true]);
        let order: Vec<u32> = u.y.segments.iter().map(|g| g.rec).collect();
        assert_eq!(order, [0, 1]);
        assert!(!u.y.segment(0).reverse && u.y.segment(1).reverse);
        assert_eq!(u.drawn(&alns).count(), 2);
    }

    #[test]
    fn reference_records_sort_along_the_query() {
        // One query record; reference r2 matches its start, r1 (reversed)
        // its end.
        let x = [rec("r1", 20_000), rec("r2", 30_000)];
        let y = [rec("q", 100_000)];
        let mut a1 = aln(0, 0, (0, 20_000), true);
        (a1.q_start, a1.q_end) = (70_000, 90_000);
        let a2 = aln(0, 1, (0, 30_000), false);
        let alns = [a1, a2];
        let s = arrange(
            &alns,
            &x,
            &y,
            &FilterConfig::default(),
            Some(SortAxis::Horizontal),
            &[],
        );
        let order: Vec<u32> = s.x.segments.iter().map(|g| g.rec).collect();
        assert_eq!(order, [1, 0]);
        assert!(!s.x.segment(1).reverse && s.x.segment(0).reverse);
        let qorder: Vec<u32> = s.y.segments.iter().map(|g| g.rec).collect();
        assert_eq!(qorder, [0]);
        // Placements still describe the vertical records.
        assert_eq!(s.placements.len(), 1);
    }
}
