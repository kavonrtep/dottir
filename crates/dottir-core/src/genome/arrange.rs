//! Filters, sorting and axis layouts for one genome plot, shared by the
//! batch CLI and the GUI.

use super::alignment::{Alignment, GenomeRecord};
use super::filter::{keep_mask, FilterConfig};
use super::layout::{AxisLayout, DEFAULT_SMALL_FRACTION};
use super::order::{sort_query, Placement};

/// What to draw and where.
#[derive(Debug, Clone, PartialEq)]
pub struct Arranged {
    /// Per alignment: passes the filters.
    pub keep: Vec<bool>,
    /// Noise length threshold, when the noise filter applied.
    pub noise: Option<f64>,
    /// Horizontal axis (reference, alignment targets), in file order.
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

/// Filter `alns` and lay out both axes. With `sort`, the vertical
/// records are ordered and oriented along the reference as D-GENIES
/// does, from the noise-filtered alignments (spec §4.7.6). `flips`
/// (per vertical record, may be empty) reverses records on top of that.
pub fn arrange(
    alns: &[Alignment],
    x_records: &[GenomeRecord],
    y_records: &[GenomeRecord],
    filter: &FilterConfig,
    sort: bool,
    flips: &[bool],
) -> Arranged {
    let (keep, noise) = keep_mask(alns, filter);
    let x = AxisLayout::new(x_records, None, &[], DEFAULT_SMALL_FRACTION);
    let flipped = |base: &[bool]| -> Vec<bool> {
        (0..y_records.len())
            .map(|r| base.get(r).copied().unwrap_or(false) ^ flips.get(r).copied().unwrap_or(false))
            .collect()
    };
    let (y, placements) = if sort {
        let (sort_keep, _) = keep_mask(
            alns,
            &FilterConfig {
                remove_noise: true,
                ..filter.clone()
            },
        );
        let used: Vec<&Alignment> = select(alns, &sort_keep).collect();
        let o = sort_query(&used, y_records, x_records, &x);
        let y = AxisLayout::new(
            y_records,
            Some(&o.order),
            &flipped(&o.flips),
            DEFAULT_SMALL_FRACTION,
        );
        (y, o.placements)
    } else {
        let used: Vec<&Alignment> = select(alns, &keep).collect();
        let placements = sort_query(&used, y_records, x_records, &x).placements;
        let y = AxisLayout::new(y_records, None, &flipped(&[]), DEFAULT_SMALL_FRACTION);
        (y, placements)
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
        let s = arrange(&alns, &x, &y, &cfg, true, &[]);
        let order: Vec<u32> = s.y.segments.iter().map(|g| g.rec).collect();
        assert_eq!(order, [1, 0]);
        assert!(s.y.segment(0).reverse && !s.y.segment(1).reverse);
        // A user flip of a undoes the sort's flip; of b adds one.
        let s = arrange(&alns, &x, &y, &cfg, true, &[true, true]);
        assert!(!s.y.segment(0).reverse && s.y.segment(1).reverse);
        let u = arrange(&alns, &x, &y, &cfg, false, &[false, true]);
        let order: Vec<u32> = u.y.segments.iter().map(|g| g.rec).collect();
        assert_eq!(order, [0, 1]);
        assert!(!u.y.segment(0).reverse && u.y.segment(1).reverse);
        assert_eq!(u.drawn(&alns).count(), 2);
    }
}
