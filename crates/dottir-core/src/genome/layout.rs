//! Axis layout: which record sits where along an axis, in which
//! orientation, and the "other" block of merged small records (D-GENIES).

use super::alignment::GenomeRecord;

/// D-GENIES merges records shorter than this share of the assembly into
/// one block at the end of the axis.
pub const DEFAULT_SMALL_FRACTION: f64 = 0.002;

/// One record placed on an axis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AxisSegment {
    pub rec: u32,
    /// Axis coordinate of the record's first drawn base.
    pub start: u64,
    pub len: u64,
    /// Drawn reverse-complemented (record end at `start`).
    pub reverse: bool,
    /// Part of the merged "other" block.
    pub small: bool,
}

/// Records of one assembly laid end to end along an axis.
#[derive(Debug, Clone, PartialEq)]
pub struct AxisLayout {
    /// In axis order.
    pub segments: Vec<AxisSegment>,
    /// Record index → index into `segments`.
    by_rec: Vec<u32>,
    pub total: u64,
    /// Axis range of the merged small records, if any.
    pub other: Option<std::ops::Range<u64>>,
}

impl AxisLayout {
    /// Lay out `records` in `order` (a permutation of record indices;
    /// `None` = file order), drawing records with `flips[rec]` reversed.
    /// Records shorter than `small_fraction` of the total go last, in the
    /// same relative order, as the "other" block — unless every record
    /// is that short.
    pub fn new(
        records: &[GenomeRecord],
        order: Option<&[u32]>,
        flips: &[bool],
        small_fraction: f64,
    ) -> Self {
        let identity: Vec<u32> = (0..records.len() as u32).collect();
        let order = order.unwrap_or(&identity);
        debug_assert_eq!(order.len(), records.len());
        let total: u64 = records.iter().map(|r| r.len).sum();
        let cutoff = (total as f64 * small_fraction) as u64;
        let is_small = |r: u32| records[r as usize].len < cutoff;
        let merge = order.iter().any(|&r| !is_small(r));
        let (big, small): (Vec<u32>, Vec<u32>) =
            order.iter().partition(|&&r| !merge || !is_small(r));

        let mut segments = Vec::with_capacity(records.len());
        let mut at = 0u64;
        let mut other_start = None;
        for (&rec, is_other) in big
            .iter()
            .map(|r| (r, false))
            .chain(small.iter().map(|r| (r, true)))
        {
            if is_other && other_start.is_none() {
                other_start = Some(at);
            }
            let len = records[rec as usize].len;
            segments.push(AxisSegment {
                rec,
                start: at,
                len,
                reverse: flips.get(rec as usize).copied().unwrap_or(false),
                small: is_other,
            });
            at += len;
        }
        let mut by_rec = vec![0u32; records.len()];
        for (i, s) in segments.iter().enumerate() {
            by_rec[s.rec as usize] = i as u32;
        }
        Self {
            segments,
            by_rec,
            total: at,
            other: other_start.map(|s| s..at),
        }
    }

    /// The segment of record `rec`.
    pub fn segment(&self, rec: u32) -> &AxisSegment {
        &self.segments[self.by_rec[rec as usize] as usize]
    }

    /// Axis coordinate of record-local position `pos` (0..=len).
    pub fn to_axis(&self, rec: u32, pos: u64) -> u64 {
        let s = self.segment(rec);
        if s.reverse {
            s.start + s.len - pos.min(s.len)
        } else {
            s.start + pos.min(s.len)
        }
    }

    /// The segment containing axis coordinate `x`, and the record-local
    /// position there.
    pub fn locate(&self, x: u64) -> Option<(&AxisSegment, u64)> {
        let i = self.segments.partition_point(|s| s.start + s.len <= x);
        let s = self.segments.get(i)?;
        let off = x - s.start;
        Some((s, if s.reverse { s.len - 1 - off } else { off }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn recs(lens: &[u64]) -> Vec<GenomeRecord> {
        lens.iter()
            .enumerate()
            .map(|(i, &len)| GenomeRecord {
                name: format!("r{i}"),
                len,
            })
            .collect()
    }

    #[test]
    fn file_order_and_coordinates() {
        let l = AxisLayout::new(&recs(&[100, 50, 30]), None, &[], 0.0);
        assert_eq!(l.total, 180);
        assert_eq!(l.to_axis(1, 10), 110);
        assert_eq!(l.locate(110).map(|(s, p)| (s.rec, p)), Some((1, 10)));
        assert_eq!(l.locate(180), None);
        assert!(l.other.is_none());
    }

    #[test]
    fn order_and_flips() {
        let l = AxisLayout::new(
            &recs(&[100, 50, 30]),
            Some(&[2, 0, 1]),
            &[true, false, false],
            0.0,
        );
        assert_eq!(l.segment(2).start, 0);
        assert_eq!(l.segment(0).start, 30);
        // Record 0 is drawn reversed: its position 0 sits at its end.
        assert_eq!(l.to_axis(0, 0), 130);
        assert_eq!(l.to_axis(0, 100), 30);
        assert_eq!(l.locate(30).map(|(s, p)| (s.rec, p)), Some((0, 99)));
    }

    #[test]
    fn small_records_merge_at_the_end() {
        // 1 + 2 bp below 0.2 % of 10_003 bp; they follow the big records.
        let l = AxisLayout::new(
            &recs(&[1, 5000, 2, 5000]),
            None,
            &[],
            DEFAULT_SMALL_FRACTION,
        );
        let order: Vec<u32> = l.segments.iter().map(|s| s.rec).collect();
        assert_eq!(order, vec![1, 3, 0, 2]);
        assert_eq!(l.other, Some(10_000..10_003));
        // All small: nothing is merged.
        let l = AxisLayout::new(&recs(&[1, 1, 1]), None, &[], 0.5);
        assert!(l.other.is_none());
    }
}
