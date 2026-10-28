//! Axis layout: which record sits where along an axis, in which
//! orientation, and which runs of small records are merged into "mix"
//! blocks (as D-GENIES `Paf.parse_index`).

use super::alignment::GenomeRecord;

/// Records shorter than this share of the assembly are small.
pub const DEFAULT_SMALL_FRACTION: f64 = 0.002;
/// A run of at least this many consecutive small records is merged into
/// one mix block; shorter runs are shown as separate records.
pub const MIN_MIX_RECORDS: usize = 5;

/// One record placed on an axis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AxisSegment {
    pub rec: u32,
    /// Axis coordinate of the record's first drawn base.
    pub start: u64,
    pub len: u64,
    /// Drawn reverse-complemented (record end at `start`).
    pub reverse: bool,
    /// Index into [`AxisLayout::mixes`] when part of a mix block.
    pub mix: Option<u32>,
}

/// Records of one assembly laid end to end along an axis.
#[derive(Debug, Clone, PartialEq)]
pub struct AxisLayout {
    /// In axis order.
    pub segments: Vec<AxisSegment>,
    /// Record index → index into `segments`.
    by_rec: Vec<u32>,
    pub total: u64,
    /// Axis ranges of the mix blocks, in axis order.
    pub mixes: Vec<std::ops::Range<u64>>,
}

impl AxisLayout {
    /// Lay out `records` in `order` (a permutation of record indices;
    /// `None` = file order), drawing records with `flips[rec]` reversed.
    /// Runs of at least [`MIN_MIX_RECORDS`] consecutive records shorter
    /// than `small_fraction` of the total become mix blocks.
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
        let cutoff = total as f64 * small_fraction;

        let mut segments: Vec<AxisSegment> = Vec::with_capacity(records.len());
        let mut at = 0u64;
        for &rec in order {
            let len = records[rec as usize].len;
            segments.push(AxisSegment {
                rec,
                start: at,
                len,
                reverse: flips.get(rec as usize).copied().unwrap_or(false),
                mix: None,
            });
            at += len;
        }
        let small = |s: &AxisSegment| (s.len as f64) < cutoff;
        let mut mixes = Vec::new();
        let mut i = 0;
        while i < segments.len() {
            if !small(&segments[i]) {
                i += 1;
                continue;
            }
            let j = i + segments[i..].iter().take_while(|s| small(s)).count();
            if j - i >= MIN_MIX_RECORDS {
                let m = mixes.len() as u32;
                mixes.push(segments[i].start..segments[j - 1].start + segments[j - 1].len);
                for s in &mut segments[i..j] {
                    s.mix = Some(m);
                }
            }
            i = j;
        }
        let mut by_rec = vec![0u32; records.len()];
        for (i, s) in segments.iter().enumerate() {
            by_rec[s.rec as usize] = i as u32;
        }
        Self {
            segments,
            by_rec,
            total: at,
            mixes,
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

    /// The record pieces under axis range `range`, in axis order: record,
    /// record-local range (forward coordinates), drawn reversed.
    pub fn regions(&self, range: std::ops::Range<u64>) -> Vec<(u32, std::ops::Range<u64>, bool)> {
        let first = self
            .segments
            .partition_point(|s| s.start + s.len <= range.start);
        self.segments[first..]
            .iter()
            .take_while(|s| s.start < range.end)
            .filter_map(|s| {
                let a = range.start.max(s.start) - s.start;
                let b = range.end.min(s.start + s.len) - s.start;
                if a >= b {
                    return None;
                }
                let local = if s.reverse {
                    s.len - b..s.len - a
                } else {
                    a..b
                };
                Some((s.rec, local, s.reverse))
            })
            .collect()
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
    fn regions_map_axis_ranges_back_to_records() {
        // r1 (50) is drawn first and reversed, then r0 (100), r2 (30).
        let l = AxisLayout::new(
            &recs(&[100, 50, 30]),
            Some(&[1, 0, 2]),
            &[false, true, false],
            0.0,
        );
        assert_eq!(l.regions(10..60), vec![(1, 0..40, true), (0, 0..10, false)]);
        assert_eq!(l.regions(0..180).len(), 3);
        assert_eq!(l.regions(155..170), vec![(2, 5..20, false)]);
        assert!(l.regions(180..200).is_empty());
        // Axis 10 is r1 position 39 (reversed); the region's last base.
        assert_eq!(l.locate(10).map(|(s, p)| (s.rec, p)), Some((1, 39)));
    }

    #[test]
    fn file_order_and_coordinates() {
        let l = AxisLayout::new(&recs(&[100, 50, 30]), None, &[], 0.0);
        assert_eq!(l.total, 180);
        assert_eq!(l.to_axis(1, 10), 110);
        assert_eq!(l.locate(110).map(|(s, p)| (s.rec, p)), Some((1, 10)));
        assert_eq!(l.locate(180), None);
        assert!(l.mixes.is_empty());
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
    fn runs_of_five_small_records_become_mix_blocks_in_place() {
        // Total 10_009; small = below 20 bp. Records 1..=5 form a run of
        // five (merged); records 7..=8 a run of two (kept separate).
        let l = AxisLayout::new(
            &recs(&[5000, 1, 1, 1, 1, 1, 4000, 2, 2, 1000]),
            None,
            &[],
            DEFAULT_SMALL_FRACTION,
        );
        assert_eq!(l.mixes, vec![5000..5005]);
        let mixed: Vec<u32> = l
            .segments
            .iter()
            .filter(|s| s.mix.is_some())
            .map(|s| s.rec)
            .collect();
        assert_eq!(mixed, vec![1, 2, 3, 4, 5]);
        assert_eq!(l.segment(6).start, 5005);
    }
}
