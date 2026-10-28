//! Record-aware sequence model (docs/REVIEW.md finding #3 / A3).
//!
//! The CLI / GUI used to call [`crate::fasta::concatenate`] right after
//! parsing and threw away the record metadata, which made spec §4.4.6
//! (breaklines) and §4.4.5 (alignment export with sane coords)
//! impossible without re-parsing. This module wraps the concatenated
//! buffer with the per-record offsets and identifiers needed by
//! downstream code.

use std::path::{Path, PathBuf};

use crate::fasta::FastaRecord;

/// One record's contribution to the concatenated buffer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecordSpan {
    /// FASTA identifier — the part after `>` up to the first whitespace.
    pub id: String,
    /// Optional remainder of the header line (description), `None` if
    /// the FASTA header has no whitespace after the ID.
    pub description: Option<String>,
    /// Range in the concatenated [`Sequence::seq`] buffer where this
    /// record's residues live. End-exclusive, like all other ranges.
    pub range: std::ops::Range<usize>,
    /// Set when the record is a region of a larger one (genome-mode
    /// drill-down), so positions can be reported in its coordinates.
    pub origin: Option<RegionOrigin>,
}

/// See [`Sequence::feature_locator`].
pub struct FeatureLocator<'a> {
    by_name: std::collections::BTreeMap<&'a str, Vec<&'a RecordSpan>>,
}

impl FeatureLocator<'_> {
    /// Concatenated-buffer ranges holding `range` (0-based, end-exclusive)
    /// of record `name`, clipped to what the sequence contains; empty when
    /// none of it is present.
    pub fn locate(&self, name: &str, range: std::ops::Range<usize>) -> Vec<std::ops::Range<usize>> {
        let Some(spans) = self.by_name.get(name) else {
            return Vec::new();
        };
        spans
            .iter()
            .filter_map(|r| {
                let (lo, hi) = match &r.origin {
                    None => (range.start.min(r.len()), range.end.min(r.len())),
                    Some(o) => {
                        let a = (range.start as u64).max(o.range.start);
                        let b = (range.end as u64).min(o.range.end);
                        if a >= b {
                            return None;
                        }
                        if o.reverse {
                            ((o.range.end - b) as usize, (o.range.end - a) as usize)
                        } else {
                            ((a - o.range.start) as usize, (b - o.range.start) as usize)
                        }
                    }
                };
                (lo < hi).then(|| r.range.start + lo..r.range.start + hi)
            })
            .collect()
    }
}

/// Where a region record came from: `range` (0-based, end-exclusive) of
/// record `record`, reverse-complemented if `reverse`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RegionOrigin {
    pub record: String,
    /// Length of the whole source record.
    pub record_len: u64,
    pub range: std::ops::Range<u64>,
    pub reverse: bool,
}

impl RegionOrigin {
    /// 0-based position in the source record of region position `pos`.
    pub fn source_pos(&self, pos: usize) -> u64 {
        if self.reverse {
            self.range.end - 1 - pos as u64
        } else {
            self.range.start + pos as u64
        }
    }
}

impl RecordSpan {
    /// Record-local interval `lo..hi` in source terms: record name, record
    /// length, forward-strand interval, and whether it was reversed.
    pub fn source_interval(&self, lo: usize, hi: usize) -> (&str, u64, u64, u64, bool) {
        match &self.origin {
            None => (&self.id, self.len() as u64, lo as u64, hi as u64, false),
            Some(o) if o.reverse => (
                &o.record,
                o.record_len,
                o.range.end - hi as u64,
                o.range.end - lo as u64,
                true,
            ),
            Some(o) => (
                &o.record,
                o.record_len,
                o.range.start + lo as u64,
                o.range.start + hi as u64,
                false,
            ),
        }
    }

    /// Length of the record (in residues).
    pub fn len(&self) -> usize {
        self.range.end - self.range.start
    }

    pub fn is_empty(&self) -> bool {
        self.range.is_empty()
    }
}

/// A loaded sequence with record-boundary metadata.
///
/// Cheap to produce from `Vec<FastaRecord>` via [`Self::from_records`].
/// The concatenated `seq` buffer is the input dottir-core's kernel
/// wants; `records` lets the GUI/CLI map positions back to original
/// FASTA records.
#[derive(Debug, Clone)]
pub struct Sequence {
    /// Concatenated residues. Uppercased, whitespace stripped (per
    /// `dottir-io::fasta`).
    pub seq: Vec<u8>,
    /// Per-record metadata, in file order.
    pub records: Vec<RecordSpan>,
    /// Soft-masked (lowercase in the input) runs, in concatenated
    /// coordinates, sorted.
    pub soft_mask: Vec<std::ops::Range<usize>>,
    /// The on-disk source path, when known. `None` for in-memory
    /// inputs (tests, future drag-and-drop from STDIN, …).
    pub source_path: Option<PathBuf>,
}

impl Sequence {
    /// Build a [`Sequence`] from a parsed `Vec<FastaRecord>`. The
    /// concatenation is in file order; per-record ranges are
    /// computed from the input lengths.
    pub fn from_records(records: Vec<FastaRecord>, source: Option<PathBuf>) -> Self {
        let total: usize = records.iter().map(|r| r.sequence.len()).sum();
        let mut seq = Vec::with_capacity(total);
        let mut spans = Vec::with_capacity(records.len());
        let mut soft_mask = Vec::new();
        let mut offset = 0;
        for r in records {
            let len = r.sequence.len();
            seq.extend_from_slice(&r.sequence);
            soft_mask.extend(
                r.soft_masked
                    .iter()
                    .map(|m| m.start + offset..m.end + offset),
            );
            spans.push(RecordSpan {
                id: r.id,
                description: r.description,
                range: offset..offset + len,
                origin: None,
            });
            offset += len;
        }
        Sequence {
            seq,
            records: spans,
            soft_mask,
            source_path: source,
        }
    }

    /// Load a FASTA file (plain or gzipped) into a [`Sequence`]. Single
    /// disk read (docs/REVIEW.md finding #6 / A4).
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, crate::fasta::FastaError> {
        let p = path.as_ref().to_path_buf();
        let records = crate::fasta::read_fasta_file(&p)?;
        Ok(Self::from_records(records, Some(p)))
    }

    /// Total residues in the concatenated buffer.
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Return a repeat-period hint from record lengths, when the FASTA
    /// contains multiple records. Equal-length repeat arrays return that
    /// length; mixed lengths return their greatest common divisor.
    pub fn record_period_hint(&self) -> Option<usize> {
        record_period_hint(&self.records)
    }

    /// Inter-record offsets in the concatenated buffer. The result has
    /// `records.len() - 1` entries — one per gap between adjacent
    /// records. Empty for single-record inputs. Used by spec §4.4.6
    /// (breakline rendering).
    ///
    /// Example: records of lengths `[10, 20, 5]` → breaks `[10, 30]`.
    pub fn breaks(&self) -> Vec<usize> {
        self.records
            .windows(2)
            .map(|pair| pair[0].range.end)
            .collect()
    }

    /// Map a concatenated-buffer coordinate to its containing record
    /// and the position within that record (0-based). `None` if
    /// `coord` is past the end of the buffer.
    ///
    /// Used by the GUI status bar to render `chr4:12345` style
    /// coordinates rather than the opaque concatenated offset.
    pub fn record_at(&self, coord: usize) -> Option<(&RecordSpan, usize)> {
        // Binary search by end; records are sorted by range.
        let i = self.records.partition_point(|r| r.range.end <= coord);
        let r = self.records.get(i)?;
        if coord >= r.range.start && coord < r.range.end {
            Some((r, coord - r.range.start))
        } else {
            None
        }
    }

    /// Borrow the concatenated bytes — what `dottir-core` actually
    /// wants. Equivalent to `&self.seq` but reads better at call
    /// sites.
    pub fn bytes(&self) -> &[u8] {
        &self.seq
    }

    /// The records in `order` (record indices), reverse-complemented
    /// where `flips[i]` is set for record `i` (DNA, IUPAC-aware). Soft
    /// masking follows the bases. Reversed records, and records that were
    /// already regions, carry a [`RegionOrigin`] so positions and
    /// annotations still map to the source records.
    pub fn reordered(&self, order: &[usize], flips: &[bool]) -> Sequence {
        let mut seq = Vec::with_capacity(self.seq.len());
        let mut records = Vec::with_capacity(order.len());
        let mut soft_mask = Vec::new();
        for &i in order {
            let r = &self.records[i];
            let flip = flips.get(i).copied().unwrap_or(false);
            let at = seq.len();
            let bases = &self.seq[r.range.clone()];
            let first = self.soft_mask.partition_point(|m| m.end <= r.range.start);
            let masks = self.soft_mask[first..]
                .iter()
                .take_while(|m| m.start < r.range.end)
                .map(|m| {
                    m.start.max(r.range.start) - r.range.start
                        ..m.end.min(r.range.end) - r.range.start
                });
            if flip {
                seq.extend(
                    bases
                        .iter()
                        .rev()
                        .map(|&b| crate::genome_export::complement(b)),
                );
                let len = r.len();
                let mut m: Vec<_> = masks
                    .map(|m| at + len - m.end..at + len - m.start)
                    .collect();
                m.reverse();
                soft_mask.extend(m);
            } else {
                seq.extend_from_slice(bases);
                soft_mask.extend(masks.map(|m| at + m.start..at + m.end));
            }
            let origin = match (&r.origin, flip) {
                (Some(o), _) => Some(RegionOrigin {
                    reverse: o.reverse ^ flip,
                    ..o.clone()
                }),
                (None, true) => Some(RegionOrigin {
                    record: r.id.clone(),
                    record_len: r.len() as u64,
                    range: 0..r.len() as u64,
                    reverse: true,
                }),
                (None, false) => None,
            };
            records.push(RecordSpan {
                id: r.id.clone(),
                description: r.description.clone(),
                range: at..seq.len(),
                origin,
            });
        }
        Sequence {
            seq,
            records,
            soft_mask,
            source_path: self.source_path.clone(),
        }
    }

    /// Maps record-local feature coordinates (GFF3 / BED) into this
    /// sequence, through region origins.
    pub fn feature_locator(&self) -> FeatureLocator<'_> {
        let mut by_name: std::collections::BTreeMap<&str, Vec<&RecordSpan>> =
            std::collections::BTreeMap::new();
        for r in &self.records {
            let name = r
                .origin
                .as_ref()
                .map_or(r.id.as_str(), |o| o.record.as_str());
            by_name.entry(name).or_default().push(r);
        }
        FeatureLocator { by_name }
    }

    /// Heuristically classify this sequence as DNA or protein.
    /// Same rule the original Dotter uses (`detectSeqType`): if the
    /// non-`ACGTNU` ratio across the first sampled residues exceeds a
    /// threshold, it's protein. Whitespace is already stripped by the
    /// FASTA reader, so every byte here is a residue.
    pub fn detect_alphabet(&self) -> DetectedAlphabet {
        detect_alphabet(&self.seq)
    }
}

/// Return a repeat-period hint from record lengths, when useful.
pub fn record_period_hint(records: &[RecordSpan]) -> Option<usize> {
    if records.len() < 2 {
        return None;
    }
    let mut period = 0usize;
    for len in records.iter().map(RecordSpan::len).filter(|&len| len > 0) {
        period = if period == 0 {
            len
        } else {
            gcd_usize(period, len)
        };
    }
    (period > 1).then_some(period)
}

fn gcd_usize(mut a: usize, mut b: usize) -> usize {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

/// Result of [`Sequence::detect_alphabet`] / [`detect_alphabet`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DetectedAlphabet {
    Dna,
    Protein,
    /// Buffer empty or too short to classify.
    Unknown,
}

/// Classify a residue buffer as DNA or protein. The rule: count the
/// fraction of bytes outside `ACGTNU` (plus `acgtnu` since the FASTA
/// reader could in theory hand back mixed case before uppercasing). If
/// that fraction is above 10% it's protein. Up to the first 4096
/// residues are sampled, which is plenty to disambiguate any real
/// FASTA without scanning a 500-Mb genome.
pub fn detect_alphabet(seq: &[u8]) -> DetectedAlphabet {
    if seq.is_empty() {
        return DetectedAlphabet::Unknown;
    }
    let sample = &seq[..seq.len().min(4096)];
    let mut non_dna: usize = 0;
    for &b in sample {
        match b {
            b'A' | b'C' | b'G' | b'T' | b'N' | b'U' => {}
            b'a' | b'c' | b'g' | b't' | b'n' | b'u' => {}
            _ if b.is_ascii_alphabetic() => non_dna += 1,
            _ => {}
        }
    }
    let ratio = non_dna as f64 / sample.len() as f64;
    if ratio > 0.10 {
        DetectedAlphabet::Protein
    } else {
        DetectedAlphabet::Dna
    }
}

#[cfg(test)]
mod reorder_tests {
    use super::*;
    use crate::fasta::parse_fasta;

    #[test]
    fn reordering_flips_bases_masks_and_feature_positions() {
        let recs = parse_fasta(">a\nAAcgT\n>b\nGGGG\n").unwrap();
        let s = Sequence::from_records(recs, None);
        let r = s.reordered(&[1, 0], &[true, false]);
        // b, then a reverse-complemented: AcgTT (mask on cg -> cg).
        assert_eq!(r.seq, b"GGGGACGTT");
        assert_eq!(r.soft_mask, vec![5..7]);
        assert_eq!(r.records[1].id, "a");
        let o = r.records[1].origin.as_ref().unwrap();
        assert!(o.reverse && o.range == (0..5));
        assert!(r.records[0].origin.is_none());
        // a[3..5) ("gT") lies at buffer 4..6 after the flip.
        let loc = r.feature_locator();
        assert_eq!(loc.locate("a", 3..5), vec![4..6]);
        assert_eq!(loc.locate("b", 1..2), vec![1..2]);
        assert!(loc.locate("c", 0..1).is_empty());
        // Flipping back restores the original.
        let back = r.reordered(&[1, 0], &[false, true]);
        assert_eq!(back.seq, s.seq);
        assert_eq!(back.soft_mask, s.soft_mask);
        assert!(!back.records[0].origin.as_ref().unwrap().reverse);
    }

    #[test]
    fn locator_clips_features_to_regions() {
        let s = Sequence::from_records(parse_fasta(">r\nACGTACGTAC\n").unwrap(), None);
        let mut region = s.clone();
        region.records[0].origin = Some(RegionOrigin {
            record: "r".into(),
            record_len: 200,
            range: 100..110,
            reverse: true,
        });
        let loc = region.feature_locator();
        // Source 95..103 overlaps the region at 100..103, which is the
        // region's last three bases when reversed.
        assert_eq!(loc.locate("r", 95..103), vec![7..10]);
        assert!(loc.locate("r", 0..50).is_empty());
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(id: &str, seq: &[u8]) -> FastaRecord {
        FastaRecord {
            id: id.to_string(),
            description: None,
            sequence: seq.to_vec(),
            soft_masked: vec![],
        }
    }

    #[test]
    fn from_records_concatenates_in_order() {
        let s = Sequence::from_records(
            vec![rec("a", b"AAAA"), rec("b", b"CCCC"), rec("c", b"GGGGGG")],
            None,
        );
        assert_eq!(s.seq, b"AAAACCCCGGGGGG");
        assert_eq!(s.records.len(), 3);
        assert_eq!(s.records[0].range, 0..4);
        assert_eq!(s.records[1].range, 4..8);
        assert_eq!(s.records[2].range, 8..14);
    }

    #[test]
    fn breaks_returns_inter_record_offsets() {
        let s = Sequence::from_records(
            vec![rec("a", b"AAAA"), rec("b", b"CCCC"), rec("c", b"GG")],
            None,
        );
        assert_eq!(s.breaks(), vec![4, 8]);
    }

    #[test]
    fn breaks_empty_for_single_record() {
        let s = Sequence::from_records(vec![rec("only", b"ACGT")], None);
        assert!(s.breaks().is_empty());
    }

    #[test]
    fn record_period_hint_uses_equal_record_length() {
        let s = Sequence::from_records(
            vec![rec("a", b"AAAA"), rec("b", b"CCCC"), rec("c", b"GGGG")],
            None,
        );
        assert_eq!(s.record_period_hint(), Some(4));
    }

    #[test]
    fn record_period_hint_uses_gcd_for_mixed_lengths() {
        let s = Sequence::from_records(vec![rec("a", b"AAAAAA"), rec("b", b"CCCCCCCCC")], None);
        assert_eq!(s.record_period_hint(), Some(3));
    }

    #[test]
    fn record_at_returns_correct_record_and_position() {
        let s = Sequence::from_records(
            vec![rec("a", b"AAAA"), rec("b", b"CCCC"), rec("c", b"GGGG")],
            None,
        );
        // First record.
        let (r, p) = s.record_at(0).unwrap();
        assert_eq!(r.id, "a");
        assert_eq!(p, 0);
        let (r, p) = s.record_at(3).unwrap();
        assert_eq!(r.id, "a");
        assert_eq!(p, 3);
        // Second record.
        let (r, p) = s.record_at(4).unwrap();
        assert_eq!(r.id, "b");
        assert_eq!(p, 0);
        let (r, p) = s.record_at(7).unwrap();
        assert_eq!(r.id, "b");
        assert_eq!(p, 3);
        // Third.
        let (r, p) = s.record_at(11).unwrap();
        assert_eq!(r.id, "c");
        assert_eq!(p, 3);
        // Past end.
        assert!(s.record_at(12).is_none());
        assert!(s.record_at(usize::MAX).is_none());
    }

    #[test]
    fn record_at_handles_empty_sequence() {
        let s = Sequence::from_records(vec![], None);
        assert!(s.record_at(0).is_none());
    }

    #[test]
    fn bytes_matches_seq() {
        let s = Sequence::from_records(vec![rec("x", b"ACGT")], None);
        assert_eq!(s.bytes(), b"ACGT");
    }

    #[test]
    fn detect_alphabet_recognises_dna_protein_and_empty() {
        assert_eq!(detect_alphabet(b""), DetectedAlphabet::Unknown);
        assert_eq!(detect_alphabet(b"ACGTACGTNACGT"), DetectedAlphabet::Dna,);
        // Real protein (P00533 prefix) — has K, T, I, etc.
        assert_eq!(
            detect_alphabet(b"MRPSGTAGAALLALLAALCPASRALEEKKVCQGTSNKLT"),
            DetectedAlphabet::Protein,
        );
        // Edge case: pure-N "DNA" classification.
        assert_eq!(detect_alphabet(b"NNNNNNNN"), DetectedAlphabet::Dna);
    }
}
