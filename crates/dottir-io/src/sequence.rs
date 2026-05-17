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
}

impl RecordSpan {
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
        let mut offset = 0;
        for r in records {
            let len = r.sequence.len();
            seq.extend_from_slice(&r.sequence);
            spans.push(RecordSpan {
                id: r.id,
                description: r.description,
                range: offset..offset + len,
            });
            offset += len;
        }
        Sequence {
            seq,
            records: spans,
            source_path: source,
        }
    }

    /// Load a FASTA file (plain or gzipped) into a [`Sequence`]. Single
    /// disk read (docs/REVIEW.md finding #6 / A4).
    pub fn load<P: AsRef<Path>>(
        path: P,
    ) -> Result<Self, crate::fasta::FastaError> {
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
        let i = self
            .records
            .partition_point(|r| r.range.end <= coord);
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
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(id: &str, seq: &[u8]) -> FastaRecord {
        FastaRecord {
            id: id.to_string(),
            description: None,
            sequence: seq.to_vec(),
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
}
