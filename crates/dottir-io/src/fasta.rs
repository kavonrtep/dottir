//! Minimal FASTA reader. Handles plain and gzipped files (`.gz`
//! transparently via `flate2`). Multi-record input is supported: records
//! are returned in file order. Sequences are uppercased; whitespace within
//! a record is dropped. IUPAC ambiguity codes are passed through unchanged
//! (the alphabet encoding in `dottir-core` filters anything not in its
//! 24- or 4-letter table).
//!
//! For Phase 1+ usage we read whole files into memory. Streaming readers
//! and `mmap` are deferred until they're a measured bottleneck (the
//! `dottir-core` kernel allocates an O(n) score vector and pixelmap, so
//! whole-file FASTA isn't the dominant cost).
//!
//! This is intentionally not a Crab-class FASTA parser — for the dottir
//! use case (a query/subject pair, plus annotation files) it covers the
//! ground and avoids a dependency that breaks MSRV-1.75.

use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use thiserror::Error;

#[derive(Debug, Clone)]
pub struct FastaRecord {
    /// Identifier — the part of the header line after `>` up to the
    /// first whitespace character.
    pub id: String,
    /// Optional remainder of the header line (after the first
    /// whitespace), `None` if there isn't any.
    pub description: Option<String>,
    /// Sequence, uppercased, whitespace stripped.
    pub sequence: Vec<u8>,
}

#[derive(Debug, Error)]
pub enum FastaError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("file does not look like FASTA: missing '>' header at start")]
    MissingHeader,
    #[error("empty record at byte offset {0}")]
    EmptyRecord(u64),
}

/// Read every record from a FASTA file. Gzip is auto-detected from the
/// `.gz` extension OR from the file's first two magic bytes (0x1f 0x8b);
/// callers don't have to special-case it.
pub fn read_fasta_file<P: AsRef<Path>>(path: P) -> Result<Vec<FastaRecord>, FastaError> {
    Ok(load_fasta_file(path)?.records)
}

/// Result of [`load_fasta_file`]: the parsed records plus the raw
/// on-disk bytes. Callers that want to hash the input file for a
/// params sidecar can use `bytes` rather than re-reading the file
/// (REVIEW.md finding #6 / A4).
#[derive(Debug, Clone)]
pub struct LoadedFasta {
    pub records: Vec<FastaRecord>,
    pub bytes: Vec<u8>,
}

/// Single-pass FASTA load: reads the file once, returns both the
/// parsed records and the raw on-disk bytes for downstream hashing
/// or re-encoding.
pub fn load_fasta_file<P: AsRef<Path>>(path: P) -> Result<LoadedFasta, FastaError> {
    let bytes = std::fs::read(path.as_ref())?;
    let is_gzipped = bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b
        || path.as_ref().extension().and_then(|s| s.to_str()) == Some("gz");
    let records = if is_gzipped {
        read_fasta(BufReader::new(MultiGzDecoder::new(&bytes[..])))?
    } else {
        parse_fasta(&String::from_utf8_lossy(&bytes))?
    };
    Ok(LoadedFasta { records, bytes })
}

/// Parse a FASTA stream into records. Generic over any `BufRead` source.
pub fn read_fasta<R: BufRead>(mut reader: R) -> Result<Vec<FastaRecord>, FastaError> {
    let mut text = String::new();
    reader.read_to_string(&mut text)?;
    parse_fasta(&text)
}

/// Parse an in-memory FASTA string. Mostly used by tests; production code
/// should prefer the streaming entry points which avoid the full buffer.
pub fn parse_fasta(text: &str) -> Result<Vec<FastaRecord>, FastaError> {
    let mut records: Vec<FastaRecord> = Vec::new();
    let mut current_header: Option<(String, Option<String>)> = None;
    let mut current_seq: Vec<u8> = Vec::new();
    let mut header_offset: u64 = 0;
    let mut offset: u64 = 0;

    for line in text.lines() {
        let bytes = line.as_bytes();
        offset += line.len() as u64 + 1; // +1 for \n
        if let Some(rest) = line.strip_prefix('>') {
            // Flush previous record.
            if let Some((id, description)) = current_header.take() {
                if current_seq.is_empty() {
                    return Err(FastaError::EmptyRecord(header_offset));
                }
                records.push(FastaRecord {
                    id,
                    description,
                    sequence: std::mem::take(&mut current_seq),
                });
            } else if !current_seq.is_empty() {
                return Err(FastaError::MissingHeader);
            }
            let trimmed = rest.trim_start();
            let mut split = trimmed.splitn(2, char::is_whitespace);
            let id = split.next().unwrap_or("").to_string();
            let description = split.next().map(|s| s.to_string());
            current_header = Some((id, description));
            header_offset = offset.saturating_sub(line.len() as u64 + 1);
        } else {
            // Sequence line: strip whitespace, uppercase, append.
            for &b in bytes {
                if !b.is_ascii_whitespace() {
                    current_seq.push(b.to_ascii_uppercase());
                }
            }
        }
    }
    // Flush trailing record.
    if let Some((id, description)) = current_header {
        if current_seq.is_empty() {
            return Err(FastaError::EmptyRecord(header_offset));
        }
        records.push(FastaRecord {
            id,
            description,
            sequence: current_seq,
        });
    } else if !current_seq.is_empty() {
        return Err(FastaError::MissingHeader);
    }
    Ok(records)
}

/// Concatenate all records into a single sequence with optional break
/// markers (used by the multi-record breakline feature). For Phase 1 the
/// CLI just uses the first record per file or concatenation without
/// markers; Phase 5+ will surface break offsets to the GUI.
pub fn concatenate(records: &[FastaRecord]) -> Vec<u8> {
    let total: usize = records.iter().map(|r| r.sequence.len()).sum();
    let mut out = Vec::with_capacity(total);
    for r in records {
        out.extend_from_slice(&r.sequence);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_single_record() {
        let text = ">seq1 description\nACGT\nACGT\n";
        let recs = parse_fasta(text).unwrap();
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].id, "seq1");
        assert_eq!(recs[0].description.as_deref(), Some("description"));
        assert_eq!(recs[0].sequence, b"ACGTACGT");
    }

    #[test]
    fn parse_multi_record() {
        let text = ">a\nAAAA\n>b some desc\nGGGG\nGGGG\n";
        let recs = parse_fasta(text).unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, "a");
        assert_eq!(recs[0].sequence, b"AAAA");
        assert_eq!(recs[1].id, "b");
        assert_eq!(recs[1].description.as_deref(), Some("some desc"));
        assert_eq!(recs[1].sequence, b"GGGGGGGG");
    }

    #[test]
    fn lowercase_uppercased() {
        let recs = parse_fasta(">a\nacgt\n").unwrap();
        assert_eq!(recs[0].sequence, b"ACGT");
    }

    #[test]
    fn missing_header_errors() {
        let recs = parse_fasta("ACGT\n");
        assert!(matches!(recs, Err(FastaError::MissingHeader)));
    }

    #[test]
    fn empty_record_errors() {
        let recs = parse_fasta(">a\n\n");
        assert!(matches!(recs, Err(FastaError::EmptyRecord(_))));
    }

    #[test]
    fn whitespace_stripped() {
        let recs = parse_fasta(">a\nA C\tG T\n\nA\n").unwrap();
        assert_eq!(recs[0].sequence, b"ACGTA");
    }

    #[test]
    fn concatenate_joins_in_order() {
        let recs = parse_fasta(">a\nAA\n>b\nCC\n>c\nGG\n").unwrap();
        let joined = concatenate(&recs);
        assert_eq!(joined, b"AACCGG");
    }
}
