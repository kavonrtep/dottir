//! Sequence input loading — FASTA, or GFF3 carrying its own sequences.
//!
//! GFF3 may embed the sequences it annotates after a `##FASTA`
//! directive (GFF3 spec, "Sequences" section). Such a file is
//! self-contained: it is both the sequence input and the annotation
//! input. This module is the single entry point the CLI and the GUI
//! use for a positional input path, so `dottir batch a.gff3 b.gff3`
//! and `dottir batch a.fasta b.fasta` take the same code path.
//!
//! Dispatch is by extension (after stripping `.gz`), with a content
//! sniff as a fallback for files named something else:
//!
//! | Input                        | Sequence                  | Annotations |
//! |------------------------------|---------------------------|-------------|
//! | `.fa/.fasta/.fna/.faa/…`     | records in file order     | none        |
//! | `.gff/.gff3` with `##FASTA`  | records in file order     | column 1-9  |
//! | `.gff/.gff3` without FASTA   | error ([`InputError::NoEmbeddedFasta`]) | — |
//!
//! Feature coordinates stay record-local, exactly as in
//! [`crate::annotation`]; binding them to the concatenated buffer is
//! the caller's job.

use std::path::{Path, PathBuf};

use thiserror::Error;

use crate::annotation::{self, AnnotError, AnnotSet, AnnotSource};
use crate::fasta::{self, FastaError};
use crate::sequence::Sequence;

/// Which on-disk shape an input path turned out to have.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InputKind {
    /// Plain FASTA. No annotations.
    Fasta,
    /// GFF3 with an embedded `##FASTA` section.
    Gff3WithFasta,
}

/// A loaded positional input: the sequence, plus the annotations that
/// came with it when the file carried both.
#[derive(Debug, Clone)]
pub struct SeqInput {
    pub sequence: Sequence,
    /// Features parsed from the same file. `None` for plain FASTA;
    /// `Some` (possibly with zero features) for GFF3.
    pub annotations: Option<AnnotSet>,
    /// Raw on-disk bytes, for the params sidecar hash. Still gzipped
    /// if the file was.
    pub bytes: Vec<u8>,
    pub path: PathBuf,
    pub kind: InputKind,
}

#[derive(Debug, Error)]
pub enum InputError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    // `path` is carried for programmatic use but kept out of Display:
    // the CLI wraps these in a "reading query <path>" context and the
    // GUI prefixes "failed to load <path>", so repeating it here reads
    // as a stutter.
    #[error("{source}")]
    Fasta {
        path: PathBuf,
        #[source]
        source: FastaError,
    },
    #[error("{source}")]
    Annot {
        path: PathBuf,
        #[source]
        source: AnnotError,
    },
    #[error(
        "{0} is a GFF3 file with no embedded sequences: either append a \
         `##FASTA` section to it, or pass the FASTA separately and give \
         this file with --gff-query/--gff-subject"
    )]
    NoEmbeddedFasta(PathBuf),
    #[error("{0} has a `##FASTA` section but no records in it")]
    EmptyEmbeddedFasta(PathBuf),
    #[error(
        "{0} is not a recognized sequence input: expected FASTA (a `>` \
         header line) or GFF3 with an embedded `##FASTA` section (a \
         `##gff-version` pragma). Plain or gzipped."
    )]
    Unrecognized(PathBuf),
}

/// Load a positional sequence input from `path`. See the module docs
/// for the dispatch table.
pub fn load_sequence_input<P: AsRef<Path>>(path: P) -> Result<SeqInput, InputError> {
    let path = path.as_ref();
    let bytes = std::fs::read(path)?;
    from_bytes(path, bytes)
}

/// [`load_sequence_input`] on bytes already in memory. Split out so
/// tests (and any future STDIN / drag-and-drop path) don't need a
/// temporary file.
pub fn from_bytes(path: &Path, bytes: Vec<u8>) -> Result<SeqInput, InputError> {
    let is_gzipped =
        (bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b) || ext_is(path, "gz");
    // Strip `.gz` before looking at the format extension, like
    // `AnnotSet::load`.
    let format_path = if ext_is(path, "gz") {
        path.with_extension("")
    } else {
        path.to_path_buf()
    };

    let plain = annotation::decompress(&bytes, is_gzipped).map_err(|source| InputError::Annot {
        path: path.to_path_buf(),
        source,
    })?;

    if !looks_like_gff3(&format_path, &plain) {
        let records =
            fasta::parse_fasta(&String::from_utf8_lossy(&plain)).map_err(
                |source| match source {
                    // No `>` anywhere and no GFF3 pragma: this is not a
                    // format we handle, which is more useful to say than
                    // "missing '>' header at start".
                    FastaError::MissingHeader => InputError::Unrecognized(path.to_path_buf()),
                    source => InputError::Fasta {
                        path: path.to_path_buf(),
                        source,
                    },
                },
            )?;
        return Ok(SeqInput {
            sequence: Sequence::from_records(records, Some(path.to_path_buf())),
            annotations: None,
            bytes,
            path: path.to_path_buf(),
            kind: InputKind::Fasta,
        });
    }

    let (feature_bytes, fasta_bytes) = annotation::split_fasta_section(&plain);
    let fasta_bytes = fasta_bytes.ok_or_else(|| InputError::NoEmbeddedFasta(path.to_path_buf()))?;
    let records = fasta::parse_fasta(&String::from_utf8_lossy(fasta_bytes)).map_err(|source| {
        InputError::Fasta {
            path: path.to_path_buf(),
            source,
        }
    })?;
    if records.is_empty() {
        return Err(InputError::EmptyEmbeddedFasta(path.to_path_buf()));
    }

    let features =
        annotation::parse_gff3_bytes(feature_bytes).map_err(|source| InputError::Annot {
            path: path.to_path_buf(),
            source,
        })?;

    Ok(SeqInput {
        sequence: Sequence::from_records(records, Some(path.to_path_buf())),
        annotations: Some(AnnotSet {
            source_path: Some(path.to_path_buf()),
            source: AnnotSource::Gff3,
            features,
        }),
        bytes,
        path: path.to_path_buf(),
        kind: InputKind::Gff3WithFasta,
    })
}

fn ext_is(path: &Path, ext: &str) -> bool {
    path.extension()
        .and_then(|s| s.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case(ext))
}

/// `.gff`/`.gff3` decides it; otherwise sniff for the mandatory
/// `##gff-version` pragma so a GFF3 named `foo.txt` still works. A
/// leading `>` always wins — that is FASTA whatever the extension.
fn looks_like_gff3(format_path: &Path, plain: &[u8]) -> bool {
    if ext_is(format_path, "gff") || ext_is(format_path, "gff3") {
        return true;
    }
    for line in plain.split(|&b| b == b'\n').take(8) {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.is_empty() {
            continue;
        }
        if line.starts_with(b">") {
            return false;
        }
        if line.to_ascii_lowercase().starts_with(b"##gff-version") {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    const GFF_WITH_FASTA: &str = "\
##gff-version 3
##sequence-region chr1 1 12
chr1\tdottir\trepeat_region\t3\t8\t.\t+\t.\tName=ALR;family=sat
chr2\tdottir\trepeat_region\t1\t4\t.\t-\t.\tName=HSAT
##FASTA
>chr1 first record
ACGTACGTACGT
>chr2
TTTTGGGG
";

    fn load(name: &str, text: &str) -> Result<SeqInput, InputError> {
        from_bytes(Path::new(name), text.as_bytes().to_vec())
    }

    #[test]
    fn gff3_with_fasta_yields_sequence_and_features() {
        let input = load("a.gff3", GFF_WITH_FASTA).unwrap();
        assert_eq!(input.kind, InputKind::Gff3WithFasta);
        assert_eq!(input.sequence.seq, b"ACGTACGTACGTTTTTGGGG");
        assert_eq!(input.sequence.records.len(), 2);
        assert_eq!(input.sequence.records[0].id, "chr1");
        assert_eq!(input.sequence.records[0].range, 0..12);
        assert_eq!(input.sequence.records[1].range, 12..20);
        assert_eq!(
            input.sequence.records[0].description.as_deref(),
            Some("first record")
        );

        let annot = input.annotations.expect("features");
        assert_eq!(annot.features.len(), 2);
        // GFF3 1-based inclusive [3,8] -> record-local 0-based [2,8).
        assert_eq!(annot.features[0].range, 2..8);
        assert_eq!(annot.features[0].record, "chr1");
        assert_eq!(
            annot.features[0].attrs.get("Name").map(String::as_str),
            Some("ALR")
        );
        assert_eq!(annot.features[1].record, "chr2");
    }

    #[test]
    fn plain_fasta_has_no_annotations() {
        let input = load("a.fasta", ">s1\nACGT\n").unwrap();
        assert_eq!(input.kind, InputKind::Fasta);
        assert!(input.annotations.is_none());
        assert_eq!(input.sequence.seq, b"ACGT");
    }

    #[test]
    fn gff3_without_fasta_section_errors() {
        let gff = "##gff-version 3\nchr1\tsrc\tgene\t1\t10\t.\t+\t.\tName=X\n";
        assert!(matches!(
            load("a.gff3", gff),
            Err(InputError::NoEmbeddedFasta(_))
        ));
    }

    #[test]
    fn gff3_with_empty_fasta_section_errors() {
        let gff = "##gff-version 3\nchr1\tsrc\tgene\t1\t10\t.\t+\t.\tName=X\n##FASTA\n";
        assert!(matches!(
            load("a.gff3", gff),
            Err(InputError::EmptyEmbeddedFasta(_))
        ));
    }

    #[test]
    fn gff3_sniffed_without_gff_extension() {
        let input = load("a.txt", GFF_WITH_FASTA).unwrap();
        assert_eq!(input.kind, InputKind::Gff3WithFasta);
    }

    #[test]
    fn neither_fasta_nor_gff3_names_the_accepted_formats() {
        let err = load("notes.txt", "just some text\nno header\n").unwrap_err();
        assert!(matches!(err, InputError::Unrecognized(_)), "got {err:?}");
        let msg = err.to_string();
        assert!(msg.contains("FASTA"), "{msg}");
        assert!(msg.contains("##FASTA"), "{msg}");
    }

    #[test]
    fn error_display_does_not_repeat_the_path() {
        // Callers prefix the path themselves; a stutter reads badly.
        let err = load("a.fasta", ">a\n\n>b\nACGT\n").unwrap_err();
        assert!(matches!(err, InputError::Fasta { .. }), "got {err:?}");
        assert!(!err.to_string().contains("a.fasta"), "{err}");
    }

    #[test]
    fn fasta_wins_over_sniff_when_file_starts_with_header() {
        let input = load("a.txt", ">s\nACGT\n").unwrap();
        assert_eq!(input.kind, InputKind::Fasta);
    }

    #[test]
    fn gzipped_gff3_with_fasta_round_trips() {
        use std::io::Write;
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(GFF_WITH_FASTA.as_bytes()).unwrap();
        let gz = enc.finish().unwrap();
        let input = from_bytes(Path::new("a.gff3.gz"), gz.clone()).unwrap();
        assert_eq!(input.kind, InputKind::Gff3WithFasta);
        assert_eq!(input.sequence.seq, b"ACGTACGTACGTTTTTGGGG");
        // The sidecar hash is over the on-disk (still compressed) bytes.
        assert_eq!(input.bytes, gz);
    }

    #[test]
    fn features_only_gff3_without_fasta_is_not_a_sequence_input() {
        // The separate-file workflow is unaffected: such a file is still
        // loadable as annotations, it just cannot serve as a sequence.
        let gff = "##gff-version 3\nchr1\tsrc\tgene\t1\t10\t.\t+\t.\tName=X\n";
        let feats = annotation::parse_gff3_bytes(gff.as_bytes()).unwrap();
        assert_eq!(feats.len(), 1);
    }
}
