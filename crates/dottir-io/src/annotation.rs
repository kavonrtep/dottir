//! Annotation interval model for the dotplot overlay (ADR 0005).
//!
//! Two on-disk formats are supported:
//!
//! * **GFF3** — parsed with `noodles-gff`, which handles the real-world
//!   warts (percent-encoded attribute values, multi-value attributes,
//!   `##` directives, embedded `##FASTA` sections). Colored by an
//!   attribute value in the GUI.
//! * **BED** — a fixed-column TSV; hand-rolled here (BED carries no rich
//!   attributes, so there is nothing for noodles to buy us). Single color
//!   per file in the GUI.
//!
//! Both formats produce the same [`Feature`]/[`AnnotSet`] shape. Coordinates
//! are normalized to **0-based, half-open** `[start, end)` *local* to each
//! record (GFF3 seqid / BED chrom), matching the byte indexing used by
//! [`crate::Sequence`]. The GUI maps the record-local range into the
//! concatenated-buffer coordinate space at bind time.

use std::collections::{BTreeMap, BTreeSet};
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use thiserror::Error;

/// Strand of a feature, normalized across formats. `Unknown` covers GFF3
/// `.`/`?` and BED features without a strand column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
    Unknown,
}

/// On-disk format an [`AnnotSet`] was loaded from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnnotSource {
    Gff3,
    Bed,
}

/// A single annotation interval on one record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Feature {
    /// Record this feature belongs to — GFF3 column 1 (seqid) or BED
    /// column 1 (chrom). Matched against [`crate::RecordSpan::id`] at bind
    /// time.
    pub record: String,
    /// 0-based, half-open `[start, end)` range *within the record*.
    pub range: std::ops::Range<usize>,
    pub strand: Strand,
    /// GFF3 column 3 (`gene`, `exon`, …). Empty for BED.
    pub feature_type: String,
    /// GFF3 column 9 key/value attributes. For BED, the optional name
    /// column is stored as `{"Name": <name>}`. Multi-value GFF3 attributes
    /// are joined with `,`.
    pub attrs: BTreeMap<String, String>,
}

impl Feature {
    /// Look up the value used to color this feature for a given key. The
    /// synthetic key [`TYPE_KEY`] returns the feature type (GFF3 column 3).
    /// Returns `None` when the key is absent.
    pub fn color_value(&self, key: &str) -> Option<&str> {
        if key == TYPE_KEY {
            (!self.feature_type.is_empty()).then_some(self.feature_type.as_str())
        } else {
            self.attrs.get(key).map(String::as_str)
        }
    }
}

/// Synthetic attribute key meaning "color by GFF3 feature type (column 3)".
/// Listed alongside the real attribute keys in the GUI dropdown.
pub const TYPE_KEY: &str = "(type)";

/// A loaded set of annotation features.
#[derive(Debug, Clone)]
pub struct AnnotSet {
    pub source_path: Option<PathBuf>,
    pub source: AnnotSource,
    pub features: Vec<Feature>,
}

#[derive(Debug, Error)]
pub enum AnnotError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("GFF3 parse error: {0}")]
    Gff(String),
    #[error("BED parse error on line {line}: {msg}")]
    Bed { line: usize, msg: String },
    #[error("unrecognized annotation extension for {0} (expected .gff/.gff3/.bed, optionally .gz)")]
    UnknownFormat(PathBuf),
}

impl AnnotSet {
    /// Load an annotation file, dispatching on extension. Gzip is detected
    /// transparently from the `.gz` suffix or the gzip magic bytes, exactly
    /// like [`crate::fasta`].
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, AnnotError> {
        let path = path.as_ref();
        let bytes = std::fs::read(path)?;
        let is_gzipped = (bytes.len() >= 2 && bytes[0] == 0x1f && bytes[1] == 0x8b)
            || ext_is(path, "gz");
        // Choose the format from the extension *after* stripping `.gz`.
        let format_path = if ext_is(path, "gz") {
            path.with_extension("")
        } else {
            path.to_path_buf()
        };
        let source = if ext_is(&format_path, "gff") || ext_is(&format_path, "gff3") {
            AnnotSource::Gff3
        } else if ext_is(&format_path, "bed") {
            AnnotSource::Bed
        } else {
            return Err(AnnotError::UnknownFormat(path.to_path_buf()));
        };

        let features = if is_gzipped {
            let r = BufReader::new(MultiGzDecoder::new(&bytes[..]));
            parse(source, r)?
        } else {
            let r = BufReader::new(&bytes[..]);
            parse(source, r)?
        };

        Ok(AnnotSet {
            source_path: Some(path.to_path_buf()),
            source,
            features,
        })
    }

    /// Union of all attribute keys present, plus the synthetic
    /// [`TYPE_KEY`] when any feature has a type. Drives the GUI
    /// "Color by" dropdown. Empty for BED (no attributes, no types).
    pub fn attribute_keys(&self) -> BTreeSet<String> {
        let mut keys = BTreeSet::new();
        let mut any_type = false;
        for f in &self.features {
            keys.extend(f.attrs.keys().cloned());
            any_type |= !f.feature_type.is_empty();
        }
        if any_type {
            keys.insert(TYPE_KEY.to_string());
        }
        keys
    }
}

fn ext_is(path: &Path, ext: &str) -> bool {
    path.extension()
        .and_then(|s| s.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case(ext))
}

fn parse<R: BufRead>(source: AnnotSource, reader: R) -> Result<Vec<Feature>, AnnotError> {
    match source {
        AnnotSource::Gff3 => parse_gff3(reader),
        AnnotSource::Bed => parse_bed(reader),
    }
}

fn parse_gff3<R: BufRead>(reader: R) -> Result<Vec<Feature>, AnnotError> {
    use noodles_gff::feature::record::Strand as GffStrand;

    let mut r = noodles_gff::io::Reader::new(reader);
    let mut out = Vec::new();
    for result in r.record_bufs() {
        let rec = result.map_err(|e| AnnotError::Gff(e.to_string()))?;
        // GFF3 is 1-based inclusive; normalize to 0-based half-open.
        let start = rec.start().get();
        let end = rec.end().get();
        let range = start.saturating_sub(1)..end;

        let strand = match rec.strand() {
            GffStrand::Forward => Strand::Forward,
            GffStrand::Reverse => Strand::Reverse,
            _ => Strand::Unknown,
        };

        let mut attrs = BTreeMap::new();
        for (tag, value) in rec.attributes().as_ref() {
            let key = tag.to_string();
            let val = if let Some(s) = value.as_string() {
                s.to_string()
            } else if let Some(arr) = value.as_array() {
                arr.iter()
                    .map(|b| b.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            } else {
                continue;
            };
            attrs.insert(key, val);
        }

        out.push(Feature {
            record: rec.reference_sequence_name().to_string(),
            range,
            strand,
            feature_type: rec.ty().to_string(),
            attrs,
        });
    }
    Ok(out)
}

fn parse_bed<R: BufRead>(reader: R) -> Result<Vec<Feature>, AnnotError> {
    let mut out = Vec::new();
    for (i, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim_end();
        if line.is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        let mut cols = line.split('\t');
        let lineno = i + 1;
        let bed_err = |msg: &str| AnnotError::Bed {
            line: lineno,
            msg: msg.to_string(),
        };

        let chrom = cols.next().ok_or_else(|| bed_err("missing chrom"))?;
        let start: usize = cols
            .next()
            .ok_or_else(|| bed_err("missing chromStart"))?
            .parse()
            .map_err(|_| bed_err("chromStart not an integer"))?;
        let end: usize = cols
            .next()
            .ok_or_else(|| bed_err("missing chromEnd"))?
            .parse()
            .map_err(|_| bed_err("chromEnd not an integer"))?;
        if end < start {
            return Err(bed_err("chromEnd < chromStart"));
        }

        // Optional columns: name (4), score (5), strand (6).
        let name = cols.next();
        let _score = cols.next();
        let strand = match cols.next() {
            Some("+") => Strand::Forward,
            Some("-") => Strand::Reverse,
            _ => Strand::Unknown,
        };

        let mut attrs = BTreeMap::new();
        if let Some(name) = name.filter(|n| !n.is_empty() && *n != ".") {
            attrs.insert("Name".to_string(), name.to_string());
        }

        out.push(Feature {
            record: chrom.to_string(),
            range: start..end, // BED is already 0-based half-open.
            strand,
            feature_type: String::new(),
            attrs,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn gff3_parses_coords_attrs_and_keys() {
        let gff = "\
##gff-version 3
chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1;Name=ALR\n\
chr1\tsrc\texon\t150\t160\t.\t-\t.\tName=HSAT;family=sat\n";
        let feats = parse_gff3(Cursor::new(gff)).unwrap();
        assert_eq!(feats.len(), 2);
        // 1-based inclusive [101,200] -> 0-based half-open [100,200).
        assert_eq!(feats[0].range, 100..200);
        assert_eq!(feats[0].record, "chr1");
        assert_eq!(feats[0].feature_type, "gene");
        assert_eq!(feats[0].strand, Strand::Forward);
        assert_eq!(feats[0].attrs.get("Name").map(String::as_str), Some("ALR"));
        assert_eq!(feats[1].range, 149..160);
        assert_eq!(feats[1].strand, Strand::Reverse);

        let set = AnnotSet {
            source_path: None,
            source: AnnotSource::Gff3,
            features: feats,
        };
        let keys = set.attribute_keys();
        assert!(keys.contains("Name"));
        assert!(keys.contains("family"));
        assert!(keys.contains("ID"));
        assert!(keys.contains(TYPE_KEY));
    }

    #[test]
    fn color_value_uses_type_key_and_attrs() {
        let f = Feature {
            record: "c".into(),
            range: 0..10,
            strand: Strand::Unknown,
            feature_type: "exon".into(),
            attrs: BTreeMap::from([("Name".to_string(), "ALR".to_string())]),
        };
        assert_eq!(f.color_value(TYPE_KEY), Some("exon"));
        assert_eq!(f.color_value("Name"), Some("ALR"));
        assert_eq!(f.color_value("missing"), None);
    }

    #[test]
    fn bed_parses_zero_based_halfopen_and_name() {
        let bed = "\
# a comment
track name=foo
chr1\t100\t200\tALR\t0\t+
chr2\t5\t8\t.\t.\t-
chr3\t0\t4
";
        let feats = parse_bed(Cursor::new(bed)).unwrap();
        assert_eq!(feats.len(), 3);
        assert_eq!(feats[0].range, 100..200); // BED already half-open.
        assert_eq!(feats[0].record, "chr1");
        assert_eq!(feats[0].strand, Strand::Forward);
        assert_eq!(feats[0].attrs.get("Name").map(String::as_str), Some("ALR"));
        assert_eq!(feats[0].feature_type, ""); // BED has no type.
        // "." name is dropped.
        assert!(feats[1].attrs.is_empty());
        assert_eq!(feats[1].strand, Strand::Reverse);
        // Minimal 3-column line.
        assert_eq!(feats[2].range, 0..4);
        assert_eq!(feats[2].strand, Strand::Unknown);

        // BED has no feature types, so no synthetic TYPE_KEY. The only key
        // is the optional Name column (the GUI ignores it and colors BED
        // sources with a single color regardless).
        let set = AnnotSet {
            source_path: None,
            source: AnnotSource::Bed,
            features: feats,
        };
        let keys = set.attribute_keys();
        assert!(!keys.contains(TYPE_KEY));
        assert_eq!(keys.into_iter().collect::<Vec<_>>(), vec!["Name".to_string()]);
    }

    #[test]
    fn bed_rejects_bad_coords() {
        let bad = "chr1\t200\t100\n";
        assert!(matches!(
            parse_bed(Cursor::new(bad)),
            Err(AnnotError::Bed { line: 1, .. })
        ));
        let nonint = "chr1\tx\t100\n";
        assert!(matches!(
            parse_bed(Cursor::new(nonint)),
            Err(AnnotError::Bed { line: 1, .. })
        ));
    }
}
