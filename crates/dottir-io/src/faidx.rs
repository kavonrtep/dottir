//! FASTA index (samtools `.fai`) builder and region reader, so genome
//! mode can work with record names and lengths and read only the regions
//! it draws (spec §4.7.8).
//!
//! Gzipped input (bgzip included) is decompressed once into the cache
//! directory by [`plain_fasta`]; this avoids a BGZF dependency at the cost
//! of disk space.

use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use dottir_core::genome::GenomeRecord;

use crate::fasta::FastaRecord;
use crate::sequence::{RegionOrigin, Sequence};
use sha2::{Digest, Sha256};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum FaidxError {
    #[error("I/O error on {path}: {source}")]
    Io {
        path: PathBuf,
        source: std::io::Error,
    },
    #[error("{path} is gzipped; pass it through plain_fasta() first")]
    Gzipped { path: PathBuf },
    #[error("{path}, line {line}: {msg}")]
    Format {
        path: PathBuf,
        line: u64,
        msg: String,
    },
    #[error("bad index {path}, line {line}")]
    BadIndex { path: PathBuf, line: u64 },
    #[error("region {name}:{start}-{end} is outside the record (length {len})")]
    Range {
        name: String,
        start: u64,
        end: u64,
        len: u64,
    },
    #[error("region {spec:?}: {msg}")]
    Spec { spec: String, msg: String },
}

fn io_err(path: &Path) -> impl FnOnce(std::io::Error) -> FaidxError + '_ {
    move |source| FaidxError::Io {
        path: path.to_path_buf(),
        source,
    }
}

/// One `.fai` line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FaiEntry {
    pub name: String,
    pub len: u64,
    /// Byte offset of the first base.
    pub offset: u64,
    pub line_bases: u64,
    /// Bytes per line, newline included.
    pub line_bytes: u64,
}

/// Scan a plain FASTA file and build its index. Like samtools, every
/// sequence line of a record except the last must have the same length.
pub fn build_fai(path: &Path) -> Result<Vec<FaiEntry>, FaidxError> {
    let mut r = BufReader::with_capacity(1 << 20, std::fs::File::open(path).map_err(io_err(path))?);
    if r.fill_buf()
        .map_err(io_err(path))?
        .starts_with(&[0x1f, 0x8b])
    {
        return Err(FaidxError::Gzipped {
            path: path.to_path_buf(),
        });
    }
    let fmt = |line: u64, msg: &str| FaidxError::Format {
        path: path.to_path_buf(),
        line,
        msg: msg.to_string(),
    };
    let mut out: Vec<FaiEntry> = Vec::new();
    let mut buf = Vec::new();
    let mut pos = 0u64;
    let mut line = 0u64;
    // State of the current record: a short (last) line or a blank line
    // has been seen, so any further sequence line is an error.
    let mut ended = false;
    loop {
        buf.clear();
        let n = r.read_until(b'\n', &mut buf).map_err(io_err(path))? as u64;
        if n == 0 {
            break;
        }
        line += 1;
        let start = pos;
        pos += n;
        let content = buf.strip_suffix(b"\n").unwrap_or(&buf);
        let content = content.strip_suffix(b"\r").unwrap_or(content);
        if let Some(header) = content.strip_prefix(b">") {
            let name = header
                .split(|b| b.is_ascii_whitespace())
                .next()
                .unwrap_or_default();
            if name.is_empty() {
                return Err(fmt(line, "empty record name"));
            }
            out.push(FaiEntry {
                name: String::from_utf8_lossy(name).into_owned(),
                len: 0,
                offset: pos,
                line_bases: 0,
                line_bytes: 0,
            });
            ended = false;
            continue;
        }
        let Some(e) = out.last_mut() else {
            if content.is_empty() {
                continue;
            }
            return Err(fmt(line, "sequence before the first '>' header"));
        };
        let bases = content.len() as u64;
        if bases == 0 {
            ended = true;
            continue;
        }
        if ended {
            return Err(fmt(line, "line lengths differ within the record"));
        }
        if e.line_bases == 0 {
            e.offset = start;
            e.line_bases = bases;
            e.line_bytes = n;
        } else if bases > e.line_bases
            // A full-width last line may lack its terminator.
            || (bases == e.line_bases && n != e.line_bytes && buf.ends_with(b"\n"))
        {
            return Err(fmt(line, "line lengths differ within the record"));
        } else if bases < e.line_bases {
            ended = true;
        }
        e.len += bases;
    }
    Ok(out)
}

/// Parse `.fai` text.
pub fn parse_fai(text: &str, path: &Path) -> Result<Vec<FaiEntry>, FaidxError> {
    text.lines()
        .enumerate()
        .filter(|(_, l)| !l.is_empty())
        .map(|(i, l)| {
            let bad = || FaidxError::BadIndex {
                path: path.to_path_buf(),
                line: i as u64 + 1,
            };
            let c: Vec<&str> = l.split('\t').collect();
            if c.len() < 5 {
                return Err(bad());
            }
            let n = |k: usize| c[k].parse::<u64>().map_err(|_| bad());
            Ok(FaiEntry {
                name: c[0].to_string(),
                len: n(1)?,
                offset: n(2)?,
                line_bases: n(3)?,
                line_bytes: n(4)?,
            })
        })
        .collect()
}

/// Write `.fai` text.
pub fn write_fai<W: Write>(mut w: W, entries: &[FaiEntry]) -> std::io::Result<()> {
    for e in entries {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}",
            e.name, e.len, e.offset, e.line_bases, e.line_bytes
        )?;
    }
    w.flush()
}

fn is_gzipped(path: &Path) -> Result<bool, FaidxError> {
    let mut magic = [0u8; 2];
    let mut f = std::fs::File::open(path).map_err(io_err(path))?;
    let n = f.read(&mut magic).map_err(io_err(path))?;
    Ok(n == 2 && magic == [0x1f, 0x8b])
}

/// Cache key for a file: canonical path, size and modification time.
fn file_key(path: &Path) -> Result<String, FaidxError> {
    let meta = std::fs::metadata(path).map_err(io_err(path))?;
    let canon = std::fs::canonicalize(path).map_err(io_err(path))?;
    let mtime = meta
        .modified()
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map_or(0, |d| d.as_nanos());
    let mut h = Sha256::new();
    h.update(canon.to_string_lossy().as_bytes());
    h.update(meta.len().to_le_bytes());
    h.update(mtime.to_le_bytes());
    Ok(h.finalize()
        .iter()
        .take(8)
        .map(|b| format!("{b:02x}"))
        .collect())
}

/// A plain FASTA path for `path`: `path` itself, or for gzipped input a
/// decompressed copy in `cache_dir` (made once, reused while the input's
/// path, size and modification time are unchanged).
pub fn plain_fasta(path: &Path, cache_dir: &Path) -> Result<PathBuf, FaidxError> {
    if !is_gzipped(path)? {
        return Ok(path.to_path_buf());
    }
    let stem = path
        .file_name()
        .map(|n| n.to_string_lossy().trim_end_matches(".gz").to_string())
        .unwrap_or_default();
    let dir = cache_dir.join("fasta");
    let out = dir.join(format!("{}-{stem}", file_key(path)?));
    if out.exists() {
        return Ok(out);
    }
    std::fs::create_dir_all(&dir).map_err(io_err(&dir))?;
    let tmp = out.with_extension(format!("tmp{}", std::process::id()));
    {
        let mut src = crate::paf::open_maybe_gz(path).map_err(io_err(path))?;
        let mut dst = std::io::BufWriter::new(std::fs::File::create(&tmp).map_err(io_err(&tmp))?);
        std::io::copy(&mut src, &mut dst).map_err(io_err(path))?;
        dst.flush().map_err(io_err(&tmp))?;
    }
    std::fs::rename(&tmp, &out).map_err(io_err(&out))?;
    Ok(out)
}

/// A region as text: `name:start-end`, 1-based inclusive like samtools,
/// with `:-` appended for the reverse complement.
pub fn region_spec(name: &str, range: &std::ops::Range<u64>, reverse: bool) -> String {
    let strand = if reverse { ":-" } else { "" };
    format!("{name}:{}-{}{strand}", range.start + 1, range.end)
}

/// The FASTA at `path` (plain or gzipped) cut to `specs` (see
/// [`region_spec`]; a bare record name means the whole record), as one
/// sequence whose records carry their [`RegionOrigin`].
pub fn load_regions(
    path: &Path,
    specs: &[String],
    cache_dir: &Path,
) -> Result<Sequence, FaidxError> {
    let plain = plain_fasta(path, cache_dir)?;
    let index = FastaIndex::open(&plain, cache_dir)?;
    let regions = specs
        .iter()
        .map(|s| index.parse_region(s))
        .collect::<Result<Vec<_>, _>>()?;
    let mut seq = index.region_sequence(&regions)?;
    seq.source_path = Some(path.to_path_buf());
    Ok(seq)
}

/// A region of an indexed record: `range` is 0-based, end-exclusive.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    pub rec: usize,
    pub range: std::ops::Range<u64>,
    pub reverse: bool,
}

/// An indexed plain FASTA file.
#[derive(Debug, Clone)]
pub struct FastaIndex {
    pub path: PathBuf,
    pub entries: Vec<FaiEntry>,
    by_name: BTreeMap<String, usize>,
}

impl FastaIndex {
    /// Index a plain FASTA file: reuse `path.fai` when it is at least as
    /// new as the FASTA, else build it and write it next to the FASTA, or
    /// to `cache_dir` when that directory is not writable.
    pub fn open(path: &Path, cache_dir: &Path) -> Result<Self, FaidxError> {
        let beside = PathBuf::from(format!("{}.fai", path.display()));
        let cached = cache_dir
            .join("fai")
            .join(format!("{}.fai", file_key(path)?));
        let mtime = |p: &Path| std::fs::metadata(p).and_then(|m| m.modified()).ok();
        let fasta_time = mtime(path);
        for idx in [&beside, &cached] {
            if mtime(idx).is_some() && mtime(idx) >= fasta_time {
                let text = std::fs::read_to_string(idx).map_err(io_err(idx))?;
                return Ok(Self::from_entries(path, parse_fai(&text, idx)?));
            }
        }
        let entries = build_fai(path)?;
        let write = |p: &Path| -> std::io::Result<()> {
            if let Some(d) = p.parent() {
                std::fs::create_dir_all(d)?;
            }
            write_fai(std::io::BufWriter::new(std::fs::File::create(p)?), &entries)
        };
        if write(&beside).is_err() {
            write(&cached).map_err(io_err(&cached))?;
        }
        Ok(Self::from_entries(path, entries))
    }

    pub fn from_entries(path: &Path, entries: Vec<FaiEntry>) -> Self {
        let by_name = entries
            .iter()
            .enumerate()
            .map(|(i, e)| (e.name.clone(), i))
            .collect();
        Self {
            path: path.to_path_buf(),
            entries,
            by_name,
        }
    }

    /// Record names and lengths, in file order.
    pub fn records(&self) -> Vec<GenomeRecord> {
        self.entries
            .iter()
            .map(|e| GenomeRecord {
                name: e.name.clone(),
                len: e.len,
            })
            .collect()
    }

    pub fn find(&self, name: &str) -> Option<usize> {
        self.by_name.get(name).copied()
    }

    /// Parse a [`region_spec`] against this index.
    pub fn parse_region(&self, spec: &str) -> Result<Region, FaidxError> {
        let bad = |msg: &str| FaidxError::Spec {
            spec: spec.to_string(),
            msg: msg.to_string(),
        };
        let (body, reverse) = match spec.strip_suffix(":-") {
            Some(b) => (b, true),
            None => (spec.strip_suffix(":+").unwrap_or(spec), false),
        };
        // A record name may itself contain ':'; try the whole text first.
        if let Some(rec) = self.find(body) {
            return Ok(Region {
                rec,
                range: 0..self.entries[rec].len,
                reverse,
            });
        }
        let (name, range) = body.rsplit_once(':').ok_or_else(|| bad("no such record"))?;
        let rec = self.find(name).ok_or_else(|| bad("no such record"))?;
        let (a, b) = range
            .replace(',', "")
            .split_once('-')
            .and_then(|(a, b)| Some((a.parse::<u64>().ok()?, b.parse::<u64>().ok()?)))
            .ok_or_else(|| bad("expected NAME:START-END"))?;
        let len = self.entries[rec].len;
        if a == 0 || a > b || b > len {
            return Err(bad(&format!("outside 1..{len} or empty")));
        }
        Ok(Region {
            rec,
            range: a - 1..b,
            reverse,
        })
    }

    /// The given regions as one multi-record [`Sequence`], reverse-
    /// complemented where asked, soft-masking kept. Records are named
    /// `name:start-end` (1-based) and carry their [`RegionOrigin`].
    pub fn region_sequence(&self, regions: &[Region]) -> Result<Sequence, FaidxError> {
        let mut records = Vec::with_capacity(regions.len());
        for r in regions {
            let mut raw = self.fetch_raw(r.rec, r.range.clone())?;
            if r.reverse {
                raw.reverse();
                for b in &mut raw {
                    *b = crate::genome_export::complement(*b);
                }
            }
            let mut soft_masked: Vec<std::ops::Range<usize>> = Vec::new();
            for (i, b) in raw.iter_mut().enumerate() {
                if b.is_ascii_lowercase() {
                    match soft_masked.last_mut() {
                        Some(m) if m.end == i => m.end += 1,
                        _ => soft_masked.push(i..i + 1),
                    }
                    b.make_ascii_uppercase();
                }
            }
            let name = &self.entries[r.rec].name;
            records.push(FastaRecord {
                id: format!("{name}:{}-{}", r.range.start + 1, r.range.end),
                description: r.reverse.then(|| "reverse complement".to_string()),
                sequence: raw,
                soft_masked,
            });
        }
        let mut seq = Sequence::from_records(records, Some(self.path.clone()));
        for (span, r) in seq.records.iter_mut().zip(regions) {
            span.origin = Some(RegionOrigin {
                record: self.entries[r.rec].name.clone(),
                record_len: self.entries[r.rec].len,
                range: r.range.clone(),
                reverse: r.reverse,
            });
        }
        Ok(seq)
    }

    /// Residues `range` of record `rec`, uppercased.
    pub fn fetch(&self, rec: usize, range: std::ops::Range<u64>) -> Result<Vec<u8>, FaidxError> {
        let mut seq = self.fetch_raw(rec, range)?;
        seq.make_ascii_uppercase();
        Ok(seq)
    }

    /// Residues `range` of record `rec` as stored (case kept).
    pub fn fetch_raw(
        &self,
        rec: usize,
        range: std::ops::Range<u64>,
    ) -> Result<Vec<u8>, FaidxError> {
        let e = &self.entries[rec];
        if range.start > range.end || range.end > e.len {
            return Err(FaidxError::Range {
                name: e.name.clone(),
                start: range.start,
                end: range.end,
                len: e.len,
            });
        }
        if range.is_empty() {
            return Ok(Vec::new());
        }
        let byte = |i: u64| e.offset + i / e.line_bases * e.line_bytes + i % e.line_bases;
        let (from, to) = (byte(range.start), byte(range.end - 1) + 1);
        let mut f = std::fs::File::open(&self.path).map_err(io_err(&self.path))?;
        f.seek(SeekFrom::Start(from)).map_err(io_err(&self.path))?;
        let mut raw = vec![0u8; (to - from) as usize];
        f.read_exact(&mut raw).map_err(io_err(&self.path))?;
        raw.retain(|b| !b.is_ascii_whitespace());
        Ok(raw)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp(name: &str) -> PathBuf {
        let d = std::env::temp_dir().join(format!("dottir-faidx-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    // samtools faidx output for this file:
    // r1 10 9 4 5 / r2 6 27 6 8
    const FA: &str = ">r1 desc\nACGT\nacgt\nAC\n>r2\r\nGGGTTT\r\n";

    #[test]
    fn region_specs_round_trip_and_load() {
        let d = tmp("spec");
        let p = d.join("a.fa");
        std::fs::write(&p, ">c:1 x\nACGTACGTAA\n>r2\nGGGTTT\n").unwrap();
        let cache = d.join("cache");
        let idx = FastaIndex::open(&p, &cache).unwrap();
        let r = idx
            .parse_region(&region_spec("c:1", &(2..6), true))
            .unwrap();
        assert_eq!(
            r,
            Region {
                rec: 0,
                range: 2..6,
                reverse: true
            }
        );
        assert_eq!(idx.parse_region("c:1").unwrap().range, 0..10);
        assert_eq!(idx.parse_region("r2:2-3").unwrap().range, 1..3);
        for bad in ["zz:1-2", "r2:0-2", "r2:3-2", "r2:1-7", "r2:a-b"] {
            assert!(idx.parse_region(bad).is_err(), "{bad}");
        }
        let seq = load_regions(&p, &["c:1:3-6:-".into(), "r2".into()], &cache).unwrap();
        assert_eq!(seq.bytes(), b"GTACGGGTTT");
        let o = seq.records[0].origin.as_ref().unwrap();
        assert_eq!(
            (o.record.as_str(), o.range.clone(), o.reverse),
            ("c:1", 2..6, true)
        );
        assert_eq!(seq.source_path.as_deref(), Some(p.as_path()));
    }

    #[test]
    fn full_last_line_without_newline() {
        // Review 2026-10-28 #5.
        let d = tmp("eof");
        let p = d.join("a.fa");
        std::fs::write(&p, ">r\nACGT\nACGT").unwrap();
        let idx = FastaIndex::open(&p, &d.join("cache")).unwrap();
        assert_eq!(idx.entries[0].len, 8);
        assert_eq!(idx.fetch(0, 2..8).unwrap(), b"GTACGT");
        // Inconsistent widths inside a record are still rejected.
        let bad = d.join("b.fa");
        std::fs::write(&bad, ">r\nACGT\nACG\nACGT\n").unwrap();
        assert!(FastaIndex::open(&bad, &d.join("cache")).is_err());
    }

    #[test]
    fn index_matches_samtools_and_regions_cross_lines() {
        let d = tmp("basic");
        let p = d.join("a.fa");
        std::fs::write(&p, FA).unwrap();
        let idx = FastaIndex::open(&p, &d.join("cache")).unwrap();
        let mut text = Vec::new();
        write_fai(&mut text, &idx.entries).unwrap();
        assert_eq!(
            String::from_utf8(text).unwrap(),
            "r1\t10\t9\t4\t5\nr2\t6\t27\t6\t8\n"
        );
        assert_eq!(idx.fetch(0, 2..9).unwrap(), b"GTACGTA");
        assert_eq!(idx.fetch(0, 0..10).unwrap(), b"ACGTACGTAC");
        assert_eq!(idx.fetch(1, 1..6).unwrap(), b"GGTTT");
        assert!(idx.fetch(1, 0..7).is_err());
        // The index was written beside the FASTA and is reused.
        assert!(d.join("a.fa.fai").exists());
        let again = FastaIndex::open(&p, &d.join("cache")).unwrap();
        assert_eq!(again.entries, idx.entries);
        assert_eq!(again.records()[0].len, 10);
    }

    #[test]
    fn region_sequence_reverse_complements_and_keeps_soft_masking() {
        let d = tmp("region");
        let p = d.join("a.fa");
        std::fs::write(&p, FA).unwrap();
        let idx = FastaIndex::open(&p, &d.join("cache")).unwrap();
        let seq = idx
            .region_sequence(&[
                Region {
                    rec: 0,
                    range: 2..7,
                    reverse: false,
                },
                Region {
                    rec: 0,
                    range: 2..7,
                    reverse: true,
                },
                Region {
                    rec: 1,
                    range: 0..3,
                    reverse: false,
                },
            ])
            .unwrap();
        // r1[2..7] = GTacg; its reverse complement = cgtAC.
        assert_eq!(seq.seq, b"GTACGCGTACGGG");
        assert_eq!(seq.soft_mask, vec![2..5, 5..8]);
        let ids: Vec<&str> = seq.records.iter().map(|r| r.id.as_str()).collect();
        assert_eq!(ids, ["r1:3-7", "r1:3-7", "r2:1-3"]);
        let o = seq.records[1].origin.as_ref().unwrap();
        assert_eq!((o.record.as_str(), o.reverse), ("r1", true));
        // Region position 0 of the reversed record is source base 6.
        assert_eq!((o.source_pos(0), o.source_pos(4)), (6, 2));
        assert_eq!(seq.records[0].origin.as_ref().unwrap().source_pos(1), 3);
    }

    #[test]
    fn ragged_lines_are_rejected() {
        let d = tmp("ragged");
        let p = d.join("a.fa");
        std::fs::write(&p, ">r\nACGT\nAC\nACGT\n").unwrap();
        assert!(matches!(
            build_fai(&p),
            Err(FaidxError::Format { line: 4, .. })
        ));
        std::fs::write(&p, ">r\nACGT\nACGTA\n").unwrap();
        assert!(matches!(
            build_fai(&p),
            Err(FaidxError::Format { line: 3, .. })
        ));
    }

    #[test]
    fn gzipped_input_is_decompressed_once_into_the_cache() {
        let d = tmp("gz");
        let p = d.join("a.fa.gz");
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        gz.write_all(FA.as_bytes()).unwrap();
        std::fs::write(&p, gz.finish().unwrap()).unwrap();
        assert!(matches!(build_fai(&p), Err(FaidxError::Gzipped { .. })));
        let plain = plain_fasta(&p, &d.join("cache")).unwrap();
        assert_ne!(plain, p);
        assert_eq!(std::fs::read_to_string(&plain).unwrap(), FA);
        assert_eq!(plain_fasta(&p, &d.join("cache")).unwrap(), plain);
        let idx = FastaIndex::open(&plain, &d.join("cache")).unwrap();
        assert_eq!(idx.fetch(1, 0..3).unwrap(), b"GGG");
    }
}
