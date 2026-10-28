//! Running an external minimap2-compatible aligner for genome mode
//! (ADR 0008 §2–3): detection, chunked query streaming, joining across
//! chunk cuts, and a cache keyed by inputs, aligner version and arguments.
//!
//! Pairwise runs cut query records into pieces of at most `chunk_size`
//! and stream them to the aligner's stdin under synthetic names, so long
//! chromosomes never reach the aligner whole (D-GENIES: this is what
//! bounds its memory). Self-comparisons send whole records under their
//! own names, because `-D` only skips the diagonal when query and target
//! names and lengths agree.

use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Mutex;

use dottir_core::genome::alignment::join_at_cuts;
use dottir_core::genome::GenomeRecord;
use sha2::{Digest, Sha256};
use thiserror::Error;

use crate::faidx::{FaidxError, FastaIndex};
use crate::paf::{read_paf, write_paf, PafError};

/// Aligners tried by `auto`, in order.
pub const AUTO_ALIGNERS: [&str; 2] = ["mm2plus", "minimap2"];
/// minimap2 cannot index a single sequence of this length or more.
pub const MAX_RECORD_LEN: u64 = 1 << 31;
/// Alignments ending this close to a chunk cut may be joined across it.
pub const JOIN_TOLERANCE: u64 = 2000;

#[derive(Debug, Error)]
pub enum AlignerError {
    #[error(
        "no aligner found on PATH (tried {tried}); install mm2plus or minimap2 \
         (e.g. from bioconda), pass --aligner PATH, or give a PAF with --paf"
    )]
    NotFound { tried: String },
    #[error("cannot run {program}: {source}")]
    Spawn {
        program: PathBuf,
        source: std::io::Error,
    },
    #[error("{program} failed ({status}):\n{stderr}")]
    Failed {
        program: PathBuf,
        status: String,
        stderr: String,
    },
    #[error("record {name} is {len} bp; minimap2 cannot index sequences of 2^31 bp or more")]
    RecordTooLong { name: String, len: u64 },
    #[error("alignment cancelled")]
    Cancelled,
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Faidx(#[from] FaidxError),
    #[error("reading aligner output: {0}")]
    Paf(#[from] PafError),
}

/// A detected aligner.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Aligner {
    pub program: PathBuf,
    /// File name without extension, e.g. `mm2plus`.
    pub name: String,
    /// `--version` output.
    pub version: String,
}

fn executable(p: &Path) -> bool {
    let Ok(meta) = std::fs::metadata(p) else {
        return false;
    };
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        meta.is_file() && meta.permissions().mode() & 0o111 != 0
    }
    #[cfg(not(unix))]
    {
        meta.is_file()
    }
}

/// Find `name` on `PATH`.
pub fn find_on_path(name: &str) -> Option<PathBuf> {
    let path = std::env::var_os("PATH")?;
    let names: Vec<String> = if cfg!(windows) {
        vec![format!("{name}.exe"), name.to_string()]
    } else {
        vec![name.to_string()]
    };
    std::env::split_paths(&path)
        .flat_map(|d| names.iter().map(move |n| d.join(n)))
        .find(|p| executable(p))
}

impl Aligner {
    /// `auto` (mm2plus, else minimap2), a name looked up on `PATH`, or a
    /// path to a minimap2-compatible binary.
    pub fn detect(choice: &str) -> Result<Self, AlignerError> {
        let program = if choice == "auto" {
            AUTO_ALIGNERS.iter().find_map(|n| find_on_path(n))
        } else if choice.contains(std::path::MAIN_SEPARATOR) || Path::new(choice).exists() {
            Some(PathBuf::from(choice)).filter(|p| executable(p))
        } else {
            find_on_path(choice)
        };
        let program = program.ok_or_else(|| AlignerError::NotFound {
            tried: if choice == "auto" {
                AUTO_ALIGNERS.join(", ")
            } else {
                choice.to_string()
            },
        })?;
        Self::at(program)
    }

    /// The aligner at `program`; asks it for its version.
    pub fn at(program: PathBuf) -> Result<Self, AlignerError> {
        let out = Command::new(&program)
            .arg("--version")
            .output()
            .map_err(|source| AlignerError::Spawn {
                program: program.clone(),
                source,
            })?;
        let last_line = |b: &[u8]| {
            String::from_utf8_lossy(b)
                .lines()
                .map(str::trim)
                .rfind(|l| !l.is_empty())
                .unwrap_or_default()
                .to_string()
        };
        let mut version = last_line(&out.stdout);
        if version.is_empty() {
            version = last_line(&out.stderr);
        }
        let name = program
            .file_stem()
            .map_or_else(String::new, |s| s.to_string_lossy().into_owned());
        Ok(Self {
            program,
            name,
            version,
        })
    }
}

/// What to align and how.
#[derive(Debug, Clone)]
pub struct AlignRequest<'a> {
    /// Horizontal axis; the aligner's target.
    pub reference: &'a FastaIndex,
    /// Vertical axis; the aligner's query. Same index for self.
    pub query: &'a FastaIndex,
    /// SHA-256 of the original input files (cache key).
    pub reference_sha256: String,
    pub query_sha256: String,
    /// `-x` preset, e.g. `asm20`.
    pub preset: String,
    pub extra_args: Vec<String>,
    pub threads: usize,
    pub chunk_size: u64,
    pub self_comparison: bool,
    /// Keep every chain (`-P`): all copies of repeats, not only the best
    /// placement of each region.
    pub all_copies: bool,
}

impl AlignRequest<'_> {
    /// Arguments that shape the output (no paths, no thread count).
    pub fn output_args(&self) -> Vec<String> {
        let mut a = vec!["-x".to_string(), self.preset.clone()];
        if self.self_comparison {
            a.extend(["-D", "--dual=no"].map(String::from));
        }
        if self.all_copies {
            a.push("-P".into());
        }
        a.extend(self.extra_args.iter().cloned());
        a
    }

    /// Cache key: inputs, aligner identity, output-shaping arguments and
    /// chunking.
    pub fn cache_key(&self, aligner: &Aligner) -> String {
        let mut h = Sha256::new();
        for part in [
            "dottir-genome-v1",
            &self.reference_sha256,
            &self.query_sha256,
            &aligner.name,
            &aligner.version,
            &self.output_args().join("\u{1f}"),
            &self.chunk_size.to_string(),
            if self.self_comparison { "self" } else { "pair" },
        ] {
            h.update(part.as_bytes());
            h.update([0]);
        }
        h.finalize().iter().map(|b| format!("{b:02x}")).collect()
    }
}

/// A piece of a query record sent to the aligner.
#[derive(Debug, Clone, PartialEq, Eq)]
struct Piece {
    rec: u32,
    offset: u64,
    len: u64,
}

/// Cut records into pieces of at most `chunk` bp (whole records when
/// `split` is false).
fn pieces(records: &[GenomeRecord], chunk: u64, split: bool) -> Vec<Piece> {
    let chunk = chunk.max(1);
    let mut out = Vec::new();
    for (rec, r) in records.iter().enumerate() {
        let step = if split { chunk } else { r.len.max(1) };
        let mut offset = 0;
        while offset < r.len {
            let len = step.min(r.len - offset);
            out.push(Piece {
                rec: rec as u32,
                offset,
                len,
            });
            offset += len;
        }
    }
    out
}

fn piece_name(i: usize) -> String {
    format!("dottir_p{i}")
}

/// Rewrite one PAF line from piece coordinates to record coordinates.
fn remap_line(line: &str, pieces: &[Piece], records: &[GenomeRecord]) -> Option<String> {
    let mut cols: Vec<&str> = line.split('\t').collect();
    let i: usize = cols.first()?.strip_prefix("dottir_p")?.parse().ok()?;
    let p = pieces.get(i)?;
    let start: u64 = cols.get(2)?.parse().ok()?;
    let end: u64 = cols.get(3)?.parse().ok()?;
    let r = &records[p.rec as usize];
    let (len, s, e) = (
        r.len.to_string(),
        (start + p.offset).to_string(),
        (end + p.offset).to_string(),
    );
    cols[0] = &r.name;
    cols[1] = &len;
    cols[2] = &s;
    cols[3] = &e;
    Some(cols.join("\t"))
}

/// The output of [`align`].
#[derive(Debug, Clone)]
pub struct AlignOutput {
    /// Joined PAF (gzipped), query = vertical axis.
    pub paf: PathBuf,
    pub cached: bool,
    /// The aligner invocations, for provenance.
    pub commands: Vec<Vec<String>>,
}

fn check_lengths(idx: &FastaIndex) -> Result<(), AlignerError> {
    match idx.entries.iter().find(|e| e.len >= MAX_RECORD_LEN) {
        Some(e) => Err(AlignerError::RecordTooLong {
            name: e.name.clone(),
            len: e.len,
        }),
        None => Ok(()),
    }
}

fn run_checked(cmd: &mut Command, program: &Path) -> Result<(), AlignerError> {
    let out = cmd.output().map_err(|source| AlignerError::Spawn {
        program: program.to_path_buf(),
        source,
    })?;
    if out.status.success() {
        return Ok(());
    }
    Err(AlignerError::Failed {
        program: program.to_path_buf(),
        status: out.status.to_string(),
        stderr: tail(&String::from_utf8_lossy(&out.stderr), 20),
    })
}

fn tail(s: &str, n: usize) -> String {
    let lines: Vec<&str> = s.lines().collect();
    lines[lines.len().saturating_sub(n)..].join("\n")
}

/// Align `req` with `aligner`, or reuse the cached result. `progress` is
/// called with (pieces mapped, pieces total) as the aligner reports them;
/// setting `cancel` stops the run.
pub fn align(
    aligner: &Aligner,
    req: &AlignRequest<'_>,
    cache_dir: &Path,
    use_cache: bool,
    progress: &(dyn Fn(u64, u64) + Sync),
    cancel: &AtomicBool,
) -> Result<AlignOutput, AlignerError> {
    let dir = cache_dir.join("genome");
    std::fs::create_dir_all(&dir)?;
    let key = req.cache_key(aligner);
    let paf = dir.join(format!("{key}.paf.gz"));
    if use_cache && paf.exists() {
        return Ok(AlignOutput {
            paf,
            cached: true,
            commands: Vec::new(),
        });
    }
    check_lengths(req.reference)?;
    check_lengths(req.query)?;
    let threads = req.threads.max(1).to_string();
    let program = &aligner.program;

    // Index the reference once.
    let mmi = dir.join(format!("{key}.mmi.tmp{}", std::process::id()));
    let mut index_cmd = vec![program.display().to_string()];
    index_cmd.extend([
        "-x".to_string(),
        req.preset.clone(),
        "-t".into(),
        threads.clone(),
    ]);
    index_cmd.extend(["-d".to_string(), mmi.display().to_string()]);
    index_cmd.push(req.reference.path.display().to_string());
    run_checked(
        Command::new(program)
            .args(&index_cmd[1..])
            .stdout(Stdio::null())
            .stderr(Stdio::piped()),
        program,
    )?;

    let q_recs = req.query.records();
    let t_recs = req.reference.records();
    let split = !req.self_comparison;
    let pcs = pieces(&q_recs, req.chunk_size, split);
    let mut map_cmd = vec![program.display().to_string()];
    map_cmd.extend(req.output_args());
    map_cmd.extend(["-t".to_string(), threads]);
    map_cmd.push(mmi.display().to_string());
    map_cmd.push("-".into());

    let result = map_pieces(program, &map_cmd, req, &pcs, &q_recs, progress, cancel);
    let _ = std::fs::remove_file(&mmi);
    let text = result?;

    let alns = read_paf(text.as_bytes(), &q_recs, &t_recs)?;
    let mut cuts: Vec<Vec<u64>> = vec![Vec::new(); q_recs.len()];
    for p in pcs.iter().filter(|p| p.offset > 0) {
        cuts[p.rec as usize].push(p.offset);
    }
    let alns = join_at_cuts(alns, &cuts, JOIN_TOLERANCE);
    let tmp = paf.with_extension(format!("tmp{}", std::process::id()));
    {
        let mut gz = flate2::write::GzEncoder::new(
            std::io::BufWriter::new(std::fs::File::create(&tmp)?),
            flate2::Compression::fast(),
        );
        write_paf(&mut gz, &alns, &q_recs, &t_recs)?;
        gz.finish()?.flush()?;
    }
    std::fs::rename(&tmp, &paf)?;
    Ok(AlignOutput {
        paf,
        cached: false,
        commands: vec![index_cmd, map_cmd],
    })
}

/// Stream the pieces into the aligner and collect its PAF, remapped to
/// record coordinates.
fn map_pieces(
    program: &Path,
    argv: &[String],
    req: &AlignRequest<'_>,
    pcs: &[Piece],
    q_recs: &[GenomeRecord],
    progress: &(dyn Fn(u64, u64) + Sync),
    cancel: &AtomicBool,
) -> Result<String, AlignerError> {
    let mut child = Command::new(program)
        .args(&argv[1..])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|source| AlignerError::Spawn {
            program: program.to_path_buf(),
            source,
        })?;
    let mut stdin = child.stdin.take().expect("piped stdin");
    let stdout = child.stdout.take().expect("piped stdout");
    let stderr = child.stderr.take().expect("piped stderr");
    let total = pcs.len() as u64;
    let stderr_tail = Mutex::new(Vec::<String>::new());
    let feed_error = Mutex::new(None::<FaidxError>);

    let text = std::thread::scope(|s| {
        s.spawn(|| {
            for (i, p) in pcs.iter().enumerate() {
                if cancel.load(Ordering::Relaxed) {
                    break;
                }
                let seq = match req
                    .query
                    .fetch_raw(p.rec as usize, p.offset..p.offset + p.len)
                {
                    Ok(seq) => seq,
                    Err(e) => {
                        *feed_error.lock().unwrap() = Some(e);
                        break;
                    }
                };
                let name = if req.self_comparison {
                    q_recs[p.rec as usize].name.clone()
                } else {
                    piece_name(i)
                };
                // A closed pipe means the aligner exited; its status says why.
                if writeln!(stdin, ">{name}")
                    .and_then(|_| stdin.write_all(&seq))
                    .and_then(|_| stdin.write_all(b"\n"))
                    .is_err()
                {
                    break;
                }
            }
            drop(stdin);
        });
        s.spawn(|| {
            let mut mapped = 0u64;
            for line in BufReader::new(stderr).lines().map_while(Result::ok) {
                // minimap2: "[M::worker_pipeline::…] mapped 123 sequences"
                if let Some(n) = line
                    .split("mapped ")
                    .nth(1)
                    .and_then(|r| r.split_whitespace().next())
                    .and_then(|n| n.parse::<u64>().ok())
                {
                    mapped += n;
                    progress(mapped.min(total), total);
                }
                let mut t = stderr_tail.lock().unwrap();
                t.push(line);
                if t.len() > 20 {
                    t.remove(0);
                }
            }
        });
        let mut text = String::new();
        for line in BufReader::new(stdout).lines() {
            let line = line?;
            if cancel.load(Ordering::Relaxed) {
                let _ = child.kill();
                return Err(AlignerError::Cancelled);
            }
            if req.self_comparison {
                text.push_str(&line);
            } else if let Some(l) = remap_line(&line, pcs, q_recs) {
                text.push_str(&l);
            } else {
                continue;
            }
            text.push('\n');
        }
        Ok(text)
    })?;
    let status = child.wait()?;
    if cancel.load(Ordering::Relaxed) {
        return Err(AlignerError::Cancelled);
    }
    if let Some(e) = feed_error.into_inner().unwrap() {
        return Err(e.into());
    }
    if !status.success() {
        return Err(AlignerError::Failed {
            program: program.to_path_buf(),
            status: status.to_string(),
            stderr: stderr_tail.into_inner().unwrap().join("\n"),
        });
    }
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn recs(v: &[(&str, u64)]) -> Vec<GenomeRecord> {
        v.iter()
            .map(|&(n, len)| GenomeRecord {
                name: n.into(),
                len,
            })
            .collect()
    }

    #[test]
    fn records_are_cut_into_pieces() {
        let r = recs(&[("a", 25), ("b", 10), ("c", 3)]);
        let p = pieces(&r, 10, true);
        let got: Vec<(u32, u64, u64)> = p.iter().map(|p| (p.rec, p.offset, p.len)).collect();
        assert_eq!(
            got,
            vec![(0, 0, 10), (0, 10, 10), (0, 20, 5), (1, 0, 10), (2, 0, 3)]
        );
        assert_eq!(pieces(&r, 10, false).len(), 3);
    }

    #[test]
    fn piece_lines_are_remapped_to_records() {
        let r = recs(&[("chr1", 25)]);
        let p = pieces(&r, 10, true);
        let line = "dottir_p2\t5\t1\t4\t+\tt\t100\t0\t3\t3\t3\t60\ttp:A:P";
        assert_eq!(
            remap_line(line, &p, &r).unwrap(),
            "chr1\t25\t21\t24\t+\tt\t100\t0\t3\t3\t3\t60\ttp:A:P"
        );
        assert!(remap_line("other\t5\t1\t4", &p, &r).is_none());
    }
}
