//! `dottir batch --engine genome` (ADR 0008): draw the alignments of two
//! assemblies — from `--paf` or an external aligner — without loading
//! their sequences.

use std::cell::OnceCell;
use std::path::{Path, PathBuf};
use std::sync::atomic::AtomicBool;

use anyhow::{bail, Context, Result};
use dottir_core::genome::alignment::sort_alignments;
use dottir_core::genome::order::Placement;
use dottir_core::genome::{
    arrange, rasterize, Alignment, Arranged, AxisLayout, FilterConfig, GenomeRecord, RasterConfig,
    GENOME_FORMAT_VERSION,
};
use dottir_io::aligner::{align, AlignRequest, Aligner};
use dottir_io::faidx::{plain_fasta, FastaIndex};
use dottir_io::genome_export::{write_association, write_fasta_records};
use dottir_io::paf::{read_paf_file, read_paf_oriented, write_paf};
use dottir_io::params::{
    sha256_file, DottirInfo, GenomeParamsInfo, HostInfo, InputInfo, ParamsSidecar, PlotParamsInfo,
};

use crate::large_batch::EngineArg;
use crate::BatchArgs;

/// Default largest plot dimension in pixels.
const DEFAULT_PLOT_PX: u32 = 4000;

/// Genome-mode options of `dottir batch`.
#[derive(clap::Args, Debug, Clone)]
pub struct GenomeArgs {
    /// Alignments to draw (PAF, plain or gzipped), e.g. from
    /// `minimap2 -x asm20 FIRST.fa SECOND.fa`, instead of running an
    /// aligner. Implies `--engine genome`. The first input is the
    /// horizontal (reference) axis.
    #[arg(long, value_name = "FILE")]
    pub paf: Option<PathBuf>,

    /// `--engine auto` uses genome mode from this many matrix cells
    /// (first × second input length) when an aligner is available.
    #[arg(long, value_name = "CELLS", default_value_t = 1e17)]
    pub genome_threshold_cells: f64,

    /// Aligner: `auto` (mm2plus, else minimap2), a name on PATH, or a
    /// path to a minimap2-compatible binary.
    #[arg(long, value_name = "NAME|PATH", default_value = "auto")]
    pub aligner: String,

    /// minimap2 preset (`-x`).
    #[arg(long, value_name = "PRESET", default_value = "asm20")]
    pub preset: String,

    /// Extra aligner arguments, space-separated (e.g. "-c" for exact
    /// identities, at several times the run time).
    #[arg(long, value_name = "ARGS", allow_hyphen_values = true)]
    pub aligner_args: Option<String>,

    /// Aligner threads (0 = all cores).
    #[arg(long, value_name = "N", default_value_t = 0)]
    pub threads: usize,

    /// Query records are sent to the aligner in pieces of at most this
    /// many bp (pairwise only), which bounds its memory.
    #[arg(long, value_name = "BP", default_value_t = 10_000_000)]
    pub chunk_size: u64,

    /// Self-comparison: keep every chain (`-P`), i.e. all copies of each
    /// repeat against each other. Can be millions of alignments.
    #[arg(long)]
    pub self_all_copies: bool,

    /// Re-run the aligner even if a cached result exists.
    #[arg(long)]
    pub no_cache: bool,

    /// Drop alignments shorter than this (bp, longer span).
    #[arg(long, value_name = "BP", default_value_t = 0)]
    pub min_aln_len: u64,

    /// Drop alignments below this identity (0–1; matches / block length,
    /// or 1 − de:f when present).
    #[arg(long, value_name = "FRACTION", default_value_t = 0.0)]
    pub min_identity: f64,

    /// Draw only primary alignments (`tp:A:P`).
    #[arg(long)]
    pub primary_only: bool,

    /// Hide short, frequent alignments (the D-GENIES noise filter;
    /// plots with more than 1000 alignments).
    #[arg(long)]
    pub remove_noise: bool,

    /// Sort and orient the vertical assembly's records along the
    /// horizontal one (D-GENIES contig sorting).
    #[arg(long)]
    pub sort_query: bool,

    /// Line thickness in pixels.
    #[arg(long, value_name = "PX", default_value_t = 2)]
    pub line_width: u32,

    /// Write the association table: each vertical-axis record, in plot
    /// order, with its best-matching horizontal-axis record (D-GENIES
    /// columns).
    #[arg(long, value_name = "FILE.tsv")]
    pub assoc_out: Option<PathBuf>,

    /// Write the vertical-axis records as FASTA in plot order,
    /// reverse-complemented where the plot shows them reversed.
    #[arg(long, value_name = "FILE.fa")]
    pub sorted_query_out: Option<PathBuf>,

    /// Write records without any alignment to `PREFIX.reference.fa`
    /// (horizontal axis) and `PREFIX.query.fa` (vertical axis).
    #[arg(long, value_name = "PREFIX")]
    pub unaligned_out: Option<PathBuf>,

    /// Write the drawn (filtered) alignments as PAF, query = vertical axis.
    #[arg(long, value_name = "FILE.paf")]
    pub paf_out: Option<PathBuf>,
}

/// Rough residue count from the file size (gzip compresses DNA ~4×).
fn estimated_len(path: &Path) -> u64 {
    let size = std::fs::metadata(path).map_or(0, |m| m.len());
    let gz = std::fs::File::open(path)
        .and_then(|mut f| {
            let mut m = [0u8; 2];
            std::io::Read::read_exact(&mut f, &mut m).map(|_| m == [0x1f, 0x8b])
        })
        .unwrap_or(false);
    if gz {
        size.saturating_mul(5)
    } else {
        size
    }
}

/// Whether `dottir batch` should take the genome path.
pub fn wanted(args: &BatchArgs) -> Result<bool> {
    let (large, genome) = (&args.large, &args.genome);
    match (large.engine, &genome.paf) {
        (EngineArg::Genome, _) | (EngineArg::Auto, Some(_)) => Ok(true),
        (EngineArg::Exact | EngineArg::Large, Some(_)) => {
            bail!("--paf is only used by --engine genome")
        }
        (EngineArg::Exact | EngineArg::Large, None) => Ok(false),
        (EngineArg::Auto, None) => auto_wants_genome(args),
    }
}

fn auto_wants_genome(args: &BatchArgs) -> Result<bool> {
    let threshold = args.genome.genome_threshold_cells;
    let second = args.subject.as_deref().unwrap_or(&args.query);
    // Cheap estimate first; only plausible candidates get indexed.
    let estimate = estimated_len(&args.query) as f64 * estimated_len(second) as f64;
    if estimate < threshold / 4.0 {
        return Ok(false);
    }
    let cache = dottir_io::cache::cache_dir();
    let len = |p: &Path| -> Result<u64> {
        Ok(open_assembly(p, &cache)?
            .records
            .iter()
            .map(|r| r.len)
            .sum())
    };
    let cells = len(&args.query)? as f64 * len(second)? as f64;
    if cells < threshold {
        return Ok(false);
    }
    match Aligner::detect(&args.genome.aligner) {
        Ok(a) => {
            tracing::info!(
                "engine auto: {cells:.2e} cells >= {threshold:.2e} and {} found, using genome mode",
                a.name
            );
            Ok(true)
        }
        Err(e) => {
            tracing::warn!("engine auto: {cells:.2e} cells, but {e}; using the large engine");
            Ok(false)
        }
    }
}

struct Assembly {
    path: PathBuf,
    index: FastaIndex,
    records: Vec<GenomeRecord>,
    sha256: OnceCell<String>,
}

impl Assembly {
    fn sha256(&self) -> Result<&str> {
        if self.sha256.get().is_none() {
            let _ = self.sha256.set(sha256_file(&self.path)?);
        }
        Ok(self.sha256.get().map(String::as_str).unwrap_or_default())
    }
}

fn open_assembly(path: &Path, cache: &Path) -> Result<Assembly> {
    let plain = plain_fasta(path, cache).with_context(|| format!("reading {}", path.display()))?;
    if plain != path {
        tracing::info!("{} is gzipped; using {}", path.display(), plain.display());
    }
    let index =
        FastaIndex::open(&plain, cache).with_context(|| format!("indexing {}", plain.display()))?;
    let records = index.records();
    tracing::info!(
        "{}: {} records, {} bp",
        path.display(),
        records.len(),
        records.iter().map(|r| r.len).sum::<u64>()
    );
    Ok(Assembly {
        path: path.to_path_buf(),
        index,
        records,
        sha256: OnceCell::new(),
    })
}

/// Read the PAF with its query on the vertical axis, swapping a PAF made
/// the other way round.
fn read_oriented(paf: &Path, x: &[GenomeRecord], y: &[GenomeRecord]) -> Result<Vec<Alignment>> {
    let (alns, swapped) =
        read_paf_oriented(paf, y, x).with_context(|| format!("reading {}", paf.display()))?;
    if swapped {
        tracing::info!("PAF query names match the first input; swapping axes");
    }
    Ok(alns)
}

fn create(path: &Path) -> Result<std::io::BufWriter<std::fs::File>> {
    Ok(std::io::BufWriter::new(
        std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?,
    ))
}

fn write_exports(
    g: &GenomeArgs,
    x: &Assembly,
    y: &Assembly,
    y_layout: &AxisLayout,
    placements: &[Option<Placement>],
    alns: &[Alignment],
    drawn: &[&Alignment],
) -> Result<()> {
    if let Some(p) = &g.assoc_out {
        write_association(create(p)?, y_layout, &y.records, &x.records, placements)?;
        tracing::info!("wrote association table {}", p.display());
    }
    if let Some(p) = &g.sorted_query_out {
        let order: Vec<u32> = y_layout.segments.iter().map(|s| s.rec).collect();
        let mut flips = vec![false; y.records.len()];
        for s in &y_layout.segments {
            flips[s.rec as usize] = s.reverse;
        }
        write_fasta_records(create(p)?, &y.index, &order, &flips)?;
        tracing::info!("wrote {}", p.display());
    }
    if let Some(prefix) = &g.unaligned_out {
        let mut x_hit = vec![false; x.records.len()];
        let mut y_hit = vec![false; y.records.len()];
        for a in alns {
            x_hit[a.t_rec as usize] = true;
            y_hit[a.q_rec as usize] = true;
        }
        for (asm, hit, side) in [(x, &x_hit, "reference"), (y, &y_hit, "query")] {
            let none: Vec<u32> = (0..hit.len() as u32)
                .filter(|&r| !hit[r as usize])
                .collect();
            let p = PathBuf::from(format!("{}.{side}.fa", prefix.display()));
            write_fasta_records(create(&p)?, &asm.index, &none, &[])?;
            tracing::info!(
                "wrote {} unaligned {side} records to {}",
                none.len(),
                p.display()
            );
        }
    }
    if let Some(p) = &g.paf_out {
        let mut out: Vec<Alignment> = drawn.iter().map(|&a| a.clone()).collect();
        sort_alignments(&mut out);
        write_paf(create(p)?, &out, &y.records, &x.records)?;
        tracing::info!("wrote {} alignments to {}", out.len(), p.display());
    }
    Ok(())
}

/// Where the alignments came from, for the sidecar.
struct Source {
    paf_path: Option<PathBuf>,
    paf_sha256: String,
    aligner: Option<Aligner>,
    aligner_args: Option<Vec<String>>,
    chunk_size: Option<u64>,
}

fn run_aligner(
    g: &GenomeArgs,
    x: &Assembly,
    y: &Assembly,
    self_comparison: bool,
) -> Result<(Vec<Alignment>, Source)> {
    let aligner = Aligner::detect(&g.aligner)?;
    tracing::info!(
        "aligner: {} ({})",
        aligner.program.display(),
        aligner.version
    );
    let threads = match g.threads {
        0 => std::thread::available_parallelism().map_or(1, |n| n.get()),
        n => n,
    };
    let req = AlignRequest {
        reference: &x.index,
        query: &y.index,
        reference_sha256: x.sha256()?.to_string(),
        query_sha256: y.sha256()?.to_string(),
        preset: g.preset.clone(),
        extra_args: g
            .aligner_args
            .as_deref()
            .map(|a| a.split_whitespace().map(String::from).collect())
            .unwrap_or_default(),
        threads,
        chunk_size: g.chunk_size,
        self_comparison,
        all_copies: g.self_all_copies,
    };
    let last = std::sync::Mutex::new(0u64);
    let progress = |done: u64, total: u64| {
        let decile = done * 10 / total.max(1);
        let mut l = last.lock().unwrap();
        if decile > *l {
            *l = decile;
            tracing::info!("aligned {done} / {total} pieces");
        }
    };
    let started = std::time::Instant::now();
    let cache = dottir_io::cache::cache_dir();
    let out = align(
        &aligner,
        &req,
        &cache,
        !g.no_cache,
        &progress,
        &AtomicBool::new(false),
    )?;
    if out.cached {
        tracing::info!("reusing cached alignments {}", out.paf.display());
    } else {
        tracing::info!(
            "aligned in {:.1?}; cached as {}",
            started.elapsed(),
            out.paf.display()
        );
    }
    let alns = read_paf_file(&out.paf, &y.records, &x.records)
        .with_context(|| format!("reading {}", out.paf.display()))?;
    let source = Source {
        paf_path: None,
        paf_sha256: sha256_file(&out.paf)?,
        aligner_args: Some(req.output_args()),
        chunk_size: (!self_comparison).then_some(g.chunk_size),
        aligner: Some(aligner),
    };
    Ok((alns, source))
}

pub fn run(args: &BatchArgs) -> Result<()> {
    let g = &args.genome;
    let self_comparison = args.subject.is_none();
    if self_comparison && g.sort_query {
        bail!("--sort-query needs two inputs");
    }
    let cache = dottir_io::cache::cache_dir();
    let x = open_assembly(&args.query, &cache)?;
    let y = match &args.subject {
        Some(p) => open_assembly(p, &cache)?,
        None => Assembly {
            path: x.path.clone(),
            index: x.index.clone(),
            records: x.records.clone(),
            sha256: x.sha256.clone(),
        },
    };

    let (alns, source) = match &g.paf {
        Some(paf) => {
            let alns = read_oriented(paf, &x.records, &y.records)?;
            let source = Source {
                paf_path: Some(paf.clone()),
                paf_sha256: sha256_file(paf)?,
                aligner: None,
                aligner_args: None,
                chunk_size: None,
            };
            (alns, source)
        }
        None => run_aligner(g, &x, &y, self_comparison)?,
    };
    let cfg = FilterConfig {
        min_len: g.min_aln_len,
        min_identity: g.min_identity,
        primary_only: g.primary_only,
        remove_noise: g.remove_noise,
    };
    let Arranged {
        keep,
        noise,
        x: x_layout,
        y: y_layout,
        placements,
    } = arrange(&alns, &x.records, &y.records, &cfg, g.sort_query, &[]);
    let drawn: Vec<&Alignment> = alns
        .iter()
        .zip(&keep)
        .filter(|(_, &k)| k)
        .map(|(a, _)| a)
        .collect();
    tracing::info!(
        "{} alignments, {} drawn{}",
        alns.len(),
        drawn.len(),
        noise.map_or(String::new(), |t| format!(" (noise below {t:.0} bp)"))
    );
    if g.sort_query {
        tracing::info!(
            "sorted {} records, {} reversed, {} without alignments",
            y_layout.segments.len(),
            y_layout.segments.iter().filter(|s| s.reverse).count(),
            placements.iter().filter(|p| p.is_none()).count()
        );
    }
    write_exports(g, &x, &y, &y_layout, &placements, &alns, &drawn)?;

    let max_px = u64::from(args.auto_zoom.unwrap_or(DEFAULT_PLOT_PX).max(1));
    let bp_per_px = x_layout.total.max(y_layout.total).div_ceil(max_px).max(1);
    let rcfg = RasterConfig {
        bp_per_px,
        line_width: g.line_width,
        self_comparison,
    };
    let raster = rasterize(drawn.iter().copied(), &x_layout, &y_layout, &rcfg);
    tracing::info!("{}×{} px, {bp_per_px} bp/px", raster.width, raster.height);

    let chunks = [
        ("dottir-engine", "genome".to_string()),
        ("dottir-version", env!("CARGO_PKG_VERSION").to_string()),
        ("dottir-genome-format", GENOME_FORMAT_VERSION.to_string()),
        ("dottir-paf-sha256", source.paf_sha256.clone()),
        ("dottir-bp-per-px", bp_per_px.to_string()),
    ];
    let chunk_refs: Vec<(&str, &str)> = chunks.iter().map(|(k, v)| (*k, v.as_str())).collect();
    dottir_io::genome_export::write_genome_png(
        &args.output,
        &raster,
        &x_layout,
        &x.records,
        &y_layout,
        &y.records,
        &chunk_refs,
    )
    .with_context(|| format!("writing {}", args.output.display()))?;
    tracing::info!("wrote {}", args.output.display());

    if !args.no_sidecar {
        let info = |a: &Assembly| -> Result<InputInfo> {
            Ok(InputInfo {
                path: a.path.display().to_string(),
                sha256: a.sha256()?.to_string(),
                size_bytes: std::fs::metadata(&a.path)?.len(),
                n_records: a.records.len(),
                total_residues: a.records.iter().map(|r| r.len).sum::<u64>() as usize,
            })
        };
        let x_info = info(&x)?;
        let y_info = if self_comparison {
            x_info.clone()
        } else {
            info(&y)?
        };
        let sidecar = ParamsSidecar {
            dottir: DottirInfo {
                version: env!("CARGO_PKG_VERSION").to_string(),
                git_sha: option_env!("GIT_SHA").map(|s| s.to_string()),
                pixelmap_format_version: GENOME_FORMAT_VERSION,
            },
            query: x_info,
            subject: y_info,
            plot: PlotParamsInfo {
                mode: "genome".to_string(),
                matrix: String::new(),
                strand: "both".to_string(),
                window_size: 0,
                zoom: u32::try_from(bp_per_px).unwrap_or(u32::MAX),
                pixel_fac: 0,
                self_comparison,
                karlin: None,
                width: raster.width,
                height: raster.height,
            },
            host: HostInfo {
                hostname: dottir_io::params::hostname(),
                os: std::env::consts::OS.to_string(),
                timestamp_utc: crate::now_iso_utc(),
            },
            large: None,
            genome: Some(GenomeParamsInfo {
                format_version: GENOME_FORMAT_VERSION,
                source: if source.aligner.is_some() {
                    "aligner"
                } else {
                    "paf"
                }
                .to_string(),
                paf_path: source.paf_path.as_ref().map(|p| p.display().to_string()),
                paf_sha256: Some(source.paf_sha256.clone()),
                aligner: source.aligner.as_ref().map(|a| a.name.clone()),
                aligner_version: source.aligner.as_ref().map(|a| a.version.clone()),
                aligner_args: source.aligner_args.clone(),
                chunk_size: source.chunk_size,
                min_aln_len: g.min_aln_len,
                min_identity: g.min_identity,
                primary_only: g.primary_only,
                remove_noise: g.remove_noise,
                noise_threshold: noise,
                sorted_query: g.sort_query,
                alignments: alns.len(),
                alignments_drawn: drawn.len(),
                bp_per_px,
                line_width: g.line_width,
            }),
        };
        let path = crate::sidecar_path(&args.output);
        std::fs::write(&path, sidecar.to_toml()?)?;
        tracing::info!("wrote sidecar {}", path.display());
    }
    Ok(())
}
