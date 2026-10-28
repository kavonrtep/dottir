//! `dottir batch` with the seed-based large-sequence engine (ADR 0007).

use std::path::PathBuf;
use std::sync::atomic::{AtomicU32, Ordering};

use anyhow::{bail, Context, Result};
use clap::ValueEnum;

use dottir_core::large::{
    seed::DEFAULT_PATTERN, ChainConfig, ClassMask, ClassThresholds, Control, DepthMethod,
    LargeConfig, LargeInput, LargeSession, RepeatClass, Shade, SpacedSeed,
};
use dottir_core::{BlastMode, DotPlot, PlotConfig};
use dottir_io::params::LargeParamsInfo;
use dottir_io::{AnnotSet, SeqInput, Sequence};

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum EngineArg {
    /// Dotter's sliding window over every cell.
    Exact,
    /// Seed-based engine for megabase inputs (DNA only).
    Large,
    /// `large` when query × subject exceeds `--large-threshold-cells`.
    Auto,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum DepthArg {
    Seed,
    Kmer,
    None,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
pub enum MaskClassArg {
    None,
    Low,
    Mid,
    High,
}

impl MaskClassArg {
    fn class(self) -> Option<RepeatClass> {
        match self {
            Self::None => None,
            Self::Low => Some(RepeatClass::Low),
            Self::Mid => Some(RepeatClass::Mid),
            Self::High => Some(RepeatClass::High),
        }
    }
}

/// Large-mode options of `dottir batch`.
#[derive(clap::Args, Debug, Clone)]
pub struct LargeArgs {
    /// Plot engine.
    #[arg(long, value_enum, default_value_t = EngineArg::Auto)]
    pub engine: EngineArg,

    /// `--engine auto` switches to `large` above this many matrix cells.
    #[arg(long, value_name = "CELLS", default_value_t = 1e11)]
    pub large_threshold_cells: f64,

    /// Spaced seed pattern (`1` = must match). Weight 6–13, span <= 32.
    #[arg(long, value_name = "PATTERN", default_value = DEFAULT_PATTERN)]
    pub seed_pattern: String,

    /// Two seed hits on one diagonal at most this far apart survive.
    #[arg(long, value_name = "BP", default_value_t = 100)]
    pub two_hit_dist: u32,

    /// Skip seeds occurring more often than this in the subject.
    #[arg(long, value_name = "N", default_value_t = 1000)]
    pub max_occ: u32,

    /// Look up one seed in N (by seed hash) for the plot; the two-hit
    /// distance scales by N. `auto` = zoom / 6000, capped at 32.
    #[arg(long, value_name = "N|auto", default_value = "auto", value_parser = parse_rate)]
    pub seed_sample: Rate,

    /// Seed sampling for `--depth seed`. `auto` = sequence length / 25 Mb,
    /// capped at 16.
    #[arg(long, value_name = "N|auto", default_value = "auto", value_parser = parse_rate)]
    pub depth_sample: Rate,

    /// Cells scored beyond each anchored run. Default: the window size.
    #[arg(long, value_name = "BP")]
    pub extend: Option<u32>,

    /// Repeat-depth estimate: `seed` counts diagonal runs in a
    /// self-comparison (sees diverged repeats; one extra seeding pass),
    /// `kmer` counts exact k-mers (fast; young repeats only), `none`
    /// makes every dot `unique` unless soft-mask/annotation floors apply.
    #[arg(long, value_enum, default_value_t = DepthArg::Seed)]
    pub depth: DepthArg,

    /// k for `--depth kmer`.
    #[arg(long, value_name = "K", default_value_t = 21)]
    pub depth_k: u8,

    /// Copy-number thresholds `LOW,MID,HIGH`: depth < LOW is unique,
    /// < MID low, < HIGH mid, otherwise high.
    #[arg(long, value_name = "L,M,H", default_value = "2,5,51")]
    pub class_thresholds: String,

    /// Repeat classes to leave out of the image (`unique,low,mid,high`).
    #[arg(long, value_delimiter = ',', value_name = "CLASSES")]
    pub hide_classes: Vec<String>,

    /// Minimum class of soft-masked (lowercase) bases.
    #[arg(long, value_enum, default_value_t = MaskClassArg::Mid)]
    pub softmask_class: MaskClassArg,

    /// GFF3/BED of repeats on the query (and on the subject in a
    /// self-comparison); features count as `high`. Defaults to the
    /// query's embedded GFF3 features when `--repeat-types` is given.
    #[arg(long, value_name = "PATH")]
    pub repeat_annot: Option<PathBuf>,

    /// GFF3/BED of repeats on the subject.
    #[arg(long, value_name = "PATH")]
    pub repeat_annot_subject: Option<PathBuf>,

    /// Feature types (GFF3 column 3) used from repeat annotations.
    /// Empty = all features.
    #[arg(long, value_delimiter = ',', value_name = "TYPES")]
    pub repeat_types: Vec<String>,

    /// Write synteny chains as PAF (no CIGAR) to this path.
    #[arg(long, value_name = "PATH")]
    pub chains_out: Option<PathBuf>,

    /// Minimum chain score (roughly bases of co-linear similarity).
    #[arg(long, value_name = "SCORE", default_value_t = 300.0)]
    pub chain_min_score: f64,

    /// Most repetitive class used for chaining.
    #[arg(long, value_name = "CLASS", default_value = "low")]
    pub chain_max_class: String,

    /// Write the depth track(s) as bedGraph: `<PREFIX>.query.bedgraph`
    /// (and `.subject.bedgraph` for a pairwise plot), binned to `--zoom`.
    #[arg(long, value_name = "PREFIX")]
    pub depth_out: Option<PathBuf>,

    /// Memory budget for the large engine (index, depth, channels).
    #[arg(long, value_name = "BYTES", default_value_t = 6 << 30)]
    pub large_memory_limit_bytes: u64,
}

impl LargeArgs {
    /// Decide whether the large engine runs.
    pub fn use_large(&self, mode: BlastMode, qlen: usize, slen: usize) -> Result<bool> {
        match self.engine {
            EngineArg::Exact => Ok(false),
            EngineArg::Large if mode != BlastMode::Blastn => {
                bail!("--engine large supports BLASTN (DNA) only")
            }
            EngineArg::Large => Ok(true),
            EngineArg::Auto => {
                let cells = qlen as f64 * slen as f64;
                let large = mode == BlastMode::Blastn && cells > self.large_threshold_cells;
                if large {
                    tracing::info!(
                        "engine auto: {cells:.2e} cells > {:.2e}, using the large engine",
                        self.large_threshold_cells
                    );
                }
                Ok(large)
            }
        }
    }
}

fn parse_class(s: &str) -> Result<RepeatClass> {
    RepeatClass::from_name(s.trim())
        .with_context(|| format!("unknown repeat class {s:?} (unique, low, mid, high)"))
}

fn parse_thresholds(s: &str) -> Result<ClassThresholds> {
    let v: Vec<u16> = s
        .split(',')
        .map(|x| x.trim().parse::<u16>())
        .collect::<Result<_, _>>()
        .with_context(|| format!("--class-thresholds {s:?}: expected three integers"))?;
    let [low, mid, high] = v[..] else {
        bail!("--class-thresholds {s:?}: expected three integers");
    };
    let t = ClassThresholds { low, mid, high };
    t.validate()?;
    Ok(t)
}

fn floor_for(
    role: &str,
    seq: &Sequence,
    embedded: Option<&AnnotSet>,
    annot_path: Option<&PathBuf>,
    args: &LargeArgs,
    thresholds: &ClassThresholds,
) -> Result<Option<Vec<u16>>> {
    let soft = args
        .softmask_class
        .class()
        .and_then(|c| dottir_io::large::soft_mask_floor(seq, thresholds.floor_of(c)));
    if let Some(f) = &soft {
        let masked = f.iter().filter(|&&x| x > 0).count();
        tracing::info!(
            "{role}: {masked} soft-masked bases raised to class {:?}",
            args.softmask_class
        );
    }
    let loaded;
    let annots = match annot_path {
        Some(p) => {
            loaded = AnnotSet::load(p).with_context(|| format!("reading {}", p.display()))?;
            Some(&loaded)
        }
        None if !args.repeat_types.is_empty() => embedded,
        None => None,
    };
    let annot = annots.map(|a| {
        let (f, skipped) = dottir_io::large::annotation_floor(
            seq,
            a,
            &args.repeat_types,
            thresholds.floor_of(RepeatClass::High),
        );
        let covered = f.iter().filter(|&&x| x > 0).count();
        tracing::info!("{role}: repeat annotation covers {covered} bases ({skipped} features on unknown records)");
        f
    });
    Ok(dottir_io::large::merge_floors(soft, annot))
}

pub struct LargeOutcome {
    pub plot: DotPlot,
    pub info: LargeParamsInfo,
}

/// Run the large engine with the resolved exact-engine `cfg` (zoom,
/// strand, window, pixel factor, triangle) plus `args`.
pub fn run(
    args: &LargeArgs,
    cfg: &PlotConfig,
    query_input: &SeqInput,
    subject_input: &SeqInput,
    is_self: bool,
) -> Result<LargeOutcome> {
    if cfg.reverse_query || cfg.reverse_subject {
        bail!("-r / -v are not supported by the large engine");
    }
    let query = &query_input.sequence;
    let subject = &subject_input.sequence;
    let thresholds = parse_thresholds(&args.class_thresholds)?;
    let hidden: Vec<RepeatClass> = args
        .hide_classes
        .iter()
        .map(|c| parse_class(c))
        .collect::<Result<_>>()?;
    let shown = hidden.iter().fold(ClassMask::ALL, |m, &c| m.with(c, false));
    let seed = SpacedSeed::parse(&args.seed_pattern)?;

    let chain = args
        .chains_out
        .as_ref()
        .map(|_| -> Result<ChainConfig> {
            Ok(ChainConfig {
                max_class: parse_class(&args.chain_max_class)?,
                min_score: args.chain_min_score,
                ..ChainConfig::default()
            })
        })
        .transpose()?;

    let lcfg = LargeConfig {
        coverage: false,
        matrix: cfg.matrix.clone(),
        window_size: cfg.window_size,
        zoom: cfg.zoom,
        pixel_fac: cfg.pixel_fac,
        strand: cfg.strand,
        self_comparison: is_self,
        triangle: cfg.triangle,
        disable_mirror: cfg.disable_mirror,
        memory_limit_bytes: args.large_memory_limit_bytes,
        seed,
        two_hit_dist: args.two_hit_dist,
        max_occ: args.max_occ,
        sample: args.seed_sample.0,
        depth_sample: args.depth_sample.0,
        extend: args.extend,
        depth: match args.depth {
            DepthArg::Seed => DepthMethod::Seed,
            DepthArg::Kmer => DepthMethod::Kmer { k: args.depth_k },
            DepthArg::None => DepthMethod::None,
        },
        classes: thresholds,
        chain,
        viewport: None,
    };

    let q_floor = floor_for(
        "query",
        query,
        query_input.annotations.as_ref(),
        args.repeat_annot.as_ref(),
        args,
        &thresholds,
    )?;
    let s_floor = if is_self {
        None
    } else {
        floor_for(
            "subject",
            subject,
            subject_input.annotations.as_ref(),
            args.repeat_annot_subject.as_ref(),
            args,
            &thresholds,
        )?
    };
    let q_breaks = query.breaks();
    let s_breaks = subject.breaks();
    let input = LargeInput {
        query: query.bytes(),
        q_breaks: &q_breaks,
        q_floor: q_floor.as_deref(),
        subject: (!is_self).then(|| subject.bytes()),
        s_breaks: if is_self { &q_breaks } else { &s_breaks },
        s_floor: s_floor.as_deref(),
    };

    let t0 = std::time::Instant::now();
    let last_decile = AtomicU32::new(0);
    let stage = std::sync::Mutex::new("prepare");
    let progress = |f: f32| {
        let d = (f * 10.0) as u32;
        if last_decile.fetch_max(d, Ordering::Relaxed) < d {
            let name = *stage.lock().unwrap_or_else(|e| e.into_inner());
            tracing::info!("{name}: {}% ({:.1}s)", d * 10, t0.elapsed().as_secs_f32());
        }
    };
    let ctl = Control {
        cancel: None,
        progress: Some(&progress),
    };
    tracing::info!(
        "large engine: seed {} (weight {}), two-hit {} bp, max-occ {}, depth {}",
        lcfg.seed.pattern(),
        lcfg.seed.weight(),
        lcfg.two_hit_dist,
        lcfg.max_occ,
        lcfg.depth.describe()
    );
    let session = LargeSession::prepare(&input, &lcfg, &ctl).context("large engine: prepare")?;
    tracing::info!(
        "prepared index ({} seeds, {:.1} MiB) in {:.1}s; W={}, pixel_fac={}",
        session.index().n_positions(),
        session.bytes() as f64 / (1 << 20) as f64,
        t0.elapsed().as_secs_f32(),
        session.window(),
        session.pixel_fac()
    );
    *stage.lock().unwrap_or_else(|e| e.into_inner()) = "compute";
    last_decile.store(0, Ordering::Relaxed);
    let t1 = std::time::Instant::now();
    let lp = session
        .compute(&lcfg, &ctl)
        .context("large engine: compute")?;
    let st = lp.stats;
    tracing::info!(
        "computed in {:.1}s (seed sampling 1/{}, depth 1/{}): {} raw hits, {} anchors, {} runs, {} cells scored, {} capped query positions",
        t1.elapsed().as_secs_f32(),
        lp.params.sample,
        lp.params.depth_sample.0.max(lp.params.depth_sample.1),
        st.raw_hits,
        st.anchors,
        st.runs,
        st.scored_cells,
        st.capped_query_positions
    );

    if let Some(path) = &args.chains_out {
        dottir_io::large::write_chains_paf(path, &lp.chains, query, subject)
            .with_context(|| format!("writing {}", path.display()))?;
        tracing::info!("wrote {} chains to {}", lp.chains.len(), path.display());
    }
    if let Some(prefix) = &args.depth_out {
        let bin = cfg.zoom as usize;
        let qp = with_suffix(prefix, "query.bedgraph");
        dottir_io::large::write_depth_bedgraph(&qp, session.query_depth(), query, bin)?;
        tracing::info!("wrote {}", qp.display());
        if !is_self {
            let sp = with_suffix(prefix, "subject.bedgraph");
            dottir_io::large::write_depth_bedgraph(&sp, session.subject_depth(), subject, bin)?;
            tracing::info!("wrote {}", sp.display());
        }
    }

    let info = LargeParamsInfo {
        format_version: lp.params.format_version,
        seed: lp.params.seed.clone(),
        two_hit_dist: lp.params.two_hit_dist,
        max_occ: lp.params.max_occ,
        seed_sample: lp.params.sample,
        depth_sample: [lp.params.depth_sample.0, lp.params.depth_sample.1],
        extend: lp.params.extend,
        depth: lp.params.depth.describe(),
        class_thresholds: [thresholds.low, thresholds.mid, thresholds.high],
        shown_classes: shown
            .classes()
            .iter()
            .map(|c| c.name().to_string())
            .collect(),
        soft_mask_class: args.softmask_class.class().map(|c| c.name().to_string()),
        repeat_annotation: args.repeat_annot.as_ref().map(|p| p.display().to_string()),
        raw_hits: st.raw_hits,
        anchors: st.anchors,
        runs: st.runs,
        scored_cells: st.scored_cells,
        capped_query_positions: st.capped_query_positions,
        chains: args.chains_out.as_ref().map(|_| lp.chains.len()),
    };
    Ok(LargeOutcome {
        plot: lp.to_dot_plot(shown, false, Shade::Score),
        info,
    })
}

fn with_suffix(prefix: &std::path::Path, suffix: &str) -> PathBuf {
    let mut s = prefix.as_os_str().to_owned();
    s.push(".");
    s.push(suffix);
    PathBuf::from(s)
}

/// A sampling rate; `None` = automatic.
#[derive(Debug, Clone, Copy)]
pub struct Rate(pub Option<u32>);

/// `auto` or a sampling rate >= 1.
fn parse_rate(s: &str) -> Result<Rate, String> {
    if s.eq_ignore_ascii_case("auto") {
        return Ok(Rate(None));
    }
    match s.parse::<u32>() {
        Ok(n) if n >= 1 => Ok(Rate(Some(n))),
        _ => Err(format!("expected `auto` or an integer >= 1, got `{s}`")),
    }
}
