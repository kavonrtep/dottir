//! Seed-based large-sequence mode (ADR 0007).
//!
//! A second engine for DNA inputs at the 100 Mb scale. Instead of scoring
//! every cell, it
//!
//! 1. indexes the subject with a spaced seed ([`seed`], [`index`]),
//! 2. streams the query on both strands and keeps only hits with a
//!    partner on the same diagonal within `two_hit_dist` (`hits`),
//! 3. scores the anchored diagonals with the *exact* sliding-window
//!    recurrence (`score`), so each written pixel is `<=` the exact
//!    engine's pixel for the same inputs,
//! 4. routes each dot into a channel by strand and repeat class, where the
//!    class comes from a per-base k-mer depth track ([`depth`]),
//! 5. optionally chains the diagonal runs into synteny segments
//!    ([`chain`]).
//!
//! [`LargeSession`] holds the expensive, viewport-independent state (index
//! and depth tracks) so the GUI can recompute viewports cheaply.
//!
//! # Example
//!
//! ```
//! use dottir_core::large::{compute_large_dotplot, DepthMethod, LargeConfig, LargeInput};
//!
//! let unit = b"ACGTTGCAAGGCTTACGATCGATCGGATCCATGCATGCAAACCCGGGTTTAGCTAGCTTAGG";
//! let seq: Vec<u8> = unit.iter().cycle().take(4000).copied().collect();
//! let mut cfg = LargeConfig::default();
//! cfg.zoom = 10;
//! cfg.self_comparison = true;
//! cfg.depth = DepthMethod::None;
//! let plot = compute_large_dotplot(&LargeInput::self_comparison(&seq, &[]), &cfg).unwrap();
//! assert_eq!(plot.width, 400);
//! // The main diagonal is drawn.
//! assert!(plot.merged(None, Default::default())[200 * 400 + 200] > 0);
//! ```

pub mod chain;
pub mod depth;
pub(crate) mod hits;
pub mod index;
pub(crate) mod score;
pub mod seed;

use std::ops::Range;
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};

pub use chain::{chain_segments, Chain, ChainConfig};
pub use depth::{compute_depth, ClassThresholds, DepthMethod, DepthTrack, RepeatClass};
pub use index::SeedIndex;
pub use seed::SpacedSeed;

use crate::alphabet::{encode, AlphabetKind};
use crate::error::DottirError;
use crate::karlin::{karlin_window_size, KarlinResult};
use crate::matrix::{BlastMode, ScoreMatrix};
use crate::pixel::{image_dimension, PixelMap};
use crate::plot::{mirror_self_comparison, DotPlot, PlotParams, Strand, Triangle};

/// Bumped whenever large-mode output for fixed inputs and parameters
/// changes. Independent of [`crate::PIXELMAP_FORMAT_VERSION`].
pub const LARGE_PIXELMAP_FORMAT_VERSION: u32 = 2;

/// Largest automatic seed-sampling rate for the plot.
pub const MAX_AUTO_SAMPLE: u32 = 32;
/// Largest automatic seed-sampling rate for the seed depth pass.
pub const MAX_AUTO_DEPTH_SAMPLE: u32 = 16;

/// Automatic plot sampling: one seed in `zoom / 6000`. Sampling loses short
/// or diverged fragments (each needs two sampled hits), but bigger pixels
/// hold more fragments, so the share of lit pixels kept stays high: on a
/// TE-rich 100 Mb benchmark at 25 kb/pixel, rate 4 kept 94 % of lit pixels
/// (rate 8: 83 %, rate 25: 48 %) and halved the compute time.
pub fn auto_sample_rate(zoom: u32) -> u32 {
    (zoom / 6000).clamp(1, MAX_AUTO_SAMPLE)
}

/// Automatic depth sampling: one seed in `len / 25 Mb`.
pub fn auto_depth_sample_rate(len: usize) -> u32 {
    ((len / 25_000_000) as u32).clamp(1, MAX_AUTO_DEPTH_SAMPLE)
}

/// Query/subject ranges, in bases, of a sub-view.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Viewport {
    pub q: Range<usize>,
    pub s: Range<usize>,
}

#[derive(Debug, Clone)]
pub struct LargeConfig {
    /// DNA score matrix. Default identity (+5/−4).
    pub matrix: ScoreMatrix,
    /// `None` = Karlin/Altschul estimate, as in the exact engine.
    pub window_size: Option<u32>,
    pub zoom: u32,
    /// `0` = auto (dotter's formula from Karlin's expected residue score).
    pub pixel_fac: u32,
    pub strand: Strand,
    pub self_comparison: bool,
    pub triangle: Triangle,
    pub disable_mirror: bool,
    /// Budget for index, depth tracks and pixel channels together.
    pub memory_limit_bytes: u64,
    pub seed: SpacedSeed,
    /// Two hits on one diagonal at most this far apart survive the filter;
    /// anchors closer than this are merged into one run.
    pub two_hit_dist: u32,
    /// Seeds occurring more often than this in the subject are skipped.
    pub max_occ: u32,
    /// Extra cells scored beyond each run end. `None` = the window size.
    pub extend: Option<u32>,
    /// How the repeat-depth track is estimated. [`DepthMethod::None`]
    /// makes every dot `Unique` unless depth floors are supplied.
    pub depth: DepthMethod,
    pub classes: ClassThresholds,
    /// Chain runs into synteny segments.
    pub chain: Option<ChainConfig>,
    /// Compute only this sub-view; `None` = whole plot.
    pub viewport: Option<Viewport>,
    /// Also count matching bases per pixel ([`LargePlot::coverage`]).
    /// Costs 8 bytes per pixel and channel while computing.
    pub coverage: bool,
    /// Seed sampling rate for the plot (1 = every seed); `None` =
    /// [`auto_sample_rate`] of the zoom. The two-hit distance scales with
    /// the rate.
    pub sample: Option<u32>,
    /// Seed sampling rate for [`DepthMethod::Seed`]; `None` =
    /// [`auto_depth_sample_rate`] of each sequence's length.
    pub depth_sample: Option<u32>,
}

impl LargeConfig {
    pub fn sample_rate(&self) -> u32 {
        self.sample
            .unwrap_or_else(|| auto_sample_rate(self.zoom))
            .max(1)
    }
    pub fn depth_sample_rate(&self, len: usize) -> u32 {
        self.depth_sample
            .unwrap_or_else(|| auto_depth_sample_rate(len))
            .max(1)
    }
}

impl Default for LargeConfig {
    fn default() -> Self {
        Self {
            matrix: ScoreMatrix::dna_identity(),
            window_size: None,
            zoom: 1,
            pixel_fac: 0,
            strand: Strand::Both,
            self_comparison: false,
            triangle: Triangle::Both,
            disable_mirror: false,
            memory_limit_bytes: 4 << 30,
            seed: SpacedSeed::default(),
            two_hit_dist: 100,
            max_occ: 1000,
            extend: None,
            depth: DepthMethod::Seed,
            classes: ClassThresholds::default(),
            chain: None,
            viewport: None,
            coverage: false,
            sample: None,
            depth_sample: None,
        }
    }
}

/// Raw inputs. `*_breaks` are record start offsets excluding 0 (as
/// returned by `dottir_io::Sequence::breaks`). `*_floor`, if given, raises
/// the per-base depth (soft-mask / annotation-derived repeats).
#[derive(Debug, Clone, Copy)]
pub struct LargeInput<'a> {
    pub query: &'a [u8],
    pub q_breaks: &'a [usize],
    pub q_floor: Option<&'a [u16]>,
    /// `None` = self-comparison against the query.
    pub subject: Option<&'a [u8]>,
    pub s_breaks: &'a [usize],
    pub s_floor: Option<&'a [u16]>,
}

impl<'a> LargeInput<'a> {
    pub fn pair(
        query: &'a [u8],
        q_breaks: &'a [usize],
        subject: &'a [u8],
        s_breaks: &'a [usize],
    ) -> Self {
        Self {
            query,
            q_breaks,
            q_floor: None,
            subject: Some(subject),
            s_breaks,
            s_floor: None,
        }
    }
    pub fn self_comparison(seq: &'a [u8], breaks: &'a [usize]) -> Self {
        Self {
            query: seq,
            q_breaks: breaks,
            q_floor: None,
            subject: None,
            s_breaks: breaks,
            s_floor: None,
        }
    }
}

/// Cancellation and progress hooks. Progress is reported as a fraction in
/// `[0, 1]` from worker threads.
#[derive(Default, Clone, Copy)]
pub struct Control<'a> {
    pub cancel: Option<&'a AtomicBool>,
    pub progress: Option<&'a (dyn Fn(f32) + Sync)>,
}

impl std::fmt::Debug for Control<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Control")
            .field("cancel", &self.cancel.map(|c| c.load(Ordering::Relaxed)))
            .field("progress", &self.progress.is_some())
            .finish()
    }
}

impl Control<'_> {
    fn cancelled(&self) -> bool {
        self.cancel.is_some_and(|c| c.load(Ordering::Relaxed))
    }
    fn check(&self) -> Result<(), DottirError> {
        if self.cancelled() {
            Err(DottirError::Cancelled)
        } else {
            Ok(())
        }
    }
    fn report(&self, f: f32) {
        if let Some(p) = self.progress {
            p(f.clamp(0.0, 1.0));
        }
    }
}

/// A diagonal run of anchors, in concatenated 0-based half-open
/// coordinates. Reverse runs pair `q_start` with `s_end - 1`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Segment {
    pub reverse: bool,
    pub q_rec: u32,
    pub s_rec: u32,
    pub q_start: u32,
    pub q_end: u32,
    pub s_start: u32,
    pub s_end: u32,
    pub n_anchors: u32,
    pub class: RepeatClass,
}

/// A set of repeat classes, for showing/hiding channels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ClassMask(u8);

impl ClassMask {
    pub const ALL: ClassMask = ClassMask(0b1111);
    pub const NONE: ClassMask = ClassMask(0);
    pub fn contains(self, c: RepeatClass) -> bool {
        self.0 & (1 << c.index()) != 0
    }
    pub fn with(self, c: RepeatClass, on: bool) -> Self {
        if on {
            ClassMask(self.0 | (1 << c.index()))
        } else {
            ClassMask(self.0 & !(1 << c.index()))
        }
    }
    pub fn from_classes(cs: &[RepeatClass]) -> Self {
        cs.iter().fold(Self::NONE, |m, &c| m.with(c, true))
    }
    pub fn classes(self) -> Vec<RepeatClass> {
        RepeatClass::ALL
            .into_iter()
            .filter(|&c| self.contains(c))
            .collect()
    }
}

impl Default for ClassMask {
    fn default() -> Self {
        Self::ALL
    }
}

#[derive(Debug, Clone)]
pub struct LargeParams {
    pub window_size: u32,
    pub pixel_fac: u32,
    pub karlin: Option<KarlinResult>,
    pub seed: String,
    pub two_hit_dist: u32,
    pub max_occ: u32,
    pub extend: u32,
    pub depth: DepthMethod,
    pub classes: ClassThresholds,
    /// Seed sampling rate used for the plot.
    pub sample: u32,
    /// Seed sampling rates used for the query / subject depth tracks
    /// (1 unless the depth method is `Seed`).
    pub depth_sample: (u32, u32),
    pub format_version: u32,
}

/// Copy-number summary of the bases under one pixel column (query) or
/// row (subject).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct DepthBin {
    /// Median copy number; 1 = unique.
    pub median: u16,
    pub max: u16,
    /// Share of bases in each [`RepeatClass`], scaled to 0–255.
    pub class_frac: [u8; RepeatClass::COUNT],
}

/// Coverage values are stored as `round(8 · log2(1 + bases))`: 8 steps per
/// doubling (±4.5 %), saturating at 255 ≈ 3.8 Gb per pixel.
pub fn coverage_code(bases: f64) -> u8 {
    if bases <= 0.0 {
        return 0;
    }
    ((1.0 + bases).log2() * 8.0).round().min(255.0) as u8
}

/// Inverse of [`coverage_code`].
pub fn coverage_bases(code: u8) -> f64 {
    (code as f64 / 8.0).exp2() - 1.0
}

/// How [`LargePlot::to_dot_plot`] fills the pixels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Shade {
    /// Best window score, as in the exact engine.
    #[default]
    Score,
    /// Matching bases per pixel, as [`coverage_code`]s. Falls back to
    /// `Score` when coverage was not computed.
    Coverage,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct LargeStats {
    pub raw_hits: u64,
    pub anchors: u64,
    pub runs: u64,
    pub scored_cells: u64,
    pub capped_query_positions: u64,
}

/// Result of a large-mode computation.
#[derive(Debug, Clone)]
pub struct LargePlot {
    pub width: u32,
    pub height: u32,
    pub zoom: u32,
    /// Base coordinates of pixel (0, 0); multiples of `zoom`.
    pub origin_q: usize,
    pub origin_s: usize,
    /// `channels[strand * RepeatClass::COUNT + class]`, strand 0 forward,
    /// 1 reverse; a strand that was not computed has empty channels.
    pub channels: Vec<Vec<u8>>,
    /// Per pixel column: query positions whose seed lookup was capped.
    pub q_capped: Vec<u32>,
    /// Per pixel row: subject positions whose seed exceeds `max_occ`.
    pub s_capped: Vec<u32>,
    /// Same layout as `channels`: [`coverage_code`] of the number of
    /// diagonal cells with a positive window score in the pixel. Empty
    /// unless [`LargeConfig::coverage`] was set.
    pub coverage: Vec<Vec<u8>>,
    /// Copy-number profile per pixel column / row.
    pub q_depth_bins: Vec<DepthBin>,
    pub s_depth_bins: Vec<DepthBin>,
    /// Diagonal runs (kept only when chaining is enabled).
    pub segments: Vec<Segment>,
    pub chains: Vec<Chain>,
    pub params: LargeParams,
    pub stats: LargeStats,
    /// Whether the plot came from a self-comparison.
    pub self_comparison: bool,
}

impl LargePlot {
    pub fn channel(&self, reverse: bool, class: RepeatClass) -> &[u8] {
        &self.channels[reverse as usize * RepeatClass::COUNT + class.index()]
    }

    /// Max-merge of the selected classes on one strand (`Some(reverse)`)
    /// or both (`None`).
    pub fn merged(&self, reverse: Option<bool>, classes: ClassMask) -> Vec<u8> {
        let n = self.width as usize * self.height as usize;
        let mut out = vec![0u8; n];
        let strands: &[bool] = match reverse {
            Some(false) => &[false],
            Some(true) => &[true],
            None => &[false, true],
        };
        for &r in strands {
            for c in classes.classes() {
                let ch = self.channel(r, c);
                if !ch.is_empty() {
                    crate::pixel::merge_max_into(&mut out, ch);
                }
            }
        }
        out
    }

    pub fn has_coverage(&self) -> bool {
        self.coverage.iter().any(|c| !c.is_empty())
    }

    /// Coverage of the selected classes on one strand or both, summed in
    /// bases and re-encoded with [`coverage_code`].
    pub fn merged_coverage(&self, reverse: Option<bool>, classes: ClassMask) -> Vec<u8> {
        let n = self.width as usize * self.height as usize;
        let strands: &[bool] = match reverse {
            Some(false) => &[false],
            Some(true) => &[true],
            None => &[false, true],
        };
        let chans: Vec<&[u8]> = strands
            .iter()
            .flat_map(|&r| {
                classes.classes().into_iter().map(move |c| {
                    self.coverage
                        .get(r as usize * RepeatClass::COUNT + c.index())
                        .map_or(&[][..], |v| v.as_slice())
                })
            })
            .filter(|c| !c.is_empty())
            .collect();
        match chans.as_slice() {
            [] => vec![0; n],
            [one] => one.to_vec(),
            _ => {
                let decode: Vec<f64> = (0..=255u8).map(coverage_bases).collect();
                (0..n)
                    .map(|i| coverage_code(chans.iter().map(|c| decode[c[i] as usize]).sum()))
                    .collect()
            }
        }
    }

    /// Present the selected classes as a [`DotPlot`], so the existing
    /// renderers and exporters can consume it.
    pub fn to_dot_plot(&self, classes: ClassMask, separate_strands: bool, shade: Shade) -> DotPlot {
        let merged = |r: Option<bool>| match shade {
            Shade::Coverage if self.has_coverage() => self.merged_coverage(r, classes),
            _ => self.merged(r, classes),
        };
        let (pixels, reverse_pixels) = if separate_strands {
            (merged(Some(false)), Some(merged(Some(true))))
        } else {
            (merged(None), None)
        };
        DotPlot {
            width: self.width,
            height: self.height,
            pixels,
            reverse_pixels,
            params: PlotParams {
                window_size: self.params.window_size,
                zoom: self.zoom,
                pixel_fac: self.params.pixel_fac,
                karlin: self.params.karlin,
            },
        }
    }
}

/// Viewport-independent state: encoded sequences, the seed index, depth
/// tracks and the resolved window/pixel factor.
#[derive(Debug, Clone)]
pub struct LargeSession {
    q_enc: Vec<u8>,
    s_enc: Option<Vec<u8>>,
    q_breaks: Vec<usize>,
    s_breaks: Vec<usize>,
    index: SeedIndex,
    q_depth: DepthTrack,
    s_depth: Option<DepthTrack>,
    karlin: Option<KarlinResult>,
    window: u32,
    pixel_fac: u32,
    matrix: ScoreMatrix,
    depth: DepthMethod,
    /// `(two_hit_dist, max_occ, depth_sample)` the seed depth was
    /// computed with.
    depth_params: (u32, u32, Option<u32>),
    /// Resolved depth sampling rates (query, subject).
    depth_rates: (u32, u32),
}

impl LargeSession {
    /// Build the index and depth tracks. Expensive: linear in input size,
    /// dominated by the depth estimate when `cfg.depth` is not `None`.
    pub fn prepare(
        input: &LargeInput<'_>,
        cfg: &LargeConfig,
        ctl: &Control<'_>,
    ) -> Result<Self, DottirError> {
        if cfg.matrix.kind != AlphabetKind::Dna {
            return Err(DottirError::InvalidConfig(
                "large mode supports DNA only".into(),
            ));
        }
        cfg.classes.validate()?;
        let query = input.query;
        let subject = input.subject.unwrap_or(query);
        if query.is_empty() || subject.is_empty() {
            return Err(DottirError::EmptySequence);
        }
        let total = query.len() as u64 + subject.len() as u64 + cfg.seed.span() as u64;
        if total >= u32::MAX as u64 {
            return Err(DottirError::InvalidConfig(format!(
                "large mode supports at most 4 Gb of query + subject, got {total} bases"
            )));
        }
        if cfg.self_comparison != input.subject.is_none() {
            return Err(DottirError::InvalidConfig(
                "self_comparison must match LargeInput::self_comparison".into(),
            ));
        }

        let pairwise_len = if input.subject.is_some() {
            subject.len() as u64
        } else {
            0
        };
        let longest = query.len().max(subject.len()) as u64;
        // Peak during prepare: subject index, encodings + depth tracks, and
        // for seed depth the coverage counters plus (pairwise) a temporary
        // query self-index.
        let depth_peak = match cfg.depth {
            DepthMethod::Seed => {
                4 * longest
                    + if input.subject.is_some() {
                        SeedIndex::estimate_bytes(&cfg.seed, query.len())
                    } else {
                        0
                    }
            }
            DepthMethod::Kmer { .. } => 16 * longest.min(1 << 25),
            DepthMethod::None => 0,
        };
        let est = SeedIndex::estimate_bytes(&cfg.seed, subject.len())
            + 3 * (query.len() as u64 + pairwise_len)
            + depth_peak;
        if est > cfg.memory_limit_bytes {
            return Err(DottirError::OutOfMemory {
                requested: est,
                per_channel: est,
                channels: 1,
                limit: cfg.memory_limit_bytes,
            });
        }

        let need_karlin = cfg.window_size.is_none() || cfg.pixel_fac == 0;
        let karlin = if need_karlin {
            Some(karlin_window_size(
                &cfg.matrix,
                query,
                subject,
                BlastMode::Blastn,
            )?)
        } else {
            None
        };
        let window = cfg
            .window_size
            .or(karlin.map(|k| k.window_size))
            .unwrap_or(1);
        if window < 1 {
            return Err(DottirError::InvalidConfig(
                "window size must be >= 1".into(),
            ));
        }
        let pixel_fac = if cfg.pixel_fac == 0 {
            let exp_res = karlin.map(|k| k.expected_residue_score).unwrap_or(0.0);
            if !(exp_res.is_finite() && exp_res > 0.0) {
                return Err(DottirError::InvalidConfig(format!(
                    "auto pixel_fac requires positive Karlin expected_residue_score, got {exp_res:.3}"
                )));
            }
            ((0.2 * 256.0 / exp_res).round() as u32).max(1)
        } else {
            cfg.pixel_fac
        };
        ctl.check()?;

        let q_enc = encode(query, AlphabetKind::Dna);
        let s_enc = input.subject.map(|s| encode(s, AlphabetKind::Dna));
        let s_ref = s_enc.as_deref().unwrap_or(&q_enc);
        ctl.report(0.05);
        let index = SeedIndex::build(&cfg.seed, s_ref, input.s_breaks)?;
        ctl.report(0.3);
        ctl.check()?;

        let depth_for = |enc: &[u8],
                         breaks: &[usize],
                         floor: Option<&[u16]>,
                         self_index: Option<&SeedIndex>|
         -> Result<DepthTrack, DottirError> {
            let mut d = match cfg.depth {
                DepthMethod::None => DepthTrack::flat(enc.len()),
                DepthMethod::Kmer { k } => compute_depth(enc, breaks, k)?,
                DepthMethod::Seed => {
                    let built;
                    let idx = match self_index {
                        Some(i) => i,
                        None => {
                            built = SeedIndex::build(&cfg.seed, enc, breaks)?;
                            &built
                        }
                    };
                    let rate = cfg.depth_sample_rate(enc.len());
                    seed_depth(idx, enc, breaks, cfg.two_hit_dist, cfg.max_occ, rate, ctl)?
                }
            };
            if let Some(f) = floor {
                d.apply_floor(f)?;
            }
            Ok(d)
        };
        // The subject index doubles as the self index of the subject (and,
        // in a self-comparison, of the query).
        let q_self_index = s_enc.is_none().then_some(&index);
        let q_depth = depth_for(&q_enc, input.q_breaks, input.q_floor, q_self_index)?;
        ctl.report(0.65);
        ctl.check()?;
        let s_depth = match &s_enc {
            Some(s) => Some(depth_for(s, input.s_breaks, input.s_floor, Some(&index))?),
            None => None,
        };
        ctl.report(1.0);
        let rate_of = |len: usize| {
            if cfg.depth == DepthMethod::Seed {
                cfg.depth_sample_rate(len)
            } else {
                1
            }
        };
        let depth_rates = (
            rate_of(q_enc.len()),
            rate_of(s_enc.as_ref().map_or(q_enc.len(), |s| s.len())),
        );

        Ok(Self {
            q_enc,
            s_enc,
            q_breaks: input.q_breaks.to_vec(),
            s_breaks: input.s_breaks.to_vec(),
            index,
            q_depth,
            s_depth,
            karlin,
            window,
            pixel_fac,
            matrix: cfg.matrix.clone(),
            depth: cfg.depth,
            depth_params: (cfg.two_hit_dist, cfg.max_occ, cfg.depth_sample),
            depth_rates,
        })
    }

    /// True if `cfg` can be computed from this session without re-running
    /// [`Self::prepare`].
    pub fn compatible(&self, cfg: &LargeConfig) -> bool {
        self.index.seed() == &cfg.seed
            && self.depth == cfg.depth
            && (self.depth != DepthMethod::Seed
                || self.depth_params == (cfg.two_hit_dist, cfg.max_occ, cfg.depth_sample))
            && self.matrix.name == cfg.matrix.name
            && self.matrix.scores == cfg.matrix.scores
            && cfg.window_size.is_none_or(|w| w == self.window)
            && (cfg.pixel_fac == 0 || cfg.pixel_fac == self.pixel_fac)
            && self.s_enc.is_none() == cfg.self_comparison
    }

    pub fn query_len(&self) -> usize {
        self.q_enc.len()
    }
    pub fn subject_len(&self) -> usize {
        self.s_enc.as_ref().map_or(self.q_enc.len(), |s| s.len())
    }
    pub fn window(&self) -> u32 {
        self.window
    }
    pub fn pixel_fac(&self) -> u32 {
        self.pixel_fac
    }
    pub fn karlin(&self) -> Option<KarlinResult> {
        self.karlin
    }
    pub fn index(&self) -> &SeedIndex {
        &self.index
    }
    pub fn query_depth(&self) -> &DepthTrack {
        &self.q_depth
    }
    pub fn subject_depth(&self) -> &DepthTrack {
        self.s_depth.as_ref().unwrap_or(&self.q_depth)
    }
    /// Approximate heap footprint.
    pub fn bytes(&self) -> u64 {
        let seq =
            self.q_enc.len() as u64 * 3 + self.s_enc.as_ref().map_or(0, |s| s.len() as u64 * 3);
        self.index.bytes() + seq
    }

    /// Compute a plot (or viewport) with `cfg`, which must be
    /// [`Self::compatible`].
    pub fn compute(&self, cfg: &LargeConfig, ctl: &Control<'_>) -> Result<LargePlot, DottirError> {
        if !self.compatible(cfg) {
            return Err(DottirError::InvalidConfig(
                "LargeConfig changed seed/depth/matrix/window; prepare a new session".into(),
            ));
        }
        if cfg.zoom == 0 {
            return Err(DottirError::InvalidConfig("zoom must be >= 1".into()));
        }
        cfg.classes.validate()?;
        let qlen = self.query_len();
        let slen = self.subject_len();
        let zoom = cfg.zoom as usize;

        let (q_view, s_view, is_viewport) = match &cfg.viewport {
            Some(v) => {
                let q0 = v.q.start.min(qlen) / zoom * zoom;
                let s0 = v.s.start.min(slen) / zoom * zoom;
                let q1 = v.q.end.clamp(q0, qlen);
                let s1 = v.s.end.clamp(s0, slen);
                if q0 >= q1 || s0 >= s1 {
                    return Err(DottirError::InvalidConfig("empty viewport".into()));
                }
                let full = q0 == 0 && s0 == 0 && q1 == qlen && s1 == slen;
                (q0..q1, s0..s1, !full)
            }
            None => (0..qlen, 0..slen, false),
        };
        let width = image_dimension(q_view.len(), cfg.zoom);
        let height = image_dimension(s_view.len(), cfg.zoom);

        let do_fwd = matches!(cfg.strand, Strand::Forward | Strand::Both);
        let do_rev = matches!(cfg.strand, Strand::Reverse | Strand::Both);
        let n_strands = do_fwd as u64 + do_rev as u64;
        let per_channel = width as u64 * height as u64;
        let channel_bytes = per_channel
            * n_strands
            * RepeatClass::COUNT as u64
            * if cfg.coverage { 1 + 8 } else { 1 };
        let requested = channel_bytes + self.bytes();
        if requested > cfg.memory_limit_bytes {
            return Err(DottirError::OutOfMemory {
                requested,
                per_channel,
                channels: (n_strands as u32) * RepeatClass::COUNT as u32,
                limit: cfg.memory_limit_bytes,
            });
        }

        // Self-comparison viewports compute both triangles directly; only
        // the full plot is computed lower-only and mirrored (as the exact
        // engine does).
        let lower_only = cfg.self_comparison && !is_viewport;
        let window = self.window as usize;
        let span = cfg.seed.span();
        let extend = cfg.extend.map_or(window, |e| e as usize);
        let sample = cfg.sample_rate();
        let dist = cfg.two_hit_dist.saturating_mul(sample) as usize;
        let margin = extend + 2 * window + span + 2 * dist;

        let new_maps = |on: bool| -> Result<Vec<PixelMap>, DottirError> {
            if !on {
                return Ok(vec![]);
            }
            (0..RepeatClass::COUNT)
                .map(|_| PixelMap::new_checked(width, height, cfg.memory_limit_bytes))
                .collect()
        };
        let mut fwd_maps = new_maps(do_fwd)?;
        let mut rev_maps = new_maps(do_rev)?;

        let s_enc = self.s_enc.as_deref().unwrap_or(&self.q_enc);
        let s_depth = self.subject_depth();
        let hit_ctx = hits::HitCtx {
            index: &self.index,
            q_enc: &self.q_enc,
            q_breaks: &self.q_breaks,
            s_breaks: &self.s_breaks,
            q_depth: &self.q_depth.depth,
            s_depth: &s_depth.depth,
            thresholds: cfg.classes,
            s_range: s_view.start.saturating_sub(margin)..(s_view.end + margin).min(slen),
            max_occ: cfg.max_occ,
            two_hit_dist: cfg.two_hit_dist.saturating_mul(sample),
            sample,
            stream_runs: true,
            forward: do_fwd,
            reverse: do_rev,
            lower_only,
        };
        let table = score::PairTable::new(&cfg.matrix);
        let score_ctx = score::ScoreCtx {
            table: &table,
            q_enc: &self.q_enc,
            s_enc,
            q_depth: &self.q_depth.depth,
            s_depth: &s_depth.depth,
            thresholds: cfg.classes,
            window,
            pixel_fac: self.pixel_fac,
            span,
            extend,
            lower_only,
        };

        let q_seed_range = q_view.start.saturating_sub(margin)..(q_view.end + margin).min(qlen);
        // Large chunks keep the stream overlap (`reach` below) small.
        const CHUNK: usize = 1 << 16;
        let chunks: Vec<Range<usize>> = q_seed_range
            .clone()
            .step_by(CHUNK)
            .map(|lo| lo..(lo + CHUNK).min(q_seed_range.end))
            .collect();
        let n_chunks = chunks.len().max(1);
        let done = AtomicUsize::new(0);
        let keep_segments = cfg.chain.is_some();
        let chain_max_class = cfg.chain.as_ref().map(|c| c.max_class);

        struct ChunkOut {
            segments: Vec<Segment>,
            capped: Vec<u32>,
            stats: LargeStats,
        }

        let cov_bufs: Vec<Vec<AtomicU64>> = if cfg.coverage {
            [do_fwd, do_rev]
                .into_iter()
                .flat_map(|on| {
                    (0..RepeatClass::COUNT).map(move |_| {
                        if on {
                            (0..per_channel).map(|_| AtomicU64::new(0)).collect()
                        } else {
                            Vec::new()
                        }
                    })
                })
                .collect()
        } else {
            Vec::new()
        };

        let (outs, s_capped) = {
            let fwd_views: Vec<_> = fwd_maps.iter_mut().map(|m| m.view_mut()).collect();
            let rev_views: Vec<_> = rev_maps.iter_mut().map(|m| m.view_mut()).collect();
            let canvas = score::Canvas {
                q0: q_view.start,
                s0: s_view.start,
                zoom,
                width: width as usize,
                height: height as usize,
                fwd: &fwd_views,
                rev: &rev_views,
                coverage: cfg.coverage.then_some(cov_bufs.as_slice()),
            };
            let work = |core: &Range<usize>| -> ChunkOut {
                if ctl.cancelled() {
                    return ChunkOut {
                        segments: vec![],
                        capped: vec![],
                        stats: LargeStats::default(),
                    };
                }
                // Wide enough that every anchor deciding whether a core
                // cell is covered is seen (see `score::score_owned`).
                let reach = extend + span + window + 2 * dist;
                let stream = core.start.saturating_sub(reach)..(core.end + reach).min(qlen);
                let ch = hits::chunk_hits(&hit_ctx, core.clone(), stream);
                let mut stats = LargeStats {
                    raw_hits: ch.stats.raw_hits,
                    anchors: ch.stats.anchors,
                    runs: ch.runs.len() as u64,
                    scored_cells: 0,
                    capped_query_positions: ch.stats.capped_positions,
                };
                stats.scored_cells +=
                    score::score_owned(&score_ctx, &canvas, &ch.stream_runs, core);
                let mut segments = Vec::new();
                for run in &ch.runs {
                    if keep_segments && chain_max_class.is_some_and(|m| run.class <= m) {
                        segments.push(run_to_segment(run, qlen, span));
                    }
                }
                let n = done.fetch_add(1, Ordering::Relaxed) + 1;
                ctl.report(n as f32 / n_chunks as f32);
                ChunkOut {
                    segments,
                    capped: ch.capped,
                    stats,
                }
            };

            #[cfg(feature = "rayon")]
            let outs: Vec<ChunkOut> = {
                use rayon::prelude::*;
                chunks.par_iter().map(work).collect()
            };
            #[cfg(not(feature = "rayon"))]
            let outs: Vec<ChunkOut> = chunks.iter().map(work).collect();

            let s_capped = self.subject_capped(cfg, &s_view, zoom, height as usize);
            (outs, s_capped)
        };
        ctl.check()?;

        let mut stats = LargeStats::default();
        let mut q_capped = vec![0u32; width as usize];
        let mut segments = Vec::new();
        for o in outs {
            stats.raw_hits += o.stats.raw_hits;
            stats.anchors += o.stats.anchors;
            stats.runs += o.stats.runs;
            stats.scored_cells += o.stats.scored_cells;
            for p in o.capped {
                let p = p as usize;
                if q_view.contains(&p) {
                    q_capped[(p - q_view.start) / zoom] += 1;
                    stats.capped_query_positions += 1;
                }
            }
            segments.extend(o.segments);
        }

        let mut channels: Vec<Vec<u8>> = Vec::with_capacity(2 * RepeatClass::COUNT);
        for maps in [fwd_maps, rev_maps] {
            if maps.is_empty() {
                channels.extend((0..RepeatClass::COUNT).map(|_| Vec::new()));
            } else {
                channels.extend(maps.into_iter().map(PixelMap::into_vec));
            }
        }
        let mut coverage: Vec<Vec<u8>> = cov_bufs
            .into_iter()
            .map(|b| {
                b.into_iter()
                    .map(|n| coverage_code(n.into_inner() as f64))
                    .collect()
            })
            .collect();
        if lower_only {
            for ch in channels
                .iter_mut()
                .chain(coverage.iter_mut())
                .filter(|c| !c.is_empty())
            {
                mirror_self_comparison(
                    ch,
                    width as usize,
                    height as usize,
                    cfg.triangle,
                    cfg.disable_mirror,
                );
            }
        }

        let chains = match &cfg.chain {
            Some(cc) => chain_segments(&segments, cc),
            None => vec![],
        };

        let q_depth_bins = depth_bins(&self.q_depth.depth, &q_view, zoom, cfg.classes);
        let s_depth_bins = depth_bins(&s_depth.depth, &s_view, zoom, cfg.classes);

        Ok(LargePlot {
            width,
            height,
            zoom: cfg.zoom,
            origin_q: q_view.start,
            origin_s: s_view.start,
            channels,
            q_capped,
            s_capped,
            coverage,
            q_depth_bins,
            s_depth_bins,
            segments,
            chains,
            params: LargeParams {
                window_size: self.window,
                pixel_fac: self.pixel_fac,
                karlin: self.karlin,
                seed: cfg.seed.pattern().to_string(),
                two_hit_dist: cfg.two_hit_dist,
                max_occ: cfg.max_occ,
                extend: extend as u32,
                depth: cfg.depth,
                classes: cfg.classes,
                sample,
                depth_sample: self.depth_rates,
                format_version: LARGE_PIXELMAP_FORMAT_VERSION,
            },
            stats,
            self_comparison: cfg.self_comparison,
        })
    }

    /// Per pixel row, subject positions in `s_view` whose seed occurs more
    /// than `max_occ` times.
    fn subject_capped(
        &self,
        cfg: &LargeConfig,
        s_view: &Range<usize>,
        zoom: usize,
        height: usize,
    ) -> Vec<u32> {
        let s_enc = self.s_enc.as_deref().unwrap_or(&self.q_enc);
        let seed = self.index.seed();
        let scan = |r: Range<usize>| {
            let mut rows = vec![];
            seed::for_each_window(seed, s_enc, &self.s_breaks, r, |p, fwd, _| {
                if self.index.count(seed.key(fwd)) > cfg.max_occ {
                    rows.push(((p - s_view.start) / zoom) as u32);
                }
            });
            rows
        };
        const CHUNK: usize = 1 << 20;
        let ranges: Vec<Range<usize>> = s_view
            .clone()
            .step_by(CHUNK)
            .map(|lo| lo..(lo + CHUNK).min(s_view.end))
            .collect();
        #[cfg(feature = "rayon")]
        let parts: Vec<Vec<u32>> = {
            use rayon::prelude::*;
            ranges.into_par_iter().map(scan).collect()
        };
        #[cfg(not(feature = "rayon"))]
        let parts: Vec<Vec<u32>> = ranges.into_iter().map(scan).collect();
        let mut out = vec![0u32; height];
        for row in parts.into_iter().flatten() {
            if let Some(c) = out.get_mut(row as usize) {
                *c += 1;
            }
        }
        out
    }
}

/// Per-base copy number from a self-comparison two-hit pass over `index`
/// (which must index `enc` itself): one plus the number of off-diagonal
/// runs covering each base. Bases whose seed exceeds `max_occ` get
/// `u16::MAX`.
#[allow(clippy::too_many_arguments)]
fn seed_depth(
    index: &SeedIndex,
    enc: &[u8],
    breaks: &[usize],
    two_hit_dist: u32,
    max_occ: u32,
    sample: u32,
    ctl: &Control<'_>,
) -> Result<DepthTrack, DottirError> {
    let two_hit_dist = two_hit_dist.saturating_mul(sample);
    use std::sync::atomic::AtomicI32;
    let n = enc.len();
    let span = index.seed().span();
    let diff: Vec<AtomicI32> = (0..=n).map(|_| AtomicI32::new(0)).collect();
    let add = |r: Range<usize>| {
        diff[r.start.min(n)].fetch_add(1, Ordering::Relaxed);
        diff[r.end.min(n)].fetch_sub(1, Ordering::Relaxed);
    };
    let ctx = hits::HitCtx {
        index,
        q_enc: enc,
        q_breaks: breaks,
        s_breaks: breaks,
        q_depth: &[],
        s_depth: &[],
        thresholds: ClassThresholds::default(),
        s_range: 0..n,
        max_occ,
        two_hit_dist,
        sample,
        stream_runs: false,
        forward: true,
        reverse: true,
        lower_only: true,
    };
    let dist = two_hit_dist as usize;
    const CHUNK: usize = 1 << 14;
    let chunks: Vec<Range<usize>> = (0..n)
        .step_by(CHUNK)
        .map(|lo| lo..(lo + CHUNK).min(n))
        .collect();
    let work = |core: &Range<usize>| -> Vec<u32> {
        if ctl.cancelled() {
            return vec![];
        }
        let stream = core.start.saturating_sub(dist)..(core.end + dist).min(n);
        let ch = hits::chunk_hits(&ctx, core.clone(), stream);
        for run in &ch.runs {
            let (first, last) = (run.first_p as usize, run.last_p as usize);
            if run.reverse {
                // s = diag + 1 - span - p, decreasing along the run.
                let c = run.diag as usize + 1 - span;
                add(first..last + span);
                add(c - last..c - first + span);
            } else {
                let off = run.diag as usize - n;
                if off == 0 {
                    continue; // the main diagonal: every base vs itself
                }
                add(first..last + span);
                add(first + off..last + off + span);
            }
        }
        ch.capped
    };
    #[cfg(feature = "rayon")]
    let capped: Vec<Vec<u32>> = {
        use rayon::prelude::*;
        chunks.par_iter().map(work).collect()
    };
    #[cfg(not(feature = "rayon"))]
    let capped: Vec<Vec<u32>> = chunks.iter().map(work).collect();
    ctl.check()?;

    let mut depth = Vec::with_capacity(n);
    let mut acc: i64 = 1;
    for d in diff.iter().take(n) {
        acc += d.load(Ordering::Relaxed) as i64;
        depth.push(acc.clamp(0, u16::MAX as i64) as u16);
    }
    for p in capped.into_iter().flatten() {
        let p = p as usize;
        depth[p..(p + span).min(n)].fill(u16::MAX);
    }
    Ok(DepthTrack {
        method: DepthMethod::Seed,
        depth,
    })
}

fn run_to_segment(run: &hits::Run, qlen: usize, span: usize) -> Segment {
    let q_start = run.first_p;
    let q_end = run.last_p + span as u32;
    let (s_start, s_end) = if run.reverse {
        // q + s = diag over the run's cells.
        let c = run.diag as i64;
        (
            (c - q_end as i64 + 1) as u32,
            (c - q_start as i64 + 1) as u32,
        )
    } else {
        let d = run.diag as i64 - qlen as i64;
        ((q_start as i64 + d) as u32, (q_end as i64 + d) as u32)
    };
    Segment {
        reverse: run.reverse,
        q_rec: run.q_rec,
        s_rec: run.s_rec,
        q_start,
        q_end,
        s_start,
        s_end,
        n_anchors: run.n_anchors,
        class: run.class,
    }
}

/// One-shot convenience: [`LargeSession::prepare`] then
/// Summarise `depth[view]` per `zoom`-base bin. Depth 0 (never seen
/// repeating) counts as one copy.
fn depth_bins(
    depth: &[u16],
    view: &Range<usize>,
    zoom: usize,
    classes: ClassThresholds,
) -> Vec<DepthBin> {
    let bin = |lo: usize| {
        let hi = (lo + zoom).min(view.end);
        let mut v: Vec<u16> = depth[lo..hi].iter().map(|&d| d.max(1)).collect();
        let mut counts = [0usize; RepeatClass::COUNT];
        for &d in &v {
            counts[classes.classify(d).index()] += 1;
        }
        let n = v.len();
        let mid = n / 2;
        let median = *v.select_nth_unstable(mid).1;
        let max = v.iter().copied().max().unwrap_or(1);
        DepthBin {
            median,
            max,
            class_frac: counts.map(|c| ((c * 255 + n / 2) / n) as u8),
        }
    };
    let starts: Vec<usize> = view.clone().step_by(zoom).collect();
    #[cfg(feature = "rayon")]
    {
        use rayon::prelude::*;
        starts.into_par_iter().map(bin).collect()
    }
    #[cfg(not(feature = "rayon"))]
    {
        starts.into_iter().map(bin).collect()
    }
}

/// [`LargeSession::compute`].
pub fn compute_large_dotplot(
    input: &LargeInput<'_>,
    cfg: &LargeConfig,
) -> Result<LargePlot, DottirError> {
    let ctl = Control::default();
    LargeSession::prepare(input, cfg, &ctl)?.compute(cfg, &ctl)
}

#[cfg(test)]
mod tests;
