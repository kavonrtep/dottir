//! Top-level egui application.

use std::path::PathBuf;
use std::sync::Arc;

use dottir_core::large::{
    coverage_bases, ChainConfig, ClassMask, ClassThresholds, DepthBin, DepthMethod, LargeConfig,
    LargePlot, RepeatClass, Shade, SpacedSeed, Viewport,
};
use dottir_core::{
    extract_ridges_from_pixels, snap_zoom_to_period_divisor, BlastMode, DotPlot, PlotConfig, Ridge,
    RidgeDirection, RidgeParams, ScoreMatrix, Strand, Triangle,
};
use dottir_io::faidx::Region;
use dottir_io::{DetectedAlphabet, Sequence};
use egui::{
    Color32, ColorImage, Context, Pos2, Rect, Sense, Slider, TextureHandle, TextureOptions, Vec2,
};

use crate::genome_view::{GenomeAction, GenomeInputs, GenomeLoad, GenomeView};
use dottir_core::genome::SortAxis;
use dottir_io::hsp::HspMode;

/// Startup-time configuration assembled from CLI args, applied once
/// by [`DottirApp::new`].
#[derive(Debug, Clone)]
pub struct StartupConfig {
    pub query: Option<PathBuf>,
    pub subject: Option<PathBuf>,
    pub mode: BlastMode,
    pub matrix_name: String,
    pub window_size: Option<u32>,
    pub zoom: u32,
    pub pixel_fac: u32,
    pub strand: Strand,
    pub self_comparison: bool,
    pub memory_limit_bytes: u64,
    /// Annotation file (GFF3/BED) to bind to the query axis on startup.
    pub annot_query: Option<PathBuf>,
    /// Annotation file (GFF3/BED) to bind to the subject axis on startup.
    pub annot_subject: Option<PathBuf>,
    /// Open a genome comparison (ADR 0008) instead of a dotplot.
    pub genome: Option<GenomeInputs>,
    /// Sort one axis's records once the sequences are loaded.
    pub sort_records: Option<SortAxis>,
    /// PAF drawn over the plot, and how.
    pub hsp_paf: Option<PathBuf>,
    pub hsp_mode: Option<HspMode>,
    /// Plot only these regions (`name:start-end[:-]`) of the inputs.
    pub query_regions: Vec<String>,
    pub subject_regions: Vec<String>,
}

impl Default for StartupConfig {
    fn default() -> Self {
        Self {
            query: None,
            subject: None,
            mode: BlastMode::Blastn,
            matrix_name: "DNA+5/-4".into(),
            window_size: None,
            zoom: 1,
            pixel_fac: 0, // auto (Karlin-derived)
            strand: Strand::Both,
            self_comparison: false,
            memory_limit_bytes: 512 * 1024 * 1024,
            annot_query: None,
            annot_subject: None,
            genome: None,
            sort_records: None,
            hsp_paf: None,
            hsp_mode: None,
            query_regions: Vec::new(),
            subject_regions: Vec::new(),
        }
    }
}

// `Sequence` itself carries the source path and record metadata, so
// the GUI just stores it directly — no wrapper needed.

/// User-tunable plot settings exposed in the Settings panel. Mirrors
/// [`PlotConfig`] with a few GUI-friendly defaults.
#[derive(Clone)]
struct Settings {
    mode: BlastMode,
    matrix_name: String,
    /// `None` = let Karlin/Altschul pick the window size.
    window_size: Option<u32>,
    zoom: u32,
    pixel_fac: u32,
    strand: Strand,
    self_comparison: bool,
    triangle: Triangle,
    /// Memory cap for the pixelmap, in bytes. Plumbed through to
    /// `PlotConfig::memory_limit_bytes` on every recompute. Default
    /// 512 MiB, matching the CLI default.
    memory_limit_bytes: u64,
    /// Pre-process: reverse-complement the query (spec §4.1.10
    /// `-r`). BLASTN only.
    reverse_query: bool,
    /// Pre-process: reverse-complement the subject (spec §4.1.10
    /// `-v`). BLASTN only.
    reverse_subject: bool,
    /// Spec §4.4.3 inverted-repeat highlighting: when on, the
    /// forward+reverse passes write into separate channels and the
    /// renderer paints reverse hits in `reverse_colour`. Only
    /// meaningful for BLASTN Strand::Both.
    inverted_repeat_colour: bool,
    /// H2 alignment-dock window size (residues, centred on the
    /// crosshair). Clamped to [20, 400] at render time.
    align_window_size: u32,
    /// Draw coherent diagonal runs from the pixmap as anti-aliased
    /// vector segments on top of the raster (default on). Hides the
    /// per-window intensity oscillation on imperfect-homology
    /// diagonals — the raster underneath is unchanged so the data
    /// view is still accessible by toggling this off.
    show_ridge_overlay: bool,
    /// Minimum continuous-lit-cells for a diagonal run to be
    /// vectorised. Smaller → more (and shorter) ridges; larger →
    /// only long runs survive.
    ridge_min_length: u32,
    /// Maximum non-lit-cells inside a ridge before it breaks into
    /// two segments. 0 = no bridging; 2 (default) bridges over
    /// occasional mismatch clusters.
    ridge_max_gap: u32,
    /// When true, the GUI picks `PlotConfig::zoom` such that the
    /// computed pixelmap is roughly *physical-canvas-sized* — Dotter's
    /// invariant of computing at display scale and rendering near 1:1.
    /// Fires on (a) sequence load once the canvas has been measured,
    /// and (b) window resizes that materially change the plot area
    /// (`>= 25%` after a 200 ms idle). When off, `settings.zoom` is
    /// respected verbatim (advanced users / reproducibility). Default
    /// on. See `docs/rendering_review.md` for the rationale.
    auto_fit_compute_zoom: bool,
    /// When true, `pixel_fac` is sent as 0 to the core (auto-derive
    /// from Karlin's expected residue score, dotter's default
    /// behaviour). When false, the slider value is used directly. The
    /// last resolved value lives on `DottirApp::resolved_pixel_fac` so
    /// the slider can seed itself sensibly when the user toggles auto
    /// off mid-session.
    auto_pixel_fac: bool,
    /// Large-sequence mode (ADR 0007).
    large: LargeSettings,
}

/// Which engine computes a view.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum EngineChoice {
    /// Large engine above `threshold_cells` visible cells, exact below.
    Auto,
    Exact,
    Large,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum DepthChoice {
    Seed,
    Kmer,
    None,
}

/// Large-mode settings. The first group needs a recompute (and, for
/// seed/depth/floors, a re-index); the display group does not.
#[derive(Clone)]
struct LargeSettings {
    engine: EngineChoice,
    /// Visible matrix cells above which `Auto` uses the large engine.
    threshold_cells: f64,
    seed_pattern: String,
    two_hit_dist: u32,
    max_occ: u32,
    depth: DepthChoice,
    depth_k: u8,
    thresholds: ClassThresholds,
    /// Soft-masked (lowercase) bases count as at least `mid`.
    softmask_as_repeat: bool,
    memory_limit_bytes: u64,
    /// Seed sampling for the plot / the seed depth pass; `None` = auto.
    seed_sample: Option<u32>,
    depth_sample: Option<u32>,
    // Display only.
    classes: ClassMask,
    colour_by_class: bool,
    show_chains: bool,
    chain_min_score: f32,
    /// What the pixel darkness means in large mode.
    shade: Shade,
    /// Coverage ramp in matching bases per pixel: white below, black at.
    cov_white_bp: f32,
    cov_black_bp: f32,
    /// Re-fit the coverage ramp to each new plot.
    cov_auto: bool,
    show_axis_tracks: bool,
}

impl Default for LargeSettings {
    fn default() -> Self {
        Self {
            engine: EngineChoice::Auto,
            threshold_cells: 2e10,
            seed_pattern: dottir_core::large::seed::DEFAULT_PATTERN.to_string(),
            two_hit_dist: 100,
            max_occ: 1000,
            depth: DepthChoice::Seed,
            depth_k: 21,
            thresholds: ClassThresholds::default(),
            softmask_as_repeat: true,
            memory_limit_bytes: 6 << 30,
            seed_sample: None,
            depth_sample: None,
            classes: ClassMask::ALL,
            colour_by_class: true,
            show_chains: false,
            chain_min_score: 1000.0,
            shade: Shade::Coverage,
            cov_white_bp: 50.0,
            cov_black_bp: 3000.0,
            cov_auto: true,
            show_axis_tracks: true,
        }
    }
}

/// Height (top) / width (left) of the dot-profile + capped-seed strip.
const AXIS_TRACK: f32 = 16.0;
/// Aligned copies per base at which the axis profile reaches full height.
const AXIS_TRACK_MAX_COPIES: f32 = 1000.0;

impl LargeSettings {
    fn to_session(&self) -> crate::session::SessionLarge {
        crate::session::SessionLarge {
            engine: match self.engine {
                EngineChoice::Auto => "auto",
                EngineChoice::Exact => "exact",
                EngineChoice::Large => "large",
            }
            .into(),
            threshold_cells: self.threshold_cells,
            seed_pattern: self.seed_pattern.clone(),
            two_hit_dist: self.two_hit_dist,
            max_occ: self.max_occ,
            depth: match self.depth {
                DepthChoice::Seed => "seed",
                DepthChoice::Kmer => "kmer",
                DepthChoice::None => "none",
            }
            .into(),
            depth_k: self.depth_k,
            class_thresholds: [
                self.thresholds.low,
                self.thresholds.mid,
                self.thresholds.high,
            ],
            softmask_as_repeat: self.softmask_as_repeat,
            memory_limit_mib: self.memory_limit_bytes >> 20,
            shown_classes: self
                .classes
                .classes()
                .iter()
                .map(|c| c.name().to_string())
                .collect(),
            colour_by_class: self.colour_by_class,
            show_chains: self.show_chains,
            chain_min_score: self.chain_min_score,
            shade: Some(
                match self.shade {
                    Shade::Coverage => "coverage",
                    Shade::Score => "score",
                }
                .into(),
            ),
            coverage_ramp_bp: Some([self.cov_white_bp, self.cov_black_bp]),
            coverage_auto: Some(self.cov_auto),
            axis_tracks: Some(self.show_axis_tracks),
            seed_sample: Some(self.seed_sample.unwrap_or(0)),
            depth_sample: Some(self.depth_sample.unwrap_or(0)),
        }
    }

    /// Unknown or invalid values fall back to the defaults.
    fn from_session(s: &crate::session::SessionLarge) -> Self {
        let d = Self::default();
        let [low, mid, high] = s.class_thresholds;
        let thresholds = ClassThresholds { low, mid, high };
        Self {
            engine: match s.engine.as_str() {
                "exact" => EngineChoice::Exact,
                "large" => EngineChoice::Large,
                _ => EngineChoice::Auto,
            },
            threshold_cells: s.threshold_cells,
            seed_pattern: if SpacedSeed::parse(&s.seed_pattern).is_ok() {
                s.seed_pattern.clone()
            } else {
                d.seed_pattern
            },
            two_hit_dist: s.two_hit_dist,
            max_occ: s.max_occ,
            depth: match s.depth.as_str() {
                "kmer" => DepthChoice::Kmer,
                "none" => DepthChoice::None,
                _ => DepthChoice::Seed,
            },
            depth_k: s.depth_k,
            thresholds: if thresholds.validate().is_ok() {
                thresholds
            } else {
                d.thresholds
            },
            softmask_as_repeat: s.softmask_as_repeat,
            memory_limit_bytes: s.memory_limit_mib << 20,
            classes: ClassMask::from_classes(
                &s.shown_classes
                    .iter()
                    .filter_map(|n| RepeatClass::from_name(n))
                    .collect::<Vec<_>>(),
            ),
            colour_by_class: s.colour_by_class,
            show_chains: s.show_chains,
            chain_min_score: s.chain_min_score,
            shade: match s.shade.as_deref() {
                Some("score") => Shade::Score,
                Some(_) => Shade::Coverage,
                None => d.shade,
            },
            cov_white_bp: s.coverage_ramp_bp.map_or(d.cov_white_bp, |r| r[0].max(1.0)),
            cov_black_bp: s
                .coverage_ramp_bp
                .map_or(d.cov_black_bp, |r| r[1].max(r[0].max(1.0) * 2.0)),
            cov_auto: s.coverage_auto.unwrap_or(d.cov_auto),
            show_axis_tracks: s.axis_tracks.unwrap_or(d.show_axis_tracks),
            seed_sample: s.seed_sample.filter(|&r| r > 0),
            depth_sample: s.depth_sample.filter(|&r| r > 0),
        }
    }

    fn depth_method(&self) -> DepthMethod {
        match self.depth {
            DepthChoice::Seed => DepthMethod::Seed,
            DepthChoice::Kmer => DepthMethod::Kmer { k: self.depth_k },
            DepthChoice::None => DepthMethod::None,
        }
    }
}

/// Colour of each repeat class when colouring by class.
const CLASS_COLOURS: [Color32; RepeatClass::COUNT] = [
    Color32::from_rgb(0, 0, 0),
    Color32::from_rgb(30, 90, 200),
    Color32::from_rgb(230, 140, 0),
    Color32::from_rgb(210, 30, 30),
];

/// Resize-retarget threshold: if the plot canvas grows or shrinks by
/// more than this fraction relative to the canvas size at last
/// successful compute, the GUI schedules a new auto-fit compute. Small
/// drags are ignored — we don't want a recompute storm during a slow
/// window-edge drag.
const RESIZE_RETARGET_THRESHOLD: f32 = 0.25;

/// Idle delay after a resize / wheel before a recompute fires. Same
/// constant the wheel zoom-settle code uses.
const RECOMPUTE_SETTLE_MS: u64 = 200;

impl Default for Settings {
    fn default() -> Self {
        Self {
            mode: BlastMode::Blastn,
            matrix_name: "DNA+5/-4".into(),
            window_size: None,
            zoom: 1,
            pixel_fac: 50, // slider seed only; ignored when auto_pixel_fac
            strand: Strand::Both,
            self_comparison: false,
            triangle: Triangle::Both,
            memory_limit_bytes: 512 * 1024 * 1024,
            reverse_query: false,
            reverse_subject: false,
            inverted_repeat_colour: false,
            align_window_size: 100,
            auto_fit_compute_zoom: true,
            auto_pixel_fac: true,
            show_ridge_overlay: false,
            ridge_min_length: 8,
            ridge_max_gap: 2,
            large: LargeSettings::default(),
        }
    }
}

impl Settings {
    fn build_matrix(&self) -> Option<ScoreMatrix> {
        if self.mode == BlastMode::Blastn {
            Some(ScoreMatrix::dna_identity())
        } else {
            ScoreMatrix::by_name(&self.matrix_name)
        }
    }
}

/// Greyramp configuration. Generates a 256-byte LUT applied to the
/// pixelmap on every redraw — the underlying [`DotPlot`] is untouched
/// (spec §4.2.1).
#[derive(Clone, Copy)]
struct Greyramp {
    /// Pixel values ≤ `white` map to white (255 in the displayed image).
    /// In C dotter terminology these are the lightest dots.
    white: u8,
    /// Pixel values ≥ `black` map to black (0). Strongest dots.
    black: u8,
    /// Swap inverts the LUT — white pixels become black and vice versa.
    swap: bool,
}

impl Default for Greyramp {
    fn default() -> Self {
        // C dotter defaults per spec §4.2.2: white=40, black=100. The
        // 60-value linear gradient between them is what gives dotter
        // its smooth-looking diagonals — marginal cells along a
        // diagonal fade through intermediate greys, which the eye
        // reads as anti-aliasing. A colocated white/black point
        // collapses to a hard threshold and produces a visibly
        // "rasterized" look.
        Self {
            white: 40,
            black: 100,
            swap: false,
        }
    }
}

impl Greyramp {
    fn lut(&self) -> [u8; 256] {
        let mut lut = [0_u8; 256];
        let (lo, hi) = if self.white <= self.black {
            (self.white as i32, self.black as i32)
        } else {
            // If user inverted the order, treat it as a swap.
            (self.black as i32, self.white as i32)
        };
        for (i, slot) in lut.iter_mut().enumerate() {
            let i = i as i32;
            // Below `lo`: displayed as white (or black if swapped).
            // Above `hi`: displayed as black (or white if swapped).
            // Between: linear ramp.
            let displayed = if i <= lo {
                if self.swap {
                    0
                } else {
                    255
                }
            } else if i >= hi {
                if self.swap {
                    255
                } else {
                    0
                }
            } else {
                let span = (hi - lo) as f32;
                let frac = (i - lo) as f32 / span;
                let v = if self.swap { frac } else { 1.0 - frac };
                (v * 255.0).round() as i32
            };
            *slot = displayed.clamp(0, 255) as u8;
        }
        lut
    }
}

pub struct DottirApp {
    query: Option<Sequence>,
    subject: Option<Sequence>,
    settings: Settings,
    plot: Option<DotPlot>,
    /// Last error from a compute / load operation. Cleared on next user action.
    last_error: Option<String>,
    greyramp: Greyramp,
    texture: Option<TextureHandle>,
    /// Whether the cached texture matches the current greyramp + pixelmap.
    texture_dirty: bool,
    /// Scroll position: top-left of the canvas in pixelmap pixel
    /// coords. The pixelmap is always drawn 1:1 (no GPU resampling);
    /// when it's bigger than the canvas, this offset chooses which
    /// portion is visible. Clamped to `[0, max(0, pixelmap_dim −
    /// canvas_dim)]` per axis.
    view_offset: Vec2,
    /// Crosshair in *absolute residue* coords `(q_seq, s_seq)`,
    /// 0-indexed into the full sequences (NOT pixelmap pixels and
    /// NOT slice-local). Keyboard nudges move by 1 residue regardless
    /// of compute zoom or slice — spec §4.2.4. Render-time conversion
    /// to slice-local pixmap pixels is `q_pix = (q_seq - slice.start) / zoom`;
    /// the crosshair hides itself when that pixel falls outside the
    /// current slice's pixmap.
    crosshair: Option<(u32, u32)>,
    show_settings: bool,
    /// True = light theme (the default — matches the greyscale plot
    /// area). False = egui's dark theme.
    light_theme: bool,
    /// Background worker that owns `compute_dotplot`. UI thread
    /// dispatches requests and polls results each frame (C2).
    worker: crate::compute_worker::ComputeWorker,
    /// Monotonic id assigned to each compute request. Results whose
    /// id doesn't match `last_dispatched_id` are discarded as stale.
    last_dispatched_id: u64,
    /// True while a worker has an in-flight request. Drives the
    /// status-bar progress indicator.
    compute_in_flight: bool,
    /// Suppress recompute during constructor — load_sequence() runs
    /// during DottirApp::new and would otherwise dispatch multiple
    /// jobs before the second sequence is loaded.
    suspend_recompute: bool,
    /// H2: whether the alignment-view dock is shown beneath the
    /// canvas. Default true; toggled via View → "Show alignment
    /// view".
    align_dock_visible: bool,
    /// Plot rectangle dimensions (logical) measured on the most recent
    /// canvas paint. `None` before the first paint. The
    /// display-matched compute zoom (Dotter invariant) needs this to
    /// choose `PlotConfig::zoom`, so the first compute is deferred
    /// until this is set.
    measured_plot_area: Option<(f32, f32)>,
    /// `true` between a sequence-load and the first canvas paint —
    /// i.e., we have sequences but no canvas measurement yet, so we
    /// can't pick a display-matched compute zoom. `update()` fires
    /// the deferred compute once `measured_plot_area` becomes `Some`.
    pending_initial_compute: bool,
    /// Timestamp of the most recent canvas-resize event significant
    /// enough to retarget auto-fit. Updated on each paint whose
    /// `plot_area` diverges by more than [`RESIZE_RETARGET_THRESHOLD`]
    /// from `last_compute_canvas_size`. Cleared once the retarget
    /// recompute fires (or once auto-fit is disabled).
    pending_resize_retarget: Option<std::time::Instant>,
    /// Canvas size (logical) used by the *most recently completed*
    /// compute. Used as the reference for resize-retarget detection
    /// — comparing against the live `measured_plot_area` would
    /// re-trigger every frame during a slow drag.
    last_compute_canvas_size: Option<(f32, f32)>,
    /// History stack for the rectangle-zoom workflow. Each entry is a
    /// view we can [`pop_history`] back to. Bounded to
    /// [`HISTORY_MAX`] entries — older ones drop off the bottom.
    history: std::collections::VecDeque<ViewSnapshot>,
    /// LRU-style cache of recently computed pixelmaps, keyed by the
    /// `(slice, compute_zoom)` pair. Hitting Back to a cached view is
    /// instant; otherwise the cache miss triggers a fresh compute.
    /// Cleared on any settings change that would invalidate the
    /// pixelmap (matrix, window size, strand, etc.).
    pixelmap_cache: std::collections::VecDeque<(ViewSlice, u32, DotPlot, Option<LargePlot>)>,
    /// View offset to apply when the next compute result lands.
    /// Set by `action_rect_zoom` / `action_back` so the new pixelmap
    /// (which has different dims, hence different coord space) gets
    /// the correct scroll position. Consumed by `poll_compute_results`.
    pending_view_offset_after_compute: Option<Vec2>,
    /// A middle-mouse drag completed during the most recent paint —
    /// we capture it before binding `&self.plot` and apply
    /// `action_rect_zoom` after the paint closure ends (where
    /// `&self.plot` is no longer alive). Carries the plot rectangle
    /// captured at the same moment, so the math uses the canvas
    /// dimensions at release time.
    pending_rect_zoom: Option<(Rect, RectSelect)>,
    /// Pending rectangle-selection state while the user is
    /// middle-mouse-dragging. `start` is the press position in screen
    /// coords; `current` is the latest mouse position. On release we
    /// compute a new zoom from the rectangle and push a history entry.
    rect_select: Option<RectSelect>,
    /// Right-drag selection: in progress (screen), then done (residues)
    /// with its action menu open.
    selection_drag: Option<RectSelect>,
    selection: Option<Selection>,
    /// Residue range covered by the most recently computed pixelmap.
    /// `None` before any compute has completed. Coordinate transforms
    /// in the paint loop add `slice.start` when converting
    /// pixmap_pixel → residue.
    current_slice: Option<ViewSlice>,
    /// Slice to use for the *next* dispatched compute. Set by
    /// `action_rect_zoom` / `action_fit` / `action_back` /
    /// `load_sequence`; consumed in `recompute`. `None` means "use the
    /// full sequence".
    target_slice: Option<ViewSlice>,
    /// Karlin-derived auto pixel_fac locked to the first compute of
    /// the current settings. Reused on subsequent computes so
    /// darkness stays consistent across slices (different sub-ranges
    /// have slightly different residue composition, which would
    /// otherwise re-derive a slightly different pixel_fac). Cleared
    /// by `invalidate_caches` when settings change in a way that
    /// affects Karlin (matrix, mode, sequences).
    locked_pixel_fac: Option<u32>,
    /// Slice that the latest dispatched compute is computing. Read
    /// by `poll_compute_results` when the matching `id` lands to
    /// update `current_slice`. Avoids plumbing slice through the
    /// worker (which only sees pre-sliced bytes).
    last_dispatched_slice: Option<ViewSlice>,
    /// Ridges extracted from the *current* plot (see
    /// `refresh_ridges`). Re-extracted whenever the plot or any
    /// parameter affecting extraction (overlay toggle, threshold,
    /// min length, gap) changes. Painted as anti-aliased line
    /// segments above the raster in the canvas draw.
    current_ridges: Vec<Ridge>,
    /// Cached `ctx.pixels_per_point()`, refreshed at the start of
    /// each `update()`. Used by `canonical_compute_zoom` (which is
    /// `&self`, so it can't access ctx) and by every texel↔screen
    /// conversion in the canvas, so the pixmap can be sized to and
    /// drawn over physical-pixel boundaries (Option C in
    /// `docs/reviews/dotter-sizing-model.md`). Default 1.0 so any
    /// path that runs before the first `update()` (load_sequence in
    /// the constructor) behaves as on a non-HiDPI display.
    pixels_per_point: f32,
    /// Annotation overlay (ADR 0005): per-axis loaded GFF3/BED bound to
    /// the query (X) / subject (Y) sequence. Drawn as full-span bands.
    query_annot: Option<crate::annotation_overlay::AxisAnnot>,
    subject_annot: Option<crate::annotation_overlay::AxisAnnot>,
    /// Master toggle for drawing annotation bands.
    annot_show: bool,
    /// Global alpha applied to every annotation band, `[0.0, 1.0]`. The
    /// single "darkness" control; overlapping bands self-darken via
    /// compositing.
    annot_alpha: f32,
    /// Large-engine inputs shared with the worker; rebuilt lazily after a
    /// sequence load or a change to the depth floors.
    large_source: Option<Arc<crate::compute_worker::LargeSource>>,
    /// Bumped whenever `large_source` is rebuilt, so the worker re-indexes.
    large_source_key: u64,
    /// Full large-engine result behind the current `plot`, if the large
    /// engine produced it.
    large_plot: Option<LargePlot>,
    /// Matching bp of the visible dots per pixel column / row and class;
    /// rebuilt with the texture.
    axis_profile: Option<AxisProfile>,
    /// Window size resolved by the large engine; reused by exact-engine
    /// detail views so both engines draw the same features.
    locked_window: Option<u32>,
    /// Genome comparison (ADR 0008), once loaded.
    genome: Option<GenomeView>,
    /// Show the genome view rather than the dotplot.
    genome_active: bool,
    /// Indexing / alignment in progress.
    genome_load: Option<GenomeLoad>,
    /// The "Open genome comparison" dialog, with its last error.
    genome_dialog: Option<(GenomeInputs, Option<String>)>,
    /// View settings to apply once a session's genome comparison loads.
    genome_pending: Option<crate::session::SessionGenome>,
    /// Keep showing the dotplot when that load finishes.
    genome_reopen_dotplot: bool,
    /// Record sorting in progress.
    record_sort_job: Option<crate::record_sort::SortJob>,
    /// Sorts still to run (a session being restored).
    record_sort_queue: std::collections::VecDeque<SortAxis>,
    /// Sorts applied to the current sequences, in order.
    record_sorts: Vec<SortAxis>,
    /// The sequences before sorting, for *File order*.
    unsorted_query: Option<Sequence>,
    unsorted_subject: Option<Sequence>,
    /// PAF alignments drawn over the plot (spec §4.2.3).
    hsps: Option<dottir_io::hsp::HspSet>,
    hsp_mode: HspMode,
}

/// Maximum number of history entries kept around for Back. Older
/// entries fall off the bottom. Also caps the pixelmap cache size.
const HISTORY_MAX: usize = 5;

/// Residue range covered by the currently displayed pixelmap. The
/// canonical (fit) view spans the full sequences; a rect-zoom narrows
/// the slice + a configurable margin (see `RECT_ZOOM_MARGIN`) so the
/// user can pan around within the cached pixelmap without
/// recomputing. Cache key includes the slice, so different sub-views
/// don't collide.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct ViewSlice {
    q_range: std::ops::Range<usize>,
    s_range: std::ops::Range<usize>,
}

impl ViewSlice {
    fn full(qlen: usize, slen: usize) -> Self {
        Self {
            q_range: 0..qlen,
            s_range: 0..slen,
        }
    }
    /// Self-comparison is the kernel's "same q axis as s axis" mode.
    /// Off-diagonal sub-rectangles aren't self-comparisons even when
    /// they come from the same underlying sequence.
    fn is_self_comparison(&self) -> bool {
        self.q_range == self.s_range
    }
}

/// Fraction of the *selection* size to add as panning margin on each
/// side of the slice during rect-zoom. With 0.5 the resulting pixelmap
/// covers the selection plus 50 % of its size in each direction
/// (clamped to the full sequence length), so the user can drag a
/// short distance without triggering a recompute.
const RECT_ZOOM_MARGIN: f64 = 0.5;

/// One slot on the history stack: enough state to restore a view.
/// `view_offset` is in pixelmap-pixel coords that only make sense at
/// this snapshot's `slice` and `compute_zoom`. `crosshair` is in
/// absolute residue coords and is slice/zoom-independent — it survives
/// a snapshot restore without translation.
#[derive(Clone, Debug)]
struct ViewSnapshot {
    slice: ViewSlice,
    compute_zoom: u32,
    view_offset: Vec2,
    crosshair: Option<(u32, u32)>,
}

/// A finished right-drag selection: absolute residue ranges, the canvas
/// it was made on, and where its menu opens.
#[derive(Clone, Debug)]
struct Selection {
    q: std::ops::Range<usize>,
    s: std::ops::Range<usize>,
    plot_area: Rect,
    menu_at: Pos2,
}

/// Middle-mouse rectangle-select drag state.
#[derive(Clone, Copy, Debug)]
struct RectSelect {
    start_screen: Pos2,
    current_screen: Pos2,
}

impl DottirApp {
    pub fn new(cc: &eframe::CreationContext<'_>, startup: StartupConfig) -> Self {
        // Default to a light theme — the plotting area is a greyscale
        // pixelmap on a near-white background, so a dark surround
        // muddles axis labels and panel text. Users who prefer dark
        // can toggle it in the View menu.
        cc.egui_ctx.set_visuals(egui::Visuals::light());

        // Worker thread that owns compute_dotplot. The repaint closure
        // wakes the UI as soon as a result arrives so we don't have to
        // spin on try_recv every frame.
        let ctx_for_repaint = cc.egui_ctx.clone();
        let worker = crate::compute_worker::ComputeWorker::spawn(move || {
            ctx_for_repaint.request_repaint();
        });

        // `startup.pixel_fac == 0` is the "auto" sentinel; map it onto
        // the auto_pixel_fac checkbox and seed the slider at 50 so it
        // has a sensible value if the user later toggles auto off.
        let startup_auto_pf = startup.pixel_fac == 0;
        let startup_slider_pf = if startup_auto_pf {
            50
        } else {
            startup.pixel_fac
        };
        let settings = Settings {
            large: LargeSettings::default(),
            mode: startup.mode,
            matrix_name: startup.matrix_name.clone(),
            window_size: startup.window_size,
            zoom: startup.zoom.max(1),
            pixel_fac: startup_slider_pf,
            strand: startup.strand,
            self_comparison: startup.self_comparison,
            triangle: Triangle::Both,
            memory_limit_bytes: startup.memory_limit_bytes,
            reverse_query: false,
            reverse_subject: false,
            inverted_repeat_colour: false,
            align_window_size: 100,
            auto_fit_compute_zoom: true,
            auto_pixel_fac: startup_auto_pf,
            show_ridge_overlay: false,
            ridge_min_length: 8,
            ridge_max_gap: 2,
        };

        let mut app = Self {
            query: None,
            subject: None,
            settings,
            plot: None,
            last_error: None,
            greyramp: Greyramp::default(),
            texture: None,
            texture_dirty: true,
            view_offset: Vec2::ZERO,
            crosshair: None,
            show_settings: false,
            light_theme: true,
            worker,
            last_dispatched_id: 0,
            compute_in_flight: false,
            // Suspend recompute during the two pre-loads — we only
            // want one job dispatched once both inputs are settled.
            suspend_recompute: true,
            align_dock_visible: true,
            measured_plot_area: None,
            pending_initial_compute: false,
            pending_resize_retarget: None,
            last_compute_canvas_size: None,
            history: std::collections::VecDeque::with_capacity(HISTORY_MAX),
            pixelmap_cache: std::collections::VecDeque::with_capacity(HISTORY_MAX),
            rect_select: None,
            selection_drag: None,
            selection: None,
            pending_view_offset_after_compute: None,
            pending_rect_zoom: None,
            current_slice: None,
            target_slice: None,
            locked_pixel_fac: None,
            last_dispatched_slice: None,
            current_ridges: Vec::new(),
            pixels_per_point: 1.0,
            query_annot: None,
            subject_annot: None,
            annot_show: true,
            annot_alpha: 0.35,
            large_source: None,
            large_source_key: 0,
            large_plot: None,
            axis_profile: None,
            locked_window: None,
            genome: None,
            genome_active: false,
            genome_load: None,
            genome_dialog: None,
            genome_pending: None,
            genome_reopen_dotplot: false,
            record_sort_job: None,
            record_sort_queue: std::collections::VecDeque::new(),
            record_sorts: Vec::new(),
            unsorted_query: None,
            unsorted_subject: None,
            hsps: None,
            hsp_mode: HspMode::Line,
        };

        // Pre-load any sequences supplied on the command line. Errors
        // are stashed in `last_error` and surfaced in the status bar
        // instead of aborting the GUI — the user can still recover via
        // File → Open.
        for (role, path, regions) in [
            (SeqRole::Query, startup.query, startup.query_regions),
            (SeqRole::Subject, startup.subject, startup.subject_regions),
        ] {
            match (path, regions.is_empty()) {
                (Some(p), true) => app.load_sequence(role, p),
                (Some(p), false) => app.load_sequence_regions(role, p, &regions),
                (None, _) => {}
            }
        }
        // Bind any startup annotation files now that the sequences are
        // loaded (set_annotation needs the target sequence present).
        if let Some(p) = startup.annot_query {
            app.load_annotation(SeqRole::Query, p);
        }
        if let Some(p) = startup.annot_subject {
            app.load_annotation(SeqRole::Subject, p);
        }
        if let Some(m) = startup.hsp_mode {
            app.hsp_mode = m;
        }
        if let Some(p) = startup.hsp_paf {
            app.load_hsps(p);
        }
        // Release the suspend. The first compute can't fire yet
        // because the canvas hasn't been measured — `pending_initial_
        // compute` is `true` and `update()` will fire it after the
        // first paint records `measured_plot_area`.
        app.suspend_recompute = false;
        if app.query.is_some() && app.subject.is_some() {
            app.pending_initial_compute = true;
        }
        if let Some(inputs) = startup.genome {
            app.start_genome_load(&cc.egui_ctx, inputs);
        }
        app.record_sort_queue.extend(startup.sort_records);
        app
    }

    fn start_genome_load(&mut self, ctx: &Context, inputs: GenomeInputs) {
        let ctx = ctx.clone();
        self.genome_dialog = None;
        self.genome_load = Some(GenomeLoad::start(inputs, move || ctx.request_repaint()));
    }

    fn poll_genome_load(&mut self) {
        let Some(result) = self.genome_load.as_ref().and_then(GenomeLoad::poll) else {
            return;
        };
        let load = self.genome_load.take().expect("polled above");
        match result {
            Ok(data) => {
                tracing::info!(
                    "genome comparison: {} alignments ({})",
                    data.alns.len(),
                    data.source
                );
                let mut view = GenomeView::new(data, load.inputs);
                if let Some(p) = self.genome_pending.take() {
                    p.apply(&mut view);
                }
                self.genome = Some(view);
                self.genome_active = !std::mem::take(&mut self.genome_reopen_dotplot);
            }
            Err(e) => {
                self.genome_pending = None;
                if e != "cancelled" {
                    self.genome_dialog = Some((load.inputs, Some(e)));
                }
            }
        }
    }

    /// Dotplot of genome regions: x regions on the query (horizontal)
    /// axis, y regions on the subject axis.
    fn drill_down(&mut self, x: Vec<Region>, y: Vec<Region>) {
        let Some(g) = &mut self.genome else {
            return;
        };
        let same = g.data.self_comparison && x == y;
        let seqs = g
            .data
            .x
            .index
            .region_sequence(&x)
            .and_then(|q| Ok((q, g.data.y.index.region_sequence(&y)?)));
        match seqs {
            Ok((q, s)) => {
                tracing::info!(
                    "genome drill-down: {} bp ({} regions) × {} bp ({} regions)",
                    q.len(),
                    q.records.len(),
                    s.len(),
                    s.records.len()
                );
                g.message = None;
                self.genome_active = false;
                self.query = Some(q);
                self.subject = Some(s);
                self.query_annot = None;
                self.subject_annot = None;
                self.settings.self_comparison = same;
                self.forget_record_sorts();
                self.history.clear();
                self.last_error = None;
                let detected = self.query.as_ref().map(Sequence::detect_alphabet);
                self.sequences_changed(detected.unwrap_or(DetectedAlphabet::Dna));
            }
            Err(e) => g.message = Some(format!("Reading regions failed: {e}")),
        }
    }

    /// Load a sequence input and bind it to an axis. The path may be a
    /// FASTA, or a GFF3 carrying its own sequences after `##FASTA` — in
    /// the latter case the features from the same file are bound as
    /// that axis's annotation overlay (an explicit
    /// `--gff-query`/`--gff-subject` or a later File-menu load
    /// replaces them).
    fn load_sequence(&mut self, role: SeqRole, path: PathBuf) {
        match dottir_io::load_sequence_input(&path) {
            Ok(input) => {
                let seq = input.sequence;
                let embedded = input.annotations.filter(|a| !a.features.is_empty());
                let detected = seq.detect_alphabet();
                tracing::info!(
                    "loaded {} ({:?}, {} residues, {} records, {} embedded features, detected {:?})",
                    path.display(),
                    input.kind,
                    seq.len(),
                    seq.records.len(),
                    embedded.as_ref().map_or(0, |a| a.features.len()),
                    detected,
                );
                match role {
                    SeqRole::Query => {
                        self.query = Some(seq);
                        // Buffer coords changed → any bound annotation is stale.
                        self.query_annot = None;
                    }
                    SeqRole::Subject => {
                        self.subject = Some(seq);
                        self.subject_annot = None;
                    }
                }
                self.forget_record_sorts();
                self.last_error = None;
                if let Some(set) = embedded {
                    self.set_annotation(role, set);
                }
                self.sequences_changed(detected);
            }
            Err(e) => {
                self.last_error = Some(format!("failed to load {}: {e}", path.display()));
            }
        }
    }

    /// Load regions of a FASTA (`--query-region` / `--subject-region`):
    /// records keep their source coordinates.
    fn load_sequence_regions(&mut self, role: SeqRole, path: PathBuf, regions: &[String]) {
        match dottir_io::faidx::load_regions(&path, regions, &dottir_io::cache::cache_dir()) {
            Ok(seq) => {
                tracing::info!(
                    "loaded {} regions of {} ({} residues)",
                    seq.records.len(),
                    path.display(),
                    seq.len()
                );
                let detected = seq.detect_alphabet();
                match role {
                    SeqRole::Query => {
                        self.query = Some(seq);
                        self.query_annot = None;
                    }
                    SeqRole::Subject => {
                        self.subject = Some(seq);
                        self.subject_annot = None;
                    }
                }
                self.forget_record_sorts();
                self.last_error = None;
                self.sequences_changed(detected);
            }
            Err(e) => {
                self.last_error = Some(format!("failed to load {}: {e}", path.display()));
            }
        }
    }

    fn forget_record_sorts(&mut self) {
        self.record_sorts.clear();
        self.unsorted_query = None;
        self.unsorted_subject = None;
        if let Some(job) = self.record_sort_job.take() {
            job.cancel();
        }
    }

    /// Why records cannot be sorted now, if they cannot.
    fn record_sort_blocker(&self) -> Option<&'static str> {
        if self.query.is_none() || self.subject.is_none() {
            Some("load two sequences first")
        } else if self.settings.self_comparison {
            Some("needs two different sequences")
        } else if self.settings.mode != BlastMode::Blastn {
            Some("DNA (BLASTN) only")
        } else if self.record_sort_job.is_some() {
            Some("sorting…")
        } else {
            None
        }
    }

    /// Order and orient one axis's records along the other (D-GENIES, from
    /// whole-plot large-engine chains), in the background.
    fn start_record_sort(&mut self, ctx: &Context, axis: SortAxis) {
        if let Some(why) = self.record_sort_blocker() {
            self.last_error = Some(format!("sort records: {why}"));
            return;
        }
        let Some(source) = self.large_source() else {
            return;
        };
        let ls = &self.settings.large;
        let seed = match SpacedSeed::parse(&ls.seed_pattern) {
            Ok(s) => s,
            Err(e) => {
                self.last_error = Some(format!("seed pattern: {e}"));
                return;
            }
        };
        let config = LargeConfig {
            matrix: ScoreMatrix::dna_identity(),
            window_size: self.settings.window_size,
            memory_limit_bytes: ls.memory_limit_bytes,
            seed,
            two_hit_dist: ls.two_hit_dist,
            max_occ: ls.max_occ,
            sample: ls.seed_sample,
            depth_sample: ls.depth_sample,
            depth: ls.depth_method(),
            classes: ls.thresholds,
            chain: Some(ChainConfig::default()),
            ..LargeConfig::default()
        };
        let ctx = ctx.clone();
        self.record_sort_job = Some(crate::record_sort::SortJob::start(
            source,
            config,
            axis,
            move || ctx.request_repaint(),
        ));
    }

    fn poll_record_sort(&mut self, ctx: &Context) {
        let Some(result) = self.record_sort_job.as_ref().and_then(|j| j.poll()) else {
            if self.record_sort_job.is_none() && self.record_sort_blocker().is_none() {
                if let Some(axis) = self.record_sort_queue.pop_front() {
                    self.start_record_sort(ctx, axis);
                }
            }
            return;
        };
        let job = self.record_sort_job.take().expect("polled above");
        match result {
            Ok(chains) => self.apply_record_sort(job.axis, &chains),
            Err(e) => {
                self.record_sort_queue.clear();
                if !e.to_lowercase().contains("cancel") {
                    self.last_error = Some(format!("sort records: {e}"));
                }
            }
        }
    }

    fn apply_record_sort(&mut self, axis: SortAxis, chains: &[dottir_core::large::Chain]) {
        let (Some(q), Some(s)) = (&self.query, &self.subject) else {
            return;
        };
        let o = crate::record_sort::order(chains, q, s, axis);
        let order: Vec<usize> = o.order.iter().map(|&r| r as usize).collect();
        let (target, unsorted, role) = match axis {
            SortAxis::Vertical => (
                &mut self.subject,
                &mut self.unsorted_subject,
                SeqRole::Subject,
            ),
            SortAxis::Horizontal => (&mut self.query, &mut self.unsorted_query, SeqRole::Query),
        };
        let Some(seq) = target.take() else {
            return;
        };
        let sorted = seq.reordered(&order, &o.flips);
        unsorted.get_or_insert(seq);
        *target = Some(sorted);
        tracing::info!(
            "sorted {} records from {} chains: {} reversed, {} without chains",
            order.len(),
            chains.len(),
            o.flips.iter().filter(|&&f| f).count(),
            o.placements.iter().filter(|p| p.is_none()).count()
        );
        self.record_sorts.push(axis);
        self.rebind_annotation(role);
        self.history.clear();
        self.sequences_changed(DetectedAlphabet::Dna);
    }

    /// Back to the records' file order.
    fn restore_record_order(&mut self) {
        let mut changed = false;
        if let Some(q) = self.unsorted_query.take() {
            self.query = Some(q);
            self.rebind_annotation(SeqRole::Query);
            changed = true;
        }
        if let Some(s) = self.unsorted_subject.take() {
            self.subject = Some(s);
            self.rebind_annotation(SeqRole::Subject);
            changed = true;
        }
        self.record_sorts.clear();
        if changed {
            self.history.clear();
            self.sequences_changed(DetectedAlphabet::Dna);
        }
    }

    /// Re-bind an axis's annotation after its sequence changed order.
    fn rebind_annotation(&mut self, role: SeqRole) {
        let (annot, seq) = match role {
            SeqRole::Query => (&mut self.query_annot, &self.query),
            SeqRole::Subject => (&mut self.subject_annot, &self.subject),
        };
        if let (Some(a), Some(seq)) = (annot.take(), seq) {
            *annot = Some(a.rebind(seq));
        }
    }

    fn draw_record_sort_window(&mut self, ctx: &Context) {
        let Some(job) = &self.record_sort_job else {
            return;
        };
        let (frac, elapsed) = job.progress();
        egui::Window::new("Sorting records")
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, Vec2::ZERO)
            .show(ctx, |ui| {
                ui.label("Synteny chains of the whole plot (large engine)…");
                ui.add(egui::ProgressBar::new(frac).show_percentage());
                ui.label(format!("{:.0} s", elapsed.as_secs_f32()));
                if ui.button("Cancel").clicked() {
                    job.cancel();
                }
            });
        ctx.request_repaint_after(std::time::Duration::from_millis(500));
    }

    fn record_sort_menu(&mut self, ui: &mut egui::Ui) {
        let blocker = self.record_sort_blocker();
        let ctx = ui.ctx().clone();
        for (axis, label) in [
            (SortAxis::Vertical, "Subject (vertical) along query"),
            (SortAxis::Horizontal, "Query (horizontal) along subject"),
        ] {
            let r = ui
                .add_enabled(blocker.is_none(), egui::Button::new(label))
                .on_hover_text(
                    "Order and reverse-complement records to follow the other \
                     sequence (the D-GENIES contig sort), from the large \
                     engine's synteny chains of the whole plot",
                );
            let r = match blocker {
                Some(why) => r.on_disabled_hover_text(why),
                None => r,
            };
            if r.clicked() {
                ui.close_menu();
                self.start_record_sort(&ctx, axis);
            }
        }
        let sorted = self.unsorted_query.is_some() || self.unsorted_subject.is_some();
        if ui
            .add_enabled(
                sorted && self.record_sort_job.is_none(),
                egui::Button::new("File order"),
            )
            .clicked()
        {
            ui.close_menu();
            self.restore_record_order();
        }
    }

    /// Fresh sequence data → reset view state and caches, then compute.
    fn sequences_changed(&mut self, detected: DetectedAlphabet) {
        // HSPs are placed on the sequences as shown; place them again.
        if let Some(path) = self.hsps.as_ref().map(|h| h.path.clone()) {
            self.load_hsps(path);
        }
        self.large_source = None;
        self.large_plot = None;
        self.view_offset = Vec2::ZERO;
        self.crosshair = None;
        self.current_slice = None;
        self.target_slice = None;
        self.invalidate_caches();
        self.maybe_switch_mode_from_alphabet(detected);
        if self.maybe_apply_auto_zoom() {
            self.pending_initial_compute = false;
            self.recompute();
        } else {
            // We don't yet know the canvas size (no paint has
            // happened) — defer the first compute until the
            // first frame measures `plot_area`. See
            // `update()`'s post-paint check.
            self.pending_initial_compute = true;
        }
    }

    /// Read a PAF and place its alignments on the plot as HSPs.
    fn load_hsps(&mut self, path: PathBuf) {
        let (Some(q), Some(s)) = (&self.query, &self.subject) else {
            self.last_error = Some("load the sequences before HSPs".into());
            return;
        };
        match dottir_io::hsp::load(&path, q, s) {
            Ok(set) => {
                tracing::info!(
                    "HSPs: {} alignments, {} placed, {} for other records",
                    set.alignments,
                    set.hsps.len(),
                    set.skipped
                );
                if set.hsps.is_empty() {
                    self.last_error = Some(format!(
                        "{}: no alignments between the plotted records",
                        path.display()
                    ));
                }
                self.hsps = Some(set);
            }
            Err(e) => {
                self.hsps = None;
                self.last_error = Some(format!("failed to load {}: {e}", path.display()));
            }
        }
        self.texture_dirty = true;
    }

    /// Load an annotation file from disk and bind it to an axis. Errors
    /// surface in the status bar (like `load_sequence`).
    fn load_annotation(&mut self, role: SeqRole, path: PathBuf) {
        match dottir_io::AnnotSet::load(&path) {
            Ok(set) => self.set_annotation(role, set),
            Err(e) => {
                self.last_error = Some(format!("failed to load {}: {e}", path.display()));
            }
        }
    }

    /// Bind a parsed annotation set to an axis and store it for the
    /// overlay. In self-comparison mode the same set is bound to both
    /// axes (each against its own — identical — sequence). The caller
    /// (File menu / CLI, Phase 5) is responsible for having loaded the
    /// target sequence first; if it isn't present the set is dropped.
    fn set_annotation(&mut self, role: SeqRole, set: dottir_io::AnnotSet) {
        use crate::annotation_overlay::AxisAnnot;
        let bind = |seq: &Option<Sequence>, set: dottir_io::AnnotSet| -> Option<AxisAnnot> {
            seq.as_ref().map(|s| AxisAnnot::new(set, s))
        };
        if self.settings.self_comparison {
            self.query_annot = bind(&self.query, set.clone());
            self.subject_annot = bind(&self.subject, set);
        } else {
            match role {
                SeqRole::Query => self.query_annot = bind(&self.query, set),
                SeqRole::Subject => self.subject_annot = bind(&self.subject, set),
            }
        }
        for (axis, ann) in [
            ("query", &self.query_annot),
            ("subject", &self.subject_annot),
        ] {
            if let Some(ax) = ann {
                tracing::info!(
                    "annotation bound to {axis}: {} features, {} skipped",
                    ax.bound.len(),
                    ax.skipped,
                );
            }
        }
        self.last_error = None;
    }

    /// Pick a *display-matched* `PlotConfig::zoom` for the current
    /// sequence pair following Dotter's invariant — the produced
    /// pixelmap is roughly the size of the *physical* canvas, so it
    /// renders close to 1:1 with no GPU resampling (no LINEAR blur,
    /// no NEAREST moiré). Per-axis target:
    /// ```text
    /// target_zoom = max(
    ///     ceil(qlen / physical_plot_width),
    ///     ceil(slen / physical_plot_height),
    ///     1,
    /// )
    /// ```
    /// Returns `true` if a target zoom could be computed (sequences
    /// loaded + canvas measured). `false` means the caller should
    /// defer the recompute until the next paint (see
    /// `pending_initial_compute`).
    /// Pure helper: per-axis display-matched compute zoom for the
    /// current sequence pair and canvas measurement. Returns `None`
    /// when sequences aren't loaded or the canvas hasn't been
    /// measured yet. Independent of `auto_fit_compute_zoom` (the
    /// setting only controls whether the result is *applied*) — also
    /// used by zoom-settle as the floor for tier-out so wheel-cycling
    /// returns to the exact same compute zoom rather than drifting
    /// through integer-truncated halves and doubles.
    fn canonical_compute_zoom(&self) -> Option<u32> {
        let q = self.query.as_ref()?;
        let s = self.subject.as_ref()?;
        let (plot_w, plot_h) = self.measured_plot_area?;
        // Option C: target the *physical* canvas dimension. The
        // pixmap then ends up with one texel per physical pixel, and
        // the draw rect uses `pw / ppp, ph / ppp` logical points
        // (see `draw_canvas`). Earlier history: this used to
        // multiply by `pixels_per_point` and overflow the canvas
        // because the draw side still used logical points; the fix
        // is to update both sides of the conversion together.
        let ppp = self.pixels_per_point as f64;
        let phys_w = (plot_w as f64).max(1.0) * ppp;
        let phys_h = (plot_h as f64).max(1.0) * ppp;
        let zoom_q = (q.len() as f64 / (phys_w * self.qx() as f64)).ceil() as u32;
        let zoom_s = (s.len() as f64 / phys_h).ceil() as u32;
        let base = zoom_q.max(zoom_s).max(1);
        Some(snap_zoom_to_period_divisor(
            base,
            &self.zoom_period_hints(),
            2.0,
        ))
    }

    fn zoom_period_hints(&self) -> Vec<usize> {
        let Some(q) = self.query.as_ref() else {
            return Vec::new();
        };
        let Some(s) = self.subject.as_ref() else {
            return Vec::new();
        };
        if self.settings.self_comparison {
            return q.record_period_hint().into_iter().collect();
        }
        match (q.record_period_hint(), s.record_period_hint()) {
            (Some(qp), Some(sp)) => vec![qp, sp],
            _ => Vec::new(),
        }
    }

    fn maybe_apply_auto_zoom(&mut self) -> bool {
        if !self.settings.auto_fit_compute_zoom {
            return true; // pass-through: respect manual settings.zoom
        }
        let Some(new_zoom) = self.canonical_compute_zoom() else {
            return false;
        };
        if new_zoom != self.settings.zoom {
            let qlen = self.query.as_ref().map(|q| q.len()).unwrap_or(0);
            let slen = self.subject.as_ref().map(|s| s.len()).unwrap_or(0);
            let (pw, ph) = self.measured_plot_area.unwrap_or((0.0, 0.0));
            tracing::info!(
                "auto-fit: zoom = {new_zoom} \
                 (qlen={qlen}, slen={slen}, plot={pw:.0}x{ph:.0} logical, periods={:?})",
                self.zoom_period_hints(),
            );
            self.settings.zoom = new_zoom;
        }
        true
    }

    /// Push the current view onto the history stack so a later `Back`
    /// can return to it. Bounded — oldest entry drops off the bottom.
    /// Re-extract ridges from the current pixmap into
    /// `current_ridges`. Called whenever the plot changes (new
    /// compute lands, cache hit, settings change) or a parameter
    /// affecting extraction is touched (greyramp white, min length,
    /// max gap, overlay toggle). Cheap — a single linear scan of
    /// the pixmap.
    fn refresh_ridges(&mut self) {
        if !self.settings.show_ridge_overlay {
            self.current_ridges.clear();
            return;
        }
        let Some(plot) = self.plot.as_ref() else {
            self.current_ridges.clear();
            return;
        };
        // Threshold tracks the greyramp's white point so the
        // overlay's "lit" definition follows the user's noise
        // floor: dropping greyramp.white hides more raster dots
        // AND prunes weak ridges in tandem.
        let params = RidgeParams {
            threshold: self.greyramp.white,
            min_length: self.settings.ridge_min_length.max(1),
            max_gap: self.settings.ridge_max_gap,
        };
        self.current_ridges =
            extract_ridges_from_pixels(&plot.pixels, plot.width, plot.height, &params);
    }

    fn push_history(&mut self) {
        let Some(plot) = &self.plot else {
            return;
        };
        let Some(slice) = self.current_slice.clone() else {
            return;
        };
        if self.history.len() >= HISTORY_MAX {
            self.history.pop_front();
        }
        self.history.push_back(ViewSnapshot {
            slice,
            compute_zoom: plot.params.zoom,
            view_offset: self.view_offset,
            crosshair: self.crosshair,
        });
    }

    /// Look up a cached pixelmap for the given `(slice, compute_zoom)`
    /// (LRU "touch on read"). Cleared whenever settings that affect
    /// the pixelmap change.
    fn cache_get(
        &mut self,
        slice: &ViewSlice,
        compute_zoom: u32,
    ) -> Option<(DotPlot, Option<LargePlot>)> {
        let pos = self
            .pixelmap_cache
            .iter()
            .position(|(s, z, _, _)| *z == compute_zoom && s == slice)?;
        let (s, z, plot, large) = self.pixelmap_cache.remove(pos)?;
        let restored = (plot.clone(), large.clone());
        self.pixelmap_cache.push_back((s, z, plot, large));
        Some(restored)
    }

    /// Insert a freshly computed pixelmap into the cache; evict
    /// oldest if at capacity.
    fn cache_insert(
        &mut self,
        slice: ViewSlice,
        compute_zoom: u32,
        plot: DotPlot,
        large: Option<LargePlot>,
    ) {
        self.pixelmap_cache
            .retain(|(s, z, _, _)| *z != compute_zoom || s != &slice);
        if self.pixelmap_cache.len() >= HISTORY_MAX {
            self.pixelmap_cache.pop_front();
        }
        self.pixelmap_cache
            .push_back((slice, compute_zoom, plot, large));
    }

    /// Drop everything pixel-dependent — call when settings that
    /// affect compute output change (matrix, window, strand, etc.).
    /// Also clears the Karlin-locked pixel_fac so the next compute
    /// re-derives it under the new settings.
    fn invalidate_caches(&mut self) {
        self.pixelmap_cache.clear();
        self.history.clear();
        self.locked_pixel_fac = None;
        self.locked_window = None;
    }

    /// Whether `slice` is computed by the large engine.
    fn use_large_engine(&self, slice: &ViewSlice) -> bool {
        if self.settings.mode != BlastMode::Blastn
            || self.settings.reverse_query
            || self.settings.reverse_subject
        {
            return false;
        }
        match self.settings.large.engine {
            EngineChoice::Exact => false,
            EngineChoice::Large => true,
            EngineChoice::Auto => {
                slice.q_range.len() as f64 * slice.s_range.len() as f64
                    > self.settings.large.threshold_cells
            }
        }
    }

    /// Inputs for the large engine, built on first use after a load.
    fn large_source(&mut self) -> Option<Arc<crate::compute_worker::LargeSource>> {
        if let Some(src) = &self.large_source {
            return Some(src.clone());
        }
        let q = self.query.as_ref()?;
        let s = self.subject.as_ref()?;
        let is_self = self.settings.self_comparison && q.bytes() == s.bytes();
        let ls = &self.settings.large;
        let floor = |seq: &Sequence| {
            ls.softmask_as_repeat
                .then(|| {
                    dottir_io::large::soft_mask_floor(seq, ls.thresholds.floor_of(RepeatClass::Mid))
                })
                .flatten()
        };
        self.large_source_key += 1;
        let src = Arc::new(crate::compute_worker::LargeSource {
            key: self.large_source_key,
            query: Arc::from(q.bytes()),
            q_breaks: q.breaks(),
            q_floor: floor(q),
            subject: (!is_self).then(|| Arc::from(s.bytes())),
            s_breaks: s.breaks(),
            s_floor: if is_self { None } else { floor(s) },
        });
        self.large_source = Some(src.clone());
        Some(src)
    }

    /// Re-derive the displayed pixelmap from the large result after a
    /// display-only change (class toggles, strand colouring).
    fn rederive_large_plot(&mut self) {
        if let Some(lp) = &self.large_plot {
            let mut plot = lp.to_dot_plot(
                self.settings.large.classes,
                self.separate_strands(),
                self.settings.large.shade,
            );
            if let Some(old) = &self.plot {
                plot.params = old.params;
            }
            self.plot = Some(plot);
            self.texture_dirty = true;
            self.refresh_ridges();
        }
    }

    /// The large plot behind the current pixelmap, if it is one.
    fn shown_large_plot(&self) -> Option<&LargePlot> {
        let (lp, plot) = (self.large_plot.as_ref()?, self.plot.as_ref()?);
        (lp.width == plot.width && lp.height == plot.height).then_some(lp)
    }

    fn coverage_shading(&self) -> bool {
        self.settings.large.shade == Shade::Coverage
            && self.shown_large_plot().is_some_and(|lp| lp.has_coverage())
    }

    /// The pixel → grey mapping on screen; image exports use the same.
    fn display_lut(&self) -> [u8; 256] {
        display_lut(
            self.coverage_shading(),
            &self.greyramp,
            self.settings.large.cov_white_bp,
            self.settings.large.cov_black_bp,
        )
    }

    /// Export metadata: the version and the sensitivity behind the greys.
    fn export_metadata(&self) -> Vec<(&'static str, String)> {
        let mut m = vec![("dottir-gui", env!("CARGO_PKG_VERSION").to_string())];
        if self.coverage_shading() {
            m.push(("shading", "coverage".into()));
            m.push((
                "coverage-white-bp",
                self.settings.large.cov_white_bp.to_string(),
            ));
            m.push((
                "coverage-black-bp",
                self.settings.large.cov_black_bp.to_string(),
            ));
        } else {
            m.push(("shading", "score".into()));
            m.push(("greyramp-white", self.greyramp.white.to_string()));
            m.push(("greyramp-black", self.greyramp.black.to_string()));
        }
        if let Some(h) = self.hsps.as_ref().filter(|_| self.hsp_mode != HspMode::Off) {
            m.push(("hsp-paf", h.path.display().to_string()));
            m.push(("hsp-mode", self.hsp_mode.name().to_string()));
        }
        m
    }

    /// Fit the coverage ramp to the shown pixels: white at the 2nd and
    /// black at the 95th percentile of the lit ones.
    fn auto_coverage_ramp(&mut self) {
        let Some(lp) = self.shown_large_plot() else {
            return;
        };
        if !lp.has_coverage() {
            return;
        }
        let mut hist = [0u64; 256];
        for c in lp.merged_coverage(None, self.settings.large.classes) {
            hist[c as usize] += 1;
        }
        let lit: u64 = hist[1..].iter().sum();
        if lit == 0 {
            return;
        }
        let pct = |p: f64| {
            let target = (lit as f64 * p).ceil() as u64;
            let mut acc = 0;
            for (code, &n) in hist.iter().enumerate().skip(1) {
                acc += n;
                if acc >= target {
                    return coverage_bases(code as u8) as f32;
                }
            }
            coverage_bases(255) as f32
        };
        // Whole bases: the sliders are integer and would otherwise report
        // the rounding as a user edit (which turns auto-fit off).
        let w = pct(0.02).max(1.0).round();
        let b = pct(0.95).max(w * 2.0).round();
        self.settings.large.cov_white_bp = w;
        self.settings.large.cov_black_bp = b;
        self.texture_dirty = true;
    }

    /// Margin reserved for the axis tracks.
    fn axis_track_width(&self) -> f32 {
        if self.settings.large.show_axis_tracks && self.large_plot.is_some() {
            AXIS_TRACK
        } else {
            0.0
        }
    }

    fn separate_strands(&self) -> bool {
        self.settings.inverted_repeat_colour
            && self.settings.mode == BlastMode::Blastn
            && matches!(self.settings.strand, Strand::Both)
    }

    /// Restore the canonical (display-matched) compute zoom over the
    /// full sequences — what you'd see right after a fresh load.
    /// Pushes the current view onto Back so the user can return.
    fn action_fit(&mut self) {
        let Some(canon) = self.canonical_compute_zoom() else {
            return;
        };
        let (Some(q), Some(s)) = (self.query.as_ref(), self.subject.as_ref()) else {
            return;
        };
        let full_slice = ViewSlice::full(q.len(), s.len());
        let already_canonical = self
            .current_slice
            .as_ref()
            .map(|cs| cs == &full_slice)
            .unwrap_or(false)
            && self.plot.as_ref().map(|p| p.params.zoom).unwrap_or(0) == canon
            && self.view_offset == Vec2::ZERO;
        if already_canonical {
            return;
        }
        self.push_history();
        self.settings.zoom = canon;
        self.target_slice = Some(full_slice);
        self.pending_view_offset_after_compute = Some(Vec2::ZERO);
        self.recompute();
    }

    /// Pop the history stack and recompute (or restore from cache)
    /// at the previous view's compute zoom + slice. No-op when
    /// history is empty.
    fn action_back(&mut self) {
        let Some(prev) = self.history.pop_back() else {
            return;
        };
        self.settings.zoom = prev.compute_zoom;
        self.crosshair = prev.crosshair;
        self.target_slice = Some(prev.slice);
        self.pending_view_offset_after_compute = Some(prev.view_offset);
        self.recompute();
    }

    /// Zoom-into-rectangle (dotter classic): the user middle-mouse-
    /// dragged a rectangle on the canvas. Compute a new compute_zoom
    /// such that the selected residue range fills the physical
    /// canvas, AND a new slice = selection ± [`RECT_ZOOM_MARGIN`] (so
    /// the user can pan a bit without recomputing). Push the current
    /// view onto Back and recompute. After the result lands
    /// `poll_compute_results` applies the pending `view_offset` so
    /// the user sees the selected region.
    fn action_rect_zoom(&mut self, plot_area: Rect, rs: RectSelect) {
        if let Some((x, y)) = self.selection_residues(plot_area, &rs) {
            self.zoom_to_residues(plot_area, x, y);
        }
    }

    /// The action menu of a finished right-drag selection.
    fn selection_menu(&mut self, ctx: &Context) {
        let Some(sel) = self.selection.clone() else {
            return;
        };
        let (Some(q), Some(s)) = (&self.query, &self.subject) else {
            self.selection = None;
            return;
        };
        let size = |r: &std::ops::Range<usize>| fmt_bp(r.len() as u64);
        let header = format!("{} × {}", size(&sel.q), size(&sel.s));
        let bedpe = selection_bedpe(q, s, &sel.q, &sel.s);
        let mut action = None;
        egui::Area::new(egui::Id::new("selection_menu"))
            .fixed_pos(sel.menu_at + Vec2::new(6.0, 6.0))
            .order(egui::Order::Foreground)
            .show(ctx, |ui| {
                egui::Frame::popup(ui.style()).show(ui, |ui| {
                    ui.label(egui::RichText::new(header).size(11.0).strong());
                    if ui.button("Zoom to selection").clicked() {
                        action = Some(0);
                    }
                    if ui
                        .button("Open in new window")
                        .on_hover_text(
                            "A new dottir window with just these regions, the same \
                             settings, HSPs and annotations",
                        )
                        .clicked()
                    {
                        action = Some(1);
                    }
                    if ui
                        .button("Copy coordinates (BEDPE)")
                        .on_hover_text(bedpe.trim_end())
                        .clicked()
                    {
                        ctx.copy_text(bedpe.clone());
                        action = Some(2);
                    }
                    if ui.button("Cancel").clicked() {
                        action = Some(2);
                    }
                });
            });
        match action {
            Some(0) => {
                self.selection = None;
                self.zoom_to_residues(sel.plot_area, sel.q, sel.s);
            }
            Some(1) => {
                self.selection = None;
                if let Err(e) = self.open_selection_window(&sel) {
                    self.last_error = Some(format!("open in new window: {e}"));
                }
            }
            Some(_) => self.selection = None,
            None => {}
        }
    }

    /// Query residues per pixelmap column at zoom 1: BLASTX plots one
    /// column per codon, so x positions are in bases but columns are 3.
    fn qx(&self) -> usize {
        if self.settings.mode == BlastMode::Blastx {
            3
        } else {
            1
        }
    }

    /// Axis span of an image export: the first residue shown on each axis
    /// and the residues covered (the zoomed-in slice, or everything).
    fn export_axes(&self) -> ((u32, u32), u32, u32) {
        let (qlen, slen) = (
            self.query.as_ref().map_or(0, |q| q.len()),
            self.subject.as_ref().map_or(0, |s| s.len()),
        );
        match &self.current_slice {
            Some(sl) => (
                (sl.q_range.start as u32, sl.s_range.start as u32),
                sl.q_range.end.min(qlen).saturating_sub(sl.q_range.start) as u32,
                sl.s_range.end.min(slen).saturating_sub(sl.s_range.start) as u32,
            ),
            None => ((0, 0), qlen as u32, slen as u32),
        }
    }

    /// Synteny chains as shown (toggle, min score, mirror), as SVG lines
    /// in the current slice's pixelmap pixels.
    fn shown_chain_lines(&self) -> Vec<dottir_io::svg_export::SvgLine> {
        let (true, Some(lp), Some(plot)) = (
            self.settings.large.show_chains,
            &self.large_plot,
            &self.plot,
        ) else {
            return Vec::new();
        };
        let min = self.settings.large.chain_min_score as f64;
        let chains: Vec<_> = lp
            .chains
            .iter()
            .filter(|c| c.score >= min)
            .cloned()
            .collect();
        let z = plot.params.zoom.max(1) as f64;
        let (q_off, s_off) = self
            .current_slice
            .as_ref()
            .map_or((0, 0), |sl| (sl.q_range.start, sl.s_range.start));
        let mirror = lp.self_comparison && lp.origin_q == lp.origin_s;
        dottir_io::large::chain_svg_lines(&chains, z as usize, mirror)
            .into_iter()
            .map(|mut l| {
                for p in [&mut l.from, &mut l.to] {
                    p.0 -= q_off as f64 / z;
                    p.1 -= s_off as f64 / z;
                }
                l
            })
            .collect()
    }

    /// Write the synteny chains (PAF, or BEDPE by extension).
    fn save_chains(&mut self) {
        let (Some(lp), Some(q), Some(s)) = (&self.large_plot, &self.query, &self.subject) else {
            return;
        };
        let Some(path) = rfd::FileDialog::new()
            .set_title("Save synteny chains")
            .set_file_name("chains.paf")
            .add_filter("PAF", &["paf"])
            .add_filter("BEDPE", &["bedpe"])
            .save_file()
        else {
            return;
        };
        let min = self.settings.large.chain_min_score as f64;
        let chains: Vec<_> = lp
            .chains
            .iter()
            .filter(|c| c.score >= min)
            .cloned()
            .collect();
        if let Err(e) = dottir_io::large::write_chains(&path, &chains, q, s) {
            self.last_error = Some(format!("writing {}: {e}", path.display()));
        }
    }

    /// Spawn another dottir-gui on the selected regions of the inputs.
    fn open_selection_window(&self, sel: &Selection) -> Result<(), String> {
        let (Some(q), Some(s)) = (&self.query, &self.subject) else {
            return Err("no sequences".into());
        };
        let set = &self.settings;
        let mut args: Vec<std::ffi::OsString> = Vec::new();
        let mut regions = Vec::new();
        for (seq, range, flipped, flag) in [
            (q, &sel.q, set.reverse_query, "--query-region"),
            (s, &sel.s, set.reverse_subject, "--subject-region"),
        ] {
            let path = seq
                .source_path
                .clone()
                .ok_or("the sequence was not read from a file")?;
            args.push(path.into());
            let specs = selection_specs(seq, range, flipped);
            for spec in &specs {
                regions.push((flag, spec.clone()));
            }
        }
        for (flag, spec) in regions {
            args.push(flag.into());
            args.push(spec.into());
        }
        let mode = match set.mode {
            BlastMode::Blastn => "blastn",
            BlastMode::Blastp => "blastp",
            BlastMode::Blastx => "blastx",
        };
        let strand = match set.strand {
            Strand::Forward => "forward",
            Strand::Reverse => "reverse",
            Strand::Both => "both",
        };
        for (k, v) in [
            ("--mode", mode),
            ("--matrix", &set.matrix_name),
            ("--strand", strand),
        ] {
            args.push(k.into());
            args.push(v.into());
        }
        if let Some(w) = set.window_size {
            args.push("-W".into());
            args.push(w.to_string().into());
        }
        if let Some(h) = &self.hsps {
            args.push("--hsp-paf".into());
            args.push(h.path.clone().into());
            args.push("--hsp-mode".into());
            args.push(self.hsp_mode.name().into());
        }
        for (annot, flag) in [
            (&self.query_annot, "--gff-query"),
            (&self.subject_annot, "--gff-subject"),
        ] {
            if let Some(p) = annot.as_ref().and_then(|a| a.set.source_path.clone()) {
                args.push(flag.into());
                args.push(p.into());
            }
        }
        let exe = std::env::current_exe().map_err(|e| e.to_string())?;
        tracing::info!("open in new window: {} {:?}", exe.display(), args);
        std::process::Command::new(exe)
            .args(&args)
            .spawn()
            .map_err(|e| e.to_string())?;
        Ok(())
    }

    /// Absolute residue ranges (query, subject) under a screen
    /// rectangle; `None` for a near-empty drag.
    fn selection_residues(
        &self,
        plot_area: Rect,
        rs: &RectSelect,
    ) -> Option<(std::ops::Range<usize>, std::ops::Range<usize>)> {
        let (Some(plot), Some(q), Some(s), Some(cur_slice)) = (
            self.plot.as_ref(),
            self.query.as_ref(),
            self.subject.as_ref(),
            self.current_slice.clone(),
        ) else {
            return None;
        };
        let cur_zoom = plot.params.zoom.max(1) as f32;
        let zoom_x = cur_zoom * self.qx() as f32;
        let ppp = self.pixels_per_point;

        // Screen → pixelmap-texel. Match the canvas's snap policy
        // (round to a whole-texel multiple in logical coords) so
        // the selection rect lands on the same texels the user
        // saw under the cursor. Logical-pt offset → texel via
        // `* ppp`.
        let snap_texel = |v: f32| (v * ppp).round() / ppp;
        let draw_off_x = snap_texel(self.view_offset.x);
        let draw_off_y = snap_texel(self.view_offset.y);
        let to_pix = |p: Pos2| -> (f32, f32) {
            let l = p - plot_area.left_top();
            ((l.x + draw_off_x) * ppp, (l.y + draw_off_y) * ppp)
        };
        let (sx, sy) = to_pix(rs.start_screen);
        let (ex, ey) = to_pix(rs.current_screen);
        let (px_lo_x, px_hi_x) = if sx <= ex { (sx, ex) } else { (ex, sx) };
        let (px_lo_y, px_hi_y) = if sy <= ey { (sy, ey) } else { (ey, sy) };
        // Drop tiny selections (likely accidental clicks).
        if (px_hi_x - px_lo_x) < 3.0 || (px_hi_y - px_lo_y) < 3.0 {
            return None;
        }
        let px_lo_x = px_lo_x.max(0.0).min(plot.width as f32);
        let px_lo_y = px_lo_y.max(0.0).min(plot.height as f32);
        let px_hi_x = px_hi_x.max(px_lo_x).min(plot.width as f32);
        let px_hi_y = px_hi_y.max(px_lo_y).min(plot.height as f32);

        // Pixelmap pixels → *absolute* residues (offset by the
        // current slice's start).
        let sel_lo_x = cur_slice.q_range.start + (px_lo_x * zoom_x).floor().max(0.0) as usize;
        let sel_lo_y = cur_slice.s_range.start + (px_lo_y * cur_zoom).floor().max(0.0) as usize;
        let sel_hi_x = (cur_slice.q_range.start + (px_hi_x * zoom_x).ceil() as usize).min(q.len());
        let sel_hi_y =
            (cur_slice.s_range.start + (px_hi_y * cur_zoom).ceil() as usize).min(s.len());
        (sel_lo_x < sel_hi_x && sel_lo_y < sel_hi_y)
            .then_some((sel_lo_x..sel_hi_x, sel_lo_y..sel_hi_y))
    }

    /// Zoom so residues `x` × `y` fill the canvas (with a panning margin).
    fn zoom_to_residues(
        &mut self,
        plot_area: Rect,
        x: std::ops::Range<usize>,
        y: std::ops::Range<usize>,
    ) {
        let (Some(plot), Some(q), Some(s)) = (
            self.plot.as_ref(),
            self.query.as_ref(),
            self.subject.as_ref(),
        ) else {
            return;
        };
        let cur_zoom_u = plot.params.zoom.max(1);
        let ppp = self.pixels_per_point;
        let (sel_lo_x, sel_hi_x, sel_lo_y, sel_hi_y) = (x.start, x.end, y.start, y.end);
        let sel_w = sel_hi_x.saturating_sub(sel_lo_x).max(1);
        let sel_h = sel_hi_y.saturating_sub(sel_lo_y).max(1);

        // New compute zoom so the *selection* (not the slice) fills
        // the *physical* canvas — matches `canonical_compute_zoom`'s
        // Option C policy (one pixmap texel = one physical pixel).
        let phys_w = (plot_area.width() as f64).max(1.0) * (ppp as f64);
        let phys_h = (plot_area.height() as f64).max(1.0) * (ppp as f64);
        let qx = self.qx();
        let zoom_x = ((sel_w as f64) / (phys_w * qx as f64)).ceil() as u32;
        let zoom_y = ((sel_h as f64) / phys_h).ceil() as u32;
        let new_zoom = zoom_x.max(zoom_y).max(1);

        // New slice = selection ± RECT_ZOOM_MARGIN of its size,
        // clamped to the full sequence. Gives a panning buffer
        // without unbounded memory.
        let margin_x = ((sel_w as f64) * RECT_ZOOM_MARGIN).ceil() as usize;
        let margin_y = ((sel_h as f64) * RECT_ZOOM_MARGIN).ceil() as usize;
        let slice_lo_x = sel_lo_x.saturating_sub(margin_x);
        let slice_lo_y = sel_lo_y.saturating_sub(margin_y);
        let slice_hi_x = sel_hi_x.saturating_add(margin_x).min(q.len());
        let slice_hi_y = sel_hi_y.saturating_add(margin_y).min(s.len());
        let new_slice = ViewSlice {
            q_range: slice_lo_x..slice_hi_x,
            s_range: slice_lo_y..slice_hi_y,
        };

        // The user's view should land on the selection inside the
        // new slice's pixelmap coords:
        //   pixel = (residue − slice.start) / new_zoom
        let new_off_x =
            ((sel_lo_x - new_slice.q_range.start) as f64 / (new_zoom as usize * qx) as f64) as f32;
        let new_off_y = ((sel_lo_y - new_slice.s_range.start) as f64 / new_zoom as f64) as f32;

        tracing::info!(
            "rect-zoom: zoom {} → {} (sel {}×{} residues, slice {}..{} × {}..{})",
            cur_zoom_u,
            new_zoom,
            sel_w,
            sel_h,
            new_slice.q_range.start,
            new_slice.q_range.end,
            new_slice.s_range.start,
            new_slice.s_range.end,
        );
        self.push_history();
        self.settings.zoom = new_zoom;
        self.target_slice = Some(new_slice);
        self.pending_view_offset_after_compute = Some(Vec2::new(new_off_x, new_off_y));
        // Crosshair is in absolute residue coords — survives slice
        // and zoom changes unchanged. If it falls outside the new
        // slice the render path hides the lines (still recoverable
        // by Back).
        self.recompute();
    }

    /// If the freshly loaded sequence's detected alphabet conflicts
    /// with the current BLAST mode, switch the mode (and the matrix
    /// to that mode's default) so the kernel produces a meaningful
    /// dotplot. Without this, opening a protein FASTA into the GUI's
    /// default BLASTN mode gives a dim, broken-looking diagonal
    /// because most residues encode to SENTINEL in the DNA alphabet.
    fn maybe_switch_mode_from_alphabet(&mut self, detected: DetectedAlphabet) {
        let want_protein = matches!(detected, DetectedAlphabet::Protein);
        let want_dna = matches!(detected, DetectedAlphabet::Dna);
        let is_protein_mode = matches!(self.settings.mode, BlastMode::Blastp | BlastMode::Blastx);
        if want_protein && !is_protein_mode {
            self.settings.mode = BlastMode::Blastp;
            self.settings.matrix_name = "BLOSUM62".into();
            self.settings.strand = Strand::Forward;
            self.last_error =
                Some("Detected protein sequence — switched to BLASTP / BLOSUM62.".into());
        } else if want_dna && is_protein_mode {
            self.settings.mode = BlastMode::Blastn;
            self.settings.matrix_name = "DNA+5/-4".into();
            self.last_error = Some("Detected DNA sequence — switched to BLASTN.".into());
        }
    }

    /// Dispatch a new compute request to the background worker.
    /// Returns immediately; the result arrives via
    /// [`Self::poll_compute_results`] on a future frame and is applied
    /// in [`Self::apply_compute_result`].
    fn recompute(&mut self) {
        if self.suspend_recompute {
            return;
        }
        if self.query.is_none() || self.subject.is_none() {
            self.plot = None;
            self.texture = None;
            self.compute_in_flight = false;
            return;
        }
        let target_zoom = self.settings.zoom.max(1);
        // Resolve the target slice: explicit (`target_slice` set by
        // action_*) wins; otherwise reuse the current slice if any;
        // otherwise default to the full sequences.
        let target_slice = self
            .target_slice
            .take()
            .or_else(|| self.current_slice.clone());
        let qlen = self.query.as_ref().map(|q| q.len()).unwrap_or(0);
        let slen = self.subject.as_ref().map(|s| s.len()).unwrap_or(0);
        let target_slice = target_slice.unwrap_or_else(|| ViewSlice::full(qlen, slen));

        // Cache hit: restore instantly, no worker dispatch.
        let use_large = self.use_large_engine(&target_slice);
        // The large engine aligns its origin to the zoom grid; snap the
        // slice so the GUI's slice-local coordinates match the pixmap.
        let target_slice = if use_large {
            let z = target_zoom as usize;
            ViewSlice {
                q_range: target_slice.q_range.start / z * z..target_slice.q_range.end,
                s_range: target_slice.s_range.start / z * z..target_slice.s_range.end,
            }
        } else {
            target_slice
        };
        if let Some((mut plot, large)) = self.cache_get(&target_slice, target_zoom) {
            if let Some(lp) = &large {
                plot = lp.to_dot_plot(
                    self.settings.large.classes,
                    self.separate_strands(),
                    self.settings.large.shade,
                );
            }
            self.large_plot = large;
            tracing::info!(
                "cache hit at zoom={target_zoom}, slice {}..{} × {}..{}",
                target_slice.q_range.start,
                target_slice.q_range.end,
                target_slice.s_range.start,
                target_slice.s_range.end,
            );
            self.current_slice = Some(target_slice);
            self.plot = Some(plot);
            self.texture_dirty = true;
            self.last_error = None;
            self.refresh_ridges();
            if self.settings.large.cov_auto {
                self.auto_coverage_ramp();
            }
            if let Some(off) = self.pending_view_offset_after_compute.take() {
                self.view_offset = off;
            }
            return;
        }
        if use_large {
            self.dispatch_large(target_slice, target_zoom);
            return;
        }
        // Re-borrow now that the cache check is done (avoids a
        // double-borrow of self).
        let (q, s) = (
            self.query.as_ref().expect("checked above"),
            self.subject.as_ref().expect("checked above"),
        );
        let matrix = match self.settings.build_matrix() {
            Some(m) => m,
            None => {
                self.last_error = Some(format!("unknown matrix '{}'", self.settings.matrix_name));
                return;
            }
        };
        // Pixel-fac resolution:
        // - Auto + first compute under current settings: pass 0 (core
        //   derives from Karlin); we'll lock the resolved value once
        //   the result lands.
        // - Auto + already locked: pass the locked value so all slices
        //   render with the same darkness scaling.
        // - Manual: pass the slider value (max(1)).
        let cfg_pixel_fac = if self.settings.auto_pixel_fac {
            self.locked_pixel_fac.unwrap_or(0)
        } else {
            self.settings.pixel_fac.max(1)
        };
        let mut cfg = PlotConfig {
            mode: self.settings.mode,
            matrix,
            window_size: self.settings.window_size.or(self.locked_window),
            zoom: self.settings.zoom.max(1),
            pixel_fac: cfg_pixel_fac,
            strand: self.settings.strand,
            // Self-comparison only applies on the diagonal — when the
            // slice ranges are equal (i.e., the user is looking at the
            // same residue range on both axes, typically the full view
            // or a square selection along the main diagonal).
            self_comparison: self.settings.self_comparison && target_slice.is_self_comparison(),
            triangle: self.settings.triangle,
            disable_mirror: false,
            memory_limit_bytes: self.settings.memory_limit_bytes,
            separate_strand_channels: self.settings.inverted_repeat_colour
                && self.settings.mode == BlastMode::Blastn
                && matches!(self.settings.strand, Strand::Both),
            reverse_query: self.settings.reverse_query,
            reverse_subject: self.settings.reverse_subject,
        };
        // BLASTP cannot use reverse strand.
        if cfg.mode == BlastMode::Blastp {
            cfg.strand = Strand::Forward;
        }

        self.last_dispatched_id = self.last_dispatched_id.wrapping_add(1);
        let id = self.last_dispatched_id;
        // Slice the sequence bytes for the compute. The core kernel
        // operates on `&[u8]` so slicing in the GUI keeps the core
        // API unchanged.
        let q_bytes = q.bytes();
        let s_bytes = s.bytes();
        let q_slice = q_bytes
            .get(target_slice.q_range.clone())
            .unwrap_or(&[])
            .to_vec();
        let s_slice = s_bytes
            .get(target_slice.s_range.clone())
            .unwrap_or(&[])
            .to_vec();
        let req = crate::compute_worker::ComputeRequest {
            id,
            // Sequences cloned per request — slice bytes only.
            query: q_slice,
            subject: s_slice,
            config: cfg,
            large: None,
        };
        tracing::info!(
            "dispatch compute id={id}, zoom={target_zoom}, slice {}..{} × {}..{}",
            target_slice.q_range.start,
            target_slice.q_range.end,
            target_slice.s_range.start,
            target_slice.s_range.end,
        );
        self.last_dispatched_slice = Some(target_slice);
        self.worker.dispatch(req);
        self.compute_in_flight = true;
    }

    fn dispatch_large(&mut self, slice: ViewSlice, zoom: u32) {
        let Some(source) = self.large_source() else {
            return;
        };
        let ls = self.settings.large.clone();
        let seed = match SpacedSeed::parse(&ls.seed_pattern) {
            Ok(s) => s,
            Err(e) => {
                self.last_error = Some(format!("seed pattern: {e}"));
                return;
            }
        };
        let qlen = source.query.len();
        let slen = source.subject.as_ref().map_or(qlen, |s| s.len());
        let full = slice.q_range == (0..qlen) && slice.s_range == (0..slen);
        let config = LargeConfig {
            matrix: ScoreMatrix::dna_identity(),
            window_size: self.settings.window_size,
            zoom,
            pixel_fac: if self.settings.auto_pixel_fac {
                0
            } else {
                self.settings.pixel_fac.max(1)
            },
            strand: self.settings.strand,
            self_comparison: source.subject.is_none(),
            triangle: self.settings.triangle,
            disable_mirror: false,
            memory_limit_bytes: ls.memory_limit_bytes,
            coverage: true,
            seed,
            two_hit_dist: ls.two_hit_dist,
            max_occ: ls.max_occ,
            sample: ls.seed_sample,
            depth_sample: ls.depth_sample,
            extend: None,
            depth: ls.depth_method(),
            classes: ls.thresholds,
            chain: Some(ChainConfig::default()),
            viewport: (!full).then(|| Viewport {
                q: slice.q_range.clone(),
                s: slice.s_range.clone(),
            }),
        };
        self.last_dispatched_id = self.last_dispatched_id.wrapping_add(1);
        let id = self.last_dispatched_id;
        tracing::info!(
            "dispatch large compute id={id}, zoom={zoom}, slice {}..{} × {}..{}",
            slice.q_range.start,
            slice.q_range.end,
            slice.s_range.start,
            slice.s_range.end,
        );
        self.last_dispatched_slice = Some(slice);
        self.worker.dispatch(crate::compute_worker::ComputeRequest {
            id,
            query: vec![],
            subject: vec![],
            config: PlotConfig::default_blastn(ScoreMatrix::dna_identity()),
            large: Some(crate::compute_worker::LargeJob {
                source,
                config,
                classes: ls.classes,
                separate_strands: self.separate_strands(),
                shade: ls.shade,
            }),
        });
        self.compute_in_flight = true;
    }

    /// Resize-retarget settle: when the canvas grew/shrunk by more
    /// than [`RESIZE_RETARGET_THRESHOLD`] relative to the canvas size
    /// at the last successful compute, recompute at a new
    /// display-matched zoom after [`RECOMPUTE_SETTLE_MS`] of resize
    /// inactivity. Only fires when `auto_fit_compute_zoom` is on.
    fn maybe_resize_retarget(&mut self, ctx: &Context) {
        let Some(t) = self.pending_resize_retarget else {
            return;
        };
        if !self.settings.auto_fit_compute_zoom {
            self.pending_resize_retarget = None;
            return;
        }
        if self.compute_in_flight {
            ctx.request_repaint_after(std::time::Duration::from_millis(RECOMPUTE_SETTLE_MS));
            return;
        }
        if t.elapsed().as_millis() < RECOMPUTE_SETTLE_MS as u128 {
            ctx.request_repaint_after(std::time::Duration::from_millis(RECOMPUTE_SETTLE_MS));
            return;
        }
        if self.maybe_apply_auto_zoom() {
            self.pending_resize_retarget = None;
            self.recompute();
        } else {
            // Couldn't compute target (sequences not loaded / no
            // measurement). Drop the pending state; the deferred-
            // initial-compute path will pick it up if/when sequences
            // arrive.
            self.pending_resize_retarget = None;
        }
    }

    /// Drain any completed worker results and apply the latest one
    /// (discarding stale results whose id is older than
    /// `last_dispatched_id`). Called once per frame.
    fn poll_compute_results(&mut self) {
        let mut latest: Option<crate::compute_worker::ComputeResult> = None;
        for r in self.worker.drain_results() {
            if r.id == self.last_dispatched_id {
                latest = Some(r);
            }
            // else: stale, drop silently.
        }
        if let Some(r) = latest {
            self.compute_in_flight = false;
            match r.plot {
                Ok(plot) => {
                    tracing::info!(
                        "computed {}×{} pixelmap (W={}, zoom={}, resolved_pixel_fac={})",
                        plot.width,
                        plot.height,
                        plot.params.window_size,
                        r.config_zoom,
                        plot.params.pixel_fac,
                    );
                    // Apply the slice that was paired with this dispatch.
                    if let Some(slice) = self.last_dispatched_slice.take() {
                        // Insert into LRU cache for instant restore on
                        // Back / Fit cycling.
                        self.cache_insert(
                            slice.clone(),
                            plot.params.zoom,
                            plot.clone(),
                            r.large.clone(),
                        );
                        self.current_slice = Some(slice);
                    }
                    if r.large.is_some() && self.settings.window_size.is_none() {
                        self.locked_window = Some(plot.params.window_size);
                    }
                    self.large_plot = r.large;
                    // Lock the resolved auto pixel_fac on the first
                    // compute after a settings change, so subsequent
                    // sliced computes use the same darkness scaling.
                    if self.settings.auto_pixel_fac && self.locked_pixel_fac.is_none() {
                        self.locked_pixel_fac = Some(plot.params.pixel_fac);
                    }
                    self.plot = Some(plot);
                    self.texture_dirty = true;
                    self.last_error = None;
                    self.refresh_ridges();
                    if self.settings.large.cov_auto {
                        self.auto_coverage_ramp();
                    }
                    // Pin the canvas size that the new compute was
                    // targeting (latest available measurement). The
                    // resize-retarget detector compares future paints
                    // against this — not against the previous frame's
                    // measurement — so a slow drag accumulates but
                    // only fires once it crosses the threshold.
                    self.last_compute_canvas_size = self.measured_plot_area;
                    // Apply any view_offset queued by the action that
                    // triggered this compute (rect-zoom / Back / Fit).
                    // The offset is in NEW pixelmap coords.
                    if let Some(off) = self.pending_view_offset_after_compute.take() {
                        self.view_offset = off;
                    } else {
                        // No queued offset → centre on origin (e.g.
                        // a fresh load or settings change).
                        self.view_offset = Vec2::ZERO;
                    }
                }
                Err(e) => {
                    self.last_error = Some(format!("compute_dotplot failed: {e}"));
                    // Failed compute leaves the previous pixelmap (if
                    // any) intact and drops the pending state.
                    self.pending_view_offset_after_compute = None;
                    self.last_dispatched_slice = None;
                }
            }
        }
    }

    /// The plot as shown: RGBA pixels (greyramp or coverage ramp, repeat-
    /// class colours, inverted repeats in magenta) and their size. Image
    /// exports write exactly this.
    ///
    /// HSPs in greyscale mode replace the dots; `hsp_lines` also draws
    /// line modes into the pixels (exports; the canvas draws vectors).
    fn display_rgba(&self, hsp_lines: bool) -> Option<(Vec<u8>, usize, usize)> {
        let plot = self.plot.as_ref()?;
        let coverage = self.coverage_shading();
        let lut = self.display_lut();
        let (pw, ph) = (plot.width as usize, plot.height as usize);
        let fwd_view = plot.pixels.as_slice();
        let rev_view = plot.reverse_pixels.as_deref();

        // Spec §4.4.3 inverted-repeat highlighting: when the plot has
        // a separate reverse channel, paint forward in grey and
        // reverse in magenta — overlapping cells take whichever
        // channel is stronger after the greyramp.
        let mut rgba = Vec::with_capacity(pw * ph * 4);
        let by_class = self.settings.large.colour_by_class
            && self
                .large_plot
                .as_ref()
                .is_some_and(|lp| lp.width == plot.width && lp.height == plot.height);
        if by_class {
            let lp = self.large_plot.as_ref().expect("checked above");
            let shown = self.settings.large.classes.classes();
            let channels: Vec<(Color32, [&[u8]; 2])> = shown
                .iter()
                .map(|&c| {
                    let ch = |r: bool| {
                        let i = r as usize * RepeatClass::COUNT + c.index();
                        if coverage {
                            lp.coverage[i].as_slice()
                        } else {
                            lp.channels[i].as_slice()
                        }
                    };
                    (CLASS_COLOURS[c.index()], [ch(false), ch(true)])
                })
                .collect();
            for i in 0..pw * ph {
                // Strongest visible class wins the pixel.
                let mut best: Option<(u8, Color32)> = None;
                for (colour, strands) in &channels {
                    for ch in strands.iter().filter(|c| !c.is_empty()) {
                        let ink = 255 - lut[ch[i] as usize];
                        if ink > 0 && best.is_none_or(|(b, _)| ink > b) {
                            best = Some((ink, *colour));
                        }
                    }
                }
                let px = match best {
                    None => [255, 255, 255, 255],
                    Some((ink, c)) => {
                        let mix = |v: u8| (255 - (ink as u16 * (255 - v as u16) / 255)) as u8;
                        [mix(c.r()), mix(c.g()), mix(c.b()), 255]
                    }
                };
                rgba.extend_from_slice(&px);
            }
        }
        match rev_view {
            _ if by_class => {}
            Some(rev) if self.settings.inverted_repeat_colour => {
                for i in 0..fwd_view.len() {
                    let f = lut[fwd_view[i] as usize];
                    let r = lut[rev[i] as usize];
                    let f_ink = 255 - f;
                    let r_ink = 255 - r;
                    let (cr, cg, cb) = if f_ink >= r_ink {
                        (f, f, f)
                    } else {
                        let ink = r_ink as u16;
                        let bg = 255_u16 - ink;
                        (
                            (bg + ink * 220 / 255) as u8,
                            (bg) as u8,
                            (bg + ink * 220 / 255) as u8,
                        )
                    };
                    rgba.extend_from_slice(&[cr, cg, cb, 255]);
                }
            }
            _ => {
                for &v in plot.combined().as_ref() {
                    let g = lut[v as usize];
                    rgba.extend_from_slice(&[g, g, g, 255]);
                }
            }
        }
        if let Some(h) = &self.hsps {
            if hsp_lines || self.hsp_mode == HspMode::Greyscale {
                let (q_off, s_off) = self
                    .current_slice
                    .as_ref()
                    .map_or((0, 0), |sl| (sl.q_range.start, sl.s_range.start));
                dottir_io::hsp::draw_rgba(
                    &mut rgba,
                    pw,
                    ph,
                    &h.hsps,
                    self.hsp_mode,
                    q_off,
                    s_off,
                    plot.params.zoom as usize,
                );
            }
        }
        Some((rgba, pw, ph))
    }

    fn ensure_texture(&mut self, ctx: &Context) {
        if !self.texture_dirty && self.texture.is_some() {
            return;
        }
        let Some((rgba, pw, ph)) = self.display_rgba(false) else {
            return;
        };
        let image = ColorImage::from_rgba_unmultiplied([pw, ph], &rgba);
        // The dotter-faithful render is always 1:1 (pixelmap pixel =
        // screen pixel). NEAREST → no GPU resampling, no sub-pixel
        // sampling phase, no moiré. The compute step is the only
        // place "zoom" happens.
        let handle = ctx.load_texture("dottir-pixelmap", image, TextureOptions::NEAREST);
        self.texture = Some(handle);
        self.texture_dirty = false;
        let (coverage, lut) = (self.coverage_shading(), self.display_lut());
        self.axis_profile = self
            .shown_large_plot()
            .and_then(|lp| AxisProfile::of(lp, self.settings.large.classes, coverage, &lut));
    }
}

/// Column (query) and row (subject) sums of the matching bp of every
/// visible dot, split by repeat class: the marginals of what the plot
/// shows. Divided by the zoom this is "aligned copies per base".
struct AxisProfile {
    q: Vec<[f32; RepeatClass::COUNT]>,
    s: Vec<[f32; RepeatClass::COUNT]>,
}

impl AxisProfile {
    /// `coverage_shading` / `lut` decide which dots are visible, exactly
    /// as the texture does. `None` when the plot has no coverage.
    fn of(
        lp: &LargePlot,
        classes: ClassMask,
        coverage_shading: bool,
        lut: &[u8; 256],
    ) -> Option<Self> {
        if !lp.has_coverage() {
            return None;
        }
        let (w, h) = (lp.width as usize, lp.height as usize);
        let decode: Vec<f32> = (0..=255u8).map(|c| coverage_bases(c) as f32).collect();
        let mut q = vec![[0f32; RepeatClass::COUNT]; w];
        let mut s = vec![[0f32; RepeatClass::COUNT]; h];
        for c in classes.classes() {
            for r in [false, true] {
                let i = r as usize * RepeatClass::COUNT + c.index();
                let cov = &lp.coverage[i];
                let shade = if coverage_shading {
                    cov
                } else {
                    &lp.channels[i]
                };
                if cov.is_empty() || shade.is_empty() {
                    continue;
                }
                for (p, (&cv, &sv)) in cov.iter().zip(shade).enumerate() {
                    if cv == 0 || lut[sv as usize] == 255 {
                        continue;
                    }
                    let bp = decode[cv as usize];
                    q[p % w][c.index()] += bp;
                    s[p / w][c.index()] += bp;
                }
            }
        }
        Some(Self { q, s })
    }
}

#[derive(Clone, Copy)]
enum SeqRole {
    Query,
    Subject,
}

impl eframe::App for DottirApp {
    fn update(&mut self, ctx: &Context, _frame: &mut eframe::Frame) {
        // Refresh cached display scale first — every coord transform
        // in the canvas closure and `canonical_compute_zoom` reads
        // this rather than calling `ctx.pixels_per_point()` directly.
        self.pixels_per_point = ctx.pixels_per_point().max(0.01);

        // Drain any completed background-compute results before any
        // panel reads `self.plot`. Stale results are discarded inside.
        self.poll_compute_results();
        self.poll_genome_load();
        self.poll_record_sort(ctx);
        self.draw_genome_windows(ctx);
        self.draw_record_sort_window(ctx);
        if self.genome_active && self.genome.is_some() {
            self.update_genome(ctx);
            return;
        }

        // Esc → Back (history pop). Backspace also; both are dotter-
        // classic for "step out".
        if ctx.input(|i| i.key_pressed(egui::Key::Escape) || i.key_pressed(egui::Key::Backspace))
            && self.selection.take().is_none()
        {
            self.action_back();
        }
        self.handle_keyboard(ctx);
        self.draw_menu(ctx);
        self.draw_greyramp_panel(ctx);
        if self.show_settings {
            self.draw_settings_window(ctx);
        }
        self.draw_status_bar(ctx);
        if self.align_dock_visible {
            self.draw_alignment_dock(ctx);
        }
        self.draw_canvas(ctx);

        // Apply a rectangle-zoom captured by the canvas paint. We do
        // this *after* the paint closure because action_rect_zoom
        // mutates self.plot via recompute, which would conflict with
        // the &self.plot borrow held during paint.
        if let Some((pa, rs)) = self.pending_rect_zoom.take() {
            self.action_rect_zoom(pa, rs);
        }
        self.selection_menu(ctx);

        // Deferred initial compute: pre-loaded sequences (CLI or
        // load_sequence before the first paint) couldn't pick a
        // display-matched zoom because no canvas size was known yet.
        // After the first paint records `measured_plot_area`, fire it.
        if self.pending_initial_compute && self.maybe_apply_auto_zoom() {
            self.pending_initial_compute = false;
            self.recompute();
        }
        // Resize-retarget: a debounced recompute fired when the canvas
        // grew/shrunk enough to invalidate the previous compute zoom.
        self.maybe_resize_retarget(ctx);
    }
}

impl DottirApp {
    fn handle_keyboard(&mut self, ctx: &Context) {
        let mods = ctx.input(|i| i.modifiers);
        // Step is in *residues* — independent of compute zoom or slice
        // (spec §4.2.4). Plain arrow = 1 residue, Shift = 10, Ctrl = 100.
        let step = if mods.ctrl {
            100_i32
        } else if mods.shift {
            10
        } else {
            1
        };
        if self.plot.is_none() {
            return;
        }
        // Clamp bounds come from the loaded sequence lengths; falling
        // back to a huge bound keeps things sensible if a sequence is
        // somehow missing (the next click will reset).
        let qlen = self
            .query
            .as_ref()
            .map(|s| s.len() as i64)
            .unwrap_or(i64::MAX);
        let slen = self
            .subject
            .as_ref()
            .map(|s| s.len() as i64)
            .unwrap_or(i64::MAX);
        let mut nudged = false;
        let mut snap = false;
        let (mut q, mut s) = self
            .crosshair
            .unwrap_or(((qlen / 2) as u32, (slen / 2) as u32));
        let mut q_i = q as i64;
        let mut s_i = s as i64;
        ctx.input(|i| {
            for ev in &i.events {
                if let egui::Event::Key {
                    key, pressed: true, ..
                } = ev
                {
                    match key {
                        // Single-axis (original GUI behaviour).
                        egui::Key::ArrowLeft => {
                            q_i -= step as i64;
                            nudged = true;
                        }
                        egui::Key::ArrowRight => {
                            q_i += step as i64;
                            nudged = true;
                        }
                        egui::Key::ArrowUp => {
                            s_i -= step as i64;
                            nudged = true;
                        }
                        egui::Key::ArrowDown => {
                            s_i += step as i64;
                            nudged = true;
                        }
                        // Main diagonal: both coords step together.
                        // Matches the original Dotter `,` / `.` keys.
                        egui::Key::Comma => {
                            q_i -= step as i64;
                            s_i -= step as i64;
                            nudged = true;
                        }
                        egui::Key::Period => {
                            q_i += step as i64;
                            s_i += step as i64;
                            nudged = true;
                        }
                        // Anti-diagonal: q and s step in opposite
                        // directions. Matches original Dotter `[`/`]`.
                        egui::Key::OpenBracket => {
                            q_i -= step as i64;
                            s_i += step as i64;
                            nudged = true;
                        }
                        egui::Key::CloseBracket => {
                            q_i += step as i64;
                            s_i -= step as i64;
                            nudged = true;
                        }
                        // Snap to the brightest pixel within a
                        // search radius — quick jump to whatever
                        // diagonal the crosshair is near.
                        egui::Key::Space => {
                            snap = true;
                        }
                        _ => {}
                    }
                }
            }
        });
        if nudged {
            q = q_i.clamp(0, qlen.saturating_sub(1).max(0)) as u32;
            s = s_i.clamp(0, slen.saturating_sub(1).max(0)) as u32;
            self.crosshair = Some((q, s));
        }
        if snap {
            self.snap_crosshair_to_line();
        }
    }

    /// Snap the crosshair to the brightest pixel within a search
    /// disc, then refine the residue position **within** that pixel
    /// by rescoring the underlying sliding window at every
    /// (q_residue, s_residue) midpoint in the pixel's `zoom × zoom`
    /// block — eliminating the off-by-up-to-(zoom-1) error that
    /// dropping back to "centre of pixel" would introduce at higher
    /// zoom tiers. Bound to **Space**.
    ///
    /// The coarse step finds the brightest pixel (current position is
    /// the tie-breaker — closest wins on equal value); the fine step
    /// recovers the exact base. For zoom = 1 the fine step is a no-op
    /// (the block has a single candidate residue).
    fn snap_crosshair_to_line(&mut self) {
        let Some(plot) = self.plot.as_ref() else {
            return;
        };
        let Some((cq_seq, cs_seq)) = self.crosshair else {
            return;
        };
        let z = plot.params.zoom.max(1) as i64;
        let zx = z * self.qx() as i64;
        let (q_off, s_off) = self
            .current_slice
            .as_ref()
            .map(|sl| (sl.q_range.start as i64, sl.s_range.start as i64))
            .unwrap_or((0, 0));
        // Convert current residue crosshair to slice-local pixel.
        // Crosshair outside the slice clamps to the nearest edge; the
        // search radius then explores inward.
        let pw = plot.width as i64;
        let ph = plot.height as i64;
        let cq_pix = ((cq_seq as i64 - q_off).max(0) / zx).clamp(0, pw - 1);
        let cs_pix = ((cs_seq as i64 - s_off).max(0) / z).clamp(0, ph - 1);
        const RADIUS: i64 = 64;
        let stride = plot.width as usize;
        let q_lo = (cq_pix - RADIUS).max(0);
        let q_hi = (cq_pix + RADIUS).min(pw - 1);
        let s_lo = (cs_pix - RADIUS).max(0);
        let s_hi = (cs_pix + RADIUS).min(ph - 1);
        // Track the best (value, -distance², q_pix, s_pix) tuple.
        // Tuple ordering gives us "max value, then min distance" for free.
        let mut best: Option<(u8, i64, i64, i64)> = None;
        for sp in s_lo..=s_hi {
            let row = sp as usize * stride;
            for qp in q_lo..=q_hi {
                let v = plot.pixels[row + qp as usize];
                if v == 0 {
                    continue;
                }
                let dq = qp - cq_pix;
                let ds = sp - cs_pix;
                let dist_sq = dq * dq + ds * ds;
                let candidate = (v, -dist_sq, qp, sp);
                match best {
                    None => best = Some(candidate),
                    Some(cur) if (candidate.0, candidate.1) > (cur.0, cur.1) => {
                        best = Some(candidate)
                    }
                    _ => {}
                }
            }
        }
        let Some((_, _, qp_best, sp_best)) = best else {
            return;
        };
        // Default: residue at the centre of the pixel's block, with
        // the slice origin added back in. The within-block refinement
        // below overrides this when we can actually score the
        // underlying window.
        let mut q_residue = q_off + qp_best * zx + zx / 2;
        let mut s_residue = s_off + sp_best * z + z / 2;
        // Fine step: walk every residue midpoint in the block,
        // compute the actual sliding-window score, pick the max.
        // Only do this when both sequences are loaded and we have a
        // score matrix — otherwise stick with the block centre.
        // (Not for BLASTX, whose windows are over translated codons.)
        if let (Some(qseq), Some(sseq), Some(matrix), false) = (
            self.query.as_ref(),
            self.subject.as_ref(),
            self.settings.build_matrix(),
            self.settings.mode == BlastMode::Blastx,
        ) {
            let w = plot.params.window_size as i64;
            let half = w / 2;
            let qbytes = qseq.bytes();
            let sbytes = sseq.bytes();
            let qlen = qbytes.len() as i64;
            let slen = sbytes.len() as i64;
            let mode = self.settings.mode;
            let strand = self.settings.strand;
            let q_mid_lo = q_off + qp_best * z;
            let q_mid_hi = (q_off + (qp_best + 1) * z).min(qlen);
            let s_mid_lo = s_off + sp_best * z;
            let s_mid_hi = (s_off + (sp_best + 1) * z).min(slen);
            let mut best_score: i64 = i64::MIN;
            for qm in q_mid_lo..q_mid_hi {
                for sm in s_mid_lo..s_mid_hi {
                    let sc = window_score(
                        qbytes, sbytes, qm, sm, half, w, qlen, slen, &matrix, mode, strand,
                    );
                    if sc > best_score {
                        best_score = sc;
                        q_residue = qm;
                        s_residue = sm;
                    }
                }
            }
        }
        let qlen_cap = self
            .query
            .as_ref()
            .map(|s| s.len() as i64)
            .unwrap_or(q_off + pw * zx);
        let slen_cap = self
            .subject
            .as_ref()
            .map(|s| s.len() as i64)
            .unwrap_or(s_off + ph * z);
        let q_residue = q_residue.clamp(0, qlen_cap.saturating_sub(1).max(0));
        let s_residue = s_residue.clamp(0, slen_cap.saturating_sub(1).max(0));
        self.crosshair = Some((q_residue as u32, s_residue as u32));
    }

    /// A frame of the genome view: the dotplot's menu, status bar and
    /// right-hand panel layout around the genome plot.
    fn update_genome(&mut self, ctx: &Context) {
        if ctx.input(|i| i.key_pressed(egui::Key::Escape) || i.key_pressed(egui::Key::Backspace)) {
            if let Some(g) = &mut self.genome {
                g.back();
            }
        }
        self.draw_menu(ctx);
        let Some(g) = &mut self.genome else {
            return;
        };
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| match &self.last_error {
                Some(err) => {
                    ui.colored_label(Color32::from_rgb(220, 90, 70), err);
                }
                None => {
                    ui.label(format!("genome mode · {}", g.status));
                }
            });
        });
        let mut action = None;
        egui::SidePanel::right("side_panel")
            .resizable(false)
            .default_width(220.0)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    action = g.panel_ui(ui);
                });
            });
        egui::CentralPanel::default().show(ctx, |ui| {
            if let Some(a) = g.canvas_ui(ui) {
                action = Some(a);
            }
        });
        match action {
            Some(GenomeAction::DrillDown { x, y }) => self.drill_down(x, y),
            Some(GenomeAction::Realign(inputs)) => {
                // Keep filters, sort and flips across the re-run.
                if let Some(g) = &self.genome {
                    let mut keep = crate::session::SessionGenome::from_view(g, true);
                    keep.inputs = inputs.clone();
                    self.genome_pending = Some(keep);
                }
                self.start_genome_load(ctx, inputs);
            }
            None => {}
        }
    }

    fn open_genome_dialog(&mut self) {
        let inputs = self
            .genome
            .as_ref()
            .map(|g| g.inputs.clone())
            .unwrap_or_default();
        self.genome_dialog = Some((inputs, None));
    }

    /// The open dialog and the progress window, in either view.
    fn draw_genome_windows(&mut self, ctx: &Context) {
        if let Some(load) = &self.genome_load {
            let (stage, frac, elapsed) = load.status();
            egui::Window::new("Genome comparison")
                .collapsible(false)
                .resizable(false)
                .anchor(egui::Align2::CENTER_CENTER, Vec2::ZERO)
                .show(ctx, |ui| {
                    ui.label(&stage);
                    match frac {
                        Some(f) => ui.add(egui::ProgressBar::new(f).show_percentage()),
                        None => ui.add(egui::ProgressBar::new(0.0).animate(true)),
                    };
                    ui.label(format!("{:.0} s", elapsed.as_secs_f32()));
                    if ui.button("Cancel").clicked() {
                        load.cancel();
                    }
                });
            // The elapsed time keeps moving.
            ctx.request_repaint_after(std::time::Duration::from_millis(500));
            return;
        }
        let Some((inputs, error)) = &mut self.genome_dialog else {
            return;
        };
        let mut open = true;
        let mut run = false;
        egui::Window::new("Open genome comparison")
            .collapsible(false)
            .open(&mut open)
            .show(ctx, |ui| {
                let pick =
                    |ui: &mut egui::Ui, label: &str, path: &mut Option<PathBuf>, fasta: bool| {
                        ui.label(label);
                        ui.horizontal(|ui| {
                            if ui.button("Choose…").clicked() {
                                let d = if fasta {
                                    rfd::FileDialog::new()
                                        .add_filter("FASTA", &["fa", "fasta", "fna", "fas", "gz"])
                                } else {
                                    rfd::FileDialog::new().add_filter("PAF", &["paf", "gz"])
                                };
                                if let Some(p) = d.add_filter("All files", &["*"]).pick_file() {
                                    *path = Some(p);
                                }
                            }
                            if path.is_some() && ui.button("Clear").clicked() {
                                *path = None;
                            }
                            let text = path
                                .as_ref()
                                .and_then(|p| p.file_name())
                                .map_or("—".to_string(), |f| f.to_string_lossy().into_owned());
                            ui.add(egui::Label::new(text).truncate());
                        });
                    };
                let mut reference =
                    (!inputs.reference.as_os_str().is_empty()).then(|| inputs.reference.clone());
                pick(ui, "Reference (horizontal axis)", &mut reference, true);
                inputs.reference = reference.unwrap_or_default();
                pick(
                    ui,
                    "Query (vertical axis; none = self-comparison)",
                    &mut inputs.query,
                    true,
                );
                ui.separator();
                let mut use_paf = inputs.paf.is_some();
                ui.horizontal(|ui| {
                    ui.radio_value(&mut use_paf, false, "Run an aligner");
                    ui.radio_value(&mut use_paf, true, "Read a PAF file");
                });
                if use_paf {
                    let mut paf = inputs.paf.clone().filter(|p| !p.as_os_str().is_empty());
                    pick(ui, "PAF (query = vertical axis)", &mut paf, false);
                    // An empty path keeps the choice until a file is picked.
                    inputs.paf = Some(paf.unwrap_or_default());
                } else {
                    inputs.paf = None;
                    egui::Grid::new("genome-aligner")
                        .num_columns(2)
                        .show(ui, |ui| {
                            ui.label("Aligner");
                            ui.text_edit_singleline(&mut inputs.aligner).on_hover_text(
                                "auto (mm2plus, else minimap2), a name on PATH, or a path",
                            );
                            ui.end_row();
                            ui.label("Preset");
                            egui::ComboBox::from_id_salt("genome-preset")
                                .selected_text(inputs.preset.clone())
                                .show_ui(ui, |ui| {
                                    for p in ["asm5", "asm10", "asm20"] {
                                        ui.selectable_value(&mut inputs.preset, p.to_string(), p);
                                    }
                                });
                            ui.end_row();
                            ui.label("Extra arguments");
                            ui.text_edit_singleline(&mut inputs.aligner_args)
                                .on_hover_text(
                                    "e.g. -c for exact identities (several times slower)",
                                );
                            ui.end_row();
                        });
                    ui.checkbox(&mut inputs.all_copies, "Keep all copies of repeats (-P)")
                        .on_hover_text(
                            "Off: each region's best placement only, so repeats \
                             are mostly hidden. On: every chain, often 10–100× \
                             more alignments.",
                        );
                    ui.weak("Results are cached; reopening the same comparison is quick.");
                }
                ui.horizontal(|ui| {
                    ui.label("Skip alignments shorter than");
                    ui.add(
                        egui::DragValue::new(&mut inputs.paf_min_len)
                            .range(0..=1_000_000)
                            .speed(100)
                            .suffix(" bp"),
                    )
                    .on_hover_text(
                        "Not read at all, to keep memory bounded for huge PAFs \
                         (e.g. -P on multi-gigabase genomes). The panel's minimum \
                         length hides alignments without dropping them.",
                    );
                });
                if let Some(e) = error {
                    ui.separator();
                    ui.colored_label(Color32::from_rgb(200, 40, 40), e.as_str());
                }
                ui.separator();
                let ready = !inputs.reference.as_os_str().is_empty()
                    && inputs
                        .paf
                        .as_ref()
                        .is_none_or(|p| !p.as_os_str().is_empty());
                if ui.add_enabled(ready, egui::Button::new("Open")).clicked() {
                    run = true;
                }
            });
        if run {
            let inputs = inputs.clone();
            self.start_genome_load(ctx, inputs);
        } else if !open {
            self.genome_dialog = None;
        }
    }

    fn draw_menu(&mut self, ctx: &Context) {
        egui::TopBottomPanel::top("menu").show(ctx, |ui| {
            // Denser menu bar: smaller text, tighter horizontal
            // padding. Doesn't affect submenu styling — those still
            // use the default body size.
            ui.style_mut().spacing.item_spacing.x = 6.0;
            ui.style_mut().spacing.button_padding = Vec2::new(4.0, 2.0);
            // Genome view: the same menu, with Fit / Back / Save PNG acting
            // on it and dotplot-only items disabled.
            let genome = self.genome_active && self.genome.is_some();
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Open genome comparison…").clicked() {
                        ui.close_menu();
                        self.open_genome_dialog();
                    }
                    if ui.button("Open query sequence (FASTA/GFF3)…").clicked() {
                        ui.close_menu();
                        self.genome_active = false;
                        pick_and_load(self, SeqRole::Query);
                    }
                    if ui.button("Open subject sequence (FASTA/GFF3)…").clicked() {
                        ui.close_menu();
                        self.genome_active = false;
                        pick_and_load(self, SeqRole::Subject);
                    }
                    ui.separator();
                    if self.settings.self_comparison {
                        if ui
                            .add_enabled(
                                !genome && self.query.is_some(),
                                egui::Button::new("Load GFF3/BED…"),
                            )
                            .clicked()
                        {
                            ui.close_menu();
                            pick_and_load_annotation(self, SeqRole::Query);
                        }
                    } else {
                        if ui
                            .add_enabled(
                                !genome && self.query.is_some(),
                                egui::Button::new("Load GFF3/BED for query…"),
                            )
                            .clicked()
                        {
                            ui.close_menu();
                            pick_and_load_annotation(self, SeqRole::Query);
                        }
                        if ui
                            .add_enabled(
                                !genome && self.subject.is_some(),
                                egui::Button::new("Load GFF3/BED for subject…"),
                            )
                            .clicked()
                        {
                            ui.close_menu();
                            pick_and_load_annotation(self, SeqRole::Subject);
                        }
                    }
                    if ui
                        .add_enabled(
                            !genome && self.query.is_some() && self.subject.is_some(),
                            egui::Button::new("Load HSPs (PAF)…"),
                        )
                        .on_hover_text(
                            "Alignments from a PAF file (minimap2, …) drawn over the \
                             plot; records not plotted are skipped",
                        )
                        .clicked()
                    {
                        ui.close_menu();
                        if let Some(path) = rfd::FileDialog::new()
                            .set_title("Load HSPs")
                            .add_filter("PAF", &["paf", "gz"])
                            .pick_file()
                        {
                            self.load_hsps(path);
                        }
                    }
                    ui.separator();
                    if ui.button("Save PNG…").clicked() {
                        ui.close_menu();
                        match (genome, &mut self.genome) {
                            (true, Some(g)) => g.export_png(),
                            _ => save_png(self),
                        }
                    }
                    if ui.button("Save SVG…").clicked() {
                        ui.close_menu();
                        match (genome, &mut self.genome) {
                            (true, Some(g)) => g.export_svg(),
                            _ => save_svg(self),
                        }
                    }
                    if ui
                        .add_enabled(!genome, egui::Button::new("Save PDF…"))
                        .on_disabled_hover_text("Not yet available for genome plots")
                        .clicked()
                    {
                        ui.close_menu();
                        save_pdf(self);
                    }
                    ui.separator();
                    if ui.button("Save session…").clicked() {
                        ui.close_menu();
                        save_session(self);
                    }
                    if ui.button("Open session…").clicked() {
                        ui.close_menu();
                        open_session(self, ctx);
                    }
                    ui.separator();
                    if ui.button("Quit").clicked() {
                        ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                    }
                });
                ui.menu_button("View", |ui| {
                    if genome {
                        if ui
                            .add_enabled(self.query.is_some(), egui::Button::new("Dot-plot view"))
                            .clicked()
                        {
                            ui.close_menu();
                            self.genome_active = false;
                        }
                    } else if self.genome.is_some() && ui.button("Genome view").clicked() {
                        ui.close_menu();
                        self.genome_active = true;
                    }
                    if ui.button("Fit (reset to canonical zoom)").clicked() {
                        ui.close_menu();
                        match (genome, &mut self.genome) {
                            (true, Some(g)) => g.fit(),
                            _ => self.action_fit(),
                        }
                    }
                    let can_back = match (genome, &self.genome) {
                        (true, Some(g)) => g.can_go_back(),
                        _ => !self.history.is_empty(),
                    };
                    if ui
                        .add_enabled(can_back, egui::Button::new("Back (Esc)"))
                        .clicked()
                    {
                        ui.close_menu();
                        match (genome, &mut self.genome) {
                            (true, Some(g)) => g.back(),
                            _ => self.action_back(),
                        }
                    }
                    ui.add_enabled_ui(!genome, |ui| {
                        ui.menu_button("Sort records", |ui| self.record_sort_menu(ui))
                    })
                    .response
                    .on_disabled_hover_text("In the genome view: Sort query along reference");
                    if ui
                        .add_enabled(!genome, egui::Button::new("Reset greyramp"))
                        .clicked()
                    {
                        self.greyramp = Greyramp::default();
                        self.texture_dirty = true;
                    }
                    ui.separator();
                    let theme_label = if self.light_theme {
                        "Switch to dark theme"
                    } else {
                        "Switch to light theme"
                    };
                    if ui.button(theme_label).clicked() {
                        self.light_theme = !self.light_theme;
                        ctx.set_visuals(if self.light_theme {
                            egui::Visuals::light()
                        } else {
                            egui::Visuals::dark()
                        });
                    }
                    ui.separator();
                    let dock_label = if self.align_dock_visible {
                        "Hide alignment view"
                    } else {
                        "Show alignment view"
                    };
                    if ui
                        .add_enabled(!genome, egui::Button::new(dock_label))
                        .clicked()
                    {
                        self.align_dock_visible = !self.align_dock_visible;
                    }
                    ui.separator();
                    // Compact shortcut crib so users can discover
                    // the keyboard map without reading the docs.
                    ui.menu_button("Keyboard shortcuts…", |ui| {
                        let style = egui::FontId::monospace(11.0);
                        let genome_keys = [
                            ("Wheel", "zoom at the pointer"),
                            ("L-drag", "pan"),
                            ("M-drag", "zoom into rectangle (or Shift-drag)"),
                            ("Dbl-click", "fit"),
                            ("R-click", "dot-plot a region, flip a record"),
                            ("Esc/Bsp", "back (pop view history)"),
                        ];
                        let dotplot_keys = [
                            (
                                "← → ↑ ↓",
                                "nudge crosshair 1 residue (Shift ×10, Ctrl ×100)",
                            ),
                            (",   .", "step along main diagonal"),
                            ("[   ]", "step along anti-diagonal"),
                            ("Space", "snap crosshair to nearest strong dot"),
                            ("L-click", "set crosshair"),
                            ("L-drag", "pan"),
                            ("M-drag", "zoom into rectangle (dotter classic)"),
                            ("Esc/Bsp", "back (pop view history)"),
                        ];
                        let keys: &[(&str, &str)] =
                            if genome { &genome_keys } else { &dotplot_keys };
                        for &(keys, desc) in keys {
                            ui.horizontal(|ui| {
                                ui.label(
                                    egui::RichText::new(format!("{keys:<10}"))
                                        .font(style.clone())
                                        .color(Color32::from_gray(70)),
                                );
                                ui.label(desc);
                            });
                        }
                    });
                    ui.separator();
                    if ui
                        .add_enabled(!genome, egui::Button::new("Settings…"))
                        .clicked()
                    {
                        self.show_settings = true;
                    }
                });
                ui.separator();
                match (genome, &self.genome) {
                    (true, Some(g)) => {
                        let (r, q) = g.input_labels();
                        ui.label(format!("reference: {r}"));
                        ui.label("·");
                        ui.label(format!("query: {q}"));
                    }
                    _ => {
                        ui.label(format_seq_label("query", self.query.as_ref()));
                        ui.label("·");
                        ui.label(format_seq_label("subject", self.subject.as_ref()));
                    }
                }
            });
        });
    }

    fn draw_greyramp_panel(&mut self, ctx: &Context) {
        egui::SidePanel::right("side_panel")
            .resizable(false)
            .default_width(220.0)
            .show(ctx, |ui| {
                // Compact, information-dense layout per user
                // feedback: smaller fonts, tighter vertical spacing.
                ui.style_mut().spacing.item_spacing = Vec2::new(4.0, 3.0);
                ui.style_mut().spacing.interact_size.y = 18.0;

                // ── View ──
                ui.label(
                    egui::RichText::new("View")
                        .strong()
                        .size(13.0)
                        .color(Color32::from_gray(50)),
                );
                ui.horizontal(|ui| {
                    if ui
                        .button("Fit")
                        .on_hover_text(
                            "Restore the canonical (display-matched) zoom showing the whole \
                             pixelmap. Pushes the current view onto the Back stack.",
                        )
                        .clicked()
                    {
                        self.action_fit();
                    }
                    let back_enabled = !self.history.is_empty();
                    if ui
                        .add_enabled(back_enabled, egui::Button::new("Back"))
                        .on_hover_text("Pop the previous view from the history stack (Esc).")
                        .clicked()
                    {
                        self.action_back();
                    }
                });
                ui.separator();

                // ── Sequences summary ──
                ui.label(
                    egui::RichText::new("Sequences")
                        .strong()
                        .size(13.0)
                        .color(Color32::from_gray(50)),
                );
                draw_seq_summary(ui, "Query  ", self.query.as_ref());
                draw_seq_summary(ui, "Subject", self.subject.as_ref());
                if let Some(p) = &self.plot {
                    ui.label(
                        egui::RichText::new(format!(
                            "Pixelmap {}×{}  W={}",
                            p.width, p.height, p.params.window_size
                        ))
                        .size(11.0)
                        .color(Color32::from_gray(80)),
                    );
                }
                ui.separator();

                if self.large_plot.is_some() {
                    self.large_panel_ui(ui);
                    ui.separator();
                }
                if self.hsps.is_some() {
                    self.hsp_panel_ui(ui);
                    ui.separator();
                }

                // ── Greyramp ── (replaced by the coverage ramp in the
                // large-mode panel when shading by coverage)
                if !self.coverage_shading() {
                    ui.label(
                        egui::RichText::new("Greyramp")
                            .strong()
                            .size(13.0)
                            .color(Color32::from_gray(50)),
                    );
                    let small = egui::RichText::new("White").size(11.0);
                    ui.label(small);
                    if ui
                        .add(
                            Slider::new(&mut self.greyramp.white, 0..=255)
                                .clamping(egui::SliderClamping::Always),
                        )
                        .changed()
                    {
                        self.texture_dirty = true;
                        // Ridge threshold tracks greyramp.white, so the
                        // overlay's "lit" definition follows the raster's
                        // noise floor.
                        self.refresh_ridges();
                    }
                    ui.label(egui::RichText::new("Black").size(11.0));
                    if ui
                        .add(
                            Slider::new(&mut self.greyramp.black, 0..=255)
                                .clamping(egui::SliderClamping::Always),
                        )
                        .changed()
                    {
                        self.texture_dirty = true;
                    }
                    ui.horizontal(|ui| {
                        if ui.small_button("Swap").clicked() {
                            self.greyramp.swap = !self.greyramp.swap;
                            self.texture_dirty = true;
                        }
                        if ui.small_button("Reset").clicked() {
                            self.greyramp = Greyramp::default();
                            self.texture_dirty = true;
                            self.refresh_ridges();
                        }
                    });
                    ui.add_space(2.0);
                    ui.label(egui::RichText::new("LUT").size(10.0));
                    let lut = self.greyramp.lut();
                    let (rect, _) = ui
                        .allocate_exact_size(Vec2::new(ui.available_width(), 20.0), Sense::hover());
                    let painter = ui.painter();
                    for (x, &g) in lut.iter().enumerate() {
                        let xp = rect.left() + (x as f32 / 256.0) * rect.width();
                        let xw = (rect.width() / 256.0).max(1.0);
                        painter.rect_filled(
                            Rect::from_min_max(
                                Pos2::new(xp, rect.top()),
                                Pos2::new(xp + xw, rect.bottom()),
                            ),
                            0.0,
                            Color32::from_gray(g),
                        );
                    }
                    painter.rect_stroke(rect, 0.0, egui::Stroke::new(1.0, Color32::from_gray(120)));
                    ui.add_space(6.0);
                    ui.separator();
                }

                // ── Vector overlays ──
                if ui
                    .checkbox(&mut self.settings.show_ridge_overlay, "Show vector ridges")
                    .changed()
                {
                    self.refresh_ridges();
                }

                // ── Annotations (display only; bands drawn live over the
                // pixelmap, so no recompute / texture invalidation needed) ──
                if self.query_annot.is_some() || self.subject_annot.is_some() {
                    ui.add_space(6.0);
                    ui.separator();
                    ui.label(
                        egui::RichText::new("Annotations")
                            .strong()
                            .size(13.0)
                            .color(Color32::from_gray(50)),
                    );
                    ui.checkbox(&mut self.annot_show, "Show bands");
                    ui.horizontal(|ui| {
                        ui.label(egui::RichText::new("Opacity").size(11.0));
                        ui.add(Slider::new(&mut self.annot_alpha, 0.0..=1.0));
                    })
                    .response
                    .on_hover_text(
                        "Dim annotation bands so the dots show through (0 = invisible).",
                    );

                    // One legend when both axes share the same set
                    // (self-comparison), per-axis otherwise.
                    let dedupe = self.settings.self_comparison
                        && self.query_annot.is_some()
                        && self.subject_annot.is_some();
                    egui::ScrollArea::vertical()
                        .max_height(220.0)
                        .auto_shrink([false, true])
                        .show(ui, |ui| {
                            if dedupe {
                                if let Some(ax) = self.query_annot.as_mut() {
                                    annotation_legend_axis(ui, ax, None);
                                }
                            } else {
                                if let Some(ax) = self.query_annot.as_mut() {
                                    annotation_legend_axis(ui, ax, Some("Query"));
                                }
                                if let Some(ax) = self.subject_annot.as_mut() {
                                    annotation_legend_axis(ui, ax, Some("Subject"));
                                }
                            }
                        });

                    // In self-comparison the single legend drives both axes:
                    // mirror the query's visibility set onto the subject so a
                    // toggle hides both the vertical (query) and horizontal
                    // (subject) bands, not just the query's.
                    if dedupe {
                        if let Some(hidden) = self.query_annot.as_ref().map(|ax| ax.hidden.clone())
                        {
                            if let Some(s) = self.subject_annot.as_mut() {
                                s.hidden = hidden;
                            }
                        }
                    }
                }

                // ── Keyboard shortcuts (reference) ──
                ui.add_space(6.0);
                ui.separator();
                ui.label(
                    egui::RichText::new("Keyboard shortcuts")
                        .strong()
                        .size(13.0)
                        .color(Color32::from_gray(50)),
                );
                let mono = egui::FontId::monospace(11.0);
                for (keys, desc) in [
                    ("← → ↑ ↓", "nudge 1 res (⇧×10, ⌃×100)"),
                    (",   .", "step along main diagonal"),
                    ("[   ]", "step along anti-diagonal"),
                    ("Space", "snap to nearest strong dot"),
                    ("L-click", "set crosshair"),
                    ("L-drag", "pan"),
                    ("M-drag", "zoom into rectangle"),
                    ("Esc/Bsp", "back (pop view history)"),
                ] {
                    ui.horizontal(|ui| {
                        ui.label(
                            egui::RichText::new(format!("{keys:<8}"))
                                .font(mono.clone())
                                .color(Color32::from_gray(70)),
                        );
                        ui.label(egui::RichText::new(desc).size(11.0));
                    });
                }
            });
    }

    fn axis_annot(&self, role: SeqRole) -> Option<&crate::annotation_overlay::AxisAnnot> {
        match role {
            SeqRole::Query => self.query_annot.as_ref(),
            SeqRole::Subject => self.subject_annot.as_ref(),
        }
    }

    fn axis_annot_mut(
        &mut self,
        role: SeqRole,
    ) -> Option<&mut crate::annotation_overlay::AxisAnnot> {
        match role {
            SeqRole::Query => self.query_annot.as_mut(),
            SeqRole::Subject => self.subject_annot.as_mut(),
        }
    }

    fn clear_annot(&mut self, role: SeqRole) {
        match role {
            SeqRole::Query => self.query_annot = None,
            SeqRole::Subject => self.subject_annot = None,
        }
    }

    /// One axis's annotation controls in the Settings panel: load/clear,
    /// a color-by dropdown (GFF3) or single swatch (BED), the value
    /// legend with per-value color + visibility, and a skipped-feature
    /// warning.
    fn annotation_axis_ui(&mut self, ui: &mut egui::Ui, role: SeqRole) {
        let (label, has_seq) = match role {
            SeqRole::Query => ("Query", self.query.is_some()),
            SeqRole::Subject => ("Subject", self.subject.is_some()),
        };
        let name = self
            .axis_annot(role)
            .and_then(|ax| ax.set.source_path.as_ref())
            .and_then(|p| p.file_name())
            .map(|f| f.to_string_lossy().into_owned())
            .unwrap_or_else(|| "—".to_string());

        ui.horizontal(|ui| {
            ui.label(format!("{label}: {name}"));
            if ui
                .add_enabled(has_seq, egui::Button::new("Load…"))
                .on_hover_text("Load a GFF3 or BED file for this axis")
                .clicked()
            {
                pick_and_load_annotation(self, role);
            }
            if self.axis_annot(role).is_some() && ui.button("Clear").clicked() {
                self.clear_annot(role);
            }
        });

        let Some(ax) = self.axis_annot_mut(role) else {
            return;
        };

        // Color-by dropdown (GFF3 only — BED has a single color).
        if !ax.is_bed() {
            let keys: Vec<String> = ax.set.attribute_keys().into_iter().collect();
            if !keys.is_empty() {
                let current = ax.color_by.clone();
                ui.horizontal(|ui| {
                    ui.label("Color by:");
                    egui::ComboBox::from_id_salt((label, "annot-colorby"))
                        .selected_text(&current)
                        .show_ui(ui, |ui| {
                            for k in &keys {
                                if ui.selectable_label(current == *k, k).clicked() {
                                    ax.set_color_by(k.clone());
                                }
                            }
                        });
                });
            }
        }

        // Legend: per-value visibility toggle + color swatch + count,
        // sorted by descending feature count.
        let counts = ax.value_counts();
        for (value, count) in &counts {
            ui.horizontal(|ui| {
                let mut visible = !ax.hidden.contains(value);
                if ui.checkbox(&mut visible, "").changed() {
                    if visible {
                        ax.hidden.remove(value);
                    } else {
                        ax.hidden.insert(value.clone());
                    }
                }
                if ax.is_bed() {
                    ui.color_edit_button_srgba(&mut ax.bed_color);
                } else if let Some(col) = ax.palette.get_mut(value) {
                    ui.color_edit_button_srgba(col);
                }
                ui.label(format!("{value}  ({count})"));
            });
        }

        if ax.skipped > 0 {
            ui.colored_label(
                Color32::from_rgb(200, 140, 40),
                format!("⚠ {} feature(s) skipped (record id not found)", ax.skipped),
            );
        }
    }

    fn draw_settings_window(&mut self, ctx: &Context) {
        let mut open = self.show_settings;
        let mut changed = false;
        egui::Window::new("Settings")
            .open(&mut open)
            .resizable(false)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.label("Mode:");
                    let mut mode_label = format!("{:?}", self.settings.mode);
                    egui::ComboBox::from_id_salt("mode")
                        .selected_text(&mode_label)
                        .show_ui(ui, |ui| {
                            for (m, label) in [
                                (BlastMode::Blastn, "Blastn"),
                                (BlastMode::Blastp, "Blastp"),
                                (BlastMode::Blastx, "Blastx"),
                            ] {
                                if ui
                                    .selectable_value(&mut self.settings.mode, m, label)
                                    .changed()
                                {
                                    changed = true;
                                    mode_label = label.to_string();
                                    // Sensible matrix default per mode.
                                    self.settings.matrix_name = match m {
                                        BlastMode::Blastn => "DNA+5/-4".into(),
                                        _ => "BLOSUM62".into(),
                                    };
                                }
                            }
                        });
                });

                if self.settings.mode != BlastMode::Blastn {
                    ui.horizontal(|ui| {
                        ui.label("Matrix:");
                        egui::ComboBox::from_id_salt("matrix")
                            .selected_text(&self.settings.matrix_name)
                            .show_ui(ui, |ui| {
                                for name in [
                                    "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90",
                                    "PAM30", "PAM70", "PAM250",
                                ] {
                                    if ui
                                        .selectable_value(
                                            &mut self.settings.matrix_name,
                                            name.to_string(),
                                            name,
                                        )
                                        .changed()
                                    {
                                        changed = true;
                                    }
                                }
                            });
                    });
                }

                ui.horizontal(|ui| {
                    ui.label("Window size:");
                    let mut auto = self.settings.window_size.is_none();
                    if ui.checkbox(&mut auto, "auto (Karlin)").changed() {
                        if auto {
                            self.settings.window_size = None;
                        } else {
                            self.settings.window_size = Some(15);
                        }
                        changed = true;
                    }
                    if let Some(w) = &mut self.settings.window_size {
                        if ui.add(Slider::new(w, 1..=200)).changed() {
                            changed = true;
                        }
                    }
                });

                ui.horizontal(|ui| {
                    ui.label("Zoom:");
                    if ui
                        .add(Slider::new(&mut self.settings.zoom, 1..=64))
                        .changed()
                    {
                        changed = true;
                    }
                    // Restore the canonical (display-matched) compute
                    // zoom — same as the View menu's Fit action.
                    if ui
                        .button("Fit")
                        .on_hover_text(
                            "Restore the canonical (display-matched) zoom showing the whole \
                             pixelmap. Pushes the current view onto the Back stack.",
                        )
                        .clicked()
                    {
                        self.action_fit();
                    }
                });
                ui.horizontal(|ui| {
                    // Auto-fit fires on sequence load (deferred until
                    // the first canvas paint measures the plot area)
                    // and on subsequent window resizes that change
                    // the canvas by > 25 %. The slider above is then
                    // driven by auto-fit; turn the checkbox off to
                    // take manual control of `Zoom`.
                    ui.checkbox(
                        &mut self.settings.auto_fit_compute_zoom,
                        "Auto-fit computation zoom",
                    )
                    .on_hover_text(
                        "Pick `Zoom` so the computed pixelmap matches the physical canvas size \
                         — dotter's invariant: compute at display scale, render near 1:1. \
                         Triggers on sequence load and when the window resizes by more than \
                         25 %. Turn off to drive `Zoom` manually (e.g. for reproducible figures).",
                    );
                });
                ui.horizontal(|ui| {
                    let prev_auto = self.settings.auto_pixel_fac;
                    let resp = ui
                        .checkbox(
                            &mut self.settings.auto_pixel_fac,
                            "Auto pixel factor (Karlin)",
                        )
                        .on_hover_text(
                            "Auto-derive the pixel factor from Karlin's expected residue \
                             score, dotter's default: positions an expected match residue \
                             at ~1/5 of the displayable range and pushes background noise \
                             toward white. Uncheck to drive the slider manually.",
                        );
                    if resp.changed() {
                        changed = true;
                        // Seed the slider from the last resolved value
                        // when the user toggles auto off — so the
                        // slider lands at "what auto picked" instead of
                        // a surprise jump.
                        if prev_auto && !self.settings.auto_pixel_fac {
                            if let Some(plot) = self.plot.as_ref() {
                                self.settings.pixel_fac = plot.params.pixel_fac.clamp(1, 255);
                            }
                        }
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("Pixel factor:");
                    let mut display_value = if self.settings.auto_pixel_fac {
                        // While auto, show the resolved value from the
                        // last compute (read-only).
                        self.plot
                            .as_ref()
                            .map(|p| p.params.pixel_fac)
                            .unwrap_or(self.settings.pixel_fac)
                    } else {
                        self.settings.pixel_fac
                    };
                    let slider = Slider::new(&mut display_value, 1..=255);
                    let resp = ui.add_enabled(!self.settings.auto_pixel_fac, slider);
                    if !self.settings.auto_pixel_fac
                        && resp.changed()
                        && display_value != self.settings.pixel_fac
                    {
                        self.settings.pixel_fac = display_value;
                        changed = true;
                    }
                });

                if self.settings.mode == BlastMode::Blastn {
                    ui.horizontal(|ui| {
                        ui.label("Strand:");
                        for (val, label) in [
                            (Strand::Forward, "Forward"),
                            (Strand::Reverse, "Reverse"),
                            (Strand::Both, "Both"),
                        ] {
                            if ui
                                .selectable_value(&mut self.settings.strand, val, label)
                                .changed()
                            {
                                changed = true;
                            }
                        }
                    });
                    // -r / -v: pre-flip the axis sequence before
                    // compute. Spec §4.1.10.
                    ui.horizontal(|ui| {
                        if ui
                            .checkbox(
                                &mut self.settings.reverse_query,
                                "Reverse-complement query (-r)",
                            )
                            .changed()
                        {
                            changed = true;
                        }
                    });
                    ui.horizontal(|ui| {
                        if ui
                            .checkbox(
                                &mut self.settings.reverse_subject,
                                "Reverse-complement subject (-v)",
                            )
                            .changed()
                        {
                            changed = true;
                        }
                    });
                    ui.horizontal(|ui| {
                        if ui
                            .checkbox(
                                &mut self.settings.inverted_repeat_colour,
                                "Highlight inverted repeats (magenta reverse strand)",
                            )
                            .changed()
                        {
                            changed = true;
                        }
                    });
                }

                if ui
                    .checkbox(&mut self.settings.self_comparison, "Self-comparison")
                    .changed()
                {
                    changed = true;
                }
                if self.settings.self_comparison {
                    ui.horizontal(|ui| {
                        ui.label("Triangle:");
                        for (val, label) in [
                            (Triangle::Both, "Both"),
                            (Triangle::Upper, "Upper"),
                            (Triangle::Lower, "Lower"),
                        ] {
                            if ui
                                .selectable_value(&mut self.settings.triangle, val, label)
                                .changed()
                            {
                                changed = true;
                            }
                        }
                    });
                }

                ui.separator();
                ui.horizontal(|ui| {
                    ui.label("Memory cap (MiB):");
                    let mut mib: u32 = (self.settings.memory_limit_bytes / (1024 * 1024))
                        .clamp(8, u32::MAX as u64) as u32;
                    if ui
                        .add(
                            Slider::new(&mut mib, 8..=16_384)
                                .logarithmic(true)
                                .clamping(egui::SliderClamping::Always),
                        )
                        .changed()
                    {
                        self.settings.memory_limit_bytes = (mib as u64) * 1024 * 1024;
                        changed = true;
                    }
                });
                ui.weak(
                    "Refuses to allocate a pixelmap larger than this. \
                     1 GiB suits ~32k × 32k at zoom 1. Halve the cap when \
                     zoom doubles.",
                );

                ui.separator();
                ui.heading("Large sequences");
                changed |= self.large_settings_ui(ui);

                // Vector ridge overlay — display-only, no recompute.
                // Tied to its own local `ridges_changed` flag so
                // settings panel adjustments don't trigger the full
                // pixelmap recompute path.
                let mut ridges_changed = false;
                ui.separator();
                ui.heading("Ridge overlay");
                ui.horizontal(|ui| {
                    if ui
                        .checkbox(&mut self.settings.show_ridge_overlay, "Show vector ridges")
                        .on_hover_text(
                            "Draw anti-aliased line segments over coherent diagonal runs \
                             detected in the raster. Hides per-window intensity oscillation \
                             on imperfect-homology lines. The raster underneath is unchanged \
                             — toggle off to see the data view.",
                        )
                        .changed()
                    {
                        ridges_changed = true;
                    }
                });
                if self.settings.show_ridge_overlay {
                    ui.horizontal(|ui| {
                        ui.label("Min length (cells):");
                        if ui
                            .add(
                                egui::DragValue::new(&mut self.settings.ridge_min_length)
                                    .range(1..=200),
                            )
                            .changed()
                        {
                            ridges_changed = true;
                        }
                    });
                    ui.horizontal(|ui| {
                        ui.label("Max gap (cells):");
                        if ui
                            .add(
                                egui::DragValue::new(&mut self.settings.ridge_max_gap)
                                    .range(0..=20),
                            )
                            .changed()
                        {
                            ridges_changed = true;
                        }
                    });
                    ui.weak(
                        "Threshold tracks the greyramp `White` slider — drop it to hide noise \
                         from both the raster AND the overlay in tandem.",
                    );
                }

                // Annotation overlay (ADR 0005) — display-only, no
                // recompute (bands are drawn straight over the existing
                // pixelmap).
                ui.separator();
                ui.heading("Annotations");
                ui.label(
                    egui::RichText::new(
                        "Show-bands toggle, opacity, and the color legend are in the side panel.",
                    )
                    .size(11.0)
                    .color(Color32::from_gray(120)),
                );
                self.annotation_axis_ui(ui, SeqRole::Query);
                self.annotation_axis_ui(ui, SeqRole::Subject);

                ui.separator();
                if ui.button("Apply").clicked() {
                    changed = true;
                }
                if ridges_changed && !changed {
                    self.refresh_ridges();
                }
            });
        self.show_settings = open;
        if changed {
            // Floors and self-comparison detection live in the source.
            self.large_source = None;
            // Any setting that affects the pixelmap invalidates the
            // cache + history (which would point at incompatible
            // pixelmaps under the new parameters).
            self.invalidate_caches();
            self.recompute();
        }
    }

    /// Large-mode parameters in the Settings window. Returns whether a
    /// recompute is needed.
    fn large_settings_ui(&mut self, ui: &mut egui::Ui) -> bool {
        let ls = &mut self.settings.large;
        let mut changed = false;
        ui.horizontal(|ui| {
            ui.label("Engine:");
            for (v, label, tip) in [
                (
                    EngineChoice::Auto,
                    "Auto",
                    "Seed-based large engine for big views, exact Dotter engine once zoomed in.",
                ),
                (
                    EngineChoice::Exact,
                    "Exact",
                    "Always score every cell (Dotter).",
                ),
                (
                    EngineChoice::Large,
                    "Large",
                    "Always use the seed-based engine (DNA only).",
                ),
            ] {
                changed |= ui
                    .selectable_value(&mut ls.engine, v, label)
                    .on_hover_text(tip)
                    .changed();
            }
        });
        if ls.engine == EngineChoice::Auto {
            ui.horizontal(|ui| {
                ui.label("Switch above (Mcells):");
                let mut mcells = ls.threshold_cells / 1e6;
                if ui
                    .add(
                        egui::DragValue::new(&mut mcells)
                            .range(1.0..=1e9)
                            .speed(100.0),
                    )
                    .on_hover_text(
                        "Views with more query × subject cells than this use the large engine.",
                    )
                    .changed()
                {
                    ls.threshold_cells = mcells * 1e6;
                    changed = true;
                }
            });
        }
        if ls.engine == EngineChoice::Exact {
            return changed;
        }
        ui.horizontal(|ui| {
            ui.label("Seed:");
            changed |= ui
                .add(egui::TextEdit::singleline(&mut ls.seed_pattern).desired_width(150.0))
                .on_hover_text("Spaced seed, 1 = must match. Weight 6–13 (default 12-of-19).")
                .lost_focus();
        });
        ui.horizontal(|ui| {
            ui.label("Two-hit distance:");
            changed |= ui
                .add(egui::DragValue::new(&mut ls.two_hit_dist).range(10..=2000))
                .changed();
            ui.label("Max occ:");
            changed |= ui
                .add(egui::DragValue::new(&mut ls.max_occ).range(10..=100_000))
                .on_hover_text(
                    "Seeds more frequent than this are skipped (and count as high-copy).",
                )
                .changed();
        });
        let rate_ui = |ui: &mut egui::Ui, label: &str, tip: &str, rate: &mut Option<u32>| {
            let mut changed = false;
            ui.horizontal(|ui| {
                ui.label(label).on_hover_text(tip);
                let mut auto = rate.is_none();
                if ui.checkbox(&mut auto, "auto").changed() {
                    *rate = if auto { None } else { Some(1) };
                    changed = true;
                }
                if let Some(r) = rate {
                    changed |= ui
                        .add(egui::DragValue::new(r).range(1..=64).prefix("1/"))
                        .changed();
                }
            });
            changed
        };
        changed |= rate_ui(
            ui,
            "Seed sampling:",
            "Look up one seed in N for the plot (the two-hit distance scales by N). \
             Auto: zoom / 6000, at most 1/32, so only overview views are sampled.",
            &mut ls.seed_sample,
        );
        if ls.depth == DepthChoice::Seed {
            changed |= rate_ui(
                ui,
                "Depth sampling:",
                "Seed sampling for the copy-number pass. Auto: length / 25 Mb, at most 1/16.",
                &mut ls.depth_sample,
            );
        }
        ui.horizontal(|ui| {
            ui.label("Repeat depth:");
            for (v, label, tip) in [
                (
                    DepthChoice::Seed,
                    "Seed",
                    "Copies found by a self-comparison seed pass; sees diverged repeats.",
                ),
                (
                    DepthChoice::Kmer,
                    "k-mer",
                    "Exact k-mer counts; fast, young repeats only.",
                ),
                (DepthChoice::None, "None", "No repeat classes."),
            ] {
                changed |= ui
                    .selectable_value(&mut ls.depth, v, label)
                    .on_hover_text(tip)
                    .changed();
            }
            if ls.depth == DepthChoice::Kmer {
                changed |= ui
                    .add(
                        egui::DragValue::new(&mut ls.depth_k)
                            .range(11..=31)
                            .prefix("k="),
                    )
                    .changed();
            }
        });
        ui.horizontal(|ui| {
            ui.label("Copies ≥ low/mid/high:");
            let t = &mut ls.thresholds;
            let before = *t;
            ui.add(egui::DragValue::new(&mut t.low).range(2..=1000));
            ui.add(egui::DragValue::new(&mut t.mid).range(3..=10_000));
            ui.add(egui::DragValue::new(&mut t.high).range(4..=60_000));
            if t.validate().is_err() {
                *t = before;
            }
            changed |= *t != before;
        });
        changed |= ui
            .checkbox(
                &mut ls.softmask_as_repeat,
                "Soft-masked (lowercase) bases are repeats",
            )
            .changed();
        ui.horizontal(|ui| {
            ui.label("Memory cap (MiB):");
            let mut mib = (ls.memory_limit_bytes >> 20) as u32;
            if ui
                .add(Slider::new(&mut mib, 256..=65_536).logarithmic(true))
                .changed()
            {
                ls.memory_limit_bytes = (mib as u64) << 20;
                changed = true;
            }
        });
        changed
    }

    /// Side-panel controls for a large-engine plot: display-only.
    fn hsp_panel_ui(&mut self, ui: &mut egui::Ui) {
        let Some(h) = &self.hsps else {
            return;
        };
        ui.label(
            egui::RichText::new("HSPs")
                .strong()
                .size(13.0)
                .color(Color32::from_gray(50)),
        );
        let name = h
            .path
            .file_name()
            .map_or_else(String::new, |f| f.to_string_lossy().into_owned());
        ui.label(
            egui::RichText::new(format!("{name}: {} shown", h.hsps.len()))
                .size(11.0)
                .color(Color32::from_gray(80)),
        )
        .on_hover_text(format!(
            "{}\n{} alignments read, {} for records not plotted",
            h.path.display(),
            h.alignments,
            h.skipped
        ));
        let mut mode = self.hsp_mode;
        ui.horizontal_wrapped(|ui| {
            for (m, label, tip) in [
                (HspMode::Off, "Off", "Hide the HSPs"),
                (HspMode::Line, "Line", "Red lines over the dots"),
                (HspMode::Score, "Score", "Lines coloured by identity"),
                (
                    HspMode::Greyscale,
                    "Grey",
                    "HSPs only, grey by identity, instead of the dots",
                ),
            ] {
                ui.selectable_value(&mut mode, m, egui::RichText::new(label).size(11.0))
                    .on_hover_text(tip);
            }
        });
        if mode != self.hsp_mode {
            self.hsp_mode = mode;
            self.texture_dirty = true;
        }
        if mode == HspMode::Score {
            ui.horizontal(|ui| {
                for (k, label) in ["<25", "<50", "<75", "≥75 %"].iter().enumerate() {
                    let [r, g, b] = dottir_io::genome_export::IDENTITY_COLOURS[k];
                    let (rect, _) = ui.allocate_exact_size(Vec2::splat(10.0), Sense::hover());
                    ui.painter()
                        .rect_filled(rect, 2.0, Color32::from_rgb(r, g, b));
                    ui.label(egui::RichText::new(*label).size(10.0));
                }
            });
        }
        if ui.small_button("Remove").clicked() {
            self.hsps = None;
            self.texture_dirty = true;
        }
    }

    fn large_panel_ui(&mut self, ui: &mut egui::Ui) {
        ui.label(
            egui::RichText::new("Repeat classes")
                .strong()
                .size(13.0)
                .color(Color32::from_gray(50)),
        );
        let mut redraw = false;
        let thresholds = self.settings.large.thresholds;
        let ranges = [
            format!("< {}", thresholds.low),
            format!("{}–{}", thresholds.low, thresholds.mid - 1),
            format!("{}–{}", thresholds.mid, thresholds.high - 1),
            format!("≥ {}", thresholds.high),
        ];
        let ls = &mut self.settings.large;
        for c in RepeatClass::ALL {
            ui.horizontal(|ui| {
                let mut on = ls.classes.contains(c);
                if ui.checkbox(&mut on, "").changed() {
                    ls.classes = ls.classes.with(c, on);
                    redraw = true;
                }
                let (rect, _) = ui.allocate_exact_size(Vec2::splat(10.0), Sense::hover());
                ui.painter()
                    .rect_filled(rect, 2.0, CLASS_COLOURS[c.index()]);
                ui.label(
                    egui::RichText::new(format!("{} ({} copies)", c.name(), ranges[c.index()]))
                        .size(11.0),
                );
            });
        }
        ui.horizontal(|ui| {
            if ui
                .small_button("Synteny")
                .on_hover_text("Unique + low only")
                .clicked()
            {
                ls.classes = ClassMask::from_classes(&[RepeatClass::Unique, RepeatClass::Low]);
                redraw = true;
            }
            if ui
                .small_button("Repeats")
                .on_hover_text("Mid + high only")
                .clicked()
            {
                ls.classes = ClassMask::from_classes(&[RepeatClass::Mid, RepeatClass::High]);
                redraw = true;
            }
            if ui.small_button("All").clicked() {
                ls.classes = ClassMask::ALL;
                redraw = true;
            }
        });
        redraw |= ui
            .checkbox(&mut ls.colour_by_class, "Colour by class")
            .changed();

        // ── Sensitivity ──
        let mut refit = false;
        let mut retexture = false;
        ui.add_space(4.0);
        ui.label(
            egui::RichText::new("Sensitivity")
                .strong()
                .size(13.0)
                .color(Color32::from_gray(50)),
        );
        ui.horizontal(|ui| {
            let before = ls.shade;
            ui.selectable_value(&mut ls.shade, Shade::Coverage, "Matching bp")
                .on_hover_text(
                    "Darkness = bases of similarity inside the pixel. \
                     Sensitive at any zoom: at overview scale the best \
                     window score is maximal almost everywhere.",
                );
            ui.selectable_value(&mut ls.shade, Shade::Score, "Window score")
                .on_hover_text("Darkness = best window score, as in Dotter (greyramp below).");
            if ls.shade != before {
                redraw = true;
                refit = ls.shade == Shade::Coverage && ls.cov_auto;
            }
        });
        if ls.shade == Shade::Coverage {
            let has_cov = self.large_plot.as_ref().is_some_and(|lp| lp.has_coverage());
            if !has_cov {
                ui.label(
                    egui::RichText::new("This view has no coverage; recompute to use it.")
                        .size(11.0)
                        .color(Color32::from_rgb(160, 80, 0)),
                );
            }
            let white = ui
                .add(
                    Slider::new(&mut ls.cov_white_bp, 1.0..=1e8)
                        .logarithmic(true)
                        .integer()
                        .text("bp: white below"),
                )
                .on_hover_text("Pixels with less matching sequence than this are hidden.");
            let black = ui
                .add(
                    Slider::new(&mut ls.cov_black_bp, 2.0..=1e9)
                        .logarithmic(true)
                        .integer()
                        .text("bp: black at"),
                )
                .on_hover_text("Pixels with at least this much matching sequence are black.");
            if white.changed() || black.changed() {
                if white.changed() && ls.cov_black_bp <= ls.cov_white_bp {
                    ls.cov_black_bp = ls.cov_white_bp * 2.0;
                }
                if black.changed() && ls.cov_white_bp >= ls.cov_black_bp {
                    ls.cov_white_bp = (ls.cov_black_bp / 2.0).max(1.0);
                }
                ls.cov_auto = false;
                retexture = true;
            }
            ui.horizontal(|ui| {
                if ui
                    .checkbox(&mut ls.cov_auto, "Auto-fit each view")
                    .on_hover_text("White at the 2nd, black at the 95th percentile of lit pixels.")
                    .changed()
                    && ls.cov_auto
                {
                    refit = true;
                }
                if ui.small_button("Fit now").clicked() {
                    refit = true;
                }
            });
        }
        ui.checkbox(&mut ls.show_axis_tracks, "Repeat tracks on axes")
            .on_hover_text(
                "Share of low/mid/high-copy bases per pixel column/row, and seeds \
                 skipped by the max-occ cap (red). Hover for copy numbers.",
            );
        ui.add_space(4.0);
        ui.checkbox(&mut ls.show_chains, "Show synteny chains")
            .on_hover_text("Co-linear chains of unique/low-copy runs (not alignments).");
        let mut save_chains = false;
        if ls.show_chains {
            ui.add(
                Slider::new(&mut ls.chain_min_score, 100.0..=100_000.0)
                    .logarithmic(true)
                    .text("min score"),
            );
            save_chains = ui
                .small_button("Save chains…")
                .on_hover_text(
                    "The chains above the min score as PAF, or BEDPE (.bedpe), in \
                     source-record coordinates",
                )
                .clicked();
        }
        if let Some(lp) = &self.large_plot {
            let shown = lp
                .chains
                .iter()
                .filter(|c| c.score >= self.settings.large.chain_min_score as f64)
                .count();
            let st = lp.stats;
            ui.label(
                egui::RichText::new(format!(
                    "{} anchors · {} runs · {} chains · seeds 1/{}{}",
                    format_thousands(st.anchors as usize),
                    format_thousands(st.runs as usize),
                    shown,
                    lp.params.sample,
                    if st.capped_query_positions > 0 {
                        format!(
                            " · {} capped",
                            format_thousands(st.capped_query_positions as usize)
                        )
                    } else {
                        String::new()
                    }
                ))
                .size(11.0)
                .color(Color32::from_gray(80)),
            );
        }
        if redraw {
            self.rederive_large_plot();
        }
        if refit {
            self.auto_coverage_ramp();
        }
        if retexture {
            self.texture_dirty = true;
        }
        if save_chains {
            self.save_chains();
        }
    }

    fn draw_status_bar(&self, ctx: &Context) {
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if self.compute_in_flight {
                    ui.spinner();
                    let (stage, f) = self.worker.progress();
                    let text = match stage {
                        crate::compute_worker::Stage::Indexing => {
                            format!("Indexing for large mode… {:.0}%", f * 100.0)
                        }
                        crate::compute_worker::Stage::Large => {
                            format!("Large engine… {:.0}%", f * 100.0)
                        }
                        _ => "Recomputing…".to_string(),
                    };
                    ui.label(text);
                    ui.separator();
                }
                if let Some(err) = &self.last_error {
                    ui.colored_label(Color32::from_rgb(220, 90, 70), err);
                    return;
                }
                if let Some(plot) = &self.plot {
                    let engine = if self.large_plot.is_some() {
                        "large"
                    } else {
                        "exact"
                    };
                    ui.label(format!(
                        "{engine} engine · pixelmap {}×{}, W={}",
                        plot.width, plot.height, plot.params.window_size
                    ));
                    if let Some((q_seq, s_seq)) = self.crosshair {
                        // Crosshair is in absolute residue coords;
                        // map to slice-local pixmap pixel to read the
                        // raster value (0 if the crosshair is outside
                        // the current slice).
                        let z = plot.params.zoom.max(1) as usize;
                        let (q_off, s_off) = self
                            .current_slice
                            .as_ref()
                            .map(|sl| (sl.q_range.start, sl.s_range.start))
                            .unwrap_or((0, 0));
                        let v = if (q_seq as usize) >= q_off && (s_seq as usize) >= s_off {
                            let qp = ((q_seq as usize) - q_off) / (z * self.qx());
                            let sp = ((s_seq as usize) - s_off) / z;
                            if qp < plot.width as usize && sp < plot.height as usize {
                                plot.pixels[sp * (plot.width as usize) + qp]
                            } else {
                                0
                            }
                        } else {
                            0
                        };
                        ui.separator();
                        ui.label(format!(
                            "q = {}, s = {}, value = {}",
                            format_coord(self.query.as_ref(), q_seq as usize),
                            format_coord(self.subject.as_ref(), s_seq as usize),
                            v,
                        ));
                    } else {
                        ui.separator();
                        ui.label(
                            "left-click = crosshair · left-drag = pan · \
                             middle-drag = zoom rectangle · Esc = Back",
                        );
                    }
                } else {
                    ui.label("load a query and subject FASTA to begin (File menu)");
                }
            });
        });
    }

    /// H2: bottom dock showing the sequence context around the
    /// crosshair as a 3-row (query / match-line / subject) monospace
    /// alignment with per-column background colour:
    ///
    /// * green — identical residues
    /// * yellow — positive-score non-identical (per the loaded matrix)
    /// * grey + `-` — out-of-bounds at slice edges
    /// * none — other (mismatch / non-positive substitution)
    fn draw_alignment_dock(&mut self, ctx: &Context) {
        egui::TopBottomPanel::bottom("alignment_dock")
            .resizable(false)
            .min_height(70.0)
            .show(ctx, |ui| {
                let Some(plot) = self.plot.as_ref() else {
                    ui.label("Alignment view: load a query and subject and click on the plot.");
                    return;
                };
                let Some((cq_seq, cs_seq)) = self.crosshair else {
                    ui.label("Alignment view: click on the plot to set the crosshair.");
                    return;
                };
                let Some(q_seq) = self.query.as_ref() else {
                    return;
                };
                let Some(s_seq) = self.subject.as_ref() else {
                    return;
                };

                // Crosshair is already in absolute residue space.
                let q_centre = cq_seq as usize;
                let s_centre = cs_seq as usize;
                let _ = plot;

                // Header row: coords + window-size spinner.
                ui.horizontal(|ui| {
                    ui.label(format!(
                        "q = {}",
                        format_coord(self.query.as_ref(), q_centre)
                    ));
                    ui.separator();
                    ui.label(format!(
                        "s = {}",
                        format_coord(self.subject.as_ref(), s_centre)
                    ));
                    ui.separator();
                    ui.label("window:");
                    ui.add(
                        egui::DragValue::new(&mut self.settings.align_window_size)
                            .range(20..=400)
                            .speed(1.0),
                    );
                });

                // Annotation overlaps at the crosshair point (F3): every
                // visible feature covering the crosshair residue on each axis,
                // so nested/overlapping annotations are all listed.
                let q_hits = self
                    .query_annot
                    .as_ref()
                    .map(|ax| ax.features_at(q_centre))
                    .unwrap_or_default();
                let s_hits = self
                    .subject_annot
                    .as_ref()
                    .map(|ax| ax.features_at(s_centre))
                    .unwrap_or_default();
                // Cell marker color: the feature color for a single hit, a
                // neutral grey when several overlap, none when there are none.
                let mark_color = |ax: Option<&crate::annotation_overlay::AxisAnnot>,
                                  hits: &[&crate::annotation_overlay::BoundFeature]|
                 -> Option<Color32> {
                    match (ax, hits.len()) {
                        (_, 0) => None,
                        (Some(ax), 1) => Some(ax.color_for(&ax.value_of(hits[0]))),
                        _ => Some(Color32::from_gray(60)),
                    }
                };
                let q_mark = mark_color(self.query_annot.as_ref(), &q_hits);
                let s_mark = mark_color(self.subject_annot.as_ref(), &s_hits);
                if self.query_annot.is_some() || self.subject_annot.is_some() {
                    draw_overlap_line(
                        ui,
                        "q",
                        self.query_annot.as_ref(),
                        &q_hits,
                        self.query.as_ref(),
                    );
                    draw_overlap_line(
                        ui,
                        "s",
                        self.subject_annot.as_ref(),
                        &s_hits,
                        self.subject.as_ref(),
                    );
                }

                let window = self.settings.align_window_size.clamp(20, 400) as usize;
                let half = window / 2;
                let q_bytes = q_seq.bytes();
                let s_bytes = s_seq.bytes();
                let matrix = self.settings.build_matrix();

                // Forward (+/+) alignment columns; BLASTX: the query
                // translated in its best-scoring frame.
                let forward_block = if self.settings.mode == BlastMode::Blastx {
                    build_blastx_columns(
                        q_bytes,
                        s_bytes,
                        q_centre,
                        s_centre,
                        half,
                        window,
                        matrix.as_ref(),
                    )
                } else {
                    build_align_columns(
                        q_bytes,
                        s_bytes,
                        q_centre,
                        s_centre,
                        half,
                        window,
                        false,
                        matrix.as_ref(),
                        self.settings.mode,
                    )
                };

                // Reverse (+/-) alignment columns. BLASTN only —
                // proteins have no reverse strand. Compares
                // query[q_c + i - half] against the complement of
                // subject[s_c + half - i] (i.e. walks the subject
                // backwards while complementing).
                let reverse_block = if self.settings.mode == BlastMode::Blastn {
                    Some(build_align_columns(
                        q_bytes,
                        s_bytes,
                        q_centre,
                        s_centre,
                        half,
                        window,
                        true,
                        matrix.as_ref(),
                        self.settings.mode,
                    ))
                } else {
                    None
                };

                // Header continuation: Copy button.
                ui.horizontal(|ui| {
                    if ui
                        .button("Copy")
                        .on_hover_text("Copy the alignment block to clipboard as plain text")
                        .clicked()
                    {
                        let text = format_alignment_clipboard(
                            &forward_block,
                            reverse_block.as_ref(),
                            q_centre,
                            s_centre,
                            self.query.as_ref(),
                            self.subject.as_ref(),
                        );
                        ctx.copy_text(text);
                    }
                });

                draw_align_block(ui, ctx, &forward_block, q_mark, s_mark);
                if let Some(rev) = &reverse_block {
                    ui.add_space(4.0);
                    draw_align_block(ui, ctx, rev, q_mark, s_mark);
                }
            });
    }

    fn draw_canvas(&mut self, ctx: &Context) {
        egui::CentralPanel::default().show(ctx, |ui| {
            self.ensure_texture(ctx);

            // Up-front canvas measurement, before any early returns
            // for the no-plot case — the display-matched compute
            // depends on this and it's required even *before* the
            // first plot exists (deferred-initial-compute path).
            let avail = ui.available_size();
            let multi_q = self
                .query
                .as_ref()
                .map(|q| q.records.len() > 1)
                .unwrap_or(false);
            let multi_s = self
                .subject
                .as_ref()
                .map(|s| s.records.len() > 1)
                .unwrap_or(false);
            let track = self.axis_track_width();
            let top_margin: f32 = if multi_q { 44.0 } else { 24.0 } + track;
            let left_margin: f32 = if multi_s { 70.0 } else { 56.0 } + track;
            let approx_plot_w = (avail.x - left_margin).max(1.0);
            let approx_plot_h = (avail.y - top_margin).max(1.0);
            let new_area = (approx_plot_w, approx_plot_h);
            // Detect resize-retarget BEFORE overwriting the previous
            // measurement: compare against the canvas size at last
            // successful compute, not the live frame-to-frame
            // measurement (continuous drag would otherwise re-trigger
            // every frame and never get below the threshold).
            if self.settings.auto_fit_compute_zoom {
                if let Some((prev_w, prev_h)) = self.last_compute_canvas_size {
                    let dw = (new_area.0 - prev_w).abs() / prev_w.max(1.0);
                    let dh = (new_area.1 - prev_h).abs() / prev_h.max(1.0);
                    if dw > RESIZE_RETARGET_THRESHOLD || dh > RESIZE_RETARGET_THRESHOLD {
                        self.pending_resize_retarget = Some(std::time::Instant::now());
                    }
                }
            }
            self.measured_plot_area = Some(new_area);

            // Allocate the canvas and resolve interactions *before*
            // binding `&self.plot` — the borrow checker can't tell
            // that the subsequent `self.action_*` mutations don't
            // touch the same field, so we capture everything we need
            // from `response` up front. The `plot` reference's
            // lifetime then doesn't overlap any self-mutating method.
            let (rect, response) = ui.allocate_exact_size(avail, Sense::click_and_drag());
            let plot_area = Rect::from_min_max(
                Pos2::new(rect.left() + left_margin, rect.top() + top_margin),
                rect.right_bottom(),
            );

            // Middle-drag = rectangle-zoom selection (dotter classic).
            // Press: start rect_select, recording the screen position.
            // Drag: update current corner. Release: stash the
            // completed selection in `pending_rect_zoom`; the action
            // fires after the paint closure ends.
            if response.drag_started_by(egui::PointerButton::Middle) {
                if let Some(p) = response.interact_pointer_pos() {
                    if plot_area.contains(p) {
                        self.rect_select = Some(RectSelect {
                            start_screen: p,
                            current_screen: p,
                        });
                    }
                }
            }
            if response.dragged_by(egui::PointerButton::Middle) {
                if let (Some(p), Some(rs)) =
                    (response.interact_pointer_pos(), self.rect_select.as_mut())
                {
                    rs.current_screen = p;
                }
            }
            if response.drag_stopped_by(egui::PointerButton::Middle) {
                if let Some(rs) = self.rect_select.take() {
                    self.pending_rect_zoom = Some((plot_area, rs));
                }
            }

            // Right-drag = selection with an action menu (zoom, copy
            // coordinates, open in a new window).
            if response.drag_started_by(egui::PointerButton::Secondary) {
                if let Some(p) = response
                    .interact_pointer_pos()
                    .filter(|p| plot_area.contains(*p))
                {
                    self.selection = None;
                    self.selection_drag = Some(RectSelect {
                        start_screen: p,
                        current_screen: p,
                    });
                }
            }
            if response.dragged_by(egui::PointerButton::Secondary) {
                if let (Some(p), Some(rs)) = (
                    response.interact_pointer_pos(),
                    self.selection_drag.as_mut(),
                ) {
                    rs.current_screen = p;
                }
            }
            if response.drag_stopped_by(egui::PointerButton::Secondary) {
                if let Some(rs) = self.selection_drag.take() {
                    self.selection =
                        self.selection_residues(plot_area, &rs)
                            .map(|(q, s)| Selection {
                                q,
                                s,
                                plot_area,
                                menu_at: rs.current_screen,
                            });
                }
            }
            if response.clicked_by(egui::PointerButton::Primary) {
                self.selection = None;
            }

            // Left-drag = pan. Plain drag, applied directly to
            // view_offset (1:1 means drag-delta is pixelmap-pixel
            // delta).
            if response.dragged_by(egui::PointerButton::Primary) {
                self.view_offset -= response.drag_delta();
            }

            let Some(plot) = &self.plot else {
                ui.centered_and_justified(|ui| {
                    ui.label("No plot. Load a query and subject FASTA (File menu).");
                });
                return;
            };
            let Some(tex) = &self.texture else {
                return;
            };

            let ppp = self.pixels_per_point;
            let pw = plot.width as f32;
            let ph = plot.height as f32;
            // Image extent in *logical points* — texture is
            // `pw × ph` texels, drawn at `pw/ppp × ph/ppp` logical
            // points so one texel lands on exactly one physical
            // pixel (Option C in
            // `docs/reviews/dotter-sizing-model.md`).
            let pw_lp = pw / ppp;
            let ph_lp = ph / ppp;
            let plot_w = plot_area.width().max(1.0);
            let plot_h = plot_area.height().max(1.0);

            // Clamp pan to keep the pixelmap on screen. `view_offset`
            // is in logical points; the image is `pw_lp × ph_lp`
            // logical points. If the image is smaller than the
            // canvas, centre it; otherwise clamp the offset to the
            // overhang.
            if pw_lp <= plot_w {
                self.view_offset.x = -(plot_w - pw_lp) / 2.0;
            } else {
                self.view_offset.x = self.view_offset.x.clamp(0.0, pw_lp - plot_w);
            }
            if ph_lp <= plot_h {
                self.view_offset.y = -(plot_h - ph_lp) / 2.0;
            } else {
                self.view_offset.y = self.view_offset.y.clamp(0.0, ph_lp - plot_h);
            }

            // Snap pan to whole-texel boundaries so the texture's
            // source texels keep landing on physical-pixel grid
            // boundaries (otherwise fractional offsets re-introduce
            // sub-pixel sampling phase). 1 texel = 1/ppp logical
            // points → snap to the nearest 1/ppp multiple.
            let snap_texel = |v: f32| (v * ppp).round() / ppp;
            let draw_offset = Vec2::new(
                snap_texel(self.view_offset.x),
                snap_texel(self.view_offset.y),
            );

            // Single click (primary, no drag): set crosshair. Clicks
            // in the margin are ignored. Use `draw_offset` so the
            // crosshair lands on the pixmap cell the user *sees*
            // under the cursor, not at a sub-pixel-shifted cell.
            // The pixmap pixel is then converted to an *absolute
            // residue* by adding the slice origin and multiplying by
            // the compute zoom (which is the residues-per-pixel
            // factor) — keeps the crosshair semantics zoom-independent.
            if response.clicked() {
                if let Some(p) = response.interact_pointer_pos() {
                    if plot_area.contains(p) {
                        let local = p - plot_area.left_top();
                        // local is in logical points; convert to
                        // texel index by multiplying by ppp (one
                        // logical pt covers `ppp` texels in this
                        // render).
                        let qp = ((local.x + draw_offset.x) * ppp).floor() as i64;
                        let sp = ((local.y + draw_offset.y) * ppp).floor() as i64;
                        if qp >= 0 && qp < plot.width as i64 && sp >= 0 && sp < plot.height as i64 {
                            let z = plot.params.zoom.max(1) as i64;
                            let (q_off, s_off) = self
                                .current_slice
                                .as_ref()
                                .map(|sl| (sl.q_range.start as i64, sl.s_range.start as i64))
                                .unwrap_or((0, 0));
                            // Centre of the pixel's residue block.
                            let zx = z * self.qx() as i64;
                            let q_seq = q_off + qp * zx + zx / 2;
                            let s_seq = s_off + sp * z + z / 2;
                            let qlen = self
                                .query
                                .as_ref()
                                .map(|s| s.len() as i64)
                                .unwrap_or(i64::MAX);
                            let slen = self
                                .subject
                                .as_ref()
                                .map(|s| s.len() as i64)
                                .unwrap_or(i64::MAX);
                            let q_seq = q_seq.clamp(0, qlen.saturating_sub(1).max(0));
                            let s_seq = s_seq.clamp(0, slen.saturating_sub(1).max(0));
                            self.crosshair = Some((q_seq as u32, s_seq as u32));
                        }
                    }
                }
            }

            // Fill the whole canvas (margins + plot area) with light
            // grey first. The margin band is left as-is; the plot
            // area gets the texture painted over it.
            ui.painter().rect_filled(rect, 0.0, Color32::from_gray(235));

            // Render the pixelmap so one texture texel lands on one
            // physical pixel. Texture is `pw × ph` texels; draw it
            // at `pw/ppp × ph/ppp` logical points. Origin snapped
            // to a physical pixel boundary; size left as the exact
            // logical extent of the texture (snapping the far edge
            // could shift it by ±1 physical pixel and re-introduce
            // sampling phase).
            let plot_screen_x = plot_area.left() - draw_offset.x;
            let plot_screen_y = plot_area.top() - draw_offset.y;
            let snap = |v: f32| ((v * ppp).round()) / ppp;
            let plot_rect = Rect::from_min_size(
                Pos2::new(snap(plot_screen_x), snap(plot_screen_y)),
                Vec2::new(pw_lp, ph_lp),
            );
            let clip_painter = ui.painter_at(plot_area);
            clip_painter.image(
                tex.id(),
                plot_rect,
                Rect::from_min_max(Pos2::ZERO, Pos2::new(1.0, 1.0)),
                Color32::WHITE,
            );
            // Boundary frame around the visible plot area.
            ui.painter().rect_stroke(
                plot_area,
                0.0,
                egui::Stroke::new(1.0, Color32::from_gray(110)),
            );

            // Rubber-band rectangle while middle-mouse drag is in
            // progress. Drawn over the pixelmap but under the
            // crosshair / labels.
            if let Some(rs) = self.rect_select {
                let r = Rect::from_two_pos(rs.start_screen, rs.current_screen).intersect(plot_area);
                if r.width() > 0.5 && r.height() > 0.5 {
                    clip_painter.rect_filled(
                        r,
                        0.0,
                        Color32::from_rgba_unmultiplied(40, 100, 200, 24),
                    );
                    clip_painter.rect_stroke(
                        r,
                        0.0,
                        egui::Stroke::new(1.0, Color32::from_rgb(40, 100, 200)),
                    );
                }
            }

            // Right-drag selection: dashed-free outline in orange, while
            // dragging and while its menu is open.
            let sel_rect = match (&self.selection_drag, &self.selection) {
                (Some(rs), _) => Some(Rect::from_two_pos(rs.start_screen, rs.current_screen)),
                (None, Some(sel)) => {
                    let z = plot.params.zoom.max(1) as f32;
                    let (q_off, s_off) = self
                        .current_slice
                        .as_ref()
                        .map_or((0, 0), |sl| (sl.q_range.start, sl.s_range.start));
                    let zx = z * self.qx() as f32;
                    let at = |q: usize, s: usize| {
                        Pos2::new(
                            plot_area.left() + (q as f32 - q_off as f32) / zx / ppp - draw_offset.x,
                            plot_area.top() + (s as f32 - s_off as f32) / z / ppp - draw_offset.y,
                        )
                    };
                    Some(Rect::from_two_pos(
                        at(sel.q.start, sel.s.start),
                        at(sel.q.end, sel.s.end),
                    ))
                }
                _ => None,
            };
            if let Some(r) = sel_rect {
                clip_painter.rect_filled(r, 0.0, Color32::from_rgba_unmultiplied(230, 120, 0, 28));
                clip_painter.rect_stroke(
                    r,
                    0.0,
                    egui::Stroke::new(1.5, Color32::from_rgb(230, 120, 0)),
                );
            }

            // Annotation bands (ADR 0005). Query features paint as
            // vertical bands spanning the plot height; subject features
            // as horizontal bands spanning the width. Drawn over the
            // pixelmap with a single global alpha so dots show through
            // and overlapping bands self-darken via compositing. Placed
            // under the breaklines / ridges / crosshair so those thin
            // marks stay legible on top.
            if self.annot_show {
                let azoom = plot.params.zoom.max(1) as usize;
                let (aq_off, as_off) = self
                    .current_slice
                    .as_ref()
                    .map(|sl| (sl.q_range.start, sl.s_range.start))
                    .unwrap_or((0, 0));
                let aa = (self.annot_alpha.clamp(0.0, 1.0) * 255.0).round() as u8;
                let with_alpha = |base: Color32| {
                    Color32::from_rgba_unmultiplied(base.r(), base.g(), base.b(), aa)
                };
                // Query → vertical bands.
                if let Some(ax) = &self.query_annot {
                    for bf in &ax.bound {
                        let value = ax.value_of(bf);
                        if ax.hidden.contains(&value) {
                            continue;
                        }
                        if let Some((x0, x1)) = crate::annotation_overlay::band_screen_span(
                            bf.range.start,
                            bf.range.end,
                            aq_off,
                            azoom * self.qx(),
                            plot_area.left(),
                            draw_offset.x,
                            ppp,
                            plot_area.left(),
                            plot_area.right(),
                        ) {
                            clip_painter.rect_filled(
                                Rect::from_min_max(
                                    Pos2::new(x0, plot_area.top()),
                                    Pos2::new(x1, plot_area.bottom()),
                                ),
                                0.0,
                                with_alpha(ax.color_for(&value)),
                            );
                        }
                    }
                }
                // Subject → horizontal bands.
                if let Some(ax) = &self.subject_annot {
                    for bf in &ax.bound {
                        let value = ax.value_of(bf);
                        if ax.hidden.contains(&value) {
                            continue;
                        }
                        if let Some((y0, y1)) = crate::annotation_overlay::band_screen_span(
                            bf.range.start,
                            bf.range.end,
                            as_off,
                            azoom,
                            plot_area.top(),
                            draw_offset.y,
                            ppp,
                            plot_area.top(),
                            plot_area.bottom(),
                        ) {
                            clip_painter.rect_filled(
                                Rect::from_min_max(
                                    Pos2::new(plot_area.left(), y0),
                                    Pos2::new(plot_area.right(), y1),
                                ),
                                0.0,
                                with_alpha(ax.color_for(&value)),
                            );
                        }
                    }
                }
            }

            // C3: breaklines for multi-record FASTA inputs. Vertical
            // lines at the query record boundaries; horizontal lines
            // at the subject record boundaries. Drawn underneath the
            // crosshair so it stays visible. The pixelmap is sliced
            // to `current_slice`, so a break at absolute residue B
            // shows at pixmap pixel (B − slice.start) / zoom (only
            // when inside the slice).
            let zoom_us = plot.params.zoom.max(1) as usize;
            let break_stroke =
                egui::Stroke::new(1.0, Color32::from_rgba_unmultiplied(80, 160, 80, 220));
            let (q_off, s_off) = self
                .current_slice
                .as_ref()
                .map(|sl| (sl.q_range.start, sl.s_range.start))
                .unwrap_or((0, 0));
            if let Some(q_seq) = &self.query {
                for &break_coord in &q_seq.breaks() {
                    if break_coord < q_off {
                        continue;
                    }
                    let pixel_x = (break_coord - q_off) / (zoom_us * self.qx());
                    if pixel_x >= plot.width as usize {
                        continue;
                    }
                    let sx = plot_area.left() + (pixel_x as f32) / ppp - draw_offset.x;
                    if sx < plot_area.left() || sx > plot_area.right() {
                        continue;
                    }
                    clip_painter.line_segment(
                        [
                            Pos2::new(sx, plot_area.top()),
                            Pos2::new(sx, plot_area.bottom()),
                        ],
                        break_stroke,
                    );
                }
            }
            if let Some(s_seq) = &self.subject {
                for &break_coord in &s_seq.breaks() {
                    if break_coord < s_off {
                        continue;
                    }
                    let pixel_y = (break_coord - s_off) / zoom_us;
                    if pixel_y >= plot.height as usize {
                        continue;
                    }
                    let sy = plot_area.top() + (pixel_y as f32) / ppp - draw_offset.y;
                    if sy < plot_area.top() || sy > plot_area.bottom() {
                        continue;
                    }
                    clip_painter.line_segment(
                        [
                            Pos2::new(plot_area.left(), sy),
                            Pos2::new(plot_area.right(), sy),
                        ],
                        break_stroke,
                    );
                }
            }

            // Vector ridge overlay: anti-aliased line segments over
            // coherent diagonal runs detected in the raster. Hides
            // the per-window intensity oscillation on imperfect-
            // homology lines without altering the pixmap. Forward
            // ridges in dark grey, reverse ridges in magenta to
            // mirror the inverted-repeat colour convention.
            if self.settings.show_ridge_overlay && !self.current_ridges.is_empty() {
                let stroke_fwd = egui::Stroke::new(0.75, Color32::from_rgb(20, 20, 20));
                let stroke_rev = egui::Stroke::new(0.75, Color32::from_rgb(170, 0, 170));
                for ridge in &self.current_ridges {
                    let p0 = Pos2::new(
                        plot_area.left() + (ridge.start.0 as f32 + 0.5) / ppp - draw_offset.x,
                        plot_area.top() + (ridge.start.1 as f32 + 0.5) / ppp - draw_offset.y,
                    );
                    let p1 = Pos2::new(
                        plot_area.left() + (ridge.end.0 as f32 + 0.5) / ppp - draw_offset.x,
                        plot_area.top() + (ridge.end.1 as f32 + 0.5) / ppp - draw_offset.y,
                    );
                    let stroke = match ridge.direction {
                        RidgeDirection::Forward => stroke_fwd,
                        RidgeDirection::Reverse => stroke_rev,
                    };
                    clip_painter.line_segment([p0, p1], stroke);
                }
            }

            // Large mode: synteny chains as vector lines, forward dark
            // green, reverse magenta. Coordinates are absolute residues.
            if let (true, Some(lp)) = (self.settings.large.show_chains, &self.large_plot) {
                let min = self.settings.large.chain_min_score as f64;
                let z = zoom_us as f32;
                let to_screen = |q: u32, s: u32| {
                    Pos2::new(
                        plot_area.left() + (q as f32 - q_off as f32) / z / ppp - draw_offset.x,
                        plot_area.top() + (s as f32 - s_off as f32) / z / ppp - draw_offset.y,
                    )
                };
                let fwd = egui::Stroke::new(2.0, Color32::from_rgba_unmultiplied(0, 140, 60, 200));
                let rev = egui::Stroke::new(2.0, Color32::from_rgba_unmultiplied(190, 0, 190, 200));
                for c in lp.chains.iter().filter(|c| c.score >= min) {
                    let (p0, p1, stroke) = if c.reverse {
                        (
                            to_screen(c.q_start, c.s_end),
                            to_screen(c.q_end, c.s_start),
                            rev,
                        )
                    } else {
                        (
                            to_screen(c.q_start, c.s_start),
                            to_screen(c.q_end, c.s_end),
                            fwd,
                        )
                    };
                    clip_painter.line_segment([p0, p1], stroke);
                    if lp.self_comparison && lp.origin_q == lp.origin_s {
                        // Mirror across the diagonal like the raster.
                        let m = |p: Pos2| {
                            let q =
                                (p.x + draw_offset.x - plot_area.left()) * ppp * z + q_off as f32;
                            let s =
                                (p.y + draw_offset.y - plot_area.top()) * ppp * z + s_off as f32;
                            to_screen(s as u32, q as u32)
                        };
                        clip_painter.line_segment([m(p0), m(p1)], stroke);
                    }
                }
            }

            // HSP overlay (line modes) as vectors; greyscale is in the
            // texture.
            if let (Some(h), true) = (
                &self.hsps,
                matches!(self.hsp_mode, HspMode::Line | HspMode::Score),
            ) {
                let z = plot.params.zoom.max(1) as f32;
                let (q_off, s_off) = self
                    .current_slice
                    .as_ref()
                    .map_or((0, 0), |sl| (sl.q_range.start, sl.s_range.start));
                let zx = z * self.qx() as f32;
                let to_screen = |x: f64, y: f64| {
                    Pos2::new(
                        plot_area.left() + (x as f32 - q_off as f32) / zx / ppp - draw_offset.x,
                        plot_area.top() + (y as f32 - s_off as f32) / z / ppp - draw_offset.y,
                    )
                };
                for hsp in &h.hsps {
                    let [r, g, b] = dottir_io::hsp::colour(self.hsp_mode, hsp.identity);
                    clip_painter.line_segment(
                        [to_screen(hsp.x0, hsp.y0), to_screen(hsp.x1, hsp.y1)],
                        egui::Stroke::new(1.5, Color32::from_rgb(r, g, b)),
                    );
                }
            }

            // C4: tick labels in the top / left margin bands —
            // outside the plot area so they don't overlap the image.
            self.draw_axis_labels(ui, rect, plot_area, plot, draw_offset, ppp);
            if let Some(tip) =
                self.draw_axis_tracks(ui, plot_area, draw_offset, ppp, response.hover_pos())
            {
                response.clone().on_hover_text_at_pointer(tip);
            }

            // Crosshair overlay + coord label — clipped to the plot
            // area so the lines never run into the axis margin. The
            // crosshair is in *absolute residue* coords; map to the
            // current slice's pixmap-pixel space (slice origin + zoom)
            // before placing on screen. Crosshairs outside the slice
            // are skipped — they'd just paint at the slice edge
            // misleadingly.
            if let Some((cq_seq, cs_seq)) = self.crosshair {
                let z = plot.params.zoom.max(1) as f32;
                let (q_off, s_off) = self
                    .current_slice
                    .as_ref()
                    .map(|sl| (sl.q_range.start as f32, sl.s_range.start as f32))
                    .unwrap_or((0.0, 0.0));
                let cq_pix = ((cq_seq as f32) - q_off) / (z * self.qx() as f32);
                let cs_pix = ((cs_seq as f32) - s_off) / z;
                // Only draw when the crosshair pixel is inside the
                // current slice's pixmap; otherwise (zoomed into a
                // region that doesn't contain the residue) hide it.
                let in_slice = cq_pix >= 0.0
                    && cq_pix < plot.width as f32
                    && cs_pix >= 0.0
                    && cs_pix < plot.height as f32;
                if in_slice {
                    let cx = plot_area.left() + (cq_pix + 0.5) / ppp - draw_offset.x;
                    let cy = plot_area.top() + (cs_pix + 0.5) / ppp - draw_offset.y;
                    let stroke = egui::Stroke::new(1.0, Color32::from_rgb(255, 80, 80));
                    clip_painter.line_segment(
                        [
                            Pos2::new(plot_area.left(), cy),
                            Pos2::new(plot_area.right(), cy),
                        ],
                        stroke,
                    );
                    clip_painter.line_segment(
                        [
                            Pos2::new(cx, plot_area.top()),
                            Pos2::new(cx, plot_area.bottom()),
                        ],
                        stroke,
                    );

                    // Coord label next to the cross — crosshair is
                    // already in absolute residue coords.
                    let q_seq = cq_seq as usize;
                    let s_seq = cs_seq as usize;
                    let label = format!(
                        "q = {}, s = {}",
                        format_coord(self.query.as_ref(), q_seq),
                        format_coord(self.subject.as_ref(), s_seq),
                    );
                    let font = egui::FontId::monospace(11.0);
                    let label_size = clip_painter
                        .layout_no_wrap(label.clone(), font.clone(), Color32::BLACK)
                        .size();
                    let pad = 4.0;
                    let mut lx = cx + 6.0;
                    let mut ly = cy + 6.0;
                    if lx + label_size.x + pad > plot_area.right() {
                        lx = cx - 6.0 - label_size.x - 2.0 * pad;
                    }
                    if ly + label_size.y + pad > plot_area.bottom() {
                        ly = cy - 6.0 - label_size.y - 2.0 * pad;
                    }
                    // 1-px white shadow in the four cardinal directions
                    // keeps the dark red label legible over a black
                    // diagonal without occluding dotplot pixels with a
                    // coloured patch.
                    let label_pos = Pos2::new(lx, ly);
                    let shadow = Color32::WHITE;
                    for &(dx, dy) in &[(-1.0, 0.0), (1.0, 0.0), (0.0, -1.0), (0.0, 1.0)] {
                        clip_painter.text(
                            Pos2::new(label_pos.x + dx, label_pos.y + dy),
                            egui::Align2::LEFT_TOP,
                            &label,
                            font.clone(),
                            shadow,
                        );
                    }
                    clip_painter.text(
                        label_pos,
                        egui::Align2::LEFT_TOP,
                        label,
                        font,
                        Color32::from_rgb(120, 0, 0),
                    );
                }
            }
        });
    }

    /// Draw tick marks and coordinate labels in *sequence* coords
    /// along the top (query) and left (subject) edges of the canvas.
    ///
    /// Tick spacing is picked so each adjacent label is at least
    /// `MIN_LABEL_SPACING_PX` apart in screen space, and the value
    /// step is a "nice" power-of-10 multiple (`1, 2, 5 × 10^k`).
    /// Large mode: strips just outside the top (query) and left (subject)
    /// plot edges. Each pixel column/row shows the matching bp of its
    /// visible dots per base of sequence ("aligned copies per base", log
    /// scale) as a bar growing away from the plot, split by the dots'
    /// repeat class, plus a red mark at the outer edge where seeds were
    /// capped. Returns the tooltip for a hovered bin.
    fn draw_axis_tracks(
        &self,
        ui: &egui::Ui,
        plot_area: Rect,
        draw_offset: Vec2,
        ppp: f32,
        hover: Option<Pos2>,
    ) -> Option<String> {
        if self.axis_track_width() == 0.0 {
            return None;
        }
        let lp = self.shown_large_plot()?;
        let prof = self.axis_profile.as_ref()?;
        let painter = ui.painter();
        let zoom = lp.zoom as usize;
        let bar = AXIS_TRACK - 5.0;
        let texel = 1.0 / ppp;
        let mut tip = None;

        // `along` = position along the axis of bin i; `away` = distance
        // from the plot edge. Returns the screen rect for (along, away).
        for query_axis in [true, false] {
            let (bins, capped, origin, profile) = if query_axis {
                (&lp.q_depth_bins, &lp.q_capped, lp.origin_q, &prof.q)
            } else {
                (&lp.s_depth_bins, &lp.s_capped, lp.origin_s, &prof.s)
            };
            let (edge, lo, hi, off) = if query_axis {
                (
                    plot_area.top(),
                    plot_area.left(),
                    plot_area.right(),
                    draw_offset.x,
                )
            } else {
                (
                    plot_area.left(),
                    plot_area.top(),
                    plot_area.bottom(),
                    draw_offset.y,
                )
            };
            let cell = |i: usize, from: f32, to: f32| {
                let a0 = lo + i as f32 * texel - off;
                let (a0, a1) = (a0.max(lo), (a0 + texel).min(hi));
                let (d0, d1) = (edge - to, edge - from);
                if query_axis {
                    Rect::from_min_max(Pos2::new(a0, d0), Pos2::new(a1, d1))
                } else {
                    Rect::from_min_max(Pos2::new(d0, a0), Pos2::new(d1, a1))
                }
            };
            let first = ((off * ppp).floor().max(0.0)) as usize;
            let last = (((off + (hi - lo)) * ppp).ceil() as usize).min(profile.len());
            for (i, bp) in profile.iter().enumerate().take(last).skip(first) {
                // Log height: 1 copy/base ≈ 10 %, AXIS_TRACK_MAX_COPIES full.
                let copies = bp.iter().sum::<f32>() / zoom as f32;
                let total_h =
                    ((1.0 + copies).ln() / (1.0 + AXIS_TRACK_MAX_COPIES).ln()).min(1.0) * bar;
                let total_bp: f32 = bp.iter().sum();
                let mut from = 1.0;
                for c in RepeatClass::ALL {
                    let h = if total_bp > 0.0 {
                        bp[c.index()] / total_bp * total_h
                    } else {
                        0.0
                    };
                    if h > 0.0 {
                        let colour = if c == RepeatClass::Unique {
                            Color32::from_gray(90)
                        } else {
                            CLASS_COLOURS[c.index()]
                        };
                        painter.rect_filled(cell(i, from, from + h), 0.0, colour);
                        from += h;
                    }
                }
                if let Some(&n) = capped.get(i).filter(|&&n| n > 0) {
                    let a = (n as f32 / zoom as f32 * 4.0).clamp(0.3, 1.0);
                    painter.rect_filled(
                        cell(i, AXIS_TRACK - 3.0, AXIS_TRACK),
                        0.0,
                        Color32::from_rgba_unmultiplied(220, 0, 0, (a * 255.0) as u8),
                    );
                }
            }
            if let Some(p) = hover {
                let band = cell(0, 0.0, AXIS_TRACK);
                let (along, inside) = if query_axis {
                    (
                        p.x,
                        p.y >= band.top() && p.y < edge && p.x >= lo && p.x < hi,
                    )
                } else {
                    (
                        p.y,
                        p.x >= band.left() && p.x < edge && p.y >= lo && p.y < hi,
                    )
                };
                let i = ((along - lo + off) * ppp).floor() as usize;
                if let (true, Some(b), Some(bp)) = (inside, bins.get(i), profile.get(i)) {
                    tip = Some(depth_bin_text(
                        if query_axis { "Query" } else { "Subject" },
                        origin + i * zoom,
                        zoom,
                        bp,
                        b,
                        capped.get(i).copied().unwrap_or(0),
                    ));
                }
            }
        }
        tip
    }

    fn draw_axis_labels(
        &self,
        ui: &mut egui::Ui,
        outer: Rect,
        plot_area: Rect,
        plot: &dottir_core::DotPlot,
        draw_offset: Vec2,
        ppp: f32,
    ) {
        const MIN_LABEL_SPACING_PX: f32 = 80.0;
        let zoom_us = plot.params.zoom.max(1) as f32;
        // Slice origin offset — labels show absolute residues, not
        // slice-local positions.
        let (q_off, s_off) = self
            .current_slice
            .as_ref()
            .map(|sl| (sl.q_range.start as u64, sl.s_range.start as u64))
            .unwrap_or((0, 0));
        // World-pixel (texel) range visible inside the plot area.
        // `draw_offset` is in logical points; one logical point
        // covers `ppp` texels, so multiply the logical extent by
        // `ppp` to get texel extent.
        let world_x_lo = draw_offset.x * ppp;
        let world_x_hi = world_x_lo + plot_area.width() * ppp;
        let world_y_lo = draw_offset.y * ppp;
        let world_y_hi = world_y_lo + plot_area.height() * ppp;
        // Convert to absolute sequence-residue range.
        // BLASTX: one column per codon, labels in bases.
        let zoom_x = zoom_us * self.qx() as f32;
        let seq_q_lo = q_off + (world_x_lo * zoom_x).max(0.0) as u64;
        let seq_q_hi = q_off + (world_x_hi * zoom_x) as u64;
        let seq_s_lo = s_off + (world_y_lo * zoom_us).max(0.0) as u64;
        let seq_s_hi = s_off + (world_y_hi * zoom_us) as u64;

        let painter = ui.painter();
        let tick_color = Color32::from_rgba_unmultiplied(80, 80, 80, 220);
        let label_color = ui.style().visuals.text_color();
        let font = egui::FontId::monospace(11.0);

        // Top axis (query) — labels above the plot, in the top
        // margin band (outer.top()..plot_area.top()).
        let top_margin_y = outer.top();
        let track = self.axis_track_width();
        let tick_baseline_y = plot_area.top() - track;
        let label_y = tick_baseline_y - 16.0; // text top
        let span_x = seq_q_hi.saturating_sub(seq_q_lo) as f32;
        // Pixels-per-residue is now expressed against *logical*
        // points: residue → texel via `/ zoom_us`, then texel →
        // logical via `/ ppp`.
        let pixels_per_residue_x = 1.0_f32 / (zoom_x * ppp);
        let step_x = nice_tick_step(span_x as f64, MIN_LABEL_SPACING_PX / pixels_per_residue_x);
        // Minor ticks (unlabeled) at step/5 intervals — gives 4
        // subdivisions between every labelled tick. Drawn first so
        // the longer labelled ticks paint on top.
        let minor_step_x = (step_x / 5.0).floor() as u64;
        if minor_step_x >= 1 {
            let mut t = (seq_q_lo / minor_step_x) * minor_step_x;
            while t < seq_q_hi.saturating_add(minor_step_x) {
                if t >= seq_q_lo && t <= seq_q_hi {
                    let sx = plot_area.left() + (t - q_off) as f32 / (zoom_x * ppp) - draw_offset.x;
                    if sx >= plot_area.left() - 1.0 && sx <= plot_area.right() + 1.0 {
                        painter.line_segment(
                            [
                                Pos2::new(sx, tick_baseline_y - 3.0),
                                Pos2::new(sx, tick_baseline_y),
                            ],
                            egui::Stroke::new(1.0, tick_color),
                        );
                    }
                }
                t = t.saturating_add(minor_step_x);
            }
        }
        let mut t = (seq_q_lo / step_x as u64) * step_x as u64;
        while t < seq_q_hi.saturating_add(step_x as u64) {
            if t >= seq_q_lo && t <= seq_q_hi {
                let sx = plot_area.left() + (t - q_off) as f32 / (zoom_x * ppp) - draw_offset.x;
                if sx < plot_area.left() - 1.0 || sx > plot_area.right() + 1.0 {
                    t = t.saturating_add(step_x as u64);
                    if step_x == 0.0 {
                        break;
                    }
                    continue;
                }
                // Tick: short line just above the plot's top edge.
                painter.line_segment(
                    [
                        Pos2::new(sx, tick_baseline_y - 5.0),
                        Pos2::new(sx, tick_baseline_y),
                    ],
                    egui::Stroke::new(1.0, tick_color),
                );
                painter.text(
                    Pos2::new(sx, label_y.max(top_margin_y + 2.0)),
                    egui::Align2::CENTER_TOP,
                    dottir_io::text_overlay::format_kb(t),
                    font.clone(),
                    label_color,
                );
            }
            t = t.saturating_add(step_x as u64);
            if step_x == 0.0 {
                break;
            }
        }

        // Left axis (subject) — labels in the left margin
        // (outer.left()..plot_area.left()).
        let tick_baseline_x = plot_area.left() - track;
        let label_x = tick_baseline_x - 8.0; // text right edge
        let span_y = seq_s_hi.saturating_sub(seq_s_lo) as f32;
        let step_y = nice_tick_step(span_y as f64, MIN_LABEL_SPACING_PX / pixels_per_residue_x);
        let minor_step_y = (step_y / 5.0).floor() as u64;
        if minor_step_y >= 1 {
            let mut t = (seq_s_lo / minor_step_y) * minor_step_y;
            while t < seq_s_hi.saturating_add(minor_step_y) {
                if t >= seq_s_lo && t <= seq_s_hi {
                    let sy = plot_area.top() + (t - s_off) as f32 / (zoom_us * ppp) - draw_offset.y;
                    if sy >= plot_area.top() - 1.0 && sy <= plot_area.bottom() + 1.0 {
                        painter.line_segment(
                            [
                                Pos2::new(tick_baseline_x - 3.0, sy),
                                Pos2::new(tick_baseline_x, sy),
                            ],
                            egui::Stroke::new(1.0, tick_color),
                        );
                    }
                }
                t = t.saturating_add(minor_step_y);
            }
        }
        let mut t = (seq_s_lo / step_y as u64) * step_y as u64;
        while t < seq_s_hi.saturating_add(step_y as u64) {
            if t >= seq_s_lo && t <= seq_s_hi {
                let sy = plot_area.top() + (t - s_off) as f32 / (zoom_us * ppp) - draw_offset.y;
                if sy < plot_area.top() - 1.0 || sy > plot_area.bottom() + 1.0 {
                    t = t.saturating_add(step_y as u64);
                    if step_y == 0.0 {
                        break;
                    }
                    continue;
                }
                painter.line_segment(
                    [
                        Pos2::new(tick_baseline_x - 5.0, sy),
                        Pos2::new(tick_baseline_x, sy),
                    ],
                    egui::Stroke::new(1.0, tick_color),
                );
                painter.text(
                    Pos2::new(label_x, sy),
                    egui::Align2::RIGHT_CENTER,
                    dottir_io::text_overlay::format_kb(t),
                    font.clone(),
                    label_color,
                );
            }
            t = t.saturating_add(step_y as u64);
            if step_y == 0.0 {
                break;
            }
        }

        // H1: record-name labels for multi-record FASTAs. Sit one
        // row above the tick labels for the query, and one column to
        // the left of the tick labels for the subject. Single-record
        // inputs render no extra labels.
        let record_font = egui::FontId::monospace(11.0);
        let record_color = ui.style().visuals.text_color();
        if let Some(q) = self.query.as_ref() {
            if q.records.len() > 1 {
                let rec_y = outer.top() + 2.0;
                for rec in &q.records {
                    // Record ranges are absolute residue positions;
                    // translate into slice-local pixmap pixels.
                    if rec.range.end as u64 <= q_off {
                        continue;
                    }
                    let r_start = (rec.range.start as u64).saturating_sub(q_off);
                    let r_end = (rec.range.end as u64).saturating_sub(q_off);
                    if r_end <= r_start {
                        continue;
                    }
                    let x0 = plot_area.left() + r_start as f32 / (zoom_us * ppp) - draw_offset.x;
                    let x1 = plot_area.left() + r_end as f32 / (zoom_us * ppp) - draw_offset.x;
                    let span = (x1 - x0).max(0.0);
                    if span < 18.0 {
                        continue;
                    }
                    let cx = ((x0 + x1) * 0.5).clamp(plot_area.left(), plot_area.right());
                    painter.text(
                        Pos2::new(cx, rec_y),
                        egui::Align2::CENTER_TOP,
                        truncate_for_span(&rec.id, span),
                        record_font.clone(),
                        record_color,
                    );
                }
            }
        }
        if let Some(s) = self.subject.as_ref() {
            if s.records.len() > 1 {
                let rec_x = outer.left() + 2.0;
                for rec in &s.records {
                    if rec.range.end as u64 <= s_off {
                        continue;
                    }
                    let r_start = (rec.range.start as u64).saturating_sub(s_off);
                    let r_end = (rec.range.end as u64).saturating_sub(s_off);
                    if r_end <= r_start {
                        continue;
                    }
                    let y0 = plot_area.top() + r_start as f32 / (zoom_us * ppp) - draw_offset.y;
                    let y1 = plot_area.top() + r_end as f32 / (zoom_us * ppp) - draw_offset.y;
                    let span = (y1 - y0).max(0.0);
                    if span < 14.0 {
                        continue;
                    }
                    let cy = ((y0 + y1) * 0.5).clamp(plot_area.top(), plot_area.bottom());
                    painter.text(
                        Pos2::new(rec_x, cy),
                        egui::Align2::LEFT_CENTER,
                        truncate_for_span(&rec.id, (plot_area.left() - rec_x).max(40.0)),
                        record_font.clone(),
                        record_color,
                    );
                }
            }
        }
    }
}

/// Greyramp over coverage codes: white below `white_bp` matching bases per
/// pixel, black from `black_bp`, linear in log(bases) between.
/// Score greyramp, or the coverage ramp in matching bases when large-mode
/// pixels hold coverage codes.
fn display_lut(coverage: bool, greyramp: &Greyramp, white_bp: f32, black_bp: f32) -> [u8; 256] {
    if coverage {
        coverage_lut(white_bp, black_bp)
    } else {
        greyramp.lut()
    }
}

fn coverage_lut(white_bp: f32, black_bp: f32) -> [u8; 256] {
    let lo = (white_bp.max(1.0) as f64).ln();
    let hi = (black_bp.max(white_bp.max(1.0) + 1.0) as f64).ln();
    let mut lut = [255u8; 256];
    for (code, slot) in lut.iter_mut().enumerate().skip(1) {
        let b = coverage_bases(code as u8).max(1.0).ln();
        if b >= lo {
            let t = ((b - lo) / (hi - lo)).clamp(0.0, 1.0);
            *slot = (255.0 * (1.0 - t)).round() as u8;
        }
    }
    lut
}

fn depth_bin_text(
    axis: &str,
    start: usize,
    zoom: usize,
    shown_bp: &[f32; RepeatClass::COUNT],
    b: &DepthBin,
    capped: u32,
) -> String {
    let total: f32 = shown_bp.iter().sum();
    let pct = |c: RepeatClass| {
        if total > 0.0 {
            shown_bp[c.index()] / total * 100.0
        } else {
            0.0
        }
    };
    let max = if b.max == u16::MAX {
        "≥ max-occ".to_string()
    } else {
        b.max.to_string()
    };
    let mut t = format!(
        "{axis} {}–{}\nvisible dots: {} bp matching = {:.1} copies per base\n\
         of which unique {:.0}% · low {:.0}% · mid {:.0}% · high {:.0}%\n\
         copy number of these bases: median {}, max {max}",
        format_thousands(start),
        format_thousands(start + zoom),
        format_thousands(total.round() as usize),
        total / zoom as f32,
        pct(RepeatClass::Unique),
        pct(RepeatClass::Low),
        pct(RepeatClass::Mid),
        pct(RepeatClass::High),
        b.median,
    );
    if capped > 0 {
        t.push_str(&format!(
            "\n{capped} positions with capped seeds (not drawn)"
        ));
    }
    t
}

/// Truncate a record name to fit roughly within `max_px` of
/// monospace 11-px font (~6 px per glyph), with a `…` marker.
/// One alignment column: query residue, subject residue (already
/// complemented for the `+/-` view if applicable), and the match
/// class for colouring.
#[derive(Debug, Clone, Copy)]
struct AlignColumn {
    q: u8,
    s: u8,
    class: MatchClass,
    /// Substitution-matrix score for this column (0 for OutOfBounds /
    /// no matrix). Used for graduated similarity shading.
    score: i32,
}

/// One block of alignment columns plus its 5'/3' coordinate
/// metadata. Carried as a struct rather than a plain `Vec` so the
/// renderer doesn't have to re-derive end coordinates.
#[derive(Debug, Clone)]
struct AlignBlock {
    cols: Vec<AlignColumn>,
    /// 1-based residue index of the first query column.
    q_start: i64,
    /// 1-based residue index of the last query column.
    q_end: i64,
    /// 1-based residue index of the first subject column. For the
    /// reverse block this is the *higher* number (the subject is
    /// drawn right-to-left to match biological 5'→3' convention).
    s_start: i64,
    s_end: i64,
    /// True iff this block is the reverse-strand view.
    reverse: bool,
    /// BLASTX: the query reading frame shown (1–3, from the sequence
    /// start); query columns are codons, `q_start`/`q_end` in bases.
    frame: Option<u8>,
    /// Index of the crosshair column within `cols` (= `half`). Used
    /// to draw the caret marker.
    crosshair_col: usize,
}

/// Walk the diagonal through the crosshair and produce one
/// `AlignColumn` per output position. `reverse = false` does the
/// `+/+` walk (`q[q_c + i - half]` vs `s[s_c + i - half]`); `true`
/// does `+/-` (`q[q_c + i - half]` vs `complement(s[s_c + half - i])`).
///
/// The returned [`AlignBlock`] also carries the 1-based residue
/// indices each row spans, so the dock can print `5' N..M 3'` labels
/// at the row ends.
#[allow(clippy::too_many_arguments)]
fn build_align_columns(
    q_bytes: &[u8],
    s_bytes: &[u8],
    q_centre: usize,
    s_centre: usize,
    half: usize,
    window: usize,
    reverse: bool,
    matrix: Option<&dottir_core::ScoreMatrix>,
    mode: BlastMode,
) -> AlignBlock {
    let mut cols = Vec::with_capacity(window);
    for i in 0..window {
        let off = i as isize - half as isize;
        let qp = q_centre as isize + off;
        let sp = if reverse {
            // Walk the subject backwards from s_centre, taking the
            // complement so the displayed strand makes biological
            // sense: at i = half we sit at s_centre; i = half + 1
            // moves to s_centre - 1 with the base complemented.
            s_centre as isize - off
        } else {
            s_centre as isize + off
        };
        let (qc, q_oob) = lookup(q_bytes, qp);
        let (raw_sc, s_oob) = lookup(s_bytes, sp);
        let sc = if reverse && !s_oob {
            dottir_core::alphabet::complement_dna_byte(raw_sc)
        } else {
            raw_sc
        };
        let mut score: i32 = 0;
        let class = if q_oob || s_oob {
            MatchClass::OutOfBounds
        } else {
            let (c, sc_val) = classify_match_with_score(qc, sc, matrix, mode);
            score = sc_val;
            c
        };
        cols.push(AlignColumn {
            q: qc,
            s: sc,
            class,
            score,
        });
    }
    // Row coordinates: 1-based first/last residue indices for the
    // query and subject rows. For +/- the subject row is *displayed*
    // reversed, so first/last are swapped accordingly.
    let q_start_0 = q_centre as isize - half as isize;
    let q_end_0 = q_start_0 + window as isize - 1;
    let s_start_0 = if reverse {
        s_centre as isize + half as isize
    } else {
        s_centre as isize - half as isize
    };
    let s_end_0 = if reverse {
        s_centre as isize - (window as isize - 1 - half as isize)
    } else {
        s_centre as isize + (window as isize - 1 - half as isize)
    };
    AlignBlock {
        cols,
        q_start: q_start_0.saturating_add(1) as i64,
        q_end: q_end_0.saturating_add(1) as i64,
        s_start: s_start_0.saturating_add(1) as i64,
        s_end: s_end_0.saturating_add(1) as i64,
        reverse,
        frame: None,
        crosshair_col: half,
    }
}

/// BLASTX alignment columns: the query translated around base `q_centre`
/// in the reading frame (of the three codons covering it) that scores
/// best against subject residues `s_centre ± half`.
fn build_blastx_columns(
    q_bytes: &[u8],
    s_bytes: &[u8],
    q_centre: usize,
    s_centre: usize,
    half: usize,
    window: usize,
    matrix: Option<&dottir_core::ScoreMatrix>,
) -> AlignBlock {
    let mut best: Option<(i64, AlignBlock)> = None;
    for k in 0..3 {
        let Some(codon) = q_centre.checked_sub(k) else {
            continue;
        };
        let mut cols = Vec::with_capacity(window);
        for i in 0..window {
            let off = i as isize - half as isize;
            let qp = codon as isize + 3 * off;
            let q_oob = qp < 0 || qp as usize + 3 > q_bytes.len();
            let qc = if q_oob {
                b'-'
            } else {
                let c = &q_bytes[qp as usize..qp as usize + 3];
                dottir_core::translation::translate_codon([c[0], c[1], c[2]])
            };
            let (sc, s_oob) = lookup(s_bytes, s_centre as isize + off);
            let (class, score) = if q_oob || s_oob {
                (MatchClass::OutOfBounds, 0)
            } else {
                classify_match_with_score(qc, sc, matrix, BlastMode::Blastp)
            };
            cols.push(AlignColumn {
                q: qc,
                s: sc,
                class,
                score,
            });
        }
        let total: i64 = cols.iter().map(|c| i64::from(c.score)).sum();
        let first = codon as i64 - 3 * half as i64;
        let block = AlignBlock {
            cols,
            q_start: first + 1,
            q_end: first + 3 * window as i64,
            s_start: s_centre as i64 - half as i64 + 1,
            s_end: s_centre as i64 + (window - half) as i64,
            reverse: false,
            frame: Some((codon % 3) as u8 + 1),
            crosshair_col: half,
        };
        if best.as_ref().is_none_or(|(b, _)| total > *b) {
            best = Some((total, block));
        }
    }
    best.map(|(_, b)| b).unwrap_or(AlignBlock {
        cols: Vec::new(),
        q_start: 0,
        q_end: 0,
        s_start: 0,
        s_end: 0,
        reverse: false,
        frame: None,
        crosshair_col: 0,
    })
}

/// Render one 3-row alignment block (query, match line, subject)
/// with per-column background colour. Prepends a strand tag (`+/+` or
/// `+/-`) and the row's 5'/3' residue numbers in the side margins,
/// and draws a caret over the crosshair column so the user can see
/// which residue the click landed on.
/// `q_mark`/`s_mark`: when set, frame the crosshair-column cell on that row in
/// the given color — the per-axis "this residue is annotated" cue (F3). The
/// detail (which / how many features) is in the dock's overlap list.
fn draw_align_block(
    ui: &mut egui::Ui,
    ctx: &Context,
    block: &AlignBlock,
    q_mark: Option<Color32>,
    s_mark: Option<Color32>,
) {
    let font = egui::FontId::monospace(12.0);
    let small_font = egui::FontId::monospace(10.0);
    let glyph_w = ctx
        .fonts(|f| f.glyph_width(&font, 'A').max(f.glyph_width(&font, 'M')))
        .max(7.0);
    let row_h = 16.0;
    // Reserve enough left/right margin to fit `5' 12345678 ` style
    // labels at the row ends. 11 glyphs of monospace 10-px ≈ 60 px
    // at our font; round up.
    let side_label_w: f32 = 78.0;
    let strand_label_w: f32 = 28.0;
    let label_w = strand_label_w + side_label_w;
    let total_w = label_w * 2.0 + glyph_w * block.cols.len() as f32;
    // 3 alignment rows + 1 caret row above.
    let caret_h = 10.0;
    let total_h = caret_h + row_h * 3.0;
    let (rect, _) = ui.allocate_exact_size(Vec2::new(total_w, total_h), Sense::hover());
    let painter = ui.painter_at(rect);

    let bg_identity = Color32::from_rgb(0x9c, 0xdc, 0x7b);
    let bg_strong = Color32::from_rgb(0xc7, 0xef, 0xb7);
    let bg_weak = Color32::from_rgb(0xff, 0xea, 0x9c);
    let bg_oob = Color32::from_gray(218);
    let text_color = ui.style().visuals.text_color();
    let label_color = Color32::from_gray(90);
    let match_color = Color32::from_gray(40);
    let caret_color = Color32::from_rgb(0xc8, 0x3e, 0x3e);

    // Strand tag in the leftmost margin column, centred on the
    // alignment block (the 3 residue rows; the caret row is above).
    let alignment_top = rect.top() + caret_h;
    let strand_label = match block.frame {
        Some(f) => format!("F{f}"),
        None if block.reverse => "+/-".to_string(),
        None => "+/+".to_string(),
    };
    painter.text(
        Pos2::new(
            rect.left() + strand_label_w * 0.5,
            alignment_top + 1.5 * row_h,
        ),
        egui::Align2::CENTER_CENTER,
        strand_label,
        egui::FontId::monospace(11.0),
        label_color,
    );

    // 5'/3' coordinate labels at the row ends. Query is always shown
    // 5'→3' left-to-right; subject is 5'→3' left-to-right for +/+
    // and 3'→5' for +/- (i.e. we render the *displayed* coordinate
    // direction).
    let cols_left = rect.left() + label_w;
    let cols_right = cols_left + glyph_w * block.cols.len() as f32;
    let left_label_x = cols_left - 4.0;
    let right_label_x = cols_right + 4.0;
    let q_row_y = alignment_top + 0.5 * row_h;
    let s_row_y = alignment_top + 2.5 * row_h;
    // Query: always reads 5'→3' as drawn.
    painter.text(
        Pos2::new(left_label_x, q_row_y),
        egui::Align2::RIGHT_CENTER,
        format!("5' {}", block.q_start),
        small_font.clone(),
        label_color,
    );
    painter.text(
        Pos2::new(right_label_x, q_row_y),
        egui::Align2::LEFT_CENTER,
        format!("{} 3'", block.q_end),
        small_font.clone(),
        label_color,
    );
    // Subject: reads 5'→3' as drawn for forward; for reverse the
    // left end is the larger residue index (displayed strand is the
    // complement), so the 5' label sits on the left.
    let (s_left_label, s_right_label) = (
        format!("5' {}", block.s_start),
        format!("{} 3'", block.s_end),
    );
    painter.text(
        Pos2::new(left_label_x, s_row_y),
        egui::Align2::RIGHT_CENTER,
        s_left_label,
        small_font.clone(),
        label_color,
    );
    painter.text(
        Pos2::new(right_label_x, s_row_y),
        egui::Align2::LEFT_CENTER,
        s_right_label,
        small_font.clone(),
        label_color,
    );

    // Caret pointing at the crosshair column.
    let caret_x = cols_left + (block.crosshair_col as f32 + 0.5) * glyph_w;
    let caret_top = rect.top() + 1.0;
    let caret_bot = rect.top() + caret_h - 1.0;
    let half_base = (glyph_w * 0.45).min(5.0);
    painter.add(egui::Shape::convex_polygon(
        vec![
            Pos2::new(caret_x - half_base, caret_top),
            Pos2::new(caret_x + half_base, caret_top),
            Pos2::new(caret_x, caret_bot),
        ],
        caret_color,
        egui::Stroke::NONE,
    ));

    for (i, col) in block.cols.iter().enumerate() {
        let x = cols_left + i as f32 * glyph_w;
        let bg = column_background(col, bg_identity, bg_strong, bg_weak, bg_oob);
        let r0 = Rect::from_min_size(Pos2::new(x, alignment_top), Vec2::new(glyph_w, row_h));
        if let Some(c) = bg {
            painter.rect_filled(r0, 0.0, c);
        }
        painter.text(
            r0.center(),
            egui::Align2::CENTER_CENTER,
            char_to_string(col.q, col.class, true),
            font.clone(),
            text_color,
        );
        let r1 = Rect::from_min_size(
            Pos2::new(x, alignment_top + row_h),
            Vec2::new(glyph_w, row_h),
        );
        let match_ch = match col.class {
            MatchClass::Identical => "|",
            MatchClass::Positive => ":",
            _ => " ",
        };
        painter.text(
            r1.center(),
            egui::Align2::CENTER_CENTER,
            match_ch,
            font.clone(),
            match_color,
        );
        let r2 = Rect::from_min_size(
            Pos2::new(x, alignment_top + 2.0 * row_h),
            Vec2::new(glyph_w, row_h),
        );
        if let Some(c) = bg {
            painter.rect_filled(r2, 0.0, c);
        }
        painter.text(
            r2.center(),
            egui::Align2::CENTER_CENTER,
            char_to_string(col.s, col.class, false),
            font.clone(),
            text_color,
        );
    }

    // Vertical highlight on the crosshair column to tie the caret to
    // the residues underneath it.
    let xh = cols_left + block.crosshair_col as f32 * glyph_w;
    painter.rect_stroke(
        Rect::from_min_size(
            Pos2::new(xh, alignment_top),
            Vec2::new(glyph_w, row_h * 3.0),
        ),
        0.0,
        egui::Stroke::new(1.0, caret_color),
    );

    // Per-axis annotation markers: frame the query / subject cell at the
    // crosshair column when that residue lies inside a visible feature.
    if let Some(c) = q_mark {
        painter.rect_stroke(
            Rect::from_min_size(Pos2::new(xh, alignment_top), Vec2::new(glyph_w, row_h)),
            0.0,
            egui::Stroke::new(2.0, c),
        );
    }
    if let Some(c) = s_mark {
        painter.rect_stroke(
            Rect::from_min_size(
                Pos2::new(xh, alignment_top + 2.0 * row_h),
                Vec2::new(glyph_w, row_h),
            ),
            0.0,
            egui::Stroke::new(2.0, c),
        );
    }
}

/// Pick a background colour for one alignment column. Identity is
/// strongest; for Positive we further split by matrix-score magnitude
/// so the user can tell a high-confidence similarity (BLOSUM62 ≥ 2)
/// from a marginal one (`= 1`). OOB columns get the muted grey so the
/// row-end gutter doesn't blend into the column area.
fn column_background(
    col: &AlignColumn,
    bg_identity: Color32,
    bg_strong: Color32,
    bg_weak: Color32,
    bg_oob: Color32,
) -> Option<Color32> {
    match col.class {
        MatchClass::Identical => Some(bg_identity),
        MatchClass::Positive => {
            if col.score >= 2 {
                Some(bg_strong)
            } else {
                Some(bg_weak)
            }
        }
        MatchClass::OutOfBounds => Some(bg_oob),
        MatchClass::Other => None,
    }
}

/// Produce a plain-text representation of the alignment block(s)
/// suitable for clipboard paste. Output shape:
///
/// ```text
/// q = chr4:1234   s = chr5:5678   window = 100
///
/// +/+
/// query    ACGT...
///          ||::...
/// subject  ACGT...
///
/// +/-
/// query    ACGT...
///          |.::...
/// subject  ACGT...
/// ```
fn format_alignment_clipboard(
    forward: &AlignBlock,
    reverse: Option<&AlignBlock>,
    q_centre: usize,
    s_centre: usize,
    q_seq: Option<&Sequence>,
    s_seq: Option<&Sequence>,
) -> String {
    let mut out = String::new();
    use std::fmt::Write as _;
    let _ = writeln!(
        out,
        "q = {}   s = {}   window = {}",
        format_coord(q_seq, q_centre),
        format_coord(s_seq, s_centre),
        forward.cols.len(),
    );
    out.push('\n');
    match forward.frame {
        Some(f) => {
            let _ = writeln!(out, "query frame {f} (translated)");
        }
        None => {
            let _ = writeln!(out, "+/+");
        }
    }
    write_align_block_text(&mut out, forward);
    if let Some(rev) = reverse {
        out.push('\n');
        let _ = writeln!(out, "+/-");
        write_align_block_text(&mut out, rev);
    }
    out
}

fn write_align_block_text(out: &mut String, block: &AlignBlock) {
    use std::fmt::Write as _;
    let q: String = block
        .cols
        .iter()
        .map(|c| char_to_string(c.q, c.class, true))
        .collect();
    let m: String = block
        .cols
        .iter()
        .map(|c| match c.class {
            MatchClass::Identical => "|",
            MatchClass::Positive => ":",
            _ => " ",
        })
        .collect();
    let s: String = block
        .cols
        .iter()
        .map(|c| char_to_string(c.s, c.class, false))
        .collect();
    let _ = writeln!(
        out,
        "query    5' {:>8} {} {:<8} 3'",
        block.q_start, q, block.q_end
    );
    let _ = writeln!(out, "                    {}", m);
    let _ = writeln!(
        out,
        "subject  5' {:>8} {} {:<8} 3'",
        block.s_start, s, block.s_end
    );
}

/// Match class used by the alignment dock to colour columns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MatchClass {
    /// Identical residues (case-insensitive for ASCII letters).
    Identical,
    /// Positive matrix score, non-identical (e.g. BLOSUM62 A/S = 1).
    Positive,
    /// Either side ran past the sequence end.
    OutOfBounds,
    /// Mismatch / zero or negative matrix score.
    Other,
}

/// Compute the sliding-window score at residue midpoint `(q_mid,
/// s_mid)` using the same window convention the kernel uses (spec
/// §4.1.4): a length-`W` window centred so that the residue range
/// is `[q_mid - half, q_mid - half + W)`. Out-of-bounds columns
/// contribute 0. Used by `snap_crosshair_to_line` to refine within
/// the brightest pixel's `zoom × zoom` block down to exact residue
/// precision.
///
/// For [`Strand::Both`] / [`Strand::Reverse`] we also evaluate the
/// reverse-complement-subject window and return whichever is
/// larger, matching the max-merged pixelmap the user is snapping
/// against.
#[allow(clippy::too_many_arguments)]
fn window_score(
    qbytes: &[u8],
    sbytes: &[u8],
    q_mid: i64,
    s_mid: i64,
    half: i64,
    w: i64,
    qlen: i64,
    slen: i64,
    matrix: &ScoreMatrix,
    mode: BlastMode,
    strand: Strand,
) -> i64 {
    let n = matrix.size();
    let encode = |b: u8| -> usize {
        let bu = b.to_ascii_uppercase();
        let idx = match mode {
            BlastMode::Blastn => dottir_core::alphabet::encode_dna(bu),
            BlastMode::Blastp | BlastMode::Blastx => dottir_core::alphabet::encode_protein(bu),
        };
        idx as usize
    };
    let do_forward = !matches!(strand, Strand::Reverse);
    let do_reverse = !matches!(strand, Strand::Forward) && mode == BlastMode::Blastn;
    let mut fwd: i64 = 0;
    let mut rev: i64 = 0;
    for k in 0..w {
        let qi = q_mid - half + k;
        if qi < 0 || qi >= qlen {
            continue;
        }
        let qb = qbytes[qi as usize];
        let qidx = encode(qb);
        if qidx >= n {
            continue;
        }
        if do_forward {
            let si = s_mid - half + k;
            if si >= 0 && si < slen {
                let sidx = encode(sbytes[si as usize]);
                if sidx < n {
                    fwd += matrix.get(qidx, sidx) as i64;
                }
            }
        }
        if do_reverse {
            // Reverse-strand pixel at (q_mid, s_mid) corresponds to
            // querying q forward against subject walked backwards
            // (and complemented). The kernel's reverse pass uses
            // `s_idx + win2` for the centre offset; mirror that here.
            let si = s_mid + half - k;
            if si >= 0 && si < slen {
                let sb = dottir_core::alphabet::complement_dna_byte(sbytes[si as usize]);
                let sidx = encode(sb);
                if sidx < n {
                    rev += matrix.get(qidx, sidx) as i64;
                }
            }
        }
    }
    match (do_forward, do_reverse) {
        (true, true) => fwd.max(rev),
        (true, false) => fwd,
        (false, true) => rev,
        (false, false) => 0,
    }
}

/// Read a residue at a possibly-out-of-bounds coordinate. Returns
/// `(b'-', true)` for out-of-bounds, `(byte, false)` otherwise.
fn lookup(seq: &[u8], pos: isize) -> (u8, bool) {
    if pos < 0 || (pos as usize) >= seq.len() {
        (b'-', true)
    } else {
        (seq[pos as usize], false)
    }
}

/// Classify a residue column and return its substitution-matrix
/// score. Identity wins; otherwise consult the score matrix;
/// otherwise "Other" with score 0. Case-insensitive on the
/// identity check so the GUI doesn't trip on a mixed-case FASTA
/// (the kernel itself uppercases via the encode tables).
fn classify_match_with_score(
    q: u8,
    s: u8,
    matrix: Option<&dottir_core::ScoreMatrix>,
    mode: BlastMode,
) -> (MatchClass, i32) {
    let qu = q.to_ascii_uppercase();
    let su = s.to_ascii_uppercase();
    // Score lookup (uppercase, valid alphabet index, else 0).
    let score = match matrix {
        None => 0,
        Some(m) => {
            let (qi, si) = match mode {
                BlastMode::Blastn => (
                    dottir_core::alphabet::encode_dna(qu),
                    dottir_core::alphabet::encode_dna(su),
                ),
                BlastMode::Blastp | BlastMode::Blastx => (
                    dottir_core::alphabet::encode_protein(qu),
                    dottir_core::alphabet::encode_protein(su),
                ),
            };
            let n = m.size();
            if (qi as usize) >= n || (si as usize) >= n {
                0
            } else {
                m.get(qi as usize, si as usize)
            }
        }
    };
    if qu == su && qu.is_ascii_alphabetic() {
        return (MatchClass::Identical, score);
    }
    if matrix.is_none() {
        return (MatchClass::Other, 0);
    }
    if score > 0 {
        (MatchClass::Positive, score)
    } else {
        (MatchClass::Other, score)
    }
}

/// Render a residue byte for display in the alignment dock. Out-of-
/// bounds columns show as `-`; everything else uses the input byte
/// as a single ASCII char. `_is_query` is currently unused but kept
/// in the signature so future logic (e.g. translated-frame hints
/// for BLASTX) can branch on which row it's drawing.
fn char_to_string(b: u8, class: MatchClass, _is_query: bool) -> String {
    if matches!(class, MatchClass::OutOfBounds) {
        "-".to_string()
    } else if b.is_ascii_graphic() {
        (b as char).to_string()
    } else {
        "?".to_string()
    }
}

fn truncate_for_span(name: &str, max_px: f32) -> String {
    let max_chars = (max_px / 6.0).floor() as usize;
    if max_chars == 0 {
        return String::new();
    }
    if name.chars().count() <= max_chars {
        return name.to_string();
    }
    if max_chars <= 1 {
        return "…".to_string();
    }
    let kept: String = name.chars().take(max_chars - 1).collect();
    format!("{kept}…")
}

/// Pick a "nice" tick step (1, 2, or 5 × 10^k) so each tick is at
/// least `min_step` residues apart when projected to the screen.
fn nice_tick_step(span: f64, min_step: f32) -> f64 {
    if span <= 0.0 || !min_step.is_finite() || min_step <= 0.0 {
        return 1.0;
    }
    let target = min_step as f64;
    let exp = target.log10().floor();
    let base = 10f64.powf(exp);
    for &m in &[1.0, 2.0, 5.0, 10.0] {
        if m * base >= target {
            return m * base;
        }
    }
    10.0 * base
}

/// Per-sequence summary row for the right-side panel: name, total
/// residue count (with thousand separators), and record count.
/// Compact side-panel color legend for one axis: small-font visibility
/// checkbox + read-only swatch + `value (count)`, sorted by descending count.
/// Color editing lives in the Settings window; this stays dense on purpose.
fn annotation_legend_axis(
    ui: &mut egui::Ui,
    ax: &mut crate::annotation_overlay::AxisAnnot,
    header: Option<&str>,
) {
    if let Some(h) = header {
        ui.label(
            egui::RichText::new(h)
                .size(11.0)
                .strong()
                .color(Color32::from_gray(80)),
        );
    }
    let counts = ax.value_counts();
    if counts.is_empty() {
        ui.label(
            egui::RichText::new("(no features)")
                .size(11.0)
                .color(Color32::from_gray(140)),
        );
        return;
    }
    for (value, count) in &counts {
        ui.horizontal(|ui| {
            let mut visible = !ax.hidden.contains(value);
            if ui
                .checkbox(&mut visible, "")
                .on_hover_text("Show/hide this value's bands")
                .changed()
            {
                if visible {
                    ax.hidden.remove(value);
                } else {
                    ax.hidden.insert(value.clone());
                }
            }
            let col = ax.color_for(value);
            let (rect, _) = ui.allocate_exact_size(Vec2::new(11.0, 11.0), Sense::hover());
            ui.painter().rect_filled(rect, 2.0, col);
            ui.painter()
                .rect_stroke(rect, 2.0, egui::Stroke::new(1.0, Color32::from_gray(120)));
            ui.label(egui::RichText::new(format!("{value}  ({count})")).size(11.0));
        });
    }
}

/// One axis's crosshair-overlap line in the alignment dock: `tag:` followed by
/// a swatch + label per covering feature (several when annotations nest), or
/// `—` when nothing covers the crosshair. Label is `<color-by value> · <Name
/// or ID>`; the hover shows feature type and local coordinate span.
fn draw_overlap_line(
    ui: &mut egui::Ui,
    tag: &str,
    ax: Option<&crate::annotation_overlay::AxisAnnot>,
    hits: &[&crate::annotation_overlay::BoundFeature],
    seq: Option<&Sequence>,
) {
    let dim = Color32::from_gray(140);
    ui.horizontal_wrapped(|ui| {
        ui.spacing_mut().item_spacing.x = 4.0;
        ui.label(
            egui::RichText::new(format!("{tag}:"))
                .size(11.0)
                .strong()
                .color(Color32::from_gray(80)),
        );
        let Some(ax) = ax else {
            ui.label(egui::RichText::new("—").size(11.0).color(dim));
            return;
        };
        if hits.is_empty() {
            ui.label(egui::RichText::new("—").size(11.0).color(dim));
            return;
        }
        for bf in hits {
            let value = ax.value_of(bf);
            let col = ax.color_for(&value);
            let (rect, _) = ui.allocate_exact_size(Vec2::new(10.0, 10.0), Sense::hover());
            ui.painter().rect_filled(rect, 2.0, col);
            ui.painter()
                .rect_stroke(rect, 2.0, egui::Stroke::new(1.0, Color32::from_gray(120)));
            let f = &ax.set.features[bf.src];
            let ident = f.attrs.get("Name").or_else(|| f.attrs.get("ID"));
            let text = match ident {
                Some(id) if *id != value => format!("{value} · {id}"),
                _ => value.clone(),
            };
            ui.label(egui::RichText::new(text).size(11.0))
                .on_hover_text(format!(
                    "{}  {}–{}",
                    f.feature_type,
                    format_coord(seq, bf.range.start),
                    format_coord(seq, bf.range.end.saturating_sub(1)),
                ));
        }
    });
}

fn draw_seq_summary(ui: &mut egui::Ui, label: &str, seq: Option<&Sequence>) {
    let body = match seq {
        None => egui::RichText::new(format!("{label}: —"))
            .size(11.0)
            .color(Color32::from_gray(120)),
        Some(s) => {
            let name = s
                .source_path
                .as_ref()
                .and_then(|p| p.file_name())
                .map(|f| f.to_string_lossy().into_owned())
                .unwrap_or_else(|| "(in-memory)".into());
            let recs = s.records.len();
            let body = if recs <= 1 {
                format!("{label}: {name}  {} bp", format_thousands(s.len()))
            } else {
                format!(
                    "{label}: {name}  {} bp · {recs} recs",
                    format_thousands(s.len())
                )
            };
            egui::RichText::new(body).size(11.0)
        }
    };
    ui.label(body);
}

/// Format a `usize` with a thin-space-style thousands grouping —
/// `123456` → `123,456`. Used by the right-panel summary so big
/// genome sizes are readable at a glance.
fn format_thousands(n: usize) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, &b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(b as char);
    }
    out
}

/// Status-bar label for the loaded query/subject. Empty case shown
/// as e.g. "query: —"; populated as
/// `"query: chr4.fa (5,123,456 bp, 12 records)"`.
fn format_seq_label(name: &str, seq: Option<&Sequence>) -> String {
    let Some(seq) = seq else {
        return format!("{name}: —");
    };
    let fname = seq
        .source_path
        .as_ref()
        .and_then(|p| p.file_name())
        .map(|f| f.to_string_lossy().into_owned())
        .unwrap_or_else(|| "(in-memory)".into());
    if let [r] = seq.records.as_slice() {
        if r.origin.is_some() {
            return format!("{name}: {fname} {} ({} bp)", r.id, seq.len());
        }
    }
    if seq.records.len() <= 1 {
        format!("{name}: {fname} ({} bp)", seq.len())
    } else {
        format!(
            "{name}: {fname} ({} bp, {} records)",
            seq.len(),
            seq.records.len()
        )
    }
}

/// Format a concatenated-buffer coord as `record_id:position` when the
/// sequence has multiple records; falls back to the bare offset
/// otherwise. Used by the status bar at the crosshair.
/// Format a 0-indexed sequence position for display. Output is 1-indexed
/// to match biology / dotter convention — so calling with `coord = 0`
/// (the first base) returns `"1"`. The multi-record branch already did
/// `pos + 1` for the per-record position; this just keeps the single-
/// record branch consistent.
fn format_coord(seq: Option<&Sequence>, coord: usize) -> String {
    let one_indexed = coord + 1;
    let Some(seq) = seq else {
        return format!("{one_indexed}");
    };
    let region = seq.records.first().is_some_and(|r| r.origin.is_some());
    if seq.records.len() <= 1 && !region {
        return format!("{one_indexed}");
    }
    match seq.record_at(coord) {
        Some((rec, pos)) => match &rec.origin {
            // Genome-mode drill-down: report the source record's position.
            Some(o) => format!(
                "{}:{}{}",
                o.record,
                o.source_pos(pos) + 1,
                if o.reverse { " (−)" } else { "" }
            ),
            None => format!("{}:{}", rec.id, pos + 1),
        },
        None => format!("{one_indexed}"),
    }
}

fn pick_and_load(app: &mut DottirApp, role: SeqRole) {
    let label = match role {
        SeqRole::Query => "Open query sequence",
        SeqRole::Subject => "Open subject sequence",
    };
    let pick = rfd::FileDialog::new()
        .set_title(label)
        .add_filter(
            "Sequences (FASTA / GFF3+FASTA)",
            &["fa", "fasta", "fna", "faa", "gff", "gff3", "gz"],
        )
        .add_filter("FASTA", &["fa", "fasta", "fna", "faa", "gz"])
        .add_filter("GFF3 with sequences", &["gff", "gff3", "gz"])
        .pick_file();
    if let Some(path) = pick {
        app.load_sequence(role, path);
    }
}

/// `n` residues as bp / kb / Mb.
fn fmt_bp(n: u64) -> String {
    match n {
        0..=9_999 => format!("{n} bp"),
        10_000..=9_999_999 => format!("{:.1} kb", n as f64 / 1e3),
        _ => format!("{:.2} Mb", n as f64 / 1e6),
    }
}

/// Pieces of plot interval `range` of `seq`, record by record: record,
/// record-local lo..hi.
fn record_pieces<'a>(
    seq: &'a Sequence,
    range: &std::ops::Range<usize>,
) -> Vec<(&'a dottir_io::sequence::RecordSpan, usize, usize)> {
    seq.records
        .iter()
        .filter_map(|r| {
            let lo = range.start.max(r.range.start);
            let hi = range.end.min(r.range.end);
            (lo < hi).then(|| (r, lo - r.range.start, hi - r.range.start))
        })
        .collect()
}

/// Region specs (`name:start-end[:-]`, source coordinates) showing plot
/// interval `range` of `seq`, in plot order. `flipped`: the whole axis
/// is shown reverse-complemented (`-r`/`-v`), so plot positions count
/// from the end.
fn selection_specs(seq: &Sequence, range: &std::ops::Range<usize>, flipped: bool) -> Vec<String> {
    let n = seq.len();
    let r = if flipped {
        n - range.end..n - range.start
    } else {
        range.clone()
    };
    let mut specs: Vec<String> = record_pieces(seq, &r)
        .into_iter()
        .map(|(rec, lo, hi)| {
            let (name, _, a, b, rev) = rec.source_interval(lo, hi);
            dottir_io::faidx::region_spec(name, &(a..b), rev ^ flipped)
        })
        .collect();
    if flipped {
        specs.reverse();
    }
    specs
}

/// BEDPE lines (query first, subject second; 0-based, end-exclusive) for
/// every record pair under a selection.
fn selection_bedpe(
    q: &Sequence,
    s: &Sequence,
    qr: &std::ops::Range<usize>,
    sr: &std::ops::Range<usize>,
) -> String {
    let mut out = String::new();
    for (qrec, qlo, qhi) in record_pieces(q, qr) {
        let (qn, _, qa, qb, qrev) = qrec.source_interval(qlo, qhi);
        for (srec, slo, shi) in record_pieces(s, sr) {
            let (sn, _, sa, sb, srev) = srec.source_interval(slo, shi);
            let strand = |r: bool| if r { '-' } else { '+' };
            out.push_str(&format!(
                "{qn}\t{qa}\t{qb}\t{sn}\t{sa}\t{sb}\tselection\t.\t{}\t{}\n",
                strand(qrev),
                strand(srev)
            ));
        }
    }
    out
}

fn pick_and_load_annotation(app: &mut DottirApp, role: SeqRole) {
    let label = match role {
        SeqRole::Query => "Load annotation for query",
        SeqRole::Subject => "Load annotation for subject",
    };
    let pick = rfd::FileDialog::new()
        .set_title(label)
        .add_filter("Annotations (GFF3/BED)", &["gff", "gff3", "bed", "gz"])
        .pick_file();
    if let Some(path) = pick {
        app.load_annotation(role, path);
    }
}

fn save_svg(app: &mut DottirApp) {
    let Some(plot) = &app.plot else {
        app.last_error = Some("nothing to save — compute a plot first".into());
        return;
    };
    let shown = app.display_rgba(true);
    let metadata = app.export_metadata();
    let chain_lines = app.shown_chain_lines();
    let pick = rfd::FileDialog::new()
        .set_title("Save SVG")
        .add_filter("SVG", &["svg"])
        .save_file();
    if let Some(path) = pick {
        let Some((rgba, _, _)) = shown else {
            return;
        };
        let pixels = ExportPixels::of(&rgba);
        let recs_x: Vec<_> = app
            .query
            .as_ref()
            .map(|q| {
                q.records
                    .iter()
                    .map(|r| {
                        dottir_io::text_overlay::AxisRecord::new(
                            r.id.clone(),
                            r.range.start as u32,
                            r.range.end as u32,
                        )
                    })
                    .collect()
            })
            .unwrap_or_default();
        let recs_y: Vec<_> = app
            .subject
            .as_ref()
            .map(|s| {
                s.records
                    .iter()
                    .map(|r| {
                        dottir_io::text_overlay::AxisRecord::new(
                            r.id.clone(),
                            r.range.start as u32,
                            r.range.end as u32,
                        )
                    })
                    .collect()
            })
            .unwrap_or_default();
        let meta: Vec<(&str, &str)> = metadata.iter().map(|(k, v)| (*k, v.as_str())).collect();
        let (w, h) = (plot.width, plot.height);
        let (origin, qlen, slen) = app.export_axes();
        let rgb = match pixels {
            ExportPixels::Grey(g) => g.iter().flat_map(|&v| [v, v, v]).collect(),
            ExportPixels::Rgb(c) => c,
        };
        let written = dottir_io::svg_export::write_rgb_svg_with_lines(
            &path,
            w,
            h,
            origin,
            qlen,
            slen,
            &rgb,
            50,
            &recs_x,
            &recs_y,
            &meta,
            &chain_lines,
        );
        match written {
            Ok(()) => {}
            Err(e) => app.last_error = Some(format!("SVG save failed: {e}")),
        }
    }
}

fn save_session(app: &mut DottirApp) {
    use crate::session::{
        codec, Session, SessionGreyramp, SessionPlot, SessionView, SESSION_VERSION,
    };
    // A genome drill-down's regions are not files of their own; the
    // genome section reopens the comparison instead.
    // Sorted sequences are saved as their files plus the sorts to redo.
    let file = |s: Option<&Sequence>| {
        s.filter(|s| s.records.iter().all(|r| r.origin.is_none()))
            .and_then(|s| s.source_path.clone())
    };
    let session = Session {
        version: SESSION_VERSION,
        query: file(app.unsorted_query.as_ref().or(app.query.as_ref())),
        subject: file(app.unsorted_subject.as_ref().or(app.subject.as_ref())),
        sorted_records: app
            .record_sorts
            .iter()
            .map(|a| match a {
                SortAxis::Vertical => "vertical".to_string(),
                SortAxis::Horizontal => "horizontal".to_string(),
            })
            .collect(),
        hsps: app.hsps.as_ref().map(|h| crate::session::SessionHsps {
            paf: h.path.clone(),
            mode: app.hsp_mode.name().to_string(),
        }),
        plot: SessionPlot {
            mode: codec::mode_to_str(app.settings.mode).to_string(),
            matrix_name: app.settings.matrix_name.clone(),
            window_size: app.settings.window_size,
            zoom: app.settings.zoom,
            // Encode auto pixel_fac as 0 in the session — matches the
            // core's sentinel convention.
            pixel_fac: if app.settings.auto_pixel_fac {
                0
            } else {
                app.settings.pixel_fac
            },
            strand: codec::strand_to_str(app.settings.strand).to_string(),
            self_comparison: app.settings.self_comparison,
            triangle: codec::triangle_to_str(app.settings.triangle).to_string(),
            memory_limit_mib: app.settings.memory_limit_bytes / (1024 * 1024),
        },
        greyramp: SessionGreyramp {
            white: app.greyramp.white,
            black: app.greyramp.black,
            swap: app.greyramp.swap,
        },
        view: SessionView {
            offset_x: app.view_offset.x,
            offset_y: app.view_offset.y,
            // Always 1.0 under the dotter-faithful 1:1 model; field
            // kept in the session schema for backward compatibility.
            display_zoom: 1.0,
            crosshair: app.crosshair.map(|(q, s)| [q, s]),
            light_theme: app.light_theme,
            align_dock_visible: app.align_dock_visible,
            align_window_size: app.settings.align_window_size,
        },
        large: Some(app.settings.large.to_session()),
        genome: app
            .genome
            .as_ref()
            .map(|g| crate::session::SessionGenome::from_view(g, app.genome_active)),
    };
    let pick = rfd::FileDialog::new()
        .set_title("Save dottir session")
        .add_filter("TOML", &["toml"])
        .set_file_name("dottir-session.toml")
        .save_file();
    if let Some(path) = pick {
        match session.save(&path) {
            Ok(()) => tracing::info!("wrote session {}", path.display()),
            Err(e) => app.last_error = Some(format!("session save failed: {e}")),
        }
    }
}

fn open_session(app: &mut DottirApp, ctx: &Context) {
    use crate::session::{codec, Session};
    let pick = rfd::FileDialog::new()
        .set_title("Open dottir session")
        .add_filter("TOML", &["toml"])
        .pick_file();
    let Some(path) = pick else { return };
    let s = match Session::load(&path) {
        Ok(s) => s,
        Err(e) => {
            app.last_error = Some(format!("session load failed: {e}"));
            return;
        }
    };
    // Apply: settings → app.settings, then re-load the FASTAs (paths
    // recorded in the session), then apply view/greyramp/crosshair.
    if let Some(m) = codec::mode_from_str(&s.plot.mode) {
        app.settings.mode = m;
    }
    app.settings.matrix_name = s.plot.matrix_name;
    app.settings.window_size = s.plot.window_size;
    app.settings.zoom = s.plot.zoom.max(1);
    // `pixel_fac = 0` in the session means "auto" (matches the core's
    // sentinel). Restore the checkbox; leave the slider at a default 50.
    if s.plot.pixel_fac == 0 {
        app.settings.auto_pixel_fac = true;
        app.settings.pixel_fac = 50;
    } else {
        app.settings.auto_pixel_fac = false;
        app.settings.pixel_fac = s.plot.pixel_fac;
    }
    if let Some(st) = codec::strand_from_str(&s.plot.strand) {
        app.settings.strand = st;
    }
    app.settings.self_comparison = s.plot.self_comparison;
    if let Some(t) = codec::triangle_from_str(&s.plot.triangle) {
        app.settings.triangle = t;
    }
    app.settings.memory_limit_bytes = s.plot.memory_limit_mib * 1024 * 1024;
    app.settings.large = s
        .large
        .as_ref()
        .map(LargeSettings::from_session)
        .unwrap_or_default();
    app.large_source = None;
    app.greyramp = Greyramp {
        white: s.greyramp.white,
        black: s.greyramp.black,
        swap: s.greyramp.swap,
    };
    app.texture_dirty = true;
    // Re-load sequences from the recorded paths. Errors don't abort
    // the load; they just leave the relevant slot empty + an error
    // in the status bar.
    if let Some(p) = s.query {
        app.load_sequence(SeqRole::Query, p);
    } else {
        app.query = None;
        app.plot = None;
    }
    if let Some(p) = s.subject {
        app.load_sequence(SeqRole::Subject, p);
    } else {
        app.subject = None;
        app.plot = None;
    }
    // Apply view state AFTER the loads (load_sequence calls recompute
    // which doesn't touch view state, so applying here is safe).
    app.view_offset = Vec2::new(s.view.offset_x, s.view.offset_y);
    // `s.view.display_zoom` is ignored under the new always-1:1
    // model — the saved field is kept for backward-compatible
    // session files but doesn't influence rendering.
    let _ = s.view.display_zoom;
    // Crosshair: schema v1 stored *full-sequence pixmap* coords;
    // schema v2+ stores absolute residue coords. Migrate v1 values
    // by remapping through the recorded compute zoom (centre of the
    // old pixel block in residue space — matches what the v1 GUI
    // showed in the status bar before the navigation rework).
    app.crosshair = s.view.crosshair.map(|[q, s_]| {
        if s.version <= 1 {
            let z = s.plot.zoom.max(1);
            (q * z + z / 2, s_ * z + z / 2)
        } else {
            (q, s_)
        }
    });
    app.light_theme = s.view.light_theme;
    app.align_dock_visible = s.view.align_dock_visible;
    app.settings.align_window_size = s.view.align_window_size.clamp(20, 400);
    if let Some(h) = &s.hsps {
        app.hsp_mode = HspMode::parse(&h.mode).unwrap_or(HspMode::Line);
        app.load_hsps(h.paf.clone());
    } else {
        app.hsps = None;
    }
    app.record_sort_queue = s
        .sorted_records
        .iter()
        .filter_map(|a| match a.as_str() {
            "vertical" => Some(SortAxis::Vertical),
            "horizontal" => Some(SortAxis::Horizontal),
            _ => None,
        })
        .collect();
    app.genome_active = false;
    if let Some(g) = s.genome {
        let inputs = g.inputs.clone();
        let active = g.active;
        app.genome_pending = Some(g);
        app.start_genome_load(ctx, inputs);
        // The dotplot stays shown when that is what was saved.
        app.genome_reopen_dotplot = !active && app.query.is_some();
    }
    tracing::info!("loaded session {}", path.display());
}

fn save_png(app: &mut DottirApp) {
    save_raster(app, false);
}

fn save_pdf(app: &mut DottirApp) {
    save_raster(app, true);
}

/// The plot as shown, with axes, to a PNG or a one-page PDF.
fn save_raster(app: &mut DottirApp, pdf: bool) {
    let Some(plot) = &app.plot else {
        app.last_error = Some("nothing to save — compute a plot first".into());
        return;
    };
    let shown = app.display_rgba(true);
    let metadata = app.export_metadata();
    let (title, kind, ext) = if pdf {
        ("Save PDF", "PDF", "pdf")
    } else {
        ("Save PNG", "PNG", "png")
    };
    let pick = rfd::FileDialog::new()
        .set_title(title)
        .add_filter(kind, &[ext])
        .save_file();
    if let Some(path) = pick {
        let Some((rgba, _, _)) = shown else {
            return;
        };
        let pixels = ExportPixels::of(&rgba);
        let recs_x: Vec<_> = app
            .query
            .as_ref()
            .map(|q| {
                q.records
                    .iter()
                    .map(|r| {
                        dottir_io::text_overlay::AxisRecord::new(
                            r.id.clone(),
                            r.range.start as u32,
                            r.range.end as u32,
                        )
                    })
                    .collect()
            })
            .unwrap_or_default();
        let recs_y: Vec<_> = app
            .subject
            .as_ref()
            .map(|s| {
                s.records
                    .iter()
                    .map(|r| {
                        dottir_io::text_overlay::AxisRecord::new(
                            r.id.clone(),
                            r.range.start as u32,
                            r.range.end as u32,
                        )
                    })
                    .collect()
            })
            .unwrap_or_default();
        // Bigger margin when we have record names so the extra label
        // row above the ticks has breathing room (matches CLI).
        let margin = if recs_x.len() > 1 || recs_y.len() > 1 {
            70
        } else {
            50
        };
        let (origin, qlen, slen) = app.export_axes();
        let meta: Vec<(&str, &str)> = metadata.iter().map(|(k, v)| (*k, v.as_str())).collect();
        let (w, h) = (plot.width, plot.height);
        use dottir_io::pdf_export::{write_pdf, PdfPixels};
        use dottir_io::text_overlay::{compose_image_with_axes_at, compose_rgb_image_with_axes_at};
        let written = match &pixels {
            ExportPixels::Grey(g) => {
                let (img, tw, th) = compose_image_with_axes_at(
                    g, w, h, origin, qlen, slen, margin, &recs_x, &recs_y,
                );
                if pdf {
                    write_pdf(&path, tw, th, PdfPixels::Grey(&img), &meta)
                        .map_err(|e| e.to_string())
                } else {
                    dottir_io::png_export::write_grayscale_png(&path, tw, th, &img, &meta)
                        .map_err(|e| e.to_string())
                }
            }
            ExportPixels::Rgb(c) => {
                let (img, tw, th) = compose_rgb_image_with_axes_at(
                    c, w, h, origin, qlen, slen, margin, &recs_x, &recs_y,
                );
                if pdf {
                    write_pdf(&path, tw, th, PdfPixels::Rgb(&img), &meta).map_err(|e| e.to_string())
                } else {
                    dottir_io::png_export::write_rgb_png(&path, tw, th, &img, &meta)
                        .map_err(|e| e.to_string())
                }
            }
        };
        if let Err(e) = written {
            app.last_error = Some(format!("{kind} save failed: {e}"));
        }
    }
}

/// Displayed pixels for an image export: greyscale when nothing on screen
/// is coloured, else RGB.
#[derive(Debug, PartialEq)]
enum ExportPixels {
    Grey(Vec<u8>),
    Rgb(Vec<u8>),
}

impl ExportPixels {
    fn of(rgba: &[u8]) -> Self {
        if rgba.chunks_exact(4).all(|p| p[0] == p[1] && p[1] == p[2]) {
            Self::Grey(rgba.chunks_exact(4).map(|p| p[0]).collect())
        } else {
            Self::Rgb(
                rgba.chunks_exact(4)
                    .flat_map(|p| [p[0], p[1], p[2]])
                    .collect(),
            )
        }
    }
}

#[cfg(test)]
mod large_display_tests {
    use super::*;
    use dottir_core::large::coverage_code;

    #[test]
    fn export_pixels_keep_colour_only_when_shown() {
        let grey = [7, 7, 7, 255, 200, 200, 200, 255];
        assert_eq!(ExportPixels::of(&grey), ExportPixels::Grey(vec![7, 200]));
        let colour = [7, 7, 7, 255, 220, 30, 220, 255];
        assert_eq!(
            ExportPixels::of(&colour),
            ExportPixels::Rgb(vec![7, 7, 7, 220, 30, 220])
        );
    }

    #[test]
    fn display_lut_follows_coverage_sensitivity() {
        // Review 2026-10-28 #4: exports used the score greyramp for
        // coverage codes, ignoring the coverage sliders.
        let g = Greyramp::default();
        let a = display_lut(true, &g, 100.0, 10_000.0);
        let b = display_lut(true, &g, 1000.0, 10_000.0);
        assert_eq!(a, coverage_lut(100.0, 10_000.0));
        assert_ne!(a, b, "coverage sensitivity changes the mapping");
        assert_ne!(a, g.lut());
        assert_eq!(display_lut(false, &g, 100.0, 10_000.0), g.lut());
    }

    #[test]
    fn coverage_lut_ramps_in_log_bases() {
        let lut = coverage_lut(100.0, 10_000.0);
        assert_eq!(lut[0], 255, "empty pixels stay white");
        assert_eq!(lut[coverage_code(50.0) as usize], 255);
        assert_eq!(lut[coverage_code(20_000.0) as usize], 0);
        // 1 kb is half-way between 100 bp and 10 kb on a log scale.
        let mid = lut[coverage_code(1000.0) as usize];
        assert!((120..=135).contains(&mid), "{mid}");
        assert!(lut[1..].windows(2).all(|w| w[0] >= w[1]), "monotone");
    }

    #[test]
    fn depth_bin_tooltip_names_copies_and_cap() {
        let b = DepthBin {
            median: 7,
            max: u16::MAX,
            class_frac: [0, 51, 153, 51],
        };
        let t = depth_bin_text("Query", 12_000, 1000, &[1000.0, 0.0, 3000.0, 0.0], &b, 3);
        assert!(t.contains("Query 12,000–13,000"), "{t}");
        assert!(t.contains("4,000 bp matching = 4.0 copies per base"), "{t}");
        assert!(t.contains("unique 25%"), "{t}");
        assert!(t.contains("mid 75%"), "{t}");
        assert!(t.contains("median 7, max ≥ max-occ"), "{t}");
        assert!(t.contains("3 positions with capped seeds"), "{t}");
    }
}

#[cfg(test)]
mod axis_profile_tests {
    use super::*;
    use dottir_core::large::{compute_large_dotplot, LargeInput};

    #[test]
    fn profile_is_the_marginal_of_visible_dots() {
        // Self-comparison of random sequence: only the main diagonal, so
        // every column carries ~zoom matching bp of unique dots.
        let mut st = 7u64;
        let q: Vec<u8> = (0..20_000)
            .map(|_| {
                st = st
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(st >> 62) as usize]
            })
            .collect();
        let cfg = LargeConfig {
            zoom: 500,
            self_comparison: true,
            coverage: true,
            depth: DepthMethod::None,
            ..LargeConfig::default()
        };
        let lp = compute_large_dotplot(&LargeInput::self_comparison(&q, &[]), &cfg).unwrap();
        let lut = coverage_lut(1.0, 1000.0);
        let p = AxisProfile::of(&lp, ClassMask::ALL, true, &lut).unwrap();
        let unique = RepeatClass::Unique.index();
        for col in &p.q[2..p.q.len() - 2] {
            let copies = col.iter().sum::<f32>() / 500.0;
            assert!((0.8..1.3).contains(&copies), "{copies}");
            assert_eq!(col.iter().sum::<f32>(), col[unique]);
        }
        // Hiding every dot (white above anything present) empties it.
        let none = AxisProfile::of(&lp, ClassMask::ALL, true, &coverage_lut(1e6, 2e6)).unwrap();
        assert!(none.q.iter().all(|c| c.iter().sum::<f32>() == 0.0));
        // Hidden classes do not count.
        let low_only = AxisProfile::of(
            &lp,
            ClassMask::from_classes(&[RepeatClass::Low]),
            true,
            &lut,
        )
        .unwrap();
        assert!(low_only.s.iter().all(|c| c.iter().sum::<f32>() == 0.0));
    }
}

#[cfg(test)]
mod selection_tests {
    use super::*;

    fn seq(text: &str) -> Sequence {
        Sequence::from_records(dottir_io::fasta::parse_fasta(text).unwrap(), None)
    }

    #[test]
    fn selections_map_to_source_regions() {
        // a: 0..10, b: 10..30 in the plot.
        let s = seq(">a\nACGTACGTAC\n>b\nACGTACGTACGTACGTACGT\n");
        assert_eq!(selection_specs(&s, &(5..15), false), ["a:6-10", "b:1-5"]);
        // Axis shown reverse-complemented: plot 0..5 is b's last 5 bp.
        assert_eq!(selection_specs(&s, &(0..5), true), ["b:16-20:-"]);
        assert_eq!(
            selection_specs(&s, &(15..25), true),
            ["b:1-5:-", "a:6-10:-"]
        );
        // A sorted (flipped) record keeps source coordinates.
        let f = s.reordered(&[1, 0], &[false, true]);
        assert_eq!(selection_specs(&f, &(0..5), false), ["b:16-20:-"]);
    }

    #[test]
    fn bedpe_lists_every_record_pair() {
        let q = seq(">a\nACGTACGTAC\n>b\nACGTACGTACGTACGTACGT\n");
        let s = seq(">s\nACGTACGTAC\n").reordered(&[0], &[true]);
        assert_eq!(
            selection_bedpe(&q, &s, &(8..12), &(0..4)),
            "a\t8\t10\ts\t6\t10\tselection\t.\t+\t-\n\
             b\t0\t2\ts\t6\t10\tselection\t.\t+\t-\n"
        );
    }
}

#[cfg(test)]
mod blastx_view_tests {
    use super::*;

    #[test]
    fn blastx_view_picks_the_frame_that_matches() {
        // MKW encoded in frame 2 (one leading base); subject MKW.
        let q = b"GATGAAATGGCC";
        let s = b"MKW";
        let m = dottir_core::ScoreMatrix::blosum62();
        // Crosshair on the second base of the K codon (base 5).
        let b = build_blastx_columns(q, s, 5, 1, 1, 3, Some(&m));
        let pep: Vec<u8> = b.cols.iter().map(|c| c.q).collect();
        assert_eq!(pep, b"MKW");
        assert_eq!(b.frame, Some(2));
        assert_eq!((b.q_start, b.q_end), (2, 10));
        assert!(b.cols.iter().all(|c| c.class == MatchClass::Identical));
        // Near the start, codons before the sequence are out of bounds.
        let b = build_blastx_columns(q, s, 1, 0, 1, 3, Some(&m));
        assert_eq!(b.cols[0].class, MatchClass::OutOfBounds);
    }
}
