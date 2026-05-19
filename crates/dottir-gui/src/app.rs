//! Top-level egui application.

use std::path::PathBuf;

use dottir_core::{
    extract_ridges_from_pixels, snap_zoom_to_period_divisor, BlastMode, DotPlot, PlotConfig, Ridge,
    RidgeDirection, RidgeParams, ScoreMatrix, Strand, Triangle,
};
use dottir_io::{DetectedAlphabet, Sequence};
use egui::{
    Color32, ColorImage, Context, Pos2, Rect, Sense, Slider, TextureHandle, TextureOptions, Vec2,
};

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
}

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
    /// Suppress recompute during constructor — load_fasta() runs
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
    pixelmap_cache: std::collections::VecDeque<(ViewSlice, u32, DotPlot)>,
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
    /// Residue range covered by the most recently computed pixelmap.
    /// `None` before any compute has completed. Coordinate transforms
    /// in the paint loop add `slice.start` when converting
    /// pixmap_pixel → residue.
    current_slice: Option<ViewSlice>,
    /// Slice to use for the *next* dispatched compute. Set by
    /// `action_rect_zoom` / `action_fit` / `action_back` /
    /// `load_fasta`; consumed in `recompute`. `None` means "use the
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
    /// path that runs before the first `update()` (load_fasta in
    /// the constructor) behaves as on a non-HiDPI display.
    pixels_per_point: f32,
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
            pending_view_offset_after_compute: None,
            pending_rect_zoom: None,
            current_slice: None,
            target_slice: None,
            locked_pixel_fac: None,
            last_dispatched_slice: None,
            current_ridges: Vec::new(),
            pixels_per_point: 1.0,
        };

        // Pre-load any sequences supplied on the command line. Errors
        // are stashed in `last_error` and surfaced in the status bar
        // instead of aborting the GUI — the user can still recover via
        // File → Open.
        if let Some(p) = startup.query {
            app.load_fasta(SeqRole::Query, p);
        }
        if let Some(p) = startup.subject {
            app.load_fasta(SeqRole::Subject, p);
        }
        // Release the suspend. The first compute can't fire yet
        // because the canvas hasn't been measured — `pending_initial_
        // compute` is `true` and `update()` will fire it after the
        // first paint records `measured_plot_area`.
        app.suspend_recompute = false;
        if app.query.is_some() && app.subject.is_some() {
            app.pending_initial_compute = true;
        }
        app
    }

    fn load_fasta(&mut self, role: SeqRole, path: PathBuf) {
        match Sequence::load(&path) {
            Ok(seq) => {
                let detected = seq.detect_alphabet();
                tracing::info!(
                    "loaded {} ({} residues, {} records, detected {:?})",
                    path.display(),
                    seq.len(),
                    seq.records.len(),
                    detected,
                );
                match role {
                    SeqRole::Query => self.query = Some(seq),
                    SeqRole::Subject => self.subject = Some(seq),
                }
                self.last_error = None;
                // Fresh sequence data → reset view state and caches.
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
            Err(e) => {
                self.last_error = Some(format!("failed to load {}: {e}", path.display()));
            }
        }
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
        let zoom_q = (q.len() as f64 / phys_w).ceil() as u32;
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
    fn cache_get(&mut self, slice: &ViewSlice, compute_zoom: u32) -> Option<DotPlot> {
        let pos = self
            .pixelmap_cache
            .iter()
            .position(|(s, z, _)| *z == compute_zoom && s == slice)?;
        let (s, z, plot) = self.pixelmap_cache.remove(pos)?;
        let restored = plot.clone();
        self.pixelmap_cache.push_back((s, z, plot));
        Some(restored)
    }

    /// Insert a freshly computed pixelmap into the cache; evict
    /// oldest if at capacity.
    fn cache_insert(&mut self, slice: ViewSlice, compute_zoom: u32, plot: DotPlot) {
        self.pixelmap_cache
            .retain(|(s, z, _)| *z != compute_zoom || s != &slice);
        if self.pixelmap_cache.len() >= HISTORY_MAX {
            self.pixelmap_cache.pop_front();
        }
        self.pixelmap_cache.push_back((slice, compute_zoom, plot));
    }

    /// Drop everything pixel-dependent — call when settings that
    /// affect compute output change (matrix, window, strand, etc.).
    /// Also clears the Karlin-locked pixel_fac so the next compute
    /// re-derives it under the new settings.
    fn invalidate_caches(&mut self) {
        self.pixelmap_cache.clear();
        self.history.clear();
        self.locked_pixel_fac = None;
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
        let (Some(plot), Some(q), Some(s), Some(cur_slice)) = (
            self.plot.as_ref(),
            self.query.as_ref(),
            self.subject.as_ref(),
            self.current_slice.clone(),
        ) else {
            return;
        };
        let cur_zoom_u = plot.params.zoom.max(1);
        let cur_zoom = cur_zoom_u as f32;
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
            return;
        }
        let px_lo_x = px_lo_x.max(0.0).min(plot.width as f32);
        let px_lo_y = px_lo_y.max(0.0).min(plot.height as f32);
        let px_hi_x = px_hi_x.max(px_lo_x).min(plot.width as f32);
        let px_hi_y = px_hi_y.max(px_lo_y).min(plot.height as f32);

        // Pixelmap pixels → *absolute* residues (offset by the
        // current slice's start).
        let sel_lo_x = cur_slice.q_range.start + (px_lo_x * cur_zoom).floor().max(0.0) as usize;
        let sel_lo_y = cur_slice.s_range.start + (px_lo_y * cur_zoom).floor().max(0.0) as usize;
        let sel_hi_x =
            (cur_slice.q_range.start + (px_hi_x * cur_zoom).ceil() as usize).min(q.len());
        let sel_hi_y =
            (cur_slice.s_range.start + (px_hi_y * cur_zoom).ceil() as usize).min(s.len());
        let sel_w = sel_hi_x.saturating_sub(sel_lo_x).max(1);
        let sel_h = sel_hi_y.saturating_sub(sel_lo_y).max(1);

        // New compute zoom so the *selection* (not the slice) fills
        // the *physical* canvas — matches `canonical_compute_zoom`'s
        // Option C policy (one pixmap texel = one physical pixel).
        let phys_w = (plot_area.width() as f64).max(1.0) * (ppp as f64);
        let phys_h = (plot_area.height() as f64).max(1.0) * (ppp as f64);
        let zoom_x = ((sel_w as f64) / phys_w).ceil() as u32;
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
        let new_off_x = ((sel_lo_x - new_slice.q_range.start) as f64 / new_zoom as f64) as f32;
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
        if let Some(plot) = self.cache_get(&target_slice, target_zoom) {
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
            if let Some(off) = self.pending_view_offset_after_compute.take() {
                self.view_offset = off;
            }
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
            window_size: self.settings.window_size,
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
                        self.cache_insert(slice.clone(), plot.params.zoom, plot.clone());
                        self.current_slice = Some(slice);
                    }
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

    fn ensure_texture(&mut self, ctx: &Context) {
        let Some(plot) = &self.plot else {
            return;
        };
        if !self.texture_dirty && self.texture.is_some() {
            return;
        }
        let lut = self.greyramp.lut();
        let (pw, ph) = (plot.width as usize, plot.height as usize);
        let fwd_view = plot.pixels.as_slice();
        let rev_view = plot.reverse_pixels.as_deref();

        // Spec §4.4.3 inverted-repeat highlighting: when the plot has
        // a separate reverse channel, paint forward in grey and
        // reverse in magenta — overlapping cells take whichever
        // channel is stronger after the greyramp.
        let mut rgba = Vec::with_capacity(pw * ph * 4);
        match rev_view {
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
        let image = ColorImage::from_rgba_unmultiplied([pw, ph], &rgba);
        // The dotter-faithful render is always 1:1 (pixelmap pixel =
        // screen pixel). NEAREST → no GPU resampling, no sub-pixel
        // sampling phase, no moiré. The compute step is the only
        // place "zoom" happens.
        let handle = ctx.load_texture("dottir-pixelmap", image, TextureOptions::NEAREST);
        self.texture = Some(handle);
        self.texture_dirty = false;
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

        // Esc → Back (history pop). Backspace also; both are dotter-
        // classic for "step out".
        if ctx.input(|i| i.key_pressed(egui::Key::Escape) || i.key_pressed(egui::Key::Backspace)) {
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

        // Deferred initial compute: pre-loaded sequences (CLI or
        // load_fasta before the first paint) couldn't pick a
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
        let cq_pix = ((cq_seq as i64 - q_off).max(0) / z).clamp(0, pw - 1);
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
        let mut q_residue = q_off + qp_best * z + z / 2;
        let mut s_residue = s_off + sp_best * z + z / 2;
        // Fine step: walk every residue midpoint in the block,
        // compute the actual sliding-window score, pick the max.
        // Only do this when both sequences are loaded and we have a
        // score matrix — otherwise stick with the block centre.
        if let (Some(qseq), Some(sseq), Some(matrix)) = (
            self.query.as_ref(),
            self.subject.as_ref(),
            self.settings.build_matrix(),
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
            .unwrap_or(q_off + pw * z);
        let slen_cap = self
            .subject
            .as_ref()
            .map(|s| s.len() as i64)
            .unwrap_or(s_off + ph * z);
        let q_residue = q_residue.clamp(0, qlen_cap.saturating_sub(1).max(0));
        let s_residue = s_residue.clamp(0, slen_cap.saturating_sub(1).max(0));
        self.crosshair = Some((q_residue as u32, s_residue as u32));
    }

    fn draw_menu(&mut self, ctx: &Context) {
        egui::TopBottomPanel::top("menu").show(ctx, |ui| {
            // Denser menu bar: smaller text, tighter horizontal
            // padding. Doesn't affect submenu styling — those still
            // use the default body size.
            ui.style_mut().spacing.item_spacing.x = 6.0;
            ui.style_mut().spacing.button_padding = Vec2::new(4.0, 2.0);
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Open query FASTA…").clicked() {
                        ui.close_menu();
                        pick_and_load(self, SeqRole::Query);
                    }
                    if ui.button("Open subject FASTA…").clicked() {
                        ui.close_menu();
                        pick_and_load(self, SeqRole::Subject);
                    }
                    ui.separator();
                    if ui.button("Save PNG…").clicked() {
                        ui.close_menu();
                        save_png(self);
                    }
                    if ui.button("Save SVG…").clicked() {
                        ui.close_menu();
                        save_svg(self);
                    }
                    ui.separator();
                    if ui.button("Save session…").clicked() {
                        ui.close_menu();
                        save_session(self);
                    }
                    if ui.button("Open session…").clicked() {
                        ui.close_menu();
                        open_session(self);
                    }
                    ui.separator();
                    if ui.button("Quit").clicked() {
                        ctx.send_viewport_cmd(egui::ViewportCommand::Close);
                    }
                });
                ui.menu_button("View", |ui| {
                    if ui.button("Fit (reset to canonical zoom)").clicked() {
                        self.action_fit();
                    }
                    if ui.button("Back (Esc)").clicked() {
                        self.action_back();
                    }
                    if ui.button("Reset greyramp").clicked() {
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
                    if ui.button(dock_label).clicked() {
                        self.align_dock_visible = !self.align_dock_visible;
                    }
                    ui.separator();
                    // Compact shortcut crib so users can discover
                    // the keyboard map without reading the docs.
                    ui.menu_button("Keyboard shortcuts…", |ui| {
                        let style = egui::FontId::monospace(11.0);
                        for (keys, desc) in [
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
                        ] {
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
                    if ui.button("Settings…").clicked() {
                        self.show_settings = true;
                    }
                });
                ui.separator();
                ui.label(format_seq_label("query", self.query.as_ref()));
                ui.label("·");
                ui.label(format_seq_label("subject", self.subject.as_ref()));
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

                // ── Greyramp ──
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
                let (rect, _) =
                    ui.allocate_exact_size(Vec2::new(ui.available_width(), 20.0), Sense::hover());
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

                // ── Vector overlays ──
                ui.add_space(6.0);
                ui.separator();
                if ui
                    .checkbox(&mut self.settings.show_ridge_overlay, "Show vector ridges")
                    .changed()
                {
                    self.refresh_ridges();
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
            // Any setting that affects the pixelmap invalidates the
            // cache + history (which would point at incompatible
            // pixelmaps under the new parameters).
            self.invalidate_caches();
            self.recompute();
        }
    }

    fn draw_status_bar(&self, ctx: &Context) {
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                if self.compute_in_flight {
                    ui.spinner();
                    ui.label("Recomputing…");
                    ui.separator();
                }
                if let Some(err) = &self.last_error {
                    ui.colored_label(Color32::from_rgb(220, 90, 70), err);
                    return;
                }
                if let Some(plot) = &self.plot {
                    ui.label(format!(
                        "pixelmap {}×{}, W={}",
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
                            let qp = ((q_seq as usize) - q_off) / z;
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
                if self.settings.mode == BlastMode::Blastx {
                    ui.label(
                        "Alignment view: BLASTX (three-frame translated) not yet supported here. \
                         Use --mode blastp or blastn for now.",
                    );
                    return;
                }
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

                let window = self.settings.align_window_size.clamp(20, 400) as usize;
                let half = window / 2;
                let q_bytes = q_seq.bytes();
                let s_bytes = s_seq.bytes();
                let matrix = self.settings.build_matrix();

                // Forward (+/+) alignment columns.
                let forward_block = build_align_columns(
                    q_bytes,
                    s_bytes,
                    q_centre,
                    s_centre,
                    half,
                    window,
                    false,
                    matrix.as_ref(),
                    self.settings.mode,
                );

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

                draw_align_block(ui, ctx, &forward_block);
                if let Some(rev) = &reverse_block {
                    ui.add_space(4.0);
                    draw_align_block(ui, ctx, rev);
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
            let top_margin: f32 = if multi_q { 44.0 } else { 24.0 };
            let left_margin: f32 = if multi_s { 70.0 } else { 56.0 };
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
                            let q_seq = q_off + qp * z + z / 2;
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
                    let pixel_x = (break_coord - q_off) / zoom_us;
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

            // C4: tick labels in the top / left margin bands —
            // outside the plot area so they don't overlap the image.
            self.draw_axis_labels(ui, rect, plot_area, plot, draw_offset, ppp);

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
                let cq_pix = ((cq_seq as f32) - q_off) / z;
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
        let seq_q_lo = q_off + (world_x_lo * zoom_us).max(0.0) as u64;
        let seq_q_hi = q_off + (world_x_hi * zoom_us) as u64;
        let seq_s_lo = s_off + (world_y_lo * zoom_us).max(0.0) as u64;
        let seq_s_hi = s_off + (world_y_hi * zoom_us) as u64;

        let painter = ui.painter();
        let tick_color = Color32::from_rgba_unmultiplied(80, 80, 80, 220);
        let label_color = ui.style().visuals.text_color();
        let font = egui::FontId::monospace(11.0);

        // Top axis (query) — labels above the plot, in the top
        // margin band (outer.top()..plot_area.top()).
        let top_margin_y = outer.top();
        let tick_baseline_y = plot_area.top();
        let label_y = tick_baseline_y - 16.0; // text top
        let span_x = seq_q_hi.saturating_sub(seq_q_lo) as f32;
        // Pixels-per-residue is now expressed against *logical*
        // points: residue → texel via `/ zoom_us`, then texel →
        // logical via `/ ppp`.
        let pixels_per_residue_x = 1.0_f32 / (zoom_us * ppp);
        let step_x = nice_tick_step(span_x as f64, MIN_LABEL_SPACING_PX / pixels_per_residue_x);
        // Minor ticks (unlabeled) at step/5 intervals — gives 4
        // subdivisions between every labelled tick. Drawn first so
        // the longer labelled ticks paint on top.
        let minor_step_x = (step_x / 5.0).floor() as u64;
        if minor_step_x >= 1 {
            let mut t = (seq_q_lo / minor_step_x) * minor_step_x;
            while t < seq_q_hi.saturating_add(minor_step_x) {
                if t >= seq_q_lo && t <= seq_q_hi {
                    let sx =
                        plot_area.left() + (t - q_off) as f32 / (zoom_us * ppp) - draw_offset.x;
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
                let sx = plot_area.left() + (t - q_off) as f32 / (zoom_us * ppp) - draw_offset.x;
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
        let tick_baseline_x = plot_area.left();
        let label_x = tick_baseline_x - 8.0; // text right edge
        let span_y = seq_s_hi.saturating_sub(seq_s_lo) as f32;
        let step_y = nice_tick_step(span_y as f64, MIN_LABEL_SPACING_PX / pixels_per_residue_x);
        let minor_step_y = (step_y / 5.0).floor() as u64;
        if minor_step_y >= 1 {
            let mut t = (seq_s_lo / minor_step_y) * minor_step_y;
            while t < seq_s_hi.saturating_add(minor_step_y) {
                if t >= seq_s_lo && t <= seq_s_hi {
                    let sy =
                        plot_area.top() + (t - s_off) as f32 / (zoom_us * ppp) - draw_offset.y;
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
        crosshair_col: half,
    }
}

/// Render one 3-row alignment block (query, match line, subject)
/// with per-column background colour. Prepends a strand tag (`+/+` or
/// `+/-`) and the row's 5'/3' residue numbers in the side margins,
/// and draws a caret over the crosshair column so the user can see
/// which residue the click landed on.
fn draw_align_block(ui: &mut egui::Ui, ctx: &Context, block: &AlignBlock) {
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
    let strand_label = if block.reverse { "+/-" } else { "+/+" };
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
    let _ = writeln!(out, "+/+");
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
    if seq.records.len() <= 1 {
        return format!("{one_indexed}");
    }
    match seq.record_at(coord) {
        Some((rec, pos)) => format!("{}:{}", rec.id, pos + 1),
        None => format!("{one_indexed}"),
    }
}

fn pick_and_load(app: &mut DottirApp, role: SeqRole) {
    let label = match role {
        SeqRole::Query => "Open query FASTA",
        SeqRole::Subject => "Open subject FASTA",
    };
    let pick = rfd::FileDialog::new()
        .set_title(label)
        .add_filter("FASTA", &["fa", "fasta", "fna", "faa", "gz"])
        .pick_file();
    if let Some(path) = pick {
        app.load_fasta(role, path);
    }
}

fn save_svg(app: &mut DottirApp) {
    let Some(plot) = &app.plot else {
        app.last_error = Some("nothing to save — compute a plot first".into());
        return;
    };
    let lut = app.greyramp.lut();
    let pick = rfd::FileDialog::new()
        .set_title("Save SVG")
        .add_filter("SVG", &["svg"])
        .save_file();
    if let Some(path) = pick {
        let mapped: Vec<u8> = plot.pixels.iter().map(|&v| lut[v as usize]).collect();
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
        match dottir_io::svg_export::write_svg(
            &path,
            plot.width,
            plot.height,
            &mapped,
            50,
            &recs_x,
            &recs_y,
            &[
                ("dottir-gui", env!("CARGO_PKG_VERSION")),
                ("greyramp-white", &app.greyramp.white.to_string()),
                ("greyramp-black", &app.greyramp.black.to_string()),
            ],
        ) {
            Ok(()) => {}
            Err(e) => app.last_error = Some(format!("SVG save failed: {e}")),
        }
    }
}

fn save_session(app: &mut DottirApp) {
    use crate::session::{
        codec, Session, SessionGreyramp, SessionPlot, SessionView, SESSION_VERSION,
    };
    let session = Session {
        version: SESSION_VERSION,
        query: app.query.as_ref().and_then(|s| s.source_path.clone()),
        subject: app.subject.as_ref().and_then(|s| s.source_path.clone()),
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

fn open_session(app: &mut DottirApp) {
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
        app.load_fasta(SeqRole::Query, p);
    } else {
        app.query = None;
        app.plot = None;
    }
    if let Some(p) = s.subject {
        app.load_fasta(SeqRole::Subject, p);
    } else {
        app.subject = None;
        app.plot = None;
    }
    // Apply view state AFTER the loads (load_fasta calls recompute
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
    tracing::info!("loaded session {}", path.display());
}

fn save_png(app: &mut DottirApp) {
    let Some(plot) = &app.plot else {
        app.last_error = Some("nothing to save — compute a plot first".into());
        return;
    };
    let lut = app.greyramp.lut();
    let pick = rfd::FileDialog::new()
        .set_title("Save PNG")
        .add_filter("PNG", &["png"])
        .save_file();
    if let Some(path) = pick {
        // Apply the greyramp LUT before saving so the on-disk image
        // matches what's on screen.
        let mapped: Vec<u8> = plot.pixels.iter().map(|&v| lut[v as usize]).collect();
        match dottir_io::png_export::write_grayscale_png(
            &path,
            plot.width,
            plot.height,
            &mapped,
            &[
                ("dottir-gui", env!("CARGO_PKG_VERSION")),
                ("greyramp-white", &app.greyramp.white.to_string()),
                ("greyramp-black", &app.greyramp.black.to_string()),
            ],
        ) {
            Ok(()) => {}
            Err(e) => app.last_error = Some(format!("PNG save failed: {e}")),
        }
    }
}
