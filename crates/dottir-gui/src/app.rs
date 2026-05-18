//! Top-level egui application.

use std::path::PathBuf;

use dottir_core::{BlastMode, DotPlot, PlotConfig, ScoreMatrix, Strand, Triangle};
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
            pixel_fac: 50,
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
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            mode: BlastMode::Blastn,
            matrix_name: "DNA+5/-4".into(),
            window_size: None,
            zoom: 1,
            pixel_fac: 50,
            strand: Strand::Both,
            self_comparison: false,
            triangle: Triangle::Both,
            memory_limit_bytes: 512 * 1024 * 1024,
            reverse_query: false,
            reverse_subject: false,
            inverted_repeat_colour: false,
            align_window_size: 100,
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
        // C dotter defaults per spec §4.2.2: white=40, black=100.
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
    /// View transform: top-left of the canvas in pixelmap coords.
    view_offset: Vec2,
    /// Pixels-per-pixelmap-pixel zoom (display zoom, separate from
    /// PlotConfig.zoom which is *computation* zoom).
    display_zoom: f32,
    /// Crosshair in pixelmap coords (q, s).
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
    /// Timestamp of the most recent scroll-wheel zoom event, used by
    /// the C2 zoom-settle recompute: when no scroll has fired for
    /// 200ms and `display_zoom` is far enough from 1.0 to warrant
    /// finer-grained computation, the kernel re-runs at a new
    /// `PlotConfig::zoom` tier.
    last_zoom_event: Option<std::time::Instant>,
    /// H2: whether the alignment-view dock is shown beneath the
    /// canvas. Default true; toggled via View → "Show alignment
    /// view".
    align_dock_visible: bool,
    /// False until the canvas has rendered its first frame for the
    /// current `plot`; on that frame the view is snapped to
    /// fit-to-canvas. Cleared whenever a new pixelmap arrives so
    /// every freshly loaded plot opens fully zoomed out.
    view_initialised: bool,
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

        let settings = Settings {
            mode: startup.mode,
            matrix_name: startup.matrix_name.clone(),
            window_size: startup.window_size,
            zoom: startup.zoom.max(1),
            pixel_fac: startup.pixel_fac.max(1),
            strand: startup.strand,
            self_comparison: startup.self_comparison,
            triangle: Triangle::Both,
            memory_limit_bytes: startup.memory_limit_bytes,
            reverse_query: false,
            reverse_subject: false,
            inverted_repeat_colour: false,
            align_window_size: 100,
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
            display_zoom: 1.0,
            crosshair: None,
            show_settings: false,
            light_theme: true,
            worker,
            last_dispatched_id: 0,
            compute_in_flight: false,
            // Suspend recompute during the two pre-loads — we only
            // want one job dispatched once both inputs are settled.
            suspend_recompute: true,
            last_zoom_event: None,
            align_dock_visible: true,
            view_initialised: false,
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
        // Release the suspend and dispatch one initial compute if
        // both inputs are present.
        app.suspend_recompute = false;
        if app.query.is_some() && app.subject.is_some() {
            app.recompute();
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
                // Fresh sequence data → reset view so the new plot
                // opens fully zoomed out on its first render.
                self.view_initialised = false;
                self.maybe_switch_mode_from_alphabet(detected);
                self.recompute();
            }
            Err(e) => {
                self.last_error = Some(format!("failed to load {}: {e}", path.display()));
            }
        }
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
        let (Some(q), Some(s)) = (&self.query, &self.subject) else {
            self.plot = None;
            self.texture = None;
            self.compute_in_flight = false;
            return;
        };
        let matrix = match self.settings.build_matrix() {
            Some(m) => m,
            None => {
                self.last_error = Some(format!("unknown matrix '{}'", self.settings.matrix_name));
                return;
            }
        };
        let mut cfg = PlotConfig {
            mode: self.settings.mode,
            matrix,
            window_size: self.settings.window_size,
            zoom: self.settings.zoom.max(1),
            pixel_fac: self.settings.pixel_fac.max(1),
            strand: self.settings.strand,
            self_comparison: self.settings.self_comparison,
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
        let req = crate::compute_worker::ComputeRequest {
            id,
            // Sequences are cloned per request — typically tens of
            // MiB at worst. The clone is dwarfed by the compute time
            // it triggers, and lets the worker outlive the borrow.
            query: q.bytes().to_vec(),
            subject: s.bytes().to_vec(),
            config: cfg,
        };
        tracing::info!("dispatch compute id={id}");
        self.worker.dispatch(req);
        self.compute_in_flight = true;
    }

    /// C2 zoom-settle recompute: when the user has stopped scrolling
    /// for [`ZOOM_SETTLE_MS`] and `display_zoom` has strayed outside
    /// `[0.5, 2.0]`, swap to a finer (or coarser) `PlotConfig::zoom`
    /// tier and recompute. Resets `display_zoom` back to 1.0 of the
    /// new pixelmap and scales `view_offset` + crosshair so the same
    /// residue stays under the same screen pixel.
    fn maybe_zoom_settle_recompute(&mut self, ctx: &Context) {
        const ZOOM_SETTLE_MS: u64 = 200;
        const ZOOM_IN_THRESHOLD: f32 = 2.0;
        const ZOOM_OUT_THRESHOLD: f32 = 0.5;

        let Some(t) = self.last_zoom_event else {
            return;
        };
        if self.compute_in_flight {
            // Wait for the in-flight job to land before scheduling
            // another tier change.
            ctx.request_repaint_after(std::time::Duration::from_millis(ZOOM_SETTLE_MS));
            return;
        }
        if t.elapsed().as_millis() < ZOOM_SETTLE_MS as u128 {
            // Schedule another frame so we re-check after the
            // settle interval, even if no other event fires.
            ctx.request_repaint_after(std::time::Duration::from_millis(ZOOM_SETTLE_MS));
            return;
        }

        let current_zoom = self.settings.zoom.max(1);
        let (new_zoom, scale) = if self.display_zoom >= ZOOM_IN_THRESHOLD && current_zoom > 1 {
            // Zoomed in past 2× — recompute at a finer tier
            // (half the computation zoom, doubling pixelmap density).
            let n = (current_zoom / 2).max(1);
            let s = current_zoom as f32 / n as f32; // > 1, e.g. 2.0
            (n, s)
        } else if self.display_zoom <= ZOOM_OUT_THRESHOLD && current_zoom < 64 {
            // Zoomed out past 0.5× — recompute at a coarser tier
            // (double the computation zoom, halving pixelmap density).
            let n = (current_zoom.saturating_mul(2)).min(64);
            let s = current_zoom as f32 / n as f32; // < 1, e.g. 0.5
            (n, s)
        } else {
            // Within [0.5, 2.0] or already at a tier extreme: clear
            // the timestamp and stop chasing.
            self.last_zoom_event = None;
            return;
        };

        tracing::info!(
            "zoom-settle: tier change {current_zoom} → {new_zoom} \
             (display_zoom {:.2} → 1.0, scale {scale})",
            self.display_zoom,
        );
        self.settings.zoom = new_zoom;
        // Rescale view state so the same residue stays under the
        // same screen position. New pixelmap dims = old × (1/scale),
        // so view_offset (in pixelmap coords) scales by 1/scale.
        let inv = 1.0 / scale;
        self.view_offset *= inv;
        self.display_zoom = 1.0;
        if let Some((cq, cs)) = self.crosshair {
            // Crosshair is in pixelmap coords; scale identically.
            let cq2 = (cq as f32 * inv).round() as u32;
            let cs2 = (cs as f32 * inv).round() as u32;
            self.crosshair = Some((cq2, cs2));
        }
        self.last_zoom_event = None;
        self.recompute();
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
                        "computed {}×{} pixelmap (W={}, zoom={})",
                        plot.width,
                        plot.height,
                        plot.params.window_size,
                        r.config_zoom,
                    );
                    self.plot = Some(plot);
                    self.texture_dirty = true;
                    self.last_error = None;
                }
                Err(e) => {
                    self.last_error = Some(format!("compute_dotplot failed: {e}"));
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
        // Spec §4.4.3 inverted-repeat highlighting: when both
        // forward+reverse channels are populated, paint forward in
        // grey and reverse in magenta (overlapping cells take
        // whichever channel is stronger after the greyramp).
        let mut rgba = Vec::with_capacity(plot.pixels.len() * 4);
        match (
            plot.forward_pixels.as_deref(),
            plot.reverse_pixels.as_deref(),
        ) {
            (Some(fwd), Some(rev)) if self.settings.inverted_repeat_colour => {
                for i in 0..plot.pixels.len() {
                    let f = lut[fwd[i] as usize]; // forward strength (0..255)
                    let r = lut[rev[i] as usize]; // reverse strength
                                                  // Forward → black on white (so darker means more
                                                  // confident); reverse channel uses magenta (255,
                                                  // 0, 255).
                                                  // Blend: take the channel with the highest "ink"
                                                  // (= 255 - lut value, since white = 255 means
                                                  // "no ink"). If forward wins, render greyscale; if
                                                  // reverse wins, render magenta-tinted.
                    let f_ink = 255 - f;
                    let r_ink = 255 - r;
                    let (cr, cg, cb) = if f_ink >= r_ink {
                        (f, f, f) // greyscale forward
                    } else {
                        // Reverse hit. Render as magenta whose
                        // intensity scales with r_ink.
                        let ink = r_ink as u16;
                        let bg = 255_u16 - ink; // white background
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
                // Plain greyscale: combined channel through the LUT.
                for &v in &plot.pixels {
                    let g = lut[v as usize];
                    rgba.extend_from_slice(&[g, g, g, 255]);
                }
            }
        }
        let image =
            ColorImage::from_rgba_unmultiplied([plot.width as usize, plot.height as usize], &rgba);
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
        // Drain any completed background-compute results before any
        // panel reads `self.plot`. Stale results are discarded inside.
        self.poll_compute_results();
        // Once the wheel has been idle for the debounce interval,
        // consider recomputing at a finer/coarser tier.
        self.maybe_zoom_settle_recompute(ctx);

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
    }
}

impl DottirApp {
    fn handle_keyboard(&mut self, ctx: &Context) {
        let mods = ctx.input(|i| i.modifiers);
        let step = if mods.ctrl {
            100_i32
        } else if mods.shift {
            10
        } else {
            1
        };
        let Some(plot) = &self.plot else {
            return;
        };
        let pw = plot.width as i64;
        let ph = plot.height as i64;
        let mut nudged = false;
        let mut snap = false;
        let (mut q, mut s) = self.crosshair.unwrap_or((plot.width / 2, plot.height / 2));
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
            q = q_i.clamp(0, pw - 1) as u32;
            s = s_i.clamp(0, ph - 1) as u32;
            self.crosshair = Some((q, s));
        }
        if snap {
            self.snap_crosshair_to_line();
        }
    }

    /// Snap the crosshair to the brightest pixel within a search
    /// disc, with current position as the tie-breaker (closest
    /// pixel wins on equal value). Bound to **Space**.
    ///
    /// This is the simplest design that actually does what users
    /// want: most strong dots on a dotplot belong to a diagonal
    /// run, so "grab the local max" lands you on the nearest line.
    /// A more sophisticated version would aggregate scores per
    /// diagonal offset, but the simple local max is robust and
    /// fast enough for arbitrarily large pixelmaps (we cap the
    /// search at ±64 pixelmap pixels around the current
    /// crosshair).
    fn snap_crosshair_to_line(&mut self) {
        let Some(plot) = self.plot.as_ref() else {
            return;
        };
        let Some((cq, cs)) = self.crosshair else {
            return;
        };
        const RADIUS: i64 = 64;
        let pw = plot.width as i64;
        let ph = plot.height as i64;
        let stride = plot.width as usize;
        let q_lo = (cq as i64 - RADIUS).max(0);
        let q_hi = (cq as i64 + RADIUS).min(pw - 1);
        let s_lo = (cs as i64 - RADIUS).max(0);
        let s_hi = (cs as i64 + RADIUS).min(ph - 1);
        // Track the best (value, -distance², q, s) tuple. Tuple
        // ordering gives us "max value, then min distance" for free.
        let mut best: Option<(u8, i64, u32, u32)> = None;
        for sp in s_lo..=s_hi {
            let row = sp as usize * stride;
            for qp in q_lo..=q_hi {
                let v = plot.pixels[row + qp as usize];
                if v == 0 {
                    continue;
                }
                let dq = qp - cq as i64;
                let ds = sp - cs as i64;
                let dist_sq = dq * dq + ds * ds;
                let candidate = (v, -dist_sq, qp as u32, sp as u32);
                match best {
                    None => best = Some(candidate),
                    Some(cur) if (candidate.0, candidate.1) > (cur.0, cur.1) => {
                        best = Some(candidate)
                    }
                    _ => {}
                }
            }
        }
        if let Some((_, _, q, s)) = best {
            self.crosshair = Some((q, s));
        }
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
                    if ui.button("Reset pan / zoom").clicked() {
                        // Defer to the canvas's first-frame logic:
                        // it'll snap to fit-zoom on the next paint.
                        self.view_initialised = false;
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
                            ("← → ↑ ↓", "nudge crosshair (Shift ×10, Ctrl ×100)"),
                            (",   .", "step along main diagonal"),
                            ("[   ]", "step along anti-diagonal"),
                            ("Space", "snap crosshair to nearest strong dot"),
                            ("scroll", "zoom on cursor"),
                            ("drag", "pan"),
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
                });
                ui.horizontal(|ui| {
                    ui.label("Pixel factor:");
                    if ui
                        .add(Slider::new(&mut self.settings.pixel_fac, 1..=255))
                        .changed()
                    {
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
                if ui.button("Apply").clicked() {
                    changed = true;
                }
            });
        self.show_settings = open;
        if changed {
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
                    if let Some((q, s)) = self.crosshair {
                        let idx = (s as usize) * (plot.width as usize) + (q as usize);
                        let v = plot.pixels.get(idx).copied().unwrap_or(0);
                        // Map pixelmap coord → sequence coord (zoom).
                        let z = plot.params.zoom as usize;
                        let q_seq = (q as usize) * z + z / 2;
                        let s_seq = (s as usize) * z + z / 2;
                        ui.separator();
                        ui.label(format!(
                            "q = {}, s = {}, value = {}",
                            format_coord(self.query.as_ref(), q_seq),
                            format_coord(self.subject.as_ref(), s_seq),
                            v,
                        ));
                    } else {
                        ui.separator();
                        ui.label("click on the plot to set the crosshair");
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
                let Some((cq_pix, cs_pix)) = self.crosshair else {
                    ui.label("Alignment view: click on the plot to set the crosshair.");
                    return;
                };
                let Some(q_seq) = self.query.as_ref() else {
                    return;
                };
                let Some(s_seq) = self.subject.as_ref() else {
                    return;
                };

                // Centre coordinates in *residue* space (translate
                // pixelmap coords back through the kernel zoom).
                let kzoom = plot.params.zoom.max(1) as usize;
                let q_centre = (cq_pix as usize) * kzoom + kzoom / 2;
                let s_centre = (cs_pix as usize) * kzoom + kzoom / 2;

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
            let Some(plot) = &self.plot else {
                ui.centered_and_justified(|ui| {
                    ui.label("No plot. Load a query and subject FASTA (File menu).");
                });
                return;
            };
            let Some(tex) = &self.texture else {
                return;
            };
            let avail = ui.available_size();
            let (rect, response) = ui.allocate_exact_size(avail, Sense::click_and_drag());

            // Reserve a margin band along the top + left for axis
            // tick labels and record-name labels. The pixelmap and
            // every world-coord overlay clip to `plot_rect`, which
            // is `rect` with these margins removed. Multi-record
            // FASTAs need extra space for the record-name label
            // strip above the tick labels.
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
            let plot_area = Rect::from_min_max(
                Pos2::new(rect.left() + left_margin, rect.top() + top_margin),
                rect.right_bottom(),
            );

            // Compute fit-to-canvas zoom first — used both as the
            // lower bound for interactive zoom and as the initial
            // zoom level for freshly loaded plots.
            let pw = plot.width as f32;
            let ph = plot.height as f32;
            let plot_w = plot_area.width().max(1.0);
            let plot_h = plot_area.height().max(1.0);
            let fit_zoom_x = plot_w / pw.max(1.0);
            let fit_zoom_y = plot_h / ph.max(1.0);
            let fit_zoom = fit_zoom_x.min(fit_zoom_y);

            // First render of a freshly loaded plot: snap to fit so
            // the user sees the whole pixelmap.
            if !self.view_initialised {
                self.display_zoom = fit_zoom;
                self.view_offset = Vec2::ZERO;
                self.view_initialised = true;
            }

            // Pan with primary drag — only inside the plot area
            // (drags into the axis margin don't pan).
            if response.dragged() {
                self.view_offset -= response.drag_delta() / self.display_zoom;
            }
            // Zoom with scroll, centered on cursor. Clamped against
            // `[fit_zoom, 100.0]` *before* the anchor math so that
            // hitting the floor doesn't re-bump `display_zoom`
            // afterwards and drift the anchor.
            let scroll = ctx.input(|i| i.smooth_scroll_delta.y);
            if scroll.abs() > 0.0 {
                if let Some(cursor) = response.hover_pos() {
                    if plot_area.contains(cursor) {
                        let factor = (scroll / 100.0).exp();
                        let new_zoom =
                            (self.display_zoom * factor).clamp(fit_zoom, 100.0);
                        let cursor_local = cursor - plot_area.left_top();
                        let world_under_cursor =
                            self.view_offset + cursor_local / self.display_zoom;
                        self.display_zoom = new_zoom;
                        self.view_offset = world_under_cursor - cursor_local / self.display_zoom;
                        // Stamp for the C2 zoom-settle recompute trigger.
                        self.last_zoom_event = Some(std::time::Instant::now());
                    }
                }
            }
            // Double-click: zoom 2× centred on the click position,
            // same anchor math as the wheel. Single-click still
            // fires before this and sets the crosshair, so a
            // double-click both moves the crosshair to the cursor
            // and zooms in around it.
            if response.double_clicked() {
                if let Some(p) = response.interact_pointer_pos() {
                    if plot_area.contains(p) {
                        let new_zoom = (self.display_zoom * 2.0).clamp(fit_zoom, 100.0);
                        let local = p - plot_area.left_top();
                        let world_under_cursor =
                            self.view_offset + local / self.display_zoom;
                        self.display_zoom = new_zoom;
                        self.view_offset = world_under_cursor - local / self.display_zoom;
                        self.last_zoom_event = Some(std::time::Instant::now());
                    }
                }
            }

            // Constrain pan so the user can't scroll the pixelmap
            // completely off-screen. With `display_zoom` already
            // clamped at `fit_zoom`, the plot is always at least as
            // big as the canvas in one axis; the other axis (if
            // aspect ratios differ) gets auto-centred.
            let cw = plot_w / self.display_zoom; // canvas width in pixelmap coords
            let ch = plot_h / self.display_zoom;
            if pw <= cw {
                self.view_offset.x = -(cw - pw) / 2.0; // centre horizontally
            } else {
                self.view_offset.x = self.view_offset.x.clamp(0.0, pw - cw);
            }
            if ph <= ch {
                self.view_offset.y = -(ch - ph) / 2.0;
            } else {
                self.view_offset.y = self.view_offset.y.clamp(0.0, ph - ch);
            }

            // Click sets the crosshair (clicks in the margin are
            // ignored).
            if response.clicked() {
                if let Some(p) = response.interact_pointer_pos() {
                    if plot_area.contains(p) {
                        let local = p - plot_area.left_top();
                        let world = self.view_offset + local / self.display_zoom;
                        let q = world.x.floor() as i64;
                        let s = world.y.floor() as i64;
                        if q >= 0 && q < plot.width as i64 && s >= 0 && s < plot.height as i64 {
                            self.crosshair = Some((q as u32, s as u32));
                        }
                    }
                }
            }

            // Fill the whole canvas (margins + plot area) with light
            // grey first. The margin band is left as-is; the plot
            // area gets the texture painted over it.
            ui.painter().rect_filled(rect, 0.0, Color32::from_gray(235));

            // Compute on-screen rect for the pixelmap, then render
            // it. Clip to `plot_area` so over-pan doesn't leak the
            // texture into the margin.
            let plot_screen_w = pw * self.display_zoom;
            let plot_screen_h = ph * self.display_zoom;
            let plot_screen_x = plot_area.left() - self.view_offset.x * self.display_zoom;
            let plot_screen_y = plot_area.top() - self.view_offset.y * self.display_zoom;
            let plot_rect = Rect::from_min_size(
                Pos2::new(plot_screen_x, plot_screen_y),
                Vec2::new(plot_screen_w, plot_screen_h),
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

            // C3: breaklines for multi-record FASTA inputs. Vertical
            // lines at the query record boundaries; horizontal lines
            // at the subject record boundaries. Drawn underneath the
            // crosshair so it stays visible.
            let zoom_us = plot.params.zoom.max(1) as usize;
            let break_stroke =
                egui::Stroke::new(1.0, Color32::from_rgba_unmultiplied(80, 160, 80, 220));
            if let Some(q_seq) = &self.query {
                for &break_coord in &q_seq.breaks() {
                    let pixel_x = break_coord / zoom_us;
                    if pixel_x >= plot.width as usize {
                        continue;
                    }
                    let sx = plot_area.left()
                        + ((pixel_x as f32) - self.view_offset.x) * self.display_zoom;
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
                    let pixel_y = break_coord / zoom_us;
                    if pixel_y >= plot.height as usize {
                        continue;
                    }
                    let sy = plot_area.top()
                        + ((pixel_y as f32) - self.view_offset.y) * self.display_zoom;
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

            // C4: tick labels in the top / left margin bands —
            // outside the plot area so they don't overlap the image.
            self.draw_axis_labels(ui, rect, plot_area, plot);

            // Crosshair overlay + coord label — clipped to the plot
            // area so the lines never run into the axis margin.
            if let Some((cq, cs)) = self.crosshair {
                let cx =
                    plot_area.left() + ((cq as f32 + 0.5) - self.view_offset.x) * self.display_zoom;
                let cy =
                    plot_area.top() + ((cs as f32 + 0.5) - self.view_offset.y) * self.display_zoom;
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

                // Coord label next to the cross.
                let z = plot.params.zoom.max(1) as usize;
                let q_seq = (cq as usize) * z + z / 2;
                let s_seq = (cs as usize) * z + z / 2;
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
    ) {
        const MIN_LABEL_SPACING_PX: f32 = 80.0;
        let zoom_us = plot.params.zoom.max(1) as f32;
        // World-pixel range visible inside the plot area.
        let world_x_lo = self.view_offset.x;
        let world_x_hi = world_x_lo + plot_area.width() / self.display_zoom;
        let world_y_lo = self.view_offset.y;
        let world_y_hi = world_y_lo + plot_area.height() / self.display_zoom;
        // Convert to sequence-residue range.
        let seq_q_lo = (world_x_lo * zoom_us).max(0.0) as u64;
        let seq_q_hi = (world_x_hi * zoom_us) as u64;
        let seq_s_lo = (world_y_lo * zoom_us).max(0.0) as u64;
        let seq_s_hi = (world_y_hi * zoom_us) as u64;

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
        let pixels_per_residue_x = self.display_zoom / zoom_us;
        let step_x = nice_tick_step(span_x as f64, MIN_LABEL_SPACING_PX / pixels_per_residue_x);
        let mut t = (seq_q_lo / step_x as u64) * step_x as u64;
        while t < seq_q_hi.saturating_add(step_x as u64) {
            if t >= seq_q_lo && t <= seq_q_hi {
                let sx = plot_area.left()
                    + (t as f32 / zoom_us - self.view_offset.x) * self.display_zoom;
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
                    format_kb(t),
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
        let mut t = (seq_s_lo / step_y as u64) * step_y as u64;
        while t < seq_s_hi.saturating_add(step_y as u64) {
            if t >= seq_s_lo && t <= seq_s_hi {
                let sy =
                    plot_area.top() + (t as f32 / zoom_us - self.view_offset.y) * self.display_zoom;
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
                    format_kb(t),
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
                    let r_start = rec.range.start as u64;
                    let r_end = rec.range.end as u64;
                    if r_end <= r_start {
                        continue;
                    }
                    let x0 = plot_area.left()
                        + (r_start as f32 / zoom_us - self.view_offset.x) * self.display_zoom;
                    let x1 = plot_area.left()
                        + (r_end as f32 / zoom_us - self.view_offset.x) * self.display_zoom;
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
                    let r_start = rec.range.start as u64;
                    let r_end = rec.range.end as u64;
                    if r_end <= r_start {
                        continue;
                    }
                    let y0 = plot_area.top()
                        + (r_start as f32 / zoom_us - self.view_offset.y) * self.display_zoom;
                    let y1 = plot_area.top()
                        + (r_end as f32 / zoom_us - self.view_offset.y) * self.display_zoom;
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

/// Format a residue coord with a `kb`/`Mb` suffix when large.
fn format_kb(n: u64) -> String {
    if n >= 1_000_000 {
        format!("{:.1}M", n as f64 / 1_000_000.0)
    } else if n >= 1_000 {
        format!("{}k", n / 1_000)
    } else {
        format!("{n}")
    }
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
fn format_coord(seq: Option<&Sequence>, coord: usize) -> String {
    let Some(seq) = seq else {
        return format!("{coord}");
    };
    if seq.records.len() <= 1 {
        return format!("{coord}");
    }
    match seq.record_at(coord) {
        Some((rec, pos)) => format!("{}:{}", rec.id, pos + 1),
        None => format!("{coord}"),
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
            pixel_fac: app.settings.pixel_fac,
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
            display_zoom: app.display_zoom,
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
    app.settings.pixel_fac = s.plot.pixel_fac.max(1);
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
    app.display_zoom = s.view.display_zoom.max(0.1);
    app.crosshair = s.view.crosshair.map(|[q, s]| (q, s));
    app.light_theme = s.view.light_theme;
    app.align_dock_visible = s.view.align_dock_visible;
    // Session restore overrides the load_fasta-driven view reset:
    // the saved transform is authoritative, don't snap-to-fit.
    app.view_initialised = true;
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
