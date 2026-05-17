//! Top-level egui application.

use std::path::PathBuf;

use dottir_core::{
    BlastMode, DotPlot, PlotConfig, ScoreMatrix, Strand, Triangle,
};
use dottir_io::Sequence;
use egui::{
    Color32, ColorImage, Context, Pos2, Rect, Sense, Slider, TextureHandle,
    TextureOptions, Vec2,
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
                if self.swap { 0 } else { 255 }
            } else if i >= hi {
                if self.swap { 255 } else { 0 }
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
                tracing::info!(
                    "loaded {} ({} residues, {} records)",
                    path.display(),
                    seq.len(),
                    seq.records.len(),
                );
                match role {
                    SeqRole::Query => self.query = Some(seq),
                    SeqRole::Subject => self.subject = Some(seq),
                }
                self.last_error = None;
                self.recompute();
            }
            Err(e) => {
                self.last_error =
                    Some(format!("failed to load {}: {e}", path.display()));
            }
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
                self.last_error = Some(format!(
                    "unknown matrix '{}'",
                    self.settings.matrix_name
                ));
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
            ctx.request_repaint_after(std::time::Duration::from_millis(
                ZOOM_SETTLE_MS,
            ));
            return;
        }
        if t.elapsed().as_millis() < ZOOM_SETTLE_MS as u128 {
            // Schedule another frame so we re-check after the
            // settle interval, even if no other event fires.
            ctx.request_repaint_after(std::time::Duration::from_millis(
                ZOOM_SETTLE_MS,
            ));
            return;
        }

        let current_zoom = self.settings.zoom.max(1);
        let (new_zoom, scale) = if self.display_zoom >= ZOOM_IN_THRESHOLD
            && current_zoom > 1
        {
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
        self.view_offset = self.view_offset * inv;
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
                    self.last_error =
                        Some(format!("compute_dotplot failed: {e}"));
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
        match (plot.forward_pixels.as_deref(), plot.reverse_pixels.as_deref()) {
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
        let image = ColorImage::from_rgba_unmultiplied(
            [plot.width as usize, plot.height as usize],
            &rgba,
        );
        let handle = ctx.load_texture(
            "dottir-pixelmap",
            image,
            TextureOptions::NEAREST,
        );
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
        let mut nudged = false;
        let (mut q, mut s) = self.crosshair.unwrap_or((plot.width / 2, plot.height / 2));
        ctx.input(|i| {
            for ev in &i.events {
                if let egui::Event::Key { key, pressed: true, .. } = ev {
                    match key {
                        egui::Key::ArrowLeft => {
                            q = q.saturating_sub(step as u32);
                            nudged = true;
                        }
                        egui::Key::ArrowRight => {
                            q = (q as i64 + step as i64)
                                .clamp(0, plot.width as i64 - 1) as u32;
                            nudged = true;
                        }
                        egui::Key::ArrowUp => {
                            s = s.saturating_sub(step as u32);
                            nudged = true;
                        }
                        egui::Key::ArrowDown => {
                            s = (s as i64 + step as i64)
                                .clamp(0, plot.height as i64 - 1) as u32;
                            nudged = true;
                        }
                        _ => {}
                    }
                }
            }
        });
        if nudged {
            self.crosshair = Some((q, s));
        }
    }

    fn draw_menu(&mut self, ctx: &Context) {
        egui::TopBottomPanel::top("menu").show(ctx, |ui| {
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
                        self.view_offset = Vec2::ZERO;
                        self.display_zoom = 1.0;
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
        egui::SidePanel::right("greyramp")
            .resizable(false)
            .default_width(220.0)
            .show(ctx, |ui| {
                ui.heading("Greyramp");
                ui.label("White point");
                if ui
                    .add(Slider::new(&mut self.greyramp.white, 0..=255).clamping(egui::SliderClamping::Always))
                    .changed()
                {
                    self.texture_dirty = true;
                }
                ui.label("Black point");
                if ui
                    .add(Slider::new(&mut self.greyramp.black, 0..=255).clamping(egui::SliderClamping::Always))
                    .changed()
                {
                    self.texture_dirty = true;
                }
                ui.horizontal(|ui| {
                    if ui.button("Swap").clicked() {
                        self.greyramp.swap = !self.greyramp.swap;
                        self.texture_dirty = true;
                    }
                    if ui.button("Reset").clicked() {
                        self.greyramp = Greyramp::default();
                        self.texture_dirty = true;
                    }
                });
                ui.separator();
                ui.label("LUT preview");
                let lut = self.greyramp.lut();
                let (rect, _) = ui.allocate_exact_size(
                    Vec2::new(ui.available_width(), 24.0),
                    Sense::hover(),
                );
                let painter = ui.painter();
                for x in 0..256 {
                    let g = lut[x];
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
                                    "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80",
                                    "BLOSUM90", "PAM30", "PAM70", "PAM250",
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
                                .selectable_value(
                                    &mut self.settings.strand,
                                    val,
                                    label,
                                )
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
                                .selectable_value(
                                    &mut self.settings.triangle,
                                    val,
                                    label,
                                )
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
                ui.weak(format!(
                    "Refuses to allocate a pixelmap larger than this. \
                     1 GiB suits ~32k × 32k at zoom 1. Halve the cap when \
                     zoom doubles."
                ));
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
                        plot.width,
                        plot.height,
                        plot.params.window_size
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

            // Pan with primary drag.
            if response.dragged() {
                self.view_offset -= response.drag_delta() / self.display_zoom;
            }
            // Zoom with scroll, centered on cursor.
            let scroll = ctx.input(|i| i.smooth_scroll_delta.y);
            if scroll.abs() > 0.0 {
                if let Some(cursor) = response.hover_pos() {
                    let factor = (scroll / 100.0).exp();
                    let cursor_local = cursor - rect.left_top();
                    let world_under_cursor =
                        self.view_offset + cursor_local / self.display_zoom;
                    self.display_zoom = (self.display_zoom * factor).clamp(0.1, 100.0);
                    self.view_offset =
                        world_under_cursor - cursor_local / self.display_zoom;
                    // Stamp for the C2 zoom-settle recompute trigger.
                    self.last_zoom_event = Some(std::time::Instant::now());
                }
            }
            // Click sets the crosshair.
            if response.clicked() {
                if let Some(p) = response.interact_pointer_pos() {
                    let local = p - rect.left_top();
                    let world = self.view_offset + local / self.display_zoom;
                    let q = world.x.floor() as i64;
                    let s = world.y.floor() as i64;
                    if q >= 0
                        && q < plot.width as i64
                        && s >= 0
                        && s < plot.height as i64
                    {
                        self.crosshair = Some((q as u32, s as u32));
                    }
                }
            }

            // Map view rect → texture UV rect.
            let world_w = avail.x / self.display_zoom;
            let world_h = avail.y / self.display_zoom;
            let u0 = self.view_offset.x / plot.width as f32;
            let v0 = self.view_offset.y / plot.height as f32;
            let u1 = (self.view_offset.x + world_w) / plot.width as f32;
            let v1 = (self.view_offset.y + world_h) / plot.height as f32;
            let uv = Rect::from_min_max(Pos2::new(u0, v0), Pos2::new(u1, v1));
            ui.painter().image(tex.id(), rect, uv, Color32::WHITE);

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
                    let sx = rect.left()
                        + ((pixel_x as f32) - self.view_offset.x) * self.display_zoom;
                    if sx < rect.left() || sx > rect.right() {
                        continue;
                    }
                    ui.painter().line_segment(
                        [Pos2::new(sx, rect.top()), Pos2::new(sx, rect.bottom())],
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
                    let sy = rect.top()
                        + ((pixel_y as f32) - self.view_offset.y) * self.display_zoom;
                    if sy < rect.top() || sy > rect.bottom() {
                        continue;
                    }
                    ui.painter().line_segment(
                        [Pos2::new(rect.left(), sy), Pos2::new(rect.right(), sy)],
                        break_stroke,
                    );
                }
            }

            // C4: tick labels along the canvas edges in sequence
            // coordinates. Anchored to sequence coords so they stay
            // crisp under viewport zoom.
            self.draw_axis_labels(ui, rect, plot);

            // Crosshair overlay + coord label.
            if let Some((cq, cs)) = self.crosshair {
                let cx = rect.left()
                    + ((cq as f32 + 0.5) - self.view_offset.x) * self.display_zoom;
                let cy = rect.top()
                    + ((cs as f32 + 0.5) - self.view_offset.y) * self.display_zoom;
                let stroke = egui::Stroke::new(1.0, Color32::from_rgb(255, 80, 80));
                ui.painter().line_segment(
                    [Pos2::new(rect.left(), cy), Pos2::new(rect.right(), cy)],
                    stroke,
                );
                ui.painter().line_segment(
                    [Pos2::new(cx, rect.top()), Pos2::new(cx, rect.bottom())],
                    stroke,
                );

                // Coord label next to the cross. Shows sequence
                // coords; for multi-record FASTAs the helper renders
                // `record_id:pos`. Placed bottom-right of the cross
                // by default, or flipped to the left/top side if
                // it'd run off the canvas edge.
                let z = plot.params.zoom.max(1) as usize;
                let q_seq = (cq as usize) * z + z / 2;
                let s_seq = (cs as usize) * z + z / 2;
                let label = format!(
                    "q = {}, s = {}",
                    format_coord(self.query.as_ref(), q_seq),
                    format_coord(self.subject.as_ref(), s_seq),
                );
                let font = egui::FontId::monospace(11.0);
                let label_size =
                    ui.painter().layout_no_wrap(label.clone(), font.clone(), Color32::BLACK).size();
                let pad = 4.0;
                let mut lx = cx + 6.0;
                let mut ly = cy + 6.0;
                if lx + label_size.x + pad > rect.right() {
                    lx = cx - 6.0 - label_size.x - 2.0 * pad;
                }
                if ly + label_size.y + pad > rect.bottom() {
                    ly = cy - 6.0 - label_size.y - 2.0 * pad;
                }
                // Backing rect for readability over dotplot ink.
                let bg = Rect::from_min_size(
                    Pos2::new(lx - pad, ly - pad),
                    Vec2::new(label_size.x + 2.0 * pad, label_size.y + 2.0 * pad),
                );
                ui.painter().rect_filled(
                    bg,
                    3.0,
                    Color32::from_rgba_unmultiplied(255, 255, 255, 220),
                );
                ui.painter().text(
                    Pos2::new(lx, ly),
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
    fn draw_axis_labels(&self, ui: &mut egui::Ui, rect: Rect, plot: &dottir_core::DotPlot) {
        const MIN_LABEL_SPACING_PX: f32 = 80.0;
        let zoom_us = plot.params.zoom.max(1) as f32;
        // World-pixel range visible in this view.
        let world_x_lo = self.view_offset.x;
        let world_x_hi = world_x_lo + rect.width() / self.display_zoom;
        let world_y_lo = self.view_offset.y;
        let world_y_hi = world_y_lo + rect.height() / self.display_zoom;
        // Convert to sequence-residue range.
        let seq_q_lo = (world_x_lo * zoom_us).max(0.0) as u64;
        let seq_q_hi = (world_x_hi * zoom_us) as u64;
        let seq_s_lo = (world_y_lo * zoom_us).max(0.0) as u64;
        let seq_s_hi = (world_y_hi * zoom_us) as u64;

        let painter = ui.painter();
        let tick_color = Color32::from_rgba_unmultiplied(80, 80, 80, 180);
        let label_color = ui.style().visuals.text_color();
        let font = egui::FontId::monospace(11.0);

        // Top axis (query).
        let span_x = seq_q_hi.saturating_sub(seq_q_lo) as f32;
        let pixels_per_residue_x = self.display_zoom / zoom_us;
        let step_x =
            nice_tick_step(span_x as f64, MIN_LABEL_SPACING_PX / pixels_per_residue_x as f32);
        let mut t = (seq_q_lo / step_x as u64) * step_x as u64;
        while t < seq_q_hi.saturating_add(step_x as u64) {
            if t >= seq_q_lo && t <= seq_q_hi {
                let sx = rect.left()
                    + (t as f32 / zoom_us - self.view_offset.x) * self.display_zoom;
                painter.line_segment(
                    [Pos2::new(sx, rect.top()), Pos2::new(sx, rect.top() + 6.0)],
                    egui::Stroke::new(1.0, tick_color),
                );
                painter.text(
                    Pos2::new(sx + 3.0, rect.top() + 2.0),
                    egui::Align2::LEFT_TOP,
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

        // Left axis (subject).
        let span_y = seq_s_hi.saturating_sub(seq_s_lo) as f32;
        let step_y = nice_tick_step(
            span_y as f64,
            MIN_LABEL_SPACING_PX / pixels_per_residue_x as f32,
        );
        let mut t = (seq_s_lo / step_y as u64) * step_y as u64;
        while t < seq_s_hi.saturating_add(step_y as u64) {
            if t >= seq_s_lo && t <= seq_s_hi {
                let sy = rect.top()
                    + (t as f32 / zoom_us - self.view_offset.y) * self.display_zoom;
                painter.line_segment(
                    [
                        Pos2::new(rect.left(), sy),
                        Pos2::new(rect.left() + 6.0, sy),
                    ],
                    egui::Stroke::new(1.0, tick_color),
                );
                painter.text(
                    Pos2::new(rect.left() + 8.0, sy + 1.0),
                    egui::Align2::LEFT_TOP,
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
    }
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
        match dottir_io::svg_export::write_svg(
            &path,
            plot.width,
            plot.height,
            &mapped,
            50,
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
    use crate::session::{codec, Session, SessionGreyramp, SessionPlot, SessionView, SESSION_VERSION};
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

