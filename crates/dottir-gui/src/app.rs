//! Top-level egui application.

use std::path::PathBuf;

use dottir_core::{
    compute_dotplot, BlastMode, DotPlot, PlotConfig, ScoreMatrix, Strand, Triangle,
};
use dottir_io::fasta;
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

/// One loaded sequence plus its source path.
struct LoadedSeq {
    path: PathBuf,
    seq: Vec<u8>,
}

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
    query: Option<LoadedSeq>,
    subject: Option<LoadedSeq>,
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
}

impl DottirApp {
    pub fn new(cc: &eframe::CreationContext<'_>, startup: StartupConfig) -> Self {
        // Default to a light theme — the plotting area is a greyscale
        // pixelmap on a near-white background, so a dark surround
        // muddles axis labels and panel text. Users who prefer dark
        // can toggle it in the View menu.
        cc.egui_ctx.set_visuals(egui::Visuals::light());

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
        app
    }

    fn load_fasta(&mut self, role: SeqRole, path: PathBuf) {
        match fasta::read_fasta_file(&path) {
            Ok(recs) => {
                let seq = fasta::concatenate(&recs);
                let n_records = recs.len();
                tracing::info!(
                    "loaded {} ({} residues, {} records)",
                    path.display(),
                    seq.len(),
                    n_records
                );
                let _ = n_records;
                let loaded = LoadedSeq { path, seq };
                match role {
                    SeqRole::Query => self.query = Some(loaded),
                    SeqRole::Subject => self.subject = Some(loaded),
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

    fn recompute(&mut self) {
        let (Some(q), Some(s)) = (&self.query, &self.subject) else {
            self.plot = None;
            self.texture = None;
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
        let cfg = PlotConfig {
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
            separate_strand_channels: false,
        };
        let qs = &q.seq;
        let ss = &s.seq;
        // BLASTP cannot use reverse strand.
        let effective_cfg = if cfg.mode == BlastMode::Blastp {
            let mut c = cfg.clone();
            c.strand = Strand::Forward;
            c
        } else {
            cfg.clone()
        };
        match compute_dotplot(qs, ss, &effective_cfg) {
            Ok(plot) => {
                tracing::info!(
                    "computed {}×{} pixelmap (W={})",
                    plot.width,
                    plot.height,
                    plot.params.window_size
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

    fn ensure_texture(&mut self, ctx: &Context) {
        let Some(plot) = &self.plot else {
            return;
        };
        if !self.texture_dirty && self.texture.is_some() {
            return;
        }
        let lut = self.greyramp.lut();
        let mut rgba = Vec::with_capacity(plot.pixels.len() * 4);
        for &v in &plot.pixels {
            let g = lut[v as usize];
            rgba.extend_from_slice(&[g, g, g, 255]);
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
                ui.label(if let Some(q) = &self.query {
                    format!(
                        "query: {} ({} bp)",
                        q.path.file_name().unwrap_or_default().to_string_lossy(),
                        q.seq.len()
                    )
                } else {
                    "query: —".into()
                });
                ui.label("·");
                ui.label(if let Some(s) = &self.subject {
                    format!(
                        "subject: {} ({} bp)",
                        s.path.file_name().unwrap_or_default().to_string_lossy(),
                        s.seq.len()
                    )
                } else {
                    "subject: —".into()
                });
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
                        let z = plot.params.zoom as u64;
                        ui.separator();
                        ui.label(format!(
                            "q = {} (≈seq {}), s = {} (≈seq {}), value = {}",
                            q,
                            (q as u64) * z + z / 2,
                            s,
                            (s as u64) * z + z / 2,
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

            // Crosshair overlay.
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
            }
        });
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

