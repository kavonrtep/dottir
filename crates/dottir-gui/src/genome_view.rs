//! Genome view (ADR 0008, G5): the alignments of two assemblies drawn as
//! identity-coloured lines, with D-GENIES sorting, display filters, record
//! flips, exports and drill-down into the dotplot engines.
//!
//! The visible window is rasterised at screen resolution whenever the
//! view changes (each pixel keeps the highest identity drawn into it, as
//! in `dottir_core::genome::raster`), so lines stay sharp at any zoom.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{mpsc, Arc, Mutex};
use std::time::{Duration, Instant};

use dottir_core::genome::raster::{identity_class, identity_code};
use dottir_core::genome::{
    arrange, Alignment, Arranged, AxisLayout, FilterConfig, GenomeRecord, RasterConfig,
};
use dottir_io::aligner::{align, AlignRequest, Aligner, AlignerError};
use dottir_io::faidx::{plain_fasta, FastaIndex, Region};
use dottir_io::genome_export::IDENTITY_COLOURS;
use dottir_io::paf::{read_paf_file, read_paf_oriented};
use egui::{Color32, Pos2, Rect, Sense, Stroke, TextureHandle, Vec2};
use serde::{Deserialize, Serialize};

/// Default largest region (per axis) handed to the dotplot engines.
pub const DEFAULT_DRILL_LIMIT_MB: f64 = 150.0;
const LEFT_MARGIN: f32 = 110.0;
const TOP_MARGIN: f32 = 28.0;
const BOTTOM_MARGIN: f32 = 6.0;
/// Hover picks the nearest alignment within this many screen pixels.
const HOVER_PX: f32 = 6.0;
const CLASS_LABELS: [&str; 4] = ["< 25 %", "25–50 %", "50–75 %", "≥ 75 %"];

/// What to compare and where the alignments come from.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GenomeInputs {
    /// Horizontal axis.
    pub reference: PathBuf,
    /// Vertical axis; `None` = self-comparison.
    pub query: Option<PathBuf>,
    /// Draw this PAF instead of running an aligner.
    pub paf: Option<PathBuf>,
    pub aligner: String,
    pub preset: String,
    pub aligner_args: String,
    /// Keep every chain (`-P`).
    #[serde(alias = "self_all_copies")]
    pub all_copies: bool,
}

impl Default for GenomeInputs {
    fn default() -> Self {
        Self {
            reference: PathBuf::new(),
            query: None,
            paf: None,
            aligner: "auto".into(),
            preset: "asm20".into(),
            aligner_args: String::new(),
            all_copies: false,
        }
    }
}

pub struct Assembly {
    /// As given (possibly gzipped).
    pub path: PathBuf,
    /// Index of the plain FASTA (decompressed into the cache if needed).
    pub index: FastaIndex,
    pub records: Vec<GenomeRecord>,
}

impl Assembly {
    fn total(&self) -> u64 {
        self.records.iter().map(|r| r.len).sum()
    }
}

pub struct GenomeData {
    pub x: Assembly,
    pub y: Assembly,
    pub self_comparison: bool,
    pub alns: Vec<Alignment>,
    /// Where the alignments came from, for the panel.
    pub source: String,
}

#[derive(Default)]
struct LoadStatus {
    stage: Mutex<String>,
    done: AtomicU64,
    total: AtomicU64,
}

impl LoadStatus {
    fn stage(&self, s: impl Into<String>) {
        *self.stage.lock().unwrap_or_else(|e| e.into_inner()) = s.into();
        self.total.store(0, Ordering::Relaxed);
    }
}

/// Indexing and alignment of a genome comparison, on its own thread.
pub struct GenomeLoad {
    rx: mpsc::Receiver<Result<GenomeData, String>>,
    status: Arc<LoadStatus>,
    cancel: Arc<AtomicBool>,
    pub inputs: GenomeInputs,
    started: Instant,
}

impl GenomeLoad {
    pub fn start(inputs: GenomeInputs, repaint: impl Fn() + Send + Sync + 'static) -> Self {
        let (tx, rx) = mpsc::channel();
        let status = Arc::new(LoadStatus::default());
        let cancel = Arc::new(AtomicBool::new(false));
        let (st, ca, inp) = (status.clone(), cancel.clone(), inputs.clone());
        let spawned = std::thread::Builder::new()
            .name("dottir-genome".into())
            .spawn(move || {
                let progress_repaint = || repaint();
                let r = load(&inp, &st, &ca, &progress_repaint);
                let _ = tx.send(r);
                repaint();
            });
        if let Err(e) = spawned {
            tracing::error!("genome loader: {e}");
        }
        Self {
            rx,
            status,
            cancel,
            inputs,
            started: Instant::now(),
        }
    }

    /// The result, once the thread has finished.
    pub fn poll(&self) -> Option<Result<GenomeData, String>> {
        match self.rx.try_recv() {
            Ok(r) => Some(r),
            Err(mpsc::TryRecvError::Empty) => None,
            Err(mpsc::TryRecvError::Disconnected) => Some(Err("genome loader stopped".into())),
        }
    }

    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::Relaxed);
    }

    /// Stage text, progress fraction when known, elapsed time.
    pub fn status(&self) -> (String, Option<f32>, Duration) {
        let stage = self
            .status
            .stage
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .clone();
        let total = self.status.total.load(Ordering::Relaxed);
        let frac =
            (total > 0).then(|| self.status.done.load(Ordering::Relaxed) as f32 / total as f32);
        (stage, frac, self.started.elapsed())
    }
}

fn load(
    inputs: &GenomeInputs,
    status: &LoadStatus,
    cancel: &AtomicBool,
    repaint: &(dyn Fn() + Sync),
) -> Result<GenomeData, String> {
    let cache = dottir_io::cache::cache_dir();
    let open = |p: &Path| -> Result<Assembly, String> {
        status.stage(format!("Indexing {}", file_name(p)));
        repaint();
        let plain = plain_fasta(p, &cache).map_err(|e| format!("{}: {e}", p.display()))?;
        let index =
            FastaIndex::open(&plain, &cache).map_err(|e| format!("{}: {e}", plain.display()))?;
        let records = index.records();
        Ok(Assembly {
            path: p.to_path_buf(),
            index,
            records,
        })
    };
    let x = open(&inputs.reference)?;
    let self_comparison = inputs.query.is_none();
    let y = match &inputs.query {
        Some(q) => open(q)?,
        None => Assembly {
            path: x.path.clone(),
            index: x.index.clone(),
            records: x.records.clone(),
        },
    };
    let (alns, source) = match &inputs.paf {
        Some(paf) => {
            status.stage(format!("Reading {}", file_name(paf)));
            repaint();
            let (alns, swapped) = read_paf_oriented(paf, &y.records, &x.records)
                .map_err(|e| format!("{}: {e}", paf.display()))?;
            let note = if swapped { ", axes swapped" } else { "" };
            (alns, format!("PAF {}{note}", file_name(paf)))
        }
        None => {
            let aligner = Aligner::detect(&inputs.aligner).map_err(|e| e.to_string())?;
            status.stage("Checksumming inputs");
            repaint();
            let sha = |a: &Assembly| {
                dottir_io::params::sha256_file(&a.path).map_err(|e| format!("{e:#}"))
            };
            let x_sha = sha(&x)?;
            let y_sha = if self_comparison {
                x_sha.clone()
            } else {
                sha(&y)?
            };
            let req = AlignRequest {
                reference: &x.index,
                query: &y.index,
                reference_sha256: x_sha,
                query_sha256: y_sha,
                preset: inputs.preset.clone(),
                extra_args: inputs
                    .aligner_args
                    .split_whitespace()
                    .map(String::from)
                    .collect(),
                threads: std::thread::available_parallelism().map_or(1, |n| n.get()),
                chunk_size: 10_000_000,
                self_comparison,
                all_copies: inputs.all_copies,
            };
            status.stage(format!(
                "Aligning with {} {}",
                aligner.name, aligner.version
            ));
            repaint();
            let progress = |done: u64, total: u64| {
                status.done.store(done, Ordering::Relaxed);
                status.total.store(total, Ordering::Relaxed);
                repaint();
            };
            let out =
                align(&aligner, &req, &cache, true, &progress, cancel).map_err(|e| match e {
                    AlignerError::Cancelled => "cancelled".to_string(),
                    e => e.to_string(),
                })?;
            status.stage("Reading alignments");
            repaint();
            let alns = read_paf_file(&out.paf, &y.records, &x.records)
                .map_err(|e| format!("{}: {e}", out.paf.display()))?;
            let cached = if out.cached { ", cached" } else { "" };
            let all = if inputs.all_copies { " -P" } else { "" };
            (
                alns,
                format!(
                    "{} {} -x {}{all}{cached}",
                    aligner.name, aligner.version, inputs.preset
                ),
            )
        }
    };
    Ok(GenomeData {
        x,
        y,
        self_comparison,
        alns,
        source,
    })
}

fn file_name(p: &Path) -> String {
    p.file_name()
        .map(|f| f.to_string_lossy().into_owned())
        .unwrap_or_else(|| p.display().to_string())
}

/// What the app should do after a genome-view frame.
pub enum GenomeAction {
    /// Compute a dotplot of these regions (x = query, y = subject).
    DrillDown { x: Vec<Region>, y: Vec<Region> },
    /// Run the aligner again with these inputs.
    Realign(GenomeInputs),
}

/// View window in axis coordinates (bp): `x0..x1` × `y0..y1`.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Window {
    x0: f64,
    x1: f64,
    y0: f64,
    y1: f64,
}

/// The hovered alignment: index, and whether its mirror image (self-
/// comparison) was hit.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Hover {
    aln: usize,
    mirrored: bool,
}

pub struct GenomeView {
    pub data: GenomeData,
    pub inputs: GenomeInputs,
    pub filter: FilterConfig,
    pub sort: bool,
    /// Per vertical record: flipped by the user (on top of sorting).
    pub flips: Vec<bool>,
    pub line_width: u32,
    pub drill_limit_mb: f64,
    arranged: Arranged,
    drawn: usize,
    view: Window,
    history: Vec<Window>,
    texture: Option<TextureHandle>,
    /// What `texture` shows; re-rasterised when this changes.
    texture_key: Option<(Window, [usize; 2], u64)>,
    /// Bumped on every change to filters, order or line width.
    version: u64,
    rect_start: Option<Pos2>,
    hover: Option<Hover>,
    cursor: Option<(f64, f64)>,
    /// Where the context menu was opened, in axis bp.
    context_at: Option<(f64, f64)>,
    pub message: Option<String>,
    /// Pointer position and counts, for the status bar.
    pub status: String,
}

impl GenomeView {
    pub fn new(data: GenomeData, inputs: GenomeInputs) -> Self {
        let filter = FilterConfig::default();
        let arranged = arrange(
            &data.alns,
            &data.x.records,
            &data.y.records,
            &filter,
            false,
            &[],
        );
        let drawn = arranged.keep.iter().filter(|&&k| k).count();
        let view = Window {
            x0: 0.0,
            x1: arranged.x.total as f64,
            y0: 0.0,
            y1: arranged.y.total as f64,
        };
        Self {
            flips: vec![false; data.y.records.len()],
            data,
            inputs,
            filter,
            sort: false,
            line_width: 2,
            drill_limit_mb: DEFAULT_DRILL_LIMIT_MB,
            arranged,
            drawn,
            view,
            history: Vec::new(),
            texture: None,
            texture_key: None,
            version: 0,
            rect_start: None,
            hover: None,
            cursor: None,
            context_at: None,
            message: None,
            status: String::new(),
        }
    }

    /// Re-apply filters, sorting and flips.
    pub fn rearrange(&mut self) {
        let d = &self.data;
        self.arranged = arrange(
            &d.alns,
            &d.x.records,
            &d.y.records,
            &self.filter,
            self.sort && !d.self_comparison,
            &self.flips,
        );
        self.drawn = self.arranged.keep.iter().filter(|&&k| k).count();
        self.hover = None;
        self.version += 1;
    }

    fn full(&self) -> Window {
        Window {
            x0: 0.0,
            x1: self.arranged.x.total.max(1) as f64,
            y0: 0.0,
            y1: self.arranged.y.total.max(1) as f64,
        }
    }

    pub fn fit(&mut self) {
        let full = self.full();
        if self.view != full {
            self.history.push(self.view);
            self.view = full;
        }
    }

    pub fn back(&mut self) {
        if let Some(v) = self.history.pop() {
            self.view = v;
        }
    }

    fn set_view(&mut self, v: Window) {
        if v != self.view {
            self.history.push(self.view);
            self.view = self.clamped(v);
        }
    }

    /// Keep the window inside the plot and no smaller than 100 bp.
    fn clamped(&self, mut v: Window) -> Window {
        let full = self.full();
        let fix = |a: &mut f64, b: &mut f64, max: f64| {
            let span = (*b - *a).clamp(100.0_f64.min(max), max);
            let mid = (*a + *b) / 2.0;
            *a = (mid - span / 2.0).clamp(0.0, max - span);
            *b = *a + span;
        };
        fix(&mut v.x0, &mut v.x1, full.x1);
        fix(&mut v.y0, &mut v.y1, full.y1);
        v
    }

    fn drill_limit(&self) -> f64 {
        self.drill_limit_mb * 1e6
    }

    fn view_fits_drill(&self) -> bool {
        let v = self.view;
        v.x1 - v.x0 <= self.drill_limit() && v.y1 - v.y0 <= self.drill_limit()
    }

    fn view_regions(&self) -> (Vec<Region>, Vec<Region>) {
        let v = self.view;
        let regions = |l: &AxisLayout, a: f64, b: f64| -> Vec<Region> {
            l.regions(a.floor() as u64..b.ceil() as u64)
                .into_iter()
                .map(|(rec, range, reverse)| Region {
                    rec: rec as usize,
                    range,
                    reverse,
                })
                .collect()
        };
        (
            regions(&self.arranged.x, v.x0, v.x1),
            regions(&self.arranged.y, v.y0, v.y1),
        )
    }

    /// Axis endpoints of alignment `a`, or of its mirror image.
    fn endpoints(&self, a: &Alignment, mirrored: bool) -> ((f64, f64), (f64, f64)) {
        let (x, y) = (&self.arranged.x, &self.arranged.y);
        let (s, e) = if mirrored {
            (
                (
                    x.to_axis(a.q_rec, a.q_start),
                    y.to_axis(a.t_rec, a.t_at_q_start()),
                ),
                (
                    x.to_axis(a.q_rec, a.q_end),
                    y.to_axis(a.t_rec, a.t_at_q_end()),
                ),
            )
        } else {
            (
                (
                    x.to_axis(a.t_rec, a.t_at_q_start()),
                    y.to_axis(a.q_rec, a.q_start),
                ),
                (
                    x.to_axis(a.t_rec, a.t_at_q_end()),
                    y.to_axis(a.q_rec, a.q_end),
                ),
            )
        };
        ((s.0 as f64, s.1 as f64), (e.0 as f64, e.1 as f64))
    }

    /// Rasterise the window into `size` pixels (RGBA, transparent where
    /// nothing is drawn).
    fn rasterize_window(&self, size: [usize; 2]) -> egui::ColorImage {
        let [w, h] = size;
        let mut codes = vec![0u8; w * h];
        let v = self.view;
        let (sx, sy) = (w as f64 / (v.x1 - v.x0), h as f64 / (v.y1 - v.y0));
        let lw = self.line_width.max(1) as i64;
        let mut line = |from: (f64, f64), to: (f64, f64), code: u8| {
            let p0 = ((from.0 - v.x0) * sx, (from.1 - v.y0) * sy);
            let p1 = ((to.0 - v.x0) * sx, (to.1 - v.y0) * sy);
            let Some((p0, p1)) = clip(p0, p1, w as f64, h as f64, lw as f64) else {
                return;
            };
            let n = (p1.0 - p0.0).abs().max((p1.1 - p0.1).abs()).ceil().max(1.0) as i64;
            let off = (lw - 1) / 2;
            for i in 0..=n {
                let t = i as f64 / n as f64;
                let cx = (p0.0 + (p1.0 - p0.0) * t).floor() as i64;
                let cy = (p0.1 + (p1.1 - p0.1) * t).floor() as i64;
                for yy in (cy - off).max(0)..(cy - off + lw).min(h as i64) {
                    for xx in (cx - off).max(0)..(cx - off + lw).min(w as i64) {
                        let p = &mut codes[yy as usize * w + xx as usize];
                        *p = (*p).max(code);
                    }
                }
            }
        };
        for a in self.arranged.drawn(&self.data.alns) {
            let code = identity_code(a.identity());
            let (s, e) = self.endpoints(a, false);
            line(s, e, code);
            if self.data.self_comparison {
                let (s, e) = self.endpoints(a, true);
                line(s, e, code);
            }
        }
        if self.data.self_comparison {
            let (x, y) = (&self.arranged.x, &self.arranged.y);
            for s in &x.segments {
                let from = (x.to_axis(s.rec, 0) as f64, y.to_axis(s.rec, 0) as f64);
                let to = (
                    x.to_axis(s.rec, s.len) as f64,
                    y.to_axis(s.rec, s.len) as f64,
                );
                line(from, to, identity_code(1.0));
            }
        }
        let pixels = codes
            .iter()
            .map(|&c| match identity_class(c) {
                None => Color32::TRANSPARENT,
                Some(k) => {
                    let [r, g, b] = IDENTITY_COLOURS[k];
                    Color32::from_rgb(r, g, b)
                }
            })
            .collect();
        egui::ColorImage { size, pixels }
    }

    /// The nearest drawn alignment to screen point `p` within
    /// [`HOVER_PX`].
    fn pick(&self, p: Pos2, to_screen: &dyn Fn((f64, f64)) -> Pos2) -> Option<Hover> {
        let mut best: Option<(f32, Hover)> = None;
        let sides: &[bool] = if self.data.self_comparison {
            &[false, true]
        } else {
            &[false]
        };
        for (i, (a, &k)) in self.data.alns.iter().zip(&self.arranged.keep).enumerate() {
            if !k {
                continue;
            }
            for &mirrored in sides {
                let (s, e) = self.endpoints(a, mirrored);
                let d = segment_distance(p, to_screen(s), to_screen(e));
                if d <= HOVER_PX && best.is_none_or(|(bd, _)| d < bd) {
                    best = Some((d, Hover { aln: i, mirrored }));
                }
            }
        }
        best.map(|(_, h)| h)
    }

    /// Record name and 1-based position at axis coordinate `v`.
    fn position(&self, layout: &AxisLayout, records: &[GenomeRecord], v: f64) -> Option<String> {
        if v < 0.0 {
            return None;
        }
        let (seg, pos) = layout.locate(v as u64)?;
        let strand = if seg.reverse { " (−)" } else { "" };
        Some(format!(
            "{}:{}{strand}",
            records[seg.rec as usize].name,
            thousands(pos + 1)
        ))
    }

    /// The canvas: plot, axes, interaction. Returns a drill-down request.
    pub fn canvas_ui(&mut self, ui: &mut egui::Ui) -> Option<GenomeAction> {
        let mut action = None;
        let (rect, response) = ui.allocate_exact_size(ui.available_size(), Sense::click_and_drag());
        let plot = Rect::from_min_max(
            Pos2::new(rect.left() + LEFT_MARGIN, rect.top() + TOP_MARGIN),
            Pos2::new(rect.right() - 8.0, rect.bottom() - BOTTOM_MARGIN),
        );
        if plot.width() < 20.0 || plot.height() < 20.0 {
            return None;
        }
        let v = self.view;
        let (kx, ky) = (
            plot.width() as f64 / (v.x1 - v.x0),
            plot.height() as f64 / (v.y1 - v.y0),
        );
        let to_bp = move |p: Pos2| {
            (
                v.x0 + f64::from(p.x - plot.left()) / kx,
                v.y0 + f64::from(p.y - plot.top()) / ky,
            )
        };

        // Interaction.
        let pointer = response.hover_pos().filter(|p| plot.contains(*p));
        self.cursor = pointer.map(to_bp);
        if let Some(p) = pointer {
            let scroll = ui.input(|i| i.smooth_scroll_delta.y);
            if scroll != 0.0 {
                let f = f64::from(-scroll * 0.003).exp();
                let (cx, cy) = to_bp(p);
                let z = Window {
                    x0: cx - (cx - v.x0) * f,
                    x1: cx + (v.x1 - cx) * f,
                    y0: cy - (cy - v.y0) * f,
                    y1: cy + (v.y1 - cy) * f,
                };
                // Wheel zooming is not recorded in the history.
                self.view = self.clamped(z);
            }
        }
        let rect_drag = |r: &egui::Response, b| {
            r.dragged_by(b) && (b == egui::PointerButton::Middle || ui.input(|i| i.modifiers.shift))
        };
        if response.drag_started()
            && (response.drag_started_by(egui::PointerButton::Middle)
                || ui.input(|i| i.modifiers.shift))
        {
            self.rect_start = response
                .interact_pointer_pos()
                .filter(|p| plot.contains(*p));
        } else if self.rect_start.is_none() && response.dragged_by(egui::PointerButton::Primary) {
            let d = response.drag_delta();
            let pan = Window {
                x0: v.x0 - f64::from(d.x) / kx,
                x1: v.x1 - f64::from(d.x) / kx,
                y0: v.y0 - f64::from(d.y) / ky,
                y1: v.y1 - f64::from(d.y) / ky,
            };
            self.view = self.clamped(pan);
        }
        let rect_sel = match (self.rect_start, response.interact_pointer_pos()) {
            (Some(a), Some(b)) => Some(Rect::from_two_pos(a, b).intersect(plot)),
            _ => None,
        };
        if response.drag_stopped() {
            if let (Some(r), Some(_)) = (rect_sel, self.rect_start.take()) {
                if r.width() > 4.0 && r.height() > 4.0 {
                    let (a, b) = (to_bp(r.min), to_bp(r.max));
                    self.set_view(Window {
                        x0: a.0,
                        x1: b.0,
                        y0: a.1,
                        y1: b.1,
                    });
                }
            }
        }
        let _ = rect_drag;
        if response.double_clicked() {
            self.fit();
        }
        if response.secondary_clicked() {
            self.context_at = self.cursor;
        }

        // Re-rasterise when the window, size or content changed.
        let ppp = ui.ctx().pixels_per_point();
        let size = [
            (plot.width() * ppp).round().max(1.0) as usize,
            (plot.height() * ppp).round().max(1.0) as usize,
        ];
        let key = (self.view, size, self.version);
        if self.texture_key != Some(key) || self.texture.is_none() {
            let img = self.rasterize_window(size);
            match &mut self.texture {
                Some(t) => t.set(img, egui::TextureOptions::NEAREST),
                None => {
                    self.texture = Some(ui.ctx().load_texture(
                        "genome-view",
                        img,
                        egui::TextureOptions::NEAREST,
                    ))
                }
            }
            self.texture_key = Some(key);
        }

        // Paint: background, small-record blocks, grid, lines, frame.
        let painter = ui.painter_at(rect);
        let v = self.view;
        let to_screen = move |p: (f64, f64)| {
            Pos2::new(
                plot.left() + ((p.0 - v.x0) * kx) as f32,
                plot.top() + ((p.1 - v.y0) * ky) as f32,
            )
        };
        let dark = ui.visuals().dark_mode;
        let (bg, grid, ink, shade) = if dark {
            (
                Color32::from_gray(250),
                Color32::from_gray(200),
                Color32::from_gray(220),
                Color32::from_gray(232),
            )
        } else {
            (
                Color32::WHITE,
                Color32::from_gray(205),
                Color32::from_gray(30),
                Color32::from_gray(238),
            )
        };
        painter.rect_filled(plot, 0.0, bg);
        let (xl, yl) = (&self.arranged.x, &self.arranged.y);
        for m in &xl.mixes {
            let a = to_screen((m.start as f64, v.y0)).x.max(plot.left());
            let b = to_screen((m.end as f64, v.y0)).x.min(plot.right());
            if b > a {
                painter.rect_filled(Rect::from_x_y_ranges(a..=b, plot.y_range()), 0.0, shade);
            }
        }
        for m in &yl.mixes {
            let a = to_screen((v.x0, m.start as f64)).y.max(plot.top());
            let b = to_screen((v.x0, m.end as f64)).y.min(plot.bottom());
            if b > a {
                painter.rect_filled(Rect::from_x_y_ranges(plot.x_range(), a..=b), 0.0, shade);
            }
        }
        // As the dotplot's axis labels.
        let font = egui::FontId::monospace(11.0);
        let small = egui::FontId::monospace(10.0);
        // Vertical grid and record names along the top.
        let mut last_label = f32::NEG_INFINITY;
        for (i, s) in xl.segments.iter().enumerate() {
            let a = to_screen((s.start as f64, 0.0)).x;
            let b = to_screen(((s.start + s.len) as f64, 0.0)).x;
            if b < plot.left() || a > plot.right() {
                continue;
            }
            if i > 0 && a >= plot.left() && s.mix.is_none_or(|m| xl.segments[i - 1].mix != Some(m))
            {
                painter.vline(a, plot.y_range(), Stroke::new(1.0, grid));
            }
            let (ca, cb) = (a.max(plot.left()), b.min(plot.right()));
            if s.mix.is_none() && cb - ca > 24.0 && ca >= last_label {
                let name = label(&self.data.x.records[s.rec as usize].name, s.reverse);
                let text = fit_text(ui, &name, cb - ca - 4.0, &font);
                painter.text(
                    Pos2::new((ca + cb) / 2.0, plot.top() - 4.0),
                    egui::Align2::CENTER_BOTTOM,
                    text,
                    font.clone(),
                    ui.visuals().text_color(),
                );
                last_label = cb;
            }
        }
        for (k, m) in xl.mixes.iter().enumerate() {
            let a = to_screen((m.start as f64, 0.0)).x.max(plot.left());
            let b = to_screen((m.end as f64, 0.0)).x.min(plot.right());
            if b - a > 40.0 {
                let n = xl
                    .segments
                    .iter()
                    .filter(|s| s.mix == Some(k as u32))
                    .count();
                painter.text(
                    Pos2::new((a + b) / 2.0, plot.top() - 4.0),
                    egui::Align2::CENTER_BOTTOM,
                    format!("{n} small"),
                    small.clone(),
                    ui.visuals().weak_text_color(),
                );
            }
        }
        // Horizontal grid and record names on the left.
        let mut last_label = f32::NEG_INFINITY;
        for (i, s) in yl.segments.iter().enumerate() {
            let a = to_screen((0.0, s.start as f64)).y;
            let b = to_screen((0.0, (s.start + s.len) as f64)).y;
            if b < plot.top() || a > plot.bottom() {
                continue;
            }
            if i > 0 && a >= plot.top() && s.mix.is_none_or(|m| yl.segments[i - 1].mix != Some(m)) {
                painter.hline(plot.x_range(), a, Stroke::new(1.0, grid));
            }
            let (ca, cb) = (a.max(plot.top()), b.min(plot.bottom()));
            if s.mix.is_none() && cb - ca > 11.0 && ca >= last_label {
                let name = label(&self.data.y.records[s.rec as usize].name, s.reverse);
                let text = fit_text(ui, &name, LEFT_MARGIN - 8.0, &font);
                painter.text(
                    Pos2::new(plot.left() - 4.0, (ca + cb) / 2.0),
                    egui::Align2::RIGHT_CENTER,
                    text,
                    font.clone(),
                    ui.visuals().text_color(),
                );
                last_label = cb;
            }
        }
        for (k, m) in yl.mixes.iter().enumerate() {
            let a = to_screen((0.0, m.start as f64)).y.max(plot.top());
            let b = to_screen((0.0, m.end as f64)).y.min(plot.bottom());
            if b - a > 14.0 {
                let n = yl
                    .segments
                    .iter()
                    .filter(|s| s.mix == Some(k as u32))
                    .count();
                painter.text(
                    Pos2::new(plot.left() - 4.0, (a + b) / 2.0),
                    egui::Align2::RIGHT_CENTER,
                    format!("{n} small"),
                    small.clone(),
                    ui.visuals().weak_text_color(),
                );
            }
        }
        if let Some(t) = &self.texture {
            painter.image(
                t.id(),
                plot,
                Rect::from_min_max(Pos2::ZERO, Pos2::new(1.0, 1.0)),
                Color32::WHITE,
            );
        }
        painter.rect_stroke(plot, 0.0, Stroke::new(1.0, ink));
        if let (Some(r), Some(_)) = (rect_sel, self.rect_start) {
            painter.rect_stroke(r, 0.0, Stroke::new(1.0, Color32::from_rgb(220, 40, 40)));
        }

        // Hover: highlight and describe the nearest alignment.
        if response.is_pointer_button_down_on() {
            self.hover = None;
        } else if let Some(p) = pointer {
            if ui.input(|i| i.pointer.delta() != Vec2::ZERO) || self.hover.is_none() {
                self.hover = self.pick(p, &to_screen);
            }
        } else {
            self.hover = None;
        }
        if let Some(h) = self.hover {
            let a = &self.data.alns[h.aln];
            let (s, e) = self.endpoints(a, h.mirrored);
            painter.line_segment(
                [to_screen(s), to_screen(e)],
                Stroke::new(3.0, Color32::from_rgb(230, 60, 30)),
            );
            let text = self.describe(a);
            response.clone().on_hover_ui_at_pointer(|ui| {
                ui.label(egui::RichText::new(text).monospace());
            });
        }

        // Status line.
        let mut status = format!(
            "{} alignments, {} drawn · view {} × {}",
            thousands(self.data.alns.len() as u64),
            thousands(self.drawn as u64),
            span_text(v.x1 - v.x0),
            span_text(v.y1 - v.y0),
        );
        if let Some((cx, cy)) = self.cursor {
            let xs = self.position(xl, &self.data.x.records, cx);
            let ys = self.position(yl, &self.data.y.records, cy);
            if let (Some(xs), Some(ys)) = (xs, ys) {
                status = format!("x {xs}   y {ys}   ·   {status}");
            }
        }
        self.status = status;

        response.context_menu(|ui| {
            if let Some(a) = self.context_menu(ui) {
                action = Some(a);
            }
        });
        action
    }

    fn describe(&self, a: &Alignment) -> String {
        let (q, t) = (
            &self.data.y.records[a.q_rec as usize].name,
            &self.data.x.records[a.t_rec as usize].name,
        );
        let exact = if a.divergence.is_some() {
            ""
        } else {
            " (approx.)"
        };
        format!(
            "{q}:{}-{}\n{t}:{}-{}\nstrand {}   {} bp   identity {:.1} %{exact}{}",
            thousands(a.q_start + 1),
            thousands(a.q_end),
            thousands(a.t_start + 1),
            thousands(a.t_end),
            if a.reverse { '−' } else { '+' },
            thousands(a.len()),
            a.identity() * 100.0,
            if a.primary { "" } else { "   secondary" },
        )
    }

    fn context_menu(&mut self, ui: &mut egui::Ui) -> Option<GenomeAction> {
        let mut action = None;
        let at = self.context_at;
        let x_seg = at.and_then(|(x, _)| {
            self.arranged
                .x
                .locate(x.max(0.0) as u64)
                .map(|(s, _)| s.clone())
        });
        let y_seg = at.and_then(|(_, y)| {
            self.arranged
                .y
                .locate(y.max(0.0) as u64)
                .map(|(s, _)| s.clone())
        });
        let limit = self.drill_limit();
        if ui
            .add_enabled(
                self.view_fits_drill(),
                egui::Button::new("Dot-plot this view"),
            )
            .on_disabled_hover_text(format!("Zoom in to at most {} per axis", span_text(limit)))
            .clicked()
        {
            let (x, y) = self.view_regions();
            action = Some(GenomeAction::DrillDown { x, y });
            ui.close_menu();
        }
        if let (Some(xs), Some(ys)) = (&x_seg, &y_seg) {
            let (xn, yn) = (
                &self.data.x.records[xs.rec as usize].name,
                &self.data.y.records[ys.rec as usize].name,
            );
            let fits = xs.len as f64 <= limit && ys.len as f64 <= limit;
            if ui
                .add_enabled(fits, egui::Button::new(format!("Dot-plot {xn} × {yn}")))
                .on_disabled_hover_text(format!("Records over {}", span_text(limit)))
                .clicked()
            {
                let whole = |s: &dottir_core::genome::AxisSegment| Region {
                    rec: s.rec as usize,
                    range: 0..s.len,
                    reverse: s.reverse,
                };
                action = Some(GenomeAction::DrillDown {
                    x: vec![whole(xs)],
                    y: vec![whole(ys)],
                });
                ui.close_menu();
            }
            ui.separator();
            if !self.data.self_comparison {
                let what = if self.flips[ys.rec as usize] {
                    "Unflip"
                } else {
                    "Flip"
                };
                if ui.button(format!("{what} {yn}")).clicked() {
                    self.flips[ys.rec as usize] ^= true;
                    self.rearrange();
                    ui.close_menu();
                }
            }
        }
        if ui.button("Fit (double-click)").clicked() {
            self.fit();
            ui.close_menu();
        }
        if ui
            .add_enabled(!self.history.is_empty(), egui::Button::new("Back (Esc)"))
            .clicked()
        {
            self.back();
            ui.close_menu();
        }
        action
    }

    /// The right-hand panel, laid out like the dotplot's: view, inputs,
    /// filters, order, legend, drill-down, exports.
    pub fn panel_ui(&mut self, ui: &mut egui::Ui) -> Option<GenomeAction> {
        let mut action = None;
        ui.style_mut().spacing.item_spacing = Vec2::new(4.0, 3.0);
        ui.style_mut().spacing.interact_size.y = 18.0;
        let small = |s: String| egui::RichText::new(s).size(11.0);
        let weak = |s: String| {
            egui::RichText::new(s)
                .size(11.0)
                .color(Color32::from_gray(110))
        };

        heading(ui, "View");
        ui.horizontal(|ui| {
            if ui
                .button("Fit")
                .on_hover_text("Show both assemblies whole (double-click).")
                .clicked()
            {
                self.fit();
            }
            if ui
                .add_enabled(!self.history.is_empty(), egui::Button::new("Back"))
                .on_hover_text("Return to the previous view (Esc).")
                .clicked()
            {
                self.back();
            }
        });
        ui.separator();

        heading(ui, "Sequences");
        let d = &self.data;
        let self_comparison = d.self_comparison;
        let approximate = d.alns.first().is_some_and(|a| a.divergence.is_none());
        let asm = |ui: &mut egui::Ui, what: &str, a: &Assembly| {
            ui.label(small(format!(
                "{what}: {}  {} · {} recs",
                file_name(&a.path),
                span_text(a.total() as f64),
                thousands(a.records.len() as u64),
            )))
            .on_hover_text(a.path.display().to_string());
        };
        asm(ui, "Reference", &d.x);
        if self_comparison {
            ui.label(small("Query: same (self-comparison)".into()));
        } else {
            asm(ui, "Query    ", &d.y);
        }
        ui.label(weak(d.source.clone()));
        ui.label(weak(format!(
            "{} alignments · {} drawn",
            thousands(d.alns.len() as u64),
            thousands(self.drawn as u64)
        )));
        ui.separator();

        if self.inputs.paf.is_none() {
            let mut all = self.inputs.all_copies;
            if ui
                .checkbox(&mut all, small("Keep all copies of repeats (-P)".into()))
                .on_hover_text(
                    "By default the aligner reports each region's best \
                     placement only, so repeats are mostly absent. -P keeps \
                     every chain: transposon and duplication blocks appear \
                     (often 10–100× more alignments). Satellite arrays stay \
                     hidden; dot-plot a region for those. Re-runs the \
                     aligner (cached).",
                )
                .changed()
            {
                let mut inputs = self.inputs.clone();
                inputs.all_copies = all;
                action = Some(GenomeAction::Realign(inputs));
            }
            if !all {
                ui.label(weak("Repeats mostly hidden (best placement only)".into()));
            }
        }
        ui.separator();

        heading(ui, "Filters");
        let mut changed = false;
        ui.horizontal(|ui| {
            ui.label(small("Min length".into()));
            changed |= ui
                .add(
                    egui::DragValue::new(&mut self.filter.min_len)
                        .speed(100.0)
                        .range(0..=10_000_000)
                        .suffix(" bp"),
                )
                .changed();
        });
        ui.label(small("Min identity".into()));
        let r =
            ui.add(egui::Slider::new(&mut self.filter.min_identity, 0.0..=1.0).fixed_decimals(2));
        // Re-filtering a million alignments per slider step stutters;
        // apply when the drag ends.
        changed |= r.drag_stopped() || (r.changed() && !r.dragged());
        changed |= ui
            .checkbox(&mut self.filter.primary_only, small("Primary only".into()))
            .changed();
        changed |= ui
            .checkbox(
                &mut self.filter.remove_noise,
                small("Remove noise (D-GENIES)".into()),
            )
            .on_hover_text(
                "With more than 1000 alignments, hide those shorter than \
                 the first sparse length bin after the most populated one",
            )
            .changed();
        if let Some(t) = self.arranged.noise {
            ui.label(weak(format!("noise below {}", span_text(t))));
        }
        ui.separator();

        heading(ui, "Order");
        changed |= ui
            .add_enabled(
                !self_comparison,
                egui::Checkbox::new(&mut self.sort, small("Sort query along reference".into())),
            )
            .changed();
        let n_flipped = self.flips.iter().filter(|&&f| f).count();
        if ui
            .add_enabled(
                n_flipped > 0,
                egui::Button::new(format!("Reset flips ({n_flipped})")).small(),
            )
            .on_hover_text("Right-click a record in the plot to flip it")
            .clicked()
        {
            self.flips.fill(false);
            changed = true;
        }
        if changed {
            self.rearrange();
        }
        ui.separator();

        heading(ui, "Identity");
        for (k, c) in IDENTITY_COLOURS.iter().enumerate() {
            ui.horizontal(|ui| {
                let (r, _) = ui.allocate_exact_size(Vec2::splat(10.0), Sense::hover());
                ui.painter()
                    .rect_filled(r, 2.0, Color32::from_rgb(c[0], c[1], c[2]));
                ui.label(small(CLASS_LABELS[k].into()));
            });
        }
        if approximate {
            ui.label(weak(
                "Approximate: no de:f tags (aligner -c for exact)".into(),
            ));
        }
        ui.horizontal(|ui| {
            ui.label(small("Line width".into()));
            if ui
                .add(egui::Slider::new(&mut self.line_width, 1..=5).suffix(" px"))
                .changed()
            {
                self.version += 1;
            }
        });
        ui.separator();

        heading(ui, "Dot-plot");
        ui.horizontal(|ui| {
            ui.label(small("Up to".into()));
            ui.add(
                egui::DragValue::new(&mut self.drill_limit_mb)
                    .speed(1.0)
                    .range(1.0..=2000.0)
                    .suffix(" Mb"),
            );
            ui.label(small("per axis".into()));
        });
        if ui
            .add_enabled(
                self.view_fits_drill(),
                egui::Button::new("Dot-plot this view"),
            )
            .on_hover_text(
                "Compute a dotplot of the visible region with the large or \
                 exact engine. Right-click a cell for a whole record pair.",
            )
            .clicked()
        {
            let (x, y) = self.view_regions();
            action = Some(GenomeAction::DrillDown { x, y });
        }
        if !self.view_fits_drill() {
            ui.label(weak(format!(
                "Zoom in to at most {} per axis",
                span_text(self.drill_limit())
            )));
        }
        ui.separator();

        heading(ui, "Export");
        ui.horizontal_wrapped(|ui| {
            if ui
                .small_button("PNG…")
                .on_hover_text("The whole plot")
                .clicked()
            {
                self.export(Export::Png);
            }
            if ui
                .small_button("Association…")
                .on_hover_text("Best reference record per query record (TSV)")
                .clicked()
            {
                self.export(Export::Association);
            }
            if ui
                .add_enabled(!self_comparison, egui::Button::new("Sorted query…").small())
                .on_hover_text("Query FASTA in plot order, reverse-complemented as shown")
                .clicked()
            {
                self.export(Export::SortedQuery);
            }
            if ui
                .small_button("PAF…")
                .on_hover_text("The alignments drawn")
                .clicked()
            {
                self.export(Export::Paf);
            }
        });
        if let Some(m) = &self.message {
            ui.label(weak(m.clone()));
        }
        action
    }

    /// Save the whole plot as PNG (also *File → Save PNG*).
    pub fn export_png(&mut self) {
        self.export(Export::Png);
    }

    pub fn export_svg(&mut self) {
        self.export(Export::Svg);
    }

    /// Menu-bar labels of the two inputs, as the dotplot shows them.
    pub fn input_labels(&self) -> (String, String) {
        let label = |a: &Assembly| {
            format!(
                "{} ({} bp, {} records)",
                file_name(&a.path),
                thousands(a.total()),
                a.records.len()
            )
        };
        let q = if self.data.self_comparison {
            "same (self-comparison)".to_string()
        } else {
            label(&self.data.y)
        };
        (label(&self.data.x), q)
    }

    pub fn can_go_back(&self) -> bool {
        !self.history.is_empty()
    }

    fn export(&mut self, what: Export) {
        let (name, ext) = match what {
            Export::Png => ("genome.png", "png"),
            Export::Svg => ("genome.svg", "svg"),
            Export::Association => ("assoc.tsv", "tsv"),
            Export::SortedQuery => ("query.sorted.fa", "fa"),
            Export::Paf => ("drawn.paf", "paf"),
        };
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(name)
            .add_filter(ext, &[ext])
            .save_file()
        else {
            return;
        };
        self.message = Some(match self.write_export(what, &path) {
            Ok(()) => format!("Wrote {}", path.display()),
            Err(e) => format!("Export failed: {e:#}"),
        });
    }

    fn write_export(&self, what: Export, path: &Path) -> anyhow::Result<()> {
        use anyhow::Context as _;
        let create = || -> anyhow::Result<std::io::BufWriter<std::fs::File>> {
            Ok(std::io::BufWriter::new(
                std::fs::File::create(path)
                    .with_context(|| format!("creating {}", path.display()))?,
            ))
        };
        let (d, a) = (&self.data, &self.arranged);
        match what {
            Export::Png => {
                let bp = a.x.total.max(a.y.total).div_ceil(4000).max(1);
                let cfg = RasterConfig {
                    bp_per_px: bp,
                    line_width: 2,
                    self_comparison: d.self_comparison,
                };
                let raster = dottir_core::genome::rasterize(a.drawn(&d.alns), &a.x, &a.y, &cfg);
                let version = dottir_core::genome::GENOME_FORMAT_VERSION.to_string();
                let bp = bp.to_string();
                dottir_io::genome_export::write_genome_png(
                    path,
                    &raster,
                    &a.x,
                    &d.x.records,
                    &a.y,
                    &d.y.records,
                    &[
                        ("dottir-engine", "genome"),
                        ("dottir-version", env!("CARGO_PKG_VERSION")),
                        ("dottir-genome-format", &version),
                        ("dottir-alignments", &d.source),
                        ("dottir-bp-per-px", &bp),
                    ],
                )?;
            }
            Export::Svg => {
                let version = dottir_core::genome::GENOME_FORMAT_VERSION.to_string();
                let lines = dottir_core::genome::plot_lines(
                    a.drawn(&d.alns),
                    &a.x,
                    &a.y,
                    d.self_comparison,
                );
                dottir_io::genome_export::write_genome_svg(
                    create()?,
                    &lines,
                    1.5,
                    &a.x,
                    &d.x.records,
                    &a.y,
                    &d.y.records,
                    &[
                        ("dottir-engine", "genome"),
                        ("dottir-version", env!("CARGO_PKG_VERSION")),
                        ("dottir-genome-format", &version),
                        ("dottir-alignments", &d.source),
                    ],
                )?;
            }
            Export::Association => dottir_io::genome_export::write_association(
                create()?,
                &a.y,
                &d.y.records,
                &d.x.records,
                &a.placements,
            )?,
            Export::SortedQuery => {
                let order: Vec<u32> = a.y.segments.iter().map(|s| s.rec).collect();
                let mut flips = vec![false; d.y.records.len()];
                for s in &a.y.segments {
                    flips[s.rec as usize] = s.reverse;
                }
                dottir_io::genome_export::write_fasta_records(
                    create()?,
                    &d.y.index,
                    &order,
                    &flips,
                )?;
            }
            Export::Paf => {
                let mut out: Vec<Alignment> = a.drawn(&d.alns).cloned().collect();
                dottir_core::genome::alignment::sort_alignments(&mut out);
                dottir_io::paf::write_paf(create()?, &out, &d.y.records, &d.x.records)?;
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
enum Export {
    Png,
    Svg,
    Association,
    SortedQuery,
    Paf,
}

fn heading(ui: &mut egui::Ui, text: &str) {
    ui.label(
        egui::RichText::new(text)
            .strong()
            .size(13.0)
            .color(Color32::from_gray(50)),
    );
}

/// Clip segment `p0`–`p1` to the box `[-m, w + m] × [-m, h + m]`
/// (Liang–Barsky).
fn clip(
    p0: (f64, f64),
    p1: (f64, f64),
    w: f64,
    h: f64,
    m: f64,
) -> Option<((f64, f64), (f64, f64))> {
    let (dx, dy) = (p1.0 - p0.0, p1.1 - p0.1);
    let (mut t0, mut t1) = (0.0_f64, 1.0_f64);
    for (p, q) in [
        (-dx, p0.0 + m),
        (dx, w + m - p0.0),
        (-dy, p0.1 + m),
        (dy, h + m - p0.1),
    ] {
        if p == 0.0 {
            if q < 0.0 {
                return None;
            }
        } else {
            let r = q / p;
            if p < 0.0 {
                t0 = t0.max(r);
            } else {
                t1 = t1.min(r);
            }
        }
    }
    (t0 <= t1).then_some((
        (p0.0 + t0 * dx, p0.1 + t0 * dy),
        (p0.0 + t1 * dx, p0.1 + t1 * dy),
    ))
}

fn segment_distance(p: Pos2, a: Pos2, b: Pos2) -> f32 {
    let ab = b - a;
    let len2 = ab.length_sq();
    let t = if len2 == 0.0 {
        0.0
    } else {
        ((p - a).dot(ab) / len2).clamp(0.0, 1.0)
    };
    (a + ab * t).distance(p)
}

fn label(name: &str, reverse: bool) -> String {
    if reverse {
        format!("{name}-")
    } else {
        name.to_string()
    }
}

/// `text`, shortened with an ellipsis to fit `max_w` points.
fn fit_text(ui: &egui::Ui, text: &str, max_w: f32, font: &egui::FontId) -> String {
    let width = |s: &str| {
        ui.fonts(|f| {
            f.layout_no_wrap(s.to_string(), font.clone(), Color32::WHITE)
                .size()
                .x
        })
    };
    if width(text) <= max_w {
        return text.to_string();
    }
    let chars: Vec<char> = text.chars().collect();
    for n in (1..chars.len()).rev() {
        let s: String = chars[..n].iter().collect::<String>() + "…";
        if width(&s) <= max_w {
            return s;
        }
    }
    String::new()
}

fn thousands(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            out.push(',');
        }
        out.push(c);
    }
    out
}

fn span_text(bp: f64) -> String {
    if bp >= 1e9 {
        format!("{:.2} Gb", bp / 1e9)
    } else if bp >= 1e6 {
        format!("{:.1} Mb", bp / 1e6)
    } else if bp >= 1e3 {
        format!("{:.1} kb", bp / 1e3)
    } else {
        format!("{bp:.0} bp")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn view(dir: &Path, sort: bool) -> GenomeView {
        let _ = std::fs::remove_dir_all(dir);
        std::fs::create_dir_all(dir).unwrap();
        let (r, q) = (dir.join("r.fa"), dir.join("q.fa"));
        // chr1 = 1000 × A then 1000 × C; ctg = reverse complement of
        // chr1[1000..1500) (GGG…), followed by 500 × T.
        std::fs::write(
            &r,
            format!(">chr1\n{}{}\n", "A".repeat(1000), "C".repeat(1000)),
        )
        .unwrap();
        std::fs::write(
            &q,
            format!(">ctg\n{}{}\n", "G".repeat(500), "T".repeat(500)),
        )
        .unwrap();
        let asm = |p: &Path| {
            let index = FastaIndex::open(p, &dir.join("cache")).unwrap();
            Assembly {
                path: p.to_path_buf(),
                records: index.records(),
                index,
            }
        };
        let aln = Alignment {
            q_rec: 0,
            q_start: 0,
            q_end: 500,
            t_rec: 0,
            t_start: 1000,
            t_end: 1500,
            reverse: true,
            matches: 500,
            block_len: 500,
            mapq: 60,
            primary: true,
            divergence: None,
        };
        let data = GenomeData {
            x: asm(&r),
            y: asm(&q),
            self_comparison: false,
            alns: vec![aln],
            source: "test".into(),
        };
        let mut v = GenomeView::new(data, GenomeInputs::default());
        v.sort = sort;
        v.rearrange();
        v
    }

    #[test]
    fn window_raster_draws_the_alignment_where_it_lies() {
        let dir = std::env::temp_dir().join(format!("dottir-gv-{}-r", std::process::id()));
        let mut v = view(&dir, false);
        v.line_width = 1;
        // Whole plot: 2000 × 1000 bp onto 200 × 100 px. The reverse
        // alignment runs from (1500, 0) to (1000, 500) bp.
        let img = v.rasterize_window([200, 100]);
        let at = |x: usize, y: usize| img.pixels[y * 200 + x] != Color32::TRANSPARENT;
        assert!(at(150, 0) && at(125, 25) && at(100, 50));
        assert!(!at(50, 25) && !at(125, 75));
        // Zoomed to x 1000..1500, y 0..500: corner to corner.
        v.view = Window {
            x0: 1000.0,
            x1: 1500.0,
            y0: 0.0,
            y1: 500.0,
        };
        let img = v.rasterize_window([100, 100]);
        let at = |x: usize, y: usize| img.pixels[y * 100 + x] != Color32::TRANSPARENT;
        assert!(at(99, 1) && at(50, 50) && at(1, 99) && !at(10, 10));
    }

    #[test]
    fn drill_down_regions_follow_the_displayed_orientation() {
        let dir = std::env::temp_dir().join(format!("dottir-gv-{}-d", std::process::id()));
        // Sorted, ctg is drawn reverse-complemented.
        let mut v = view(&dir, true);
        assert!(v.arranged.y.segment(0).reverse);
        v.view = Window {
            x0: 900.0,
            x1: 1600.0,
            y0: 0.0,
            y1: 600.0,
        };
        assert!(v.view_fits_drill());
        let (x, y) = v.view_regions();
        assert_eq!(
            x,
            vec![Region {
                rec: 0,
                range: 900..1600,
                reverse: false
            }]
        );
        // Axis 0..600 of the reversed 1000 bp record is record 400..1000.
        assert_eq!(
            y,
            vec![Region {
                rec: 0,
                range: 400..1000,
                reverse: true
            }]
        );
        let q = v.data.x.index.region_sequence(&x).unwrap();
        let s = v.data.y.index.region_sequence(&y).unwrap();
        assert_eq!(&q.seq[..100], "A".repeat(100).as_bytes());
        // Reverse complement of G×100 T×500 = A×500 C×100: the aligned
        // part now runs forward, as on screen.
        assert_eq!(s.seq, [vec![b'A'; 500], vec![b'C'; 100]].concat());
        assert_eq!(s.records[0].id, "ctg:401-1000");
        v.drill_limit_mb = 0.0005;
        assert!(!v.view_fits_drill());
    }

    #[test]
    fn clipping_keeps_the_visible_part() {
        let c = clip((-10.0, 5.0), (20.0, 5.0), 10.0, 10.0, 0.0).unwrap();
        assert_eq!(c, ((0.0, 5.0), (10.0, 5.0)));
        assert!(clip((-10.0, -5.0), (-1.0, -1.0), 10.0, 10.0, 0.0).is_none());
        let d = clip((-5.0, -5.0), (15.0, 15.0), 10.0, 10.0, 1.0).unwrap();
        assert_eq!(d, ((-1.0, -1.0), (11.0, 11.0)));
    }

    #[test]
    fn distance_to_segment() {
        let (a, b) = (Pos2::new(0.0, 0.0), Pos2::new(10.0, 0.0));
        assert_eq!(segment_distance(Pos2::new(5.0, 3.0), a, b), 3.0);
        assert_eq!(segment_distance(Pos2::new(13.0, 4.0), a, b), 5.0);
        assert_eq!(segment_distance(Pos2::new(1.0, 1.0), a, a), 2f32.sqrt());
    }

    #[test]
    fn number_formats() {
        assert_eq!(thousands(1_234_567), "1,234,567");
        assert_eq!(thousands(12), "12");
        assert_eq!(span_text(2.5e6), "2.5 Mb");
        assert_eq!(span_text(1.2e9), "1.20 Gb");
    }
}
