//! Background compute worker (C2 part 1).
//!
//! Long-lived `std::thread` that owns the heavy dot-plot computation
//! and communicates with the egui app via `mpsc` channels. The UI
//! thread dispatches a `ComputeRequest` whenever a setting changes;
//! the worker computes (using rayon internally for parallelism) and
//! sends a `ComputeResult` back. The worker wakes the UI via the
//! supplied `repaint` closure so we don't have to spin on
//! `try_recv` every frame.
//!
//! Stale results — those whose `id` is older than the most recent
//! dispatched one — are discarded by the receiver. The worker itself
//! processes requests serially; the dispatcher in `app.rs` only ever
//! has at most one request "in flight" plus at most one "pending"
//! request buffered.

use std::sync::atomic::{AtomicBool, AtomicU32, AtomicU8, Ordering};
use std::sync::{mpsc, Arc};
use std::thread;

use dottir_core::large::{ClassMask, Control, LargeConfig, LargeInput, LargePlot, LargeSession};
use dottir_core::{compute_dotplot, DotPlot, DottirError, PlotConfig};

/// Inputs of the large engine, shared with the worker without copying
/// per request. `key` identifies the inputs: the worker re-prepares its
/// cached [`LargeSession`] when it changes.
pub struct LargeSource {
    pub key: u64,
    pub query: Arc<[u8]>,
    pub q_breaks: Vec<usize>,
    pub q_floor: Option<Vec<u16>>,
    /// `None` = self-comparison.
    pub subject: Option<Arc<[u8]>>,
    pub s_breaks: Vec<usize>,
    pub s_floor: Option<Vec<u16>>,
}

/// A large-engine job (ADR 0007). `config.viewport` carries the slice.
pub struct LargeJob {
    pub source: Arc<LargeSource>,
    pub config: LargeConfig,
    pub classes: ClassMask,
    pub separate_strands: bool,
}

/// A single computation job.
pub struct ComputeRequest {
    /// Monotonically increasing id assigned by the dispatcher.
    pub id: u64,
    pub query: Vec<u8>,
    pub subject: Vec<u8>,
    pub config: PlotConfig,
    /// When set, the large engine runs instead and `query`/`subject`/
    /// `config` are ignored.
    pub large: Option<LargeJob>,
}

/// Output of a single job.
pub struct ComputeResult {
    /// Echoes the request id so the receiver can discard stale results.
    pub id: u64,
    /// The `PlotConfig::zoom` that was used to compute this plot.
    /// Cached separately so the receiver can index its multi-resolution
    /// cache without going through `DotPlot::params`.
    pub config_zoom: u32,
    pub plot: Result<DotPlot, DottirError>,
    /// The full large-engine result (all class channels, chains), kept so
    /// class toggles can re-derive `plot` without recomputing.
    pub large: Option<LargePlot>,
}

/// What the worker is doing, for the status bar.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stage {
    Idle,
    Exact,
    Indexing,
    Large,
}

#[derive(Default)]
struct Progress {
    stage: AtomicU8,
    permille: AtomicU32,
}

impl Progress {
    fn set(&self, stage: Stage) {
        self.stage.store(stage as u8, Ordering::Relaxed);
        self.permille.store(0, Ordering::Relaxed);
    }
}

/// Handle to the worker thread held by the egui app.
pub struct ComputeWorker {
    request_tx: mpsc::Sender<ComputeRequest>,
    result_rx: mpsc::Receiver<ComputeResult>,
    progress: Arc<Progress>,
    cancel: Arc<AtomicBool>,
    /// Drop guard — JoinHandle keeps the OS thread alive for the
    /// app's lifetime. The thread exits when the request channel
    /// closes (i.e. when `ComputeWorker` is dropped).
    _handle: thread::JoinHandle<()>,
}

impl ComputeWorker {
    /// Spawn a fresh worker thread. `repaint` is called whenever a
    /// result is sent so the UI wakes up promptly; typically this is
    /// `let ctx = ctx.clone(); move || ctx.request_repaint()`.
    pub fn spawn<F>(repaint: F) -> Self
    where
        F: Fn() + Send + Sync + 'static,
    {
        let (request_tx, request_rx) = mpsc::channel::<ComputeRequest>();
        let (result_tx, result_rx) = mpsc::channel::<ComputeResult>();
        let progress = Arc::new(Progress::default());
        let cancel = Arc::new(AtomicBool::new(false));
        let (wp, wc) = (progress.clone(), cancel.clone());
        let handle = thread::Builder::new()
            .name("dottir-compute".into())
            .spawn(move || {
                let mut session: Option<(u64, LargeSession)> = None;
                while let Ok(mut req) = request_rx.recv() {
                    // Only the newest queued request matters; older ones
                    // would be discarded as stale by the receiver anyway.
                    while let Ok(newer) = request_rx.try_recv() {
                        req = newer;
                    }
                    wc.store(false, Ordering::Relaxed);
                    let zoom = req
                        .large
                        .as_ref()
                        .map_or(req.config.zoom, |j| j.config.zoom);
                    let (plot, large) = match req.large.take() {
                        Some(job) => run_large(job, &mut session, &wp, &wc, &repaint),
                        None => {
                            wp.set(Stage::Exact);
                            (compute_dotplot(&req.query, &req.subject, &req.config), None)
                        }
                    };
                    wp.set(Stage::Idle);
                    if result_tx
                        .send(ComputeResult {
                            id: req.id,
                            config_zoom: zoom,
                            plot,
                            large,
                        })
                        .is_err()
                    {
                        // Receiver dropped — exit cleanly.
                        break;
                    }
                    repaint();
                }
            })
            .expect("spawn compute worker thread");
        Self {
            request_tx,
            result_rx,
            progress,
            cancel,
            _handle: handle,
        }
    }

    /// Dispatch a new request. Never blocks. An in-flight large job is
    /// cancelled, since its result would be stale.
    pub fn dispatch(&self, req: ComputeRequest) {
        self.cancel.store(true, Ordering::Relaxed);
        // `send` only fails if the worker has exited, which we don't
        // expect during normal use. Treat as fatal-ish: log and drop.
        if let Err(e) = self.request_tx.send(req) {
            tracing::error!("compute worker dispatch failed: {e}");
        }
    }

    /// Current stage and progress fraction of the running job.
    pub fn progress(&self) -> (Stage, f32) {
        let stage = match self.progress.stage.load(Ordering::Relaxed) {
            1 => Stage::Exact,
            2 => Stage::Indexing,
            3 => Stage::Large,
            _ => Stage::Idle,
        };
        let f = self.progress.permille.load(Ordering::Relaxed) as f32 / 1000.0;
        (stage, f)
    }

    /// Non-blocking poll for completed results. Returns all that have
    /// arrived since the last call.
    pub fn drain_results(&self) -> Vec<ComputeResult> {
        let mut out = Vec::new();
        while let Ok(r) = self.result_rx.try_recv() {
            out.push(r);
        }
        out
    }
}

type LargeOutcome = (Result<DotPlot, DottirError>, Option<LargePlot>);

fn run_large(
    job: LargeJob,
    session: &mut Option<(u64, LargeSession)>,
    progress: &Progress,
    cancel: &AtomicBool,
    repaint: &(dyn Fn() + Sync),
) -> LargeOutcome {
    let report = |f: f32| {
        let pm = (f * 1000.0) as u32;
        if progress.permille.swap(pm, Ordering::Relaxed) / 50 != pm / 50 {
            repaint();
        }
    };
    let src = &job.source;
    let stale = match session {
        Some((key, s)) => *key != src.key || !s.compatible(&job.config),
        None => true,
    };
    if stale {
        *session = None;
        progress.set(Stage::Indexing);
        let input = LargeInput {
            query: &src.query,
            q_breaks: &src.q_breaks,
            q_floor: src.q_floor.as_deref(),
            subject: src.subject.as_deref(),
            s_breaks: if src.subject.is_some() {
                &src.s_breaks
            } else {
                &src.q_breaks
            },
            s_floor: src.s_floor.as_deref(),
        };
        // Preparing is not cancelled: the session survives the request
        // that triggered it and serves the next viewport.
        let ctl = Control {
            cancel: None,
            progress: Some(&report),
        };
        match LargeSession::prepare(&input, &job.config, &ctl) {
            Ok(s) => *session = Some((src.key, s)),
            Err(e) => return (Err(e), None),
        }
    }
    let (_, s) = session.as_ref().expect("prepared above");
    progress.set(Stage::Large);
    let ctl = Control {
        cancel: Some(cancel),
        progress: Some(&report),
    };
    match s.compute(&job.config, &ctl) {
        Ok(lp) => {
            let plot = lp.to_dot_plot(job.classes, job.separate_strands);
            (Ok(plot), Some(lp))
        }
        Err(e) => (Err(e), None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use dottir_core::{PlotConfig, ScoreMatrix, Strand};
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Arc;
    use std::time::Duration;

    /// The worker actually runs jobs and returns results. The repaint
    /// closure is called once per result.
    #[test]
    fn worker_runs_and_repaints() {
        let repaint_count = Arc::new(AtomicUsize::new(0));
        let rc = repaint_count.clone();
        let worker = ComputeWorker::spawn(move || {
            rc.fetch_add(1, Ordering::SeqCst);
        });

        let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
        cfg.strand = Strand::Forward;
        cfg.window_size = Some(5);
        cfg.zoom = 1;

        worker.dispatch(ComputeRequest {
            id: 1,
            query: b"ACGTACGTACGTACGTACGT".to_vec(),
            subject: b"ACGTACGTACGTACGTACGT".to_vec(),
            config: cfg.clone(),
            large: None,
        });

        // Poll for the result with a small timeout. Worker is single-
        // threaded and the job is trivial, so this completes quickly.
        let mut got = None;
        for _ in 0..200 {
            let mut r = worker.drain_results();
            if let Some(res) = r.pop() {
                got = Some(res);
                break;
            }
            thread::sleep(Duration::from_millis(5));
        }
        let r = got.expect("worker never produced a result");
        assert_eq!(r.id, 1);
        assert!(r.plot.is_ok());
        assert!(repaint_count.load(Ordering::SeqCst) >= 1);
    }

    /// Requests run in dispatch order; queued requests superseded by a
    /// newer one may be skipped, but the newest always runs.
    #[test]
    fn requests_processed_in_dispatch_order() {
        let worker = ComputeWorker::spawn(|| {});
        let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
        cfg.strand = Strand::Forward;
        cfg.window_size = Some(5);
        for id in 1..=3 {
            worker.dispatch(ComputeRequest {
                id,
                query: b"ACGTACGTACGTACGTACGT".to_vec(),
                subject: b"ACGTACGTACGTACGTACGT".to_vec(),
                config: cfg.clone(),
                large: None,
            });
        }
        let mut received = Vec::new();
        for _ in 0..200 {
            received.extend(worker.drain_results());
            if received.last().is_some_and(|r: &ComputeResult| r.id == 3) {
                break;
            }
            thread::sleep(Duration::from_millis(5));
        }
        let ids: Vec<u64> = received.iter().map(|r| r.id).collect();
        assert_eq!(ids.last(), Some(&3));
        assert!(ids.windows(2).all(|w| w[0] < w[1]), "{ids:?}");
    }

    /// A large job prepares a session once and reuses it for viewports.
    #[test]
    fn large_job_runs_and_reuses_session() {
        use dottir_core::large::{DepthMethod, Viewport};
        let worker = ComputeWorker::spawn(|| {});
        let unit = b"ACGTTGCAAGGCTTACGATCGATCGGATCCATGCATGCAAACCCGGGTTTAGCTAGCTTAGG";
        let seq: Arc<[u8]> = unit.iter().cycle().take(6000).copied().collect();
        let source = Arc::new(LargeSource {
            key: 7,
            query: seq,
            q_breaks: vec![],
            q_floor: None,
            subject: None,
            s_breaks: vec![],
            s_floor: None,
        });
        let mut config = LargeConfig {
            zoom: 10,
            self_comparison: true,
            depth: DepthMethod::None,
            ..LargeConfig::default()
        };
        for id in 1..=2 {
            if id == 2 {
                config.zoom = 2;
                config.viewport = Some(Viewport {
                    q: 1000..2000,
                    s: 1000..2000,
                });
            }
            worker.dispatch(ComputeRequest {
                id,
                query: vec![],
                subject: vec![],
                config: PlotConfig::default_blastn(ScoreMatrix::dna_identity()),
                large: Some(LargeJob {
                    source: source.clone(),
                    config: config.clone(),
                    classes: ClassMask::ALL,
                    separate_strands: false,
                }),
            });
            let mut got = None;
            for _ in 0..400 {
                if let Some(r) = worker.drain_results().pop() {
                    got = Some(r);
                    break;
                }
                thread::sleep(Duration::from_millis(5));
            }
            let r = got.expect("no result");
            let plot = r.plot.expect("large compute failed");
            assert!(r.large.is_some());
            let side = if id == 1 { 600 } else { 500 };
            assert_eq!((plot.width, plot.height), (side, side));
        }
    }
}
