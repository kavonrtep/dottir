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

use std::sync::mpsc;
use std::thread;

use dottir_core::{compute_dotplot, DotPlot, DottirError, PlotConfig};

/// A single computation job.
pub struct ComputeRequest {
    /// Monotonically increasing id assigned by the dispatcher.
    pub id: u64,
    pub query: Vec<u8>,
    pub subject: Vec<u8>,
    pub config: PlotConfig,
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
}

/// Handle to the worker thread held by the egui app.
pub struct ComputeWorker {
    request_tx: mpsc::Sender<ComputeRequest>,
    result_rx: mpsc::Receiver<ComputeResult>,
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
        F: Fn() + Send + 'static,
    {
        let (request_tx, request_rx) = mpsc::channel::<ComputeRequest>();
        let (result_tx, result_rx) = mpsc::channel::<ComputeResult>();
        let handle = thread::Builder::new()
            .name("dottir-compute".into())
            .spawn(move || {
                while let Ok(req) = request_rx.recv() {
                    let plot = compute_dotplot(&req.query, &req.subject, &req.config);
                    let zoom = req.config.zoom;
                    if result_tx
                        .send(ComputeResult {
                            id: req.id,
                            config_zoom: zoom,
                            plot,
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
        Self { request_tx, result_rx, _handle: handle }
    }

    /// Dispatch a new request. Returns the id assigned. Never blocks.
    pub fn dispatch(&self, req: ComputeRequest) {
        // `send` only fails if the worker has exited, which we don't
        // expect during normal use. Treat as fatal-ish: log and drop.
        if let Err(e) = self.request_tx.send(req) {
            tracing::error!("compute worker dispatch failed: {e}");
        }
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

    /// Multiple requests run in order; ids are echoed correctly.
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
            });
        }
        let mut received = Vec::new();
        for _ in 0..200 {
            received.extend(worker.drain_results());
            if received.len() >= 3 {
                break;
            }
            thread::sleep(Duration::from_millis(5));
        }
        assert_eq!(received.len(), 3);
        let ids: Vec<u64> = received.iter().map(|r| r.id).collect();
        assert_eq!(ids, vec![1, 2, 3]);
    }
}
