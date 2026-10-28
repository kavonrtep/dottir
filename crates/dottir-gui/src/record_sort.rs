//! Background record sorting (the D-GENIES order of genome mode) for
//! ordinary dotplots: whole-plot large-engine chains on a thread of its
//! own, then `dottir_core::genome::sort_by_chains`.

use std::sync::atomic::{AtomicBool, AtomicU32, Ordering};
use std::sync::{mpsc, Arc};
use std::time::{Duration, Instant};

use dottir_core::genome::{sort_by_chains, InputRecord, QueryOrder, SortAxis};
use dottir_core::large::{whole_plot_chains, Chain, Control, LargeConfig, LargeInput};
use dottir_io::Sequence;

use crate::compute_worker::LargeSource;

pub struct SortJob {
    pub axis: SortAxis,
    rx: mpsc::Receiver<Result<Vec<Chain>, String>>,
    cancel: Arc<AtomicBool>,
    permille: Arc<AtomicU32>,
    started: Instant,
}

impl SortJob {
    pub fn start(
        source: Arc<LargeSource>,
        config: LargeConfig,
        axis: SortAxis,
        repaint: impl Fn() + Send + Sync + 'static,
    ) -> Self {
        let (tx, rx) = mpsc::channel();
        let cancel = Arc::new(AtomicBool::new(false));
        let permille = Arc::new(AtomicU32::new(0));
        let (c, p) = (cancel.clone(), permille.clone());
        let spawned = std::thread::Builder::new()
            .name("dottir-sort".into())
            .spawn(move || {
                let report = |f: f32| {
                    let pm = (f * 1000.0) as u32;
                    if p.swap(pm, Ordering::Relaxed) / 20 != pm / 20 {
                        repaint();
                    }
                };
                let input = LargeInput {
                    query: &source.query,
                    q_breaks: &source.q_breaks,
                    q_floor: source.q_floor.as_deref(),
                    subject: source.subject.as_deref(),
                    s_breaks: &source.s_breaks,
                    s_floor: source.s_floor.as_deref(),
                };
                let ctl = Control {
                    cancel: Some(&c),
                    progress: Some(&report),
                };
                let r = whole_plot_chains(&input, &config, &ctl).map_err(|e| e.to_string());
                let _ = tx.send(r);
                repaint();
            });
        if let Err(e) = spawned {
            tracing::error!("record sort: {e}");
        }
        Self {
            axis,
            rx,
            cancel,
            permille,
            started: Instant::now(),
        }
    }

    pub fn poll(&self) -> Option<Result<Vec<Chain>, String>> {
        match self.rx.try_recv() {
            Ok(r) => Some(r),
            Err(mpsc::TryRecvError::Empty) => None,
            Err(mpsc::TryRecvError::Disconnected) => Some(Err("record sort stopped".into())),
        }
    }

    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::Relaxed);
    }

    pub fn progress(&self) -> (f32, Duration) {
        (
            self.permille.load(Ordering::Relaxed) as f32 / 1000.0,
            self.started.elapsed(),
        )
    }
}

/// The order of `axis`'s records from `chains` (query = horizontal).
pub fn order(chains: &[Chain], query: &Sequence, subject: &Sequence, axis: SortAxis) -> QueryOrder {
    let recs = |s: &Sequence| -> Vec<InputRecord> {
        s.records
            .iter()
            .map(|r| InputRecord {
                name: r.id.clone(),
                range: r.range.clone(),
            })
            .collect()
    };
    sort_by_chains(chains, &recs(query), &recs(subject), axis)
}
