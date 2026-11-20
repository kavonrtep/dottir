//! `compute_dotplot_with`: identical pixels without cancellation,
//! `Cancelled` with it, monotone progress ending at 1.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Mutex;

use dottir_core::{
    compute_dotplot, compute_dotplot_with, BlastMode, Control, DottirError, PlotConfig,
    ScoreMatrix, Strand,
};

fn dna(n: usize, mut state: u64) -> Vec<u8> {
    (0..n)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 62) as usize]
        })
        .collect()
}

fn configs() -> Vec<(PlotConfig, Vec<u8>, Vec<u8>)> {
    let q = dna(6000, 1);
    let mut s = q[1000..5000].to_vec();
    s.extend(dna(3000, 2));
    let mut out = Vec::new();
    for strand in [Strand::Forward, Strand::Reverse, Strand::Both] {
        let mut c = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
        c.strand = strand;
        c.zoom = 4;
        out.push((c, q.clone(), s.clone()));
    }
    let mut c = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    c.self_comparison = true;
    c.strand = Strand::Both;
    c.zoom = 3;
    out.push((c, q.clone(), q.clone()));
    let mut c = PlotConfig::default_blastp(ScoreMatrix::blosum62());
    c.mode = BlastMode::Blastx;
    c.zoom = 2;
    let pep: Vec<u8> = b"MKVLAAGIVGLLLAAAQPRSTWYHEDNC".repeat(40);
    out.push((c, q[..3000].to_vec(), pep));
    out
}

#[test]
fn uncancelled_run_matches_compute_dotplot() {
    for (cfg, q, s) in configs() {
        let flag = AtomicBool::new(false);
        let seen = Mutex::new(Vec::new());
        let progress = |f: f32| seen.lock().unwrap().push(f);
        let ctl = Control {
            cancel: Some(&flag),
            progress: Some(&progress),
        };
        let a = compute_dotplot(&q, &s, &cfg).unwrap();
        let b = compute_dotplot_with(&q, &s, &cfg, &ctl).unwrap();
        assert_eq!(a.pixels, b.pixels, "{:?} {:?}", cfg.mode, cfg.strand);
        assert_eq!(a.reverse_pixels, b.reverse_pixels);
        let seen = seen.into_inner().unwrap();
        let last = *seen.last().expect("progress reported");
        assert!((last - 1.0).abs() < 0.05, "final progress {last}");
        assert!(seen.iter().all(|f| (0.0..=1.0).contains(f)));
    }
}

#[test]
fn cancel_before_start_returns_cancelled() {
    let (cfg, q, s) = configs().remove(0);
    let flag = AtomicBool::new(true);
    let ctl = Control {
        cancel: Some(&flag),
        progress: None,
    };
    assert!(matches!(
        compute_dotplot_with(&q, &s, &cfg, &ctl),
        Err(DottirError::Cancelled)
    ));
}

#[test]
fn cancel_mid_run_returns_cancelled() {
    // Large enough for several progress ticks (2^24 cells each).
    let q = dna(20_000, 3);
    let s = dna(20_000, 4);
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Both;
    cfg.zoom = 10;
    let flag = AtomicBool::new(false);
    let progress = |f: f32| {
        if f > 0.1 {
            flag.store(true, Ordering::Relaxed);
        }
    };
    let ctl = Control {
        cancel: Some(&flag),
        progress: Some(&progress),
    };
    assert!(matches!(
        compute_dotplot_with(&q, &s, &cfg, &ctl),
        Err(DottirError::Cancelled)
    ));
}
