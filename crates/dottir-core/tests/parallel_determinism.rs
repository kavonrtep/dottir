//! Phase 3: byte-identical output across thread counts.
//!
//! Per spec §4.1.11 and CLAUDE.md, the rayon-chunked driver MUST produce
//! pixel-identical output regardless of how many threads do the work.
//! These tests configure rayon to 1, 2, 4, and 8 threads and assert all
//! produce the same pixelmap.

use dottir_core::{compute_dotplot, BlastMode, PlotConfig, ScoreMatrix, Strand};

fn run_with_thread_pool<F, R>(n: usize, f: F) -> R
where
    F: FnOnce() -> R + Send,
    R: Send,
{
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .expect("build thread pool");
    pool.install(f)
}

fn make_input() -> (Vec<u8>, Vec<u8>) {
    // Large enough to actually trigger chunking (the driver only splits
    // when slen >= window*64 and there are >1 threads).
    let q = b"AAACCCGGGTAACTGAACCTTAGGCAAATTTGGCCAAGGTTACAACTGAACCTTAGGCAAATTTGGCC"
        .repeat(64); // ~4 kb
    let s = b"TAACTGAACCTTAGGCAAATTTGGCCAAGGTTAAACCCGGGTAACTGAACCTTAGGCAAATTTGGCCAA"
        .repeat(64);
    (q, s)
}

#[test]
fn blastn_forward_byte_identical_across_thread_counts() {
    let (q, s) = make_input();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(10);
    cfg.zoom = 1;

    let baseline = run_with_thread_pool(1, || compute_dotplot(&q, &s, &cfg).unwrap());
    for n in [2_usize, 4, 8] {
        let p = run_with_thread_pool(n, || compute_dotplot(&q, &s, &cfg).unwrap());
        assert_eq!(
            p.pixels, baseline.pixels,
            "Strand::Forward differs at n_threads={n}"
        );
    }
}

#[test]
fn blastn_both_strands_byte_identical_across_thread_counts() {
    let (q, s) = make_input();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Both;
    cfg.window_size = Some(10);
    cfg.zoom = 1;
    cfg.separate_strand_channels = true;

    let baseline = run_with_thread_pool(1, || compute_dotplot(&q, &s, &cfg).unwrap());
    for n in [2_usize, 4, 8] {
        let p = run_with_thread_pool(n, || compute_dotplot(&q, &s, &cfg).unwrap());
        assert_eq!(p.pixels, baseline.pixels, "combined differs at n={n}");
        assert_eq!(p.forward_pixels, baseline.forward_pixels, "fwd differs at n={n}");
        assert_eq!(p.reverse_pixels, baseline.reverse_pixels, "rev differs at n={n}");
    }
}

#[test]
fn blastp_byte_identical_across_thread_counts() {
    let q: Vec<u8> =
        b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKAL\
          RSALEFNAHVDEMVRLRREVGNQLEELQNRLREYIQRDHRGHEALQQYRVKQVHLDQEEIA"
            .repeat(20);
    let s = q.clone();
    let mut cfg = PlotConfig::default_blastp(ScoreMatrix::blosum62());
    cfg.mode = BlastMode::Blastp;
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(8);
    cfg.zoom = 2;

    let baseline = run_with_thread_pool(1, || compute_dotplot(&q, &s, &cfg).unwrap());
    for n in [2_usize, 4, 8] {
        let p = run_with_thread_pool(n, || compute_dotplot(&q, &s, &cfg).unwrap());
        assert_eq!(p.pixels, baseline.pixels, "BLASTP differs at n_threads={n}");
    }
}

/// Self-comparison + Both strand + separate channels — the full
/// machinery — must still be deterministic across thread counts.
#[test]
fn full_machinery_byte_identical_across_thread_counts() {
    use dottir_core::Triangle;
    let q = b"AAACCCGGGTAACTGAACCTTAGGCAAATTTGGCCAAGGTTACAACTGAACCTTAGGCAAATTTGGCC"
        .repeat(32);
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Both;
    cfg.window_size = Some(8);
    cfg.zoom = 1;
    cfg.self_comparison = true;
    cfg.triangle = Triangle::Both;
    cfg.separate_strand_channels = true;

    let baseline = run_with_thread_pool(1, || compute_dotplot(&q, &q, &cfg).unwrap());
    for n in [2_usize, 4, 8] {
        let p = run_with_thread_pool(n, || compute_dotplot(&q, &q, &cfg).unwrap());
        assert_eq!(p.pixels, baseline.pixels, "differs at n_threads={n}");
    }
}
