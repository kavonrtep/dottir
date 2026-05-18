//! Phase A1 regression: parallel compute must not multiply the
//! pixelmap allocation by the chunk count.
//!
//! Before A1, `run_pass` allocated `n_chunks + 1` full pixelmaps and
//! collected them before merging. With `target_chunks = 4 × n_threads`
//! that's up to ~33 maps for 8 threads — a hard violation of
//! `memory_limit_bytes`. After A1 there's exactly one pixelmap per
//! strand for the whole pass.
//!
//! These tests don't measure RSS — they exercise behaviour that's only
//! possible under the new design:
//!
//! 1. A `memory_limit_bytes` sized to fit *exactly one* pixelmap must
//!    succeed at any thread count. The old design would
//!    `OutOfMemory` for every per-chunk allocation past the first.
//! 2. Sweeping thread counts at a tight memory cap that's well below
//!    `2 × W × H` still produces byte-identical output.

use dottir_core::{compute_dotplot, BlastMode, DottirError, PlotConfig, ScoreMatrix, Strand};

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

/// A non-trivial input that's large enough to trigger chunking
/// (the driver only splits when slen >= window*64).
fn make_input() -> (Vec<u8>, Vec<u8>) {
    let q = b"AAACCCGGGTAACTGAACCTTAGGCAAATTTGGCCAAGGTTACAACTGAACCTTAGGCAAATTTGGCC".repeat(48); // ~3.2 kb
    let s = b"TAACTGAACCTTAGGCAAATTTGGCCAAGGTTAAACCCGGGTAACTGAACCTTAGGCAAATTTGGCCAA".repeat(48);
    (q, s)
}

/// A `memory_limit_bytes` set to *exactly* one pixelmap worth must
/// succeed at every thread count. If A1 regresses and per-chunk
/// allocation comes back, this fires `OutOfMemory` on n ≥ 2.
#[test]
fn one_pixelmap_budget_succeeds_at_every_thread_count() {
    let (q, s) = make_input();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(8);
    cfg.zoom = 1;
    // For this input width=qlen=3264, height=slen=3264, so the pixelmap
    // is ~10.4 MiB. Allocate exactly that (no slack).
    let pixelmap_bytes = (q.len() as u64) * (s.len() as u64);
    cfg.memory_limit_bytes = pixelmap_bytes;

    for n in [1_usize, 2, 4, 8] {
        let result = run_with_thread_pool(n, || compute_dotplot(&q, &s, &cfg));
        assert!(
            result.is_ok(),
            "n_threads={n} hit the cap even though one pixelmap fits: {:?}",
            result.err()
        );
    }
}

/// At the same tight cap, the output is byte-identical across thread
/// counts — proves the shared-pixelmap atomic max-merge converges to
/// the serial value.
#[test]
fn tight_budget_still_byte_identical() {
    let (q, s) = make_input();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(8);
    cfg.zoom = 1;
    cfg.memory_limit_bytes = (q.len() as u64) * (s.len() as u64);

    let baseline = run_with_thread_pool(1, || compute_dotplot(&q, &s, &cfg).unwrap());
    for n in [2_usize, 4, 8] {
        let p = run_with_thread_pool(n, || compute_dotplot(&q, &s, &cfg).unwrap());
        assert_eq!(
            p.pixels, baseline.pixels,
            "n_threads={n} produced a different pixelmap under the tight cap"
        );
    }
}

/// Both-strand BLASTN with `separate_strand_channels` allocates *two*
/// pixelmaps (forward + reverse). A cap of exactly 2 × W × H must
/// suffice at every thread count.
#[test]
fn both_strands_fit_in_two_pixelmaps() {
    let (q, s) = make_input();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Both;
    cfg.window_size = Some(8);
    cfg.zoom = 1;
    cfg.separate_strand_channels = true;
    let pixelmap_bytes = (q.len() as u64) * (s.len() as u64);
    cfg.memory_limit_bytes = 2 * pixelmap_bytes;

    for n in [1_usize, 2, 4, 8] {
        let result = run_with_thread_pool(n, || compute_dotplot(&q, &s, &cfg));
        assert!(
            result.is_ok(),
            "n_threads={n} hit the 2× cap with separate channels: {:?}",
            result.err()
        );
    }
}

/// Halving the cap to one pixelmap should still let single-channel
/// both-strand work (forward+reverse max-merge into one pixelmap).
#[test]
fn shared_channel_both_strands_fits_in_one_pixelmap() {
    let (q, s) = make_input();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Both;
    cfg.window_size = Some(8);
    cfg.zoom = 1;
    cfg.separate_strand_channels = false;
    cfg.memory_limit_bytes = (q.len() as u64) * (s.len() as u64);

    for n in [1_usize, 2, 4, 8] {
        let result = run_with_thread_pool(n, || compute_dotplot(&q, &s, &cfg));
        assert!(
            result.is_ok(),
            "n_threads={n} hit the cap with shared-channel both-strand: {:?}",
            result.err()
        );
    }
}

/// And the inverse: a cap that's *smaller* than one pixelmap should
/// fail uniformly across thread counts. (Sanity check on
/// `new_checked`.)
#[test]
fn undersized_budget_fails_uniformly() {
    let (q, s) = make_input();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(8);
    cfg.zoom = 1;
    let pixelmap_bytes = (q.len() as u64) * (s.len() as u64);
    cfg.memory_limit_bytes = pixelmap_bytes - 1;

    for n in [1_usize, 2, 4, 8] {
        let result = run_with_thread_pool(n, || compute_dotplot(&q, &s, &cfg));
        assert!(
            result.is_err(),
            "n_threads={n} succeeded with a cap one byte too small",
        );
    }

    // BLASTP also flows through the same allocator path.
    let mut cfg2 = PlotConfig::default_blastp(ScoreMatrix::blosum62());
    cfg2.mode = BlastMode::Blastp;
    cfg2.window_size = Some(5);
    cfg2.memory_limit_bytes = 10; // way too small
    let r = run_with_thread_pool(4, || {
        compute_dotplot(b"MKTAYIAKQRQI", b"MAATKRIIRQRY", &cfg2)
    });
    assert!(r.is_err());
}

/// D: the upfront budget check accounts for *all* retained channels.
/// `separate_strand_channels=true` with both strands needs 2 × W × H,
/// so a cap that's enough for one pixelmap but less than two must
/// reject at the entry, with an error that reports `channels = 2`.
#[test]
fn upfront_check_rejects_separate_channels_at_one_pixelmap_cap() {
    let (q, s) = make_input();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Both;
    cfg.window_size = Some(8);
    cfg.zoom = 1;
    cfg.separate_strand_channels = true;
    let pixelmap_bytes = (q.len() as u64) * (s.len() as u64);
    cfg.memory_limit_bytes = pixelmap_bytes; // enough for ONE channel, not two

    let err = compute_dotplot(&q, &s, &cfg).expect_err("should reject upfront");
    match err {
        DottirError::OutOfMemory {
            requested,
            per_channel,
            channels,
            limit,
        } => {
            assert_eq!(per_channel, pixelmap_bytes);
            assert_eq!(channels, 2);
            assert_eq!(requested, 2 * pixelmap_bytes);
            assert_eq!(limit, pixelmap_bytes);
        }
        other => panic!("unexpected error variant: {other:?}"),
    }
}
