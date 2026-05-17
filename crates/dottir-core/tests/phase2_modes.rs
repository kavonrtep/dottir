//! Phase 2 integration tests: dual strand, BLASTP, self-comparison.

use dottir_core::{
    compute_dotplot, reverse_complement, BlastMode, PlotConfig, ScoreMatrix, Strand, Triangle,
};

/// Pick a sequence that is *not* self-reverse-complementary so the
/// forward vs reverse passes produce visibly different plots. ACGT
/// repeats are palindromic under reverse-complement, which makes them
/// useless for these tests.
const NONPAL_DNA: &[u8] = b"AAACCCGGGTAACTGAACCTTAGGCAAATTTGGCCAAGGTTACAACTGAACCTTAGGCAAATTTGGCC";

/// Both-strand BLASTN: query vs reverse_complement(query) should light
/// up the anti-diagonal (a reverse-strand self-match).
#[test]
fn both_strand_picks_up_reverse_match() {
    let q = NONPAL_DNA.to_vec();
    let s = reverse_complement(&q);
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Both;
    cfg.window_size = Some(8);
    cfg.zoom = 1;
    cfg.pixel_fac = 80;

    let p = compute_dotplot(&q, &s, &cfg).unwrap();
    let n = q.len();
    let w = cfg.window_size.unwrap() as usize;

    let mut anti_lit = 0;
    for i in w..n - w {
        if p.pixels[(n - 1 - i) * n + i] > 0 {
            anti_lit += 1;
        }
    }
    assert!(
        anti_lit > (n - 2 * w) / 3,
        "reverse pass didn't light up the anti-diagonal: {anti_lit}/{}",
        n - 2 * w
    );
}

/// Strand::Reverse on a non-palindromic self pair: the main diagonal
/// should mostly stay dark because reverse-strand matches of `q` against
/// itself only fire where `q` contains palindromic substrings, which
/// NONPAL_DNA avoids.
#[test]
fn reverse_only_skips_forward_diagonal() {
    let q = NONPAL_DNA.to_vec();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Reverse;
    cfg.window_size = Some(8);
    cfg.zoom = 1;

    let p = compute_dotplot(&q, &q, &cfg).unwrap();
    let n = q.len();
    let w = cfg.window_size.unwrap() as usize;

    let mut diag_lit = 0;
    for i in w..n {
        if p.pixels[i * n + i] > 0 {
            diag_lit += 1;
        }
    }
    assert!(
        diag_lit < (n - w) / 4,
        "reverse-only lit too many forward-diagonal pixels: {diag_lit}/{}",
        n - w
    );
}

/// BLASTP on a protein pair (identical inputs) lights up the main
/// diagonal at scores ≈ matrix diagonal sums.
#[test]
fn blastp_self_comparison_diagonal() {
    let p = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKAL".to_vec();
    let mut cfg = PlotConfig::default_blastp(ScoreMatrix::blosum62());
    cfg.mode = BlastMode::Blastp;
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(5);
    cfg.zoom = 1;
    cfg.pixel_fac = 80;

    let plot = compute_dotplot(&p, &p, &cfg).unwrap();
    let n = p.len();
    let w = cfg.window_size.unwrap() as usize;

    let mut diag_hits = 0;
    for i in w..n {
        if plot.pixels[i * n + i] > 0 {
            diag_hits += 1;
        }
    }
    // BLOSUM62 self-pair scores are heterogeneous; some windows will sum
    // strongly positive, others not. We just check we get a meaningful
    // signal — at least half the diagonal pixels are non-zero.
    assert!(
        diag_hits > (n - w) / 2,
        "BLASTP self diagonal too dark: {diag_hits}/{}",
        n - w
    );
}

/// Self-comparison with Triangle::Both mirrors the (kernel-filled) lower
/// triangle into the upper one. The resulting plot is fully symmetric.
#[test]
fn self_comparison_both_mirror_is_symmetric() {
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(6);
    cfg.zoom = 1;
    cfg.pixel_fac = 80;
    cfg.self_comparison = true;
    cfg.triangle = Triangle::Both;

    let p = compute_dotplot(&seq, &seq, &cfg).unwrap();
    let n = seq.len();
    let mut diffs = 0;
    for s in 0..n {
        for q in 0..n {
            if p.pixels[s * n + q] != p.pixels[q * n + s] {
                diffs += 1;
            }
        }
    }
    assert_eq!(diffs, 0, "Triangle::Both should produce a symmetric plot");
}

/// Triangle::Upper copies lower→upper then zeros the lower triangle, so
/// strictly-lower-triangle pixels (row > col, i.e. q < s) are all zero
/// after post-processing.
#[test]
fn self_comparison_upper_zeros_lower_triangle() {
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(6);
    cfg.zoom = 1;
    cfg.self_comparison = true;
    cfg.triangle = Triangle::Upper;

    let p = compute_dotplot(&seq, &seq, &cfg).unwrap();
    let n = seq.len();
    for s in 1..n {
        for q in 0..s {
            assert_eq!(
                p.pixels[s * n + q],
                0,
                "lower-triangle pixel (s={s}, q={q}) should be zero after Triangle::Upper"
            );
        }
    }
}

/// Triangle::Lower (the kernel's natural output) leaves the upper
/// triangle untouched (zero). Verify the upper triangle is all-zero.
#[test]
fn self_comparison_lower_leaves_upper_zero() {
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(6);
    cfg.zoom = 1;
    cfg.self_comparison = true;
    cfg.triangle = Triangle::Lower;

    let p = compute_dotplot(&seq, &seq, &cfg).unwrap();
    let n = seq.len();
    for s in 0..n {
        for q in (s + 1)..n {
            assert_eq!(
                p.pixels[s * n + q],
                0,
                "upper-triangle pixel (s={s}, q={q}) should be zero with Triangle::Lower"
            );
        }
    }
}

/// disable_mirror short-circuits all post-processing, regardless of
/// `triangle`. With Triangle::Both + disable_mirror, the upper triangle
/// remains as the kernel left it — which for self-comparison is all
/// zeros (the kernel caps qmax at s+1).
#[test]
fn self_comparison_disable_mirror_skips_postprocess() {
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(6);
    cfg.zoom = 1;
    cfg.self_comparison = true;
    cfg.triangle = Triangle::Both;
    cfg.disable_mirror = true;

    let p = compute_dotplot(&seq, &seq, &cfg).unwrap();
    let n = seq.len();

    // Kernel only fills q <= s. With mirror disabled, the upper
    // triangle must remain zero even though Triangle::Both was set.
    for s in 0..n {
        for q in (s + 1)..n {
            assert_eq!(p.pixels[s * n + q], 0, "(s={s}, q={q}) should be zero");
        }
    }
}

/// Separate strand channels: forward_pixels and reverse_pixels populated,
/// and their element-wise max equals pixels.
#[test]
fn separate_strand_channels_split_and_recombine() {
    let q = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let s = reverse_complement(&q);
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Both;
    cfg.separate_strand_channels = true;
    cfg.window_size = Some(6);
    cfg.zoom = 1;

    let p = compute_dotplot(&q, &s, &cfg).unwrap();
    let fwd = p.forward_pixels.as_ref().expect("forward_pixels populated");
    let rev = p.reverse_pixels.as_ref().expect("reverse_pixels populated");
    assert_eq!(fwd.len(), p.pixels.len());
    assert_eq!(rev.len(), p.pixels.len());

    for (i, &v) in p.pixels.iter().enumerate() {
        assert_eq!(v, fwd[i].max(rev[i]), "combined != max(fwd, rev) at {i}");
    }
}

/// Self-comparison size mismatch is rejected before computation.
#[test]
fn self_comparison_with_unequal_lengths_errors() {
    let q = b"ACGTACGT".to_vec();
    let s = b"ACGTACG".to_vec();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.self_comparison = true;
    cfg.window_size = Some(4);
    let err = compute_dotplot(&q, &s, &cfg).unwrap_err();
    let msg = format!("{err}");
    assert!(msg.contains("self_comparison"), "got: {msg}");
}

/// F2: `reverse_query = true` is equivalent to passing
/// `reverse_complement(query)` and the default config — it pre-flips
/// the query axis before computation. So Forward+RC(q) vs s should
/// match Forward+q vs s where q has been externally reverse-
/// complemented.
#[test]
fn reverse_query_flag_matches_external_revcomp() {
    let q = b"AAACCCGGGTAACTGAACCTTAGGCAAATTTGGCCAAGGTTACAACTGAACC".to_vec();
    let s = q.clone();

    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(8);
    cfg.zoom = 1;

    // Path A: flag.
    let mut cfg_a = cfg.clone();
    cfg_a.reverse_query = true;
    let plot_a = compute_dotplot(&q, &s, &cfg_a).unwrap();

    // Path B: external revcomp.
    let q_rc = reverse_complement(&q);
    let plot_b = compute_dotplot(&q_rc, &s, &cfg).unwrap();

    assert_eq!(plot_a.pixels, plot_b.pixels);
}

/// F2: BLASTP ignores reverse_query/reverse_subject (proteins have no
/// reverse strand). Setting either flag must not change the result.
#[test]
fn rev_flags_are_noop_for_blastp() {
    let q = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQ".to_vec();
    let s = q.clone();
    use dottir_core::BlastMode;
    let mut cfg = PlotConfig::default_blastp(ScoreMatrix::blosum62());
    cfg.mode = BlastMode::Blastp;
    cfg.window_size = Some(5);

    let plot_off = compute_dotplot(&q, &s, &cfg).unwrap();
    cfg.reverse_query = true;
    cfg.reverse_subject = true;
    let plot_on = compute_dotplot(&q, &s, &cfg).unwrap();
    assert_eq!(plot_off.pixels, plot_on.pixels);
}

/// Determinism across thread counts placeholder — for Phase 3 we'll
/// re-run this with rayon configured. For now, single-threaded
/// repetition is enough.
#[test]
fn determinism_extended_modes() {
    let q = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let s = reverse_complement(&q);
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Both;
    cfg.window_size = Some(7);
    cfg.zoom = 2;
    cfg.separate_strand_channels = true;
    let p1 = compute_dotplot(&q, &s, &cfg).unwrap();
    let p2 = compute_dotplot(&q, &s, &cfg).unwrap();
    assert_eq!(p1.pixels, p2.pixels);
    assert_eq!(p1.forward_pixels, p2.forward_pixels);
    assert_eq!(p1.reverse_pixels, p2.reverse_pixels);
}
