//! Phase 1 integration tests for the BLASTN forward kernel.

use dottir_core::{compute_dotplot, PlotConfig, ScoreMatrix, Strand};

/// Self-comparison of a 40 bp sequence at zoom 1 must produce a maximal
/// hit on the main diagonal: every (i, i) pixel inside the valid range
/// [W, len) must hit 255 (since each window position sums to exactly W·5
/// from the +5 diagonal, and pixel_fac=50 / W=8 → 50·8·5/8 = 250 < 255
/// for some, but with a high pixel_fac it saturates).
#[test]
fn self_comparison_main_diagonal_lights_up() {
    let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(8);
    cfg.zoom = 1;
    cfg.pixel_fac = 80; // each W·5 window → 80·40/8 = 400 → clamped to 255

    let p = compute_dotplot(&seq, &seq, &cfg).unwrap();
    assert_eq!(p.width, seq.len() as u32);
    assert_eq!(p.height, seq.len() as u32);

    // Inspect the steady-state diagonal pixels (i = W..len).
    let w = cfg.window_size.unwrap() as usize;
    let n = seq.len();
    let mut diag_lit = 0;
    for i in w..n {
        let v = p.pixels[i * n + i];
        if v == 255 {
            diag_lit += 1;
        }
    }
    assert!(
        diag_lit > (n - w) / 2,
        "expected most main-diagonal pixels to saturate; got {diag_lit}/{}",
        n - w
    );
}

/// Completely unrelated random-ish sequences should produce mostly-dark
/// pixels (occasional spurious matches OK).
#[test]
fn unrelated_sequences_mostly_dark() {
    let q = b"AAAAAAAAAACCCCCCCCCC".to_vec(); // 20 nt
    let s = b"GGGGGGGGGGTTTTTTTTTT".to_vec(); // 20 nt — no shared bases
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(5);
    cfg.zoom = 1;

    let p = compute_dotplot(&q, &s, &cfg).unwrap();
    // No 5-mer of A-stretch or C-stretch matches any 5-mer of G-stretch
    // or T-stretch, so every pixel should be 0.
    assert!(
        p.pixels.iter().all(|&v| v == 0),
        "expected all-dark pixelmap; got nonzero at {:?}",
        p.pixels.iter().position(|&v| v != 0)
    );
}

/// Determinism: identical inputs and parameters MUST produce identical
/// pixelmaps across repeated invocations (spec §4.1.11).
#[test]
fn determinism_same_inputs_same_pixelmap() {
    let q = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT".repeat(3);
    let s = b"GTACACGTACGTGTACACGTACGTGTACACGTACGTGTACACGTACGT".repeat(3);
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(7);
    cfg.zoom = 4;

    let p1 = compute_dotplot(&q, &s, &cfg).unwrap();
    let p2 = compute_dotplot(&q, &s, &cfg).unwrap();
    assert_eq!(p1.pixels, p2.pixels);
}

/// Zoom > 1 reduces the output dimensions consistently.
#[test]
fn zoom_factor_reduces_dimensions() {
    let q = vec![b'A'; 100];
    let s = vec![b'A'; 80];
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(5);
    cfg.zoom = 10;

    let p = compute_dotplot(&q, &s, &cfg).unwrap();
    assert_eq!(p.width, 10);
    assert_eq!(p.height, 8);
    assert_eq!(p.pixels.len(), 80);
}

/// Empty inputs return an explicit error rather than panicking.
#[test]
fn empty_sequence_errors() {
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(5);
    assert!(compute_dotplot(b"", b"ACGT", &cfg).is_err());
    assert!(compute_dotplot(b"ACGT", b"", &cfg).is_err());
}

/// Memory limit refuses a huge pixelmap rather than allocating.
#[test]
fn memory_limit_enforced() {
    let q = vec![b'A'; 10_000];
    let s = vec![b'A'; 10_000];
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.strand = Strand::Forward;
    cfg.window_size = Some(5);
    cfg.zoom = 1;
    cfg.memory_limit_bytes = 1024; // way too small
    let err = compute_dotplot(&q, &s, &cfg).unwrap_err();
    let s = format!("{err}");
    assert!(s.contains("memory_limit"), "got: {s}");
}

/// BLASTX is still NotImplemented (Phase 2-extra).
#[test]
fn blastx_returns_not_implemented_error() {
    use dottir_core::BlastMode;
    let q = b"ACGTACGTACGT";
    let s = b"MKTAYIAKQRQI";
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.mode = BlastMode::Blastx;
    cfg.window_size = Some(3);
    assert!(compute_dotplot(q, s, &cfg).is_err());
}

/// BLASTP + reverse strand is meaningless and returns InvalidConfig.
#[test]
fn blastp_rejects_reverse_strand() {
    use dottir_core::BlastMode;
    let q = b"MKTAYIAKQRQI";
    let s = b"MAATKRIIRQRY";
    let mut cfg = PlotConfig::default_blastp(ScoreMatrix::blosum62());
    cfg.mode = BlastMode::Blastp;
    cfg.window_size = Some(3);
    cfg.strand = Strand::Reverse;
    assert!(compute_dotplot(q, s, &cfg).is_err());
}
