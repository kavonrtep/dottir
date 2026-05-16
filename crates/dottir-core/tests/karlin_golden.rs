//! Golden test: λ/K/H/window numerically identical to the C reference
//! `tests/golden_gen/karlin_ref.c` (a verbatim extraction of `karlin()` and
//! `winsizeFromlambdak()` from `dotterKarlin.c`).
//!
//! Tolerance is zero on the integer outputs (predicted, window) and `1e-15`
//! relative on the f64 outputs. The doubles will only differ if a deviation
//! from the C control flow was introduced; if a real algorithmic change is
//! made, regenerate the goldens AND bump `PIXELMAP_FORMAT_VERSION` (see
//! CLAUDE.md).

use dottir_core::karlin::{karlin_window_size, KarlinResult};
use dottir_core::matrix::{BlastMode, ScoreMatrix};

#[derive(Debug)]
struct Golden {
    name: &'static str,
    lambda: f64,
    k: f64,
    h: f64,
    exp_res: f64,
    exp_msp: f64,
    predicted: u32,
    window: u32,
}

/// Output captured from `tests/golden_gen/karlin_ref` on 2026-05-16.
/// See `tests/golden/karlin/README.md` for the regeneration command.
const GOLDENS: &[Golden] = &[
    Golden {
        name: "dna_uniform",
        lambda: 0.19152927398681641,
        k: 0.17334580950735284,
        h: 0.35672379416682848,
        exp_res: 1.8625027221237367,
        exp_msp: 38.938557203279061,
        predicted: 21,
        window: 21,
    },
    Golden {
        name: "dna_at_rich_vs_gc",
        lambda: 0.19408679008483887,
        k: 0.17628106821340828,
        h: 0.36590167746510743,
        exp_res: 1.8852477147216717,
        exp_msp: 38.5119708014956,
        predicted: 20,
        window: 20,
    },
    Golden {
        name: "dna_self_repeat",
        lambda: 0.19152927398681641,
        k: 0.17334580950735284,
        h: 0.35672379416682848,
        exp_res: 1.8625027221237369,
        exp_msp: 38.938557203279061,
        predicted: 21,
        window: 21,
    },
    Golden {
        name: "prot_uniform",
        lambda: 0.20744946599006653,
        k: 0.078393203331640751,
        h: 0.14503534878044433,
        exp_res: 0.69913580200485692,
        exp_msp: 32.125039669724075,
        predicted: 46,
        window: 46,
    },
];

fn dna_q1_long() -> Vec<u8> {
    let seed = b"ACGTACGTACGTACGTACGT";
    seed.repeat(40)
}

fn dna_s1_long() -> Vec<u8> {
    let seed = b"GTGTACGAGCATCGTCTACT";
    seed.repeat(40)
}

const DNA_AT_VS_GC_Q: &[u8] =
    b"AAAACCCGGGTTTAACAGCTAGCTACGATCGATCGATCGTAGCTAGCTAGCT";
const DNA_AT_VS_GC_S: &[u8] =
    b"TTTTGGGCCCAATTGCTAGCTACGATCGATCGATCGTAGCTAGCTAGCTACG";

const DNA_SELF: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

const PROT_Q1: &[u8] = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKAL\
                        RSALEFNAHVDEMVRLRREVGNQLEELQNRLREYIQRDHRGHEALQQYRVKQVHLDQEEIA";
const PROT_S1: &[u8] = b"MAATKRIIRQRYTIKHYVTRLREHIDHEEQVRKDLDEHKHRADRMLEELAGAILAAEHRLRD\
                        AREAFEQLLDKLEEHLRYAEELQEKFAKLERELAEHRLEEIEGRLAQAEEEFVEQHRRLENEL";

fn run_case(name: &str) -> KarlinResult {
    let (matrix, mode, q, s): (ScoreMatrix, BlastMode, Vec<u8>, Vec<u8>) = match name {
        "dna_uniform" => (
            ScoreMatrix::dna_identity(),
            BlastMode::Blastn,
            dna_q1_long(),
            dna_s1_long(),
        ),
        "dna_at_rich_vs_gc" => (
            ScoreMatrix::dna_identity(),
            BlastMode::Blastn,
            DNA_AT_VS_GC_Q.to_vec(),
            DNA_AT_VS_GC_S.to_vec(),
        ),
        "dna_self_repeat" => (
            ScoreMatrix::dna_identity(),
            BlastMode::Blastn,
            DNA_SELF.to_vec(),
            DNA_SELF.to_vec(),
        ),
        "prot_uniform" => (
            ScoreMatrix::blosum62(),
            BlastMode::Blastp,
            PROT_Q1.to_vec(),
            PROT_S1.to_vec(),
        ),
        other => panic!("unknown fixture {other}"),
    };
    karlin_window_size(&matrix, &q, &s, mode).expect("karlin failed")
}

/// Strict bit-identical check between two f64s, with a friendly message.
#[track_caller]
fn assert_f64_eq(actual: f64, expected: f64, name: &str, field: &str) {
    if actual.to_bits() != expected.to_bits() {
        let rel = ((actual - expected) / expected).abs();
        panic!(
            "[{name}] {field} differs: rust={actual:.17e}, c-ref={expected:.17e}, \
             relative diff = {rel:e}"
        );
    }
}

#[test]
fn matches_c_reference_bit_identical() {
    for g in GOLDENS {
        let r = run_case(g.name);
        assert_f64_eq(r.lambda, g.lambda, g.name, "lambda");
        assert_f64_eq(r.k, g.k, g.name, "k");
        assert_f64_eq(r.h, g.h, g.name, "h");
        assert_f64_eq(r.expected_residue_score, g.exp_res, g.name, "exp_res");
        assert_f64_eq(r.expected_msp_score, g.exp_msp, g.name, "exp_msp");
        assert_eq!(r.predicted_msp_length, g.predicted, "{}: predicted", g.name);
        assert_eq!(r.window_size, g.window, "{}: window", g.name);
    }
}
