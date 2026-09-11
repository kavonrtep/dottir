//! End-to-end `dottir batch` runs with GFF3-with-sequence positional
//! inputs. Corpus: `tests/corpora/gff3_with_sequence/`.

use std::path::{Path, PathBuf};
use std::process::Command;

fn corpus(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/corpora")
        .join(rel)
}

fn tmpdir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dottir_gff3_input_{}_{tag}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn batch(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_dottir"))
        .arg("batch")
        .args(args)
        .output()
        .expect("running dottir batch")
}

#[test]
fn self_comparison_from_one_gff3() {
    let dir = tmpdir("self");
    let out = dir.join("self.png");
    let result = batch(&[
        corpus("gff3_with_sequence/pair_a.gff3").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(out.exists(), "no PNG written");

    // The sidecar must describe the GFF3 file itself, and count the
    // records that came out of its embedded FASTA.
    let sidecar = std::fs::read_to_string(dir.join("self.png.params.toml")).unwrap();
    assert!(sidecar.contains("pair_a.gff3"), "{sidecar}");
    assert!(sidecar.contains("n_records = 2"), "{sidecar}");
    assert!(sidecar.contains("total_residues = 1200"), "{sidecar}");
    assert!(sidecar.contains("self_comparison = true"), "{sidecar}");

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn pairwise_from_two_gff3_files() {
    let dir = tmpdir("pair");
    let out = dir.join("pair.png");
    let result = batch(&[
        corpus("gff3_with_sequence/pair_a.gff3").to_str().unwrap(),
        corpus("gff3_with_sequence/pair_b.gff3").to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(out.exists(), "no PNG written");

    let sidecar = std::fs::read_to_string(dir.join("pair.png.params.toml")).unwrap();
    assert!(sidecar.contains("pair_a.gff3"), "{sidecar}");
    assert!(sidecar.contains("pair_b.gff3"), "{sidecar}");
    assert!(sidecar.contains("self_comparison = false"), "{sidecar}");

    std::fs::remove_dir_all(&dir).ok();
}

/// A GFF3 input and the equivalent separate FASTA must produce the
/// same pixelmap — the sequence path is the only thing that changed.
#[test]
fn gff3_input_matches_the_equivalent_fasta() {
    let dir = tmpdir("equiv");
    let from_gff = dir.join("gff.png");
    let from_fasta = dir.join("fasta.png");

    for (input, out) in [
        ("gff3_with_sequence/tir_simple.with_seq.gff3", &from_gff),
        ("annotation_overlay/tir_simple.fasta", &from_fasta),
    ] {
        let result = batch(&[
            corpus(input).to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--zoom",
            "8",
            "--width",
            "0",
            "--no-sidecar",
        ]);
        assert!(
            result.status.success(),
            "{input} stderr: {}",
            String::from_utf8_lossy(&result.stderr)
        );
    }

    assert_eq!(
        std::fs::read(&from_gff).unwrap(),
        std::fs::read(&from_fasta).unwrap(),
        "GFF3-sourced and FASTA-sourced pixelmaps differ"
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn gff3_without_sequences_fails_with_an_actionable_message() {
    let dir = tmpdir("nofasta");
    let result = batch(&[
        corpus("gff3_with_sequence/no_sequence.gff3")
            .to_str()
            .unwrap(),
        "-o",
        dir.join("x.png").to_str().unwrap(),
    ]);
    assert!(!result.status.success(), "expected a non-zero exit");
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(stderr.contains("no embedded sequences"), "{stderr}");
    assert!(stderr.contains("--gff-query"), "{stderr}");

    std::fs::remove_dir_all(&dir).ok();
}
