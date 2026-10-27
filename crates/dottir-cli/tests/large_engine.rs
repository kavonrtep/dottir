//! End-to-end `dottir batch --engine large` runs on generated inputs.

use std::fmt::Write as _;
use std::path::PathBuf;
use std::process::Command;

fn tmpdir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dottir_large_{}_{tag}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn dna(n: usize, mut state: u64) -> String {
    (0..n)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 62) as usize] as char
        })
        .collect()
}

fn batch(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_dottir"))
        .arg("batch")
        .args(args)
        .output()
        .expect("running dottir batch")
}

#[test]
fn pairwise_large_writes_png_sidecar_and_chains() {
    let dir = tmpdir("pair");
    // Two records per genome sharing a 20 kb block; the subject's copy of
    // a 3 kb element is soft-masked.
    let block = dna(20_000, 1);
    let te = dna(3000, 2);
    let mut q = String::new();
    writeln!(q, ">q1\n{}{}{}", dna(5000, 3), block, dna(5000, 4)).unwrap();
    writeln!(q, ">q2\n{}", dna(30_000, 5)).unwrap();
    let mut s = String::new();
    writeln!(s, ">s1\n{}{}", dna(8000, 6), block).unwrap();
    writeln!(s, ">s2\n{}{}", te.to_lowercase(), dna(10_000, 7)).unwrap();
    let (qp, sp) = (dir.join("q.fa"), dir.join("s.fa"));
    std::fs::write(&qp, q).unwrap();
    std::fs::write(&sp, s).unwrap();
    let out = dir.join("pair.png");
    let paf = dir.join("pair.paf");
    let result = batch(&[
        qp.to_str().unwrap(),
        sp.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--engine",
        "large",
        "--zoom",
        "50",
        "--hide-classes",
        "high",
        "--chains-out",
        paf.to_str().unwrap(),
    ]);
    assert!(
        result.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(out.exists());

    let sidecar = std::fs::read_to_string(dir.join("pair.png.params.toml")).unwrap();
    assert!(sidecar.contains("[large]"), "{sidecar}");
    assert!(
        sidecar.contains("seed = \"1110100110010101111\""),
        "{sidecar}"
    );
    assert!(sidecar.contains("depth = \"seed\""), "{sidecar}");
    assert!(sidecar.contains("soft_mask_class = \"mid\""), "{sidecar}");

    // The shared block: q1:5000-25000 vs s1:8000-28000, forward.
    let paf = std::fs::read_to_string(&paf).unwrap();
    let line = paf
        .lines()
        .find(|l| l.starts_with("q1\t"))
        .unwrap_or_else(|| panic!("no q1 chain in:\n{paf}"));
    let f: Vec<&str> = line.split('\t').collect();
    assert_eq!((f[4], f[5]), ("+", "s1"), "{line}");
    let (qs, qe): (i64, i64) = (f[2].parse().unwrap(), f[3].parse().unwrap());
    let (ss, se): (i64, i64) = (f[7].parse().unwrap(), f[8].parse().unwrap());
    assert!(
        (qs - 5000).abs() < 100 && (qe - 25_000).abs() < 100,
        "{line}"
    );
    assert!(
        (ss - 8000).abs() < 100 && (se - 28_000).abs() < 100,
        "{line}"
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn large_engine_rejects_protein() {
    let dir = tmpdir("prot");
    let p = dir.join("p.fa");
    std::fs::write(&p, ">p\nMKVLAAGIVGLLLAAAQ\n").unwrap();
    let result = batch(&[
        p.to_str().unwrap(),
        "-o",
        dir.join("p.png").to_str().unwrap(),
        "--mode",
        "blastp",
        "--engine",
        "large",
    ]);
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("BLASTN"));
    std::fs::remove_dir_all(&dir).ok();
}
