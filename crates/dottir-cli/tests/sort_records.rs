//! `dottir batch --sort-records` (D-GENIES order from large-engine
//! chains) on a generated reference and shuffled, partly reversed
//! contigs.

use std::path::{Path, PathBuf};
use std::process::Command;

fn random_dna(n: usize, mut state: u64) -> Vec<u8> {
    (0..n)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 62) as usize]
        })
        .collect()
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

fn fasta(path: &Path, records: &[(&str, &[u8])]) {
    let mut s = String::new();
    for (name, seq) in records {
        s.push_str(&format!(">{name}\n"));
        for line in seq.chunks(60) {
            s.push_str(std::str::from_utf8(line).unwrap());
            s.push('\n');
        }
    }
    std::fs::write(path, s).unwrap();
}

fn run(dir: &Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_dottir"))
        .arg("batch")
        .args(args)
        .current_dir(dir)
        .output()
        .unwrap()
}

#[test]
fn contigs_are_sorted_and_oriented_along_the_reference() {
    let dir: PathBuf = std::env::temp_dir().join(format!("dottir_sortrec_{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    let chr = random_dna(60_000, 7);
    fasta(&dir.join("ref.fa"), &[("chr", &chr)]);
    // File order c, a, b, n; a and c reversed; n unrelated.
    let a = revcomp(&chr[0..20_000]);
    let b = chr[20_000..40_000].to_vec();
    let c = revcomp(&chr[40_000..60_000]);
    let n = random_dna(5000, 99);
    fasta(
        &dir.join("ctg.fa"),
        &[("c", &c), ("a", &a), ("b", &b), ("n", &n)],
    );
    let out = run(
        &dir,
        &[
            "ref.fa",
            "ctg.fa",
            "--engine",
            "large",
            "--sort-records",
            "vertical",
            "-o",
            "s.png",
        ],
    );
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let sidecar = std::fs::read_to_string(dir.join("s.png.params.toml")).unwrap();
    assert!(sidecar.contains("[record_sort]"), "{sidecar}");
    assert!(
        sidecar.contains(r#"order = ["a-", "b", "c-", "n"]"#)
            || sidecar.contains("order = [\n    \"a-\",\n    \"b\",\n    \"c-\",\n    \"n\",\n]"),
        "{sidecar}"
    );

    // The other axis, and the refusals.
    let out = run(
        &dir,
        &[
            "ctg.fa",
            "ref.fa",
            "--engine",
            "large",
            "--sort-records",
            "horizontal",
            "-o",
            "h.png",
        ],
    );
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let sidecar = std::fs::read_to_string(dir.join("h.png.params.toml")).unwrap();
    assert!(sidecar.contains("axis = \"horizontal\""));
    assert!(sidecar.contains("\"a-\""));
    let out = run(
        &dir,
        &["ref.fa", "--sort-records", "vertical", "-o", "x.png"],
    );
    assert!(String::from_utf8_lossy(&out.stderr).contains("needs two inputs"));
}
