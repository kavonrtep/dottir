//! End-to-end `dottir batch --paf` (genome mode, ADR 0008) on generated
//! inputs; no aligner needed.

use std::fmt::Write as _;
use std::path::{Path, PathBuf};
use std::process::Command;

fn tmpdir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dottir_genome_{}_{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn batch(dir: &Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_dottir"))
        .arg("batch")
        .args(args)
        .env("DOTTIR_CACHE_DIR", dir.join("cache"))
        .output()
        .expect("running dottir batch")
}

fn fasta(records: &[(&str, usize)]) -> String {
    let mut s = String::new();
    for (name, len) in records {
        writeln!(s, ">{name}").unwrap();
        for chunk in "ACGT".repeat(len / 4 + 1).as_bytes()[..*len].chunks(60) {
            writeln!(s, "{}", std::str::from_utf8(chunk).unwrap()).unwrap();
        }
    }
    s
}

/// Reference chr1/chr2; query contigs: c1 = chr2 start (forward), c2 =
/// chr1 end (reverse), c3 = chr1 start (forward), c4 unaligned.
const PAF: &str = "\
c1\t30000\t0\t30000\t+\tchr2\t80000\t0\t30000\t27000\t30000\t60\ttp:A:P
c2\t20000\t0\t20000\t-\tchr1\t100000\t80000\t100000\t18000\t20000\t60\ttp:A:P
c3\t25000\t0\t25000\t+\tchr1\t100000\t0\t25000\t20000\t25000\t60\ttp:A:P
c3\t25000\t1000\t3000\t+\tchr2\t80000\t50000\t52000\t1500\t2000\t0\ttp:A:S
";

fn inputs(dir: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let (r, q, p) = (dir.join("ref.fa"), dir.join("qry.fa"), dir.join("aln.paf"));
    std::fs::write(&r, fasta(&[("chr1", 100_000), ("chr2", 80_000)])).unwrap();
    std::fs::write(
        &q,
        fasta(&[("c1", 30_000), ("c2", 20_000), ("c3", 25_000), ("c4", 5000)]),
    )
    .unwrap();
    std::fs::write(&p, PAF).unwrap();
    (r, q, p)
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

#[test]
fn paf_plot_sorting_exports_and_sidecar() {
    let dir = tmpdir("sort");
    let (r, q, paf) = inputs(&dir);
    let (out, assoc, sorted) = (dir.join("g.png"), dir.join("a.tsv"), dir.join("s.fa"));
    let res = batch(
        &dir,
        &[
            s(&r),
            s(&q),
            "--paf",
            s(&paf),
            "--sort-query",
            "--assoc-out",
            s(&assoc),
            "--sorted-query-out",
            s(&sorted),
            "--auto-zoom",
            "500",
            "-o",
            s(&out),
        ],
    );
    assert!(
        res.status.success(),
        "{}",
        String::from_utf8_lossy(&res.stderr)
    );
    assert!(out.exists());
    let table = std::fs::read_to_string(&assoc).unwrap();
    let rows: Vec<Vec<&str>> = table
        .lines()
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect();
    let order: Vec<(&str, &str, &str)> = rows.iter().map(|r| (r[0], r[1], r[2])).collect();
    assert_eq!(
        order,
        vec![
            ("c3", "chr1", "+"),
            ("c2", "chr1", "-"),
            ("c1", "chr2", "+"),
            ("c4", "None", "+"),
        ]
    );
    let fa = std::fs::read_to_string(&sorted).unwrap();
    let headers: Vec<&str> = fa.lines().filter(|l| l.starts_with('>')).collect();
    assert_eq!(headers, vec![">c3", ">c2 dottir:revcomp", ">c1", ">c4"]);
    let sidecar = std::fs::read_to_string(dir.join("g.png.params.toml")).unwrap();
    assert!(sidecar.contains("[genome]"), "{sidecar}");
    assert!(sidecar.contains("sorted_query = true"));
    assert!(sidecar.contains("alignments = 4"));
}

#[test]
fn output_ignores_paf_order_and_orientation() {
    let dir = tmpdir("order");
    let (r, q, paf) = inputs(&dir);
    let run = |paf: &Path, first: &Path, second: &Path, name: &str| -> Vec<u8> {
        let out = dir.join(name);
        let res = batch(
            &dir,
            &[
                s(first),
                s(second),
                "--paf",
                s(paf),
                "--no-sidecar",
                "-o",
                s(&out),
            ],
        );
        assert!(
            res.status.success(),
            "{}",
            String::from_utf8_lossy(&res.stderr)
        );
        // Pixels only: the PNG text chunks carry the PAF's SHA-256.
        let mut reader = png::Decoder::new(std::fs::File::open(out).unwrap())
            .read_info()
            .unwrap();
        let mut buf = vec![0; reader.output_buffer_size()];
        let info = reader.next_frame(&mut buf).unwrap();
        buf.truncate(info.buffer_size());
        buf
    };
    let base = run(&paf, &r, &q, "a.png");
    let shuffled = dir.join("shuffled.paf");
    let mut lines: Vec<&str> = PAF.lines().collect();
    lines.reverse();
    std::fs::write(&shuffled, lines.join("\n") + "\n").unwrap();
    assert_eq!(base, run(&shuffled, &r, &q, "b.png"));
    // The same PAF with the inputs given the other way round is swapped
    // on reading: a different plot, but no error.
    let swapped = run(&paf, &q, &r, "c.png");
    assert_ne!(base, swapped);
}

#[test]
fn paf_from_other_sequences_is_rejected() {
    let dir = tmpdir("bad");
    let (r, q, _) = inputs(&dir);
    let paf = dir.join("bad.paf");
    std::fs::write(&paf, PAF.replace("chr2\t80000", "chr2\t80001")).unwrap();
    let res = batch(
        &dir,
        &[s(&r), s(&q), "--paf", s(&paf), "-o", s(&dir.join("x.png"))],
    );
    assert!(!res.status.success());
    assert!(String::from_utf8_lossy(&res.stderr).contains("different sequences"));
    let res = batch(
        &dir,
        &[
            s(&r),
            s(&q),
            "--paf",
            s(&paf),
            "--engine",
            "large",
            "-o",
            s(&dir.join("x.png")),
        ],
    );
    assert!(String::from_utf8_lossy(&res.stderr).contains("--paf is only used by --engine genome"));
}
