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

/// A minimap2 stand-in: `--version`; `-d OUT` writes an index; otherwise
/// it maps every stdin record to `chr1` at `[0, len)` and logs each run
/// to `runs.log`. `FAKE_FAIL` makes mapping fail.
#[cfg(unix)]
fn fake_aligner(dir: &Path) -> PathBuf {
    use std::os::unix::fs::PermissionsExt;
    let p = dir.join("fake-minimap2");
    let script = format!(
        r#"#!/bin/sh
if [ "$1" = "--version" ]; then echo fake-1.0; exit 0; fi
for a in "$@"; do
  if [ "$prev" = "-d" ]; then : > "$a"; exit 0; fi
  prev="$a"
done
if [ -n "$FAKE_FAIL" ]; then echo "[ERROR] fake failure" >&2; exit 3; fi
echo run >> "{log}"
awk '/^>/ {{ if (n) print n "\t" l "\t0\t" l "\t+\tchr1\t100000\t0\t" l "\t" l "\t" l "\t60\ttp:A:P"; n = substr($1, 2); l = 0; next }}
     {{ l += length($0) }}
     END {{ if (n) print n "\t" l "\t0\t" l "\t+\tchr1\t100000\t0\t" l "\t" l "\t" l "\t60\ttp:A:P"; print "[M::worker_pipeline] mapped " NR / 2 " sequences" > "/dev/stderr" }}'
"#,
        log = dir.join("runs.log").display()
    );
    std::fs::write(&p, script).unwrap();
    std::fs::set_permissions(&p, std::fs::Permissions::from_mode(0o755)).unwrap();
    p
}

#[cfg(unix)]
#[test]
fn aligner_runs_are_chunked_remapped_and_cached() {
    let dir = tmpdir("aligner");
    let (r, q, _) = inputs(&dir);
    let fake = fake_aligner(&dir);
    let run = |extra: &[&str]| {
        let out = dir.join("a.png");
        let mut args = vec![
            s(&r),
            s(&q),
            "--engine",
            "genome",
            "--aligner",
            s(&fake),
            "--chunk-size",
            "10000",
            "--paf-out",
        ];
        let paf_out = dir.join("out.paf");
        args.push(s(&paf_out));
        args.extend(["-o", s(&out)]);
        args.extend(extra);
        let res = batch(&dir, &args);
        assert!(
            res.status.success(),
            "{}",
            String::from_utf8_lossy(&res.stderr)
        );
        std::fs::read_to_string(&paf_out).unwrap()
    };
    let paf = run(&[]);
    // c1 (30 kb) is sent as three 10 kb pieces; each maps to chr1 at
    // [0, 10 kb), so the pieces are not co-linear and stay apart.
    let c1: Vec<(u64, u64)> = paf
        .lines()
        .filter(|l| l.starts_with("c1\t"))
        .map(|l| {
            let c: Vec<&str> = l.split('\t').collect();
            (c[2].parse().unwrap(), c[3].parse().unwrap())
        })
        .collect();
    assert_eq!(c1, vec![(0, 10_000), (10_000, 20_000), (20_000, 30_000)]);
    // Record length, not piece length, in column 2.
    assert!(paf
        .lines()
        .filter(|l| l.starts_with("c1\t"))
        .all(|l| l.split('\t').nth(1) == Some("30000")));
    let sidecar = std::fs::read_to_string(dir.join("a.png.params.toml")).unwrap();
    assert!(sidecar.contains("source = \"aligner\""), "{sidecar}");
    assert!(sidecar.contains("aligner_version = \"fake-1.0\""));
    assert!(sidecar.contains("chunk_size = 10000"));

    // The second run reuses the cache; --no-cache maps again.
    run(&[]);
    let runs = || {
        std::fs::read_to_string(dir.join("runs.log"))
            .unwrap()
            .lines()
            .count()
    };
    assert_eq!(runs(), 1);
    run(&["--no-cache"]);
    assert_eq!(runs(), 2);
}

#[cfg(unix)]
#[test]
fn aligner_failure_is_reported_with_its_stderr() {
    let dir = tmpdir("fail");
    let (r, q, _) = inputs(&dir);
    let fake = fake_aligner(&dir);
    let res = Command::new(env!("CARGO_BIN_EXE_dottir"))
        .args([
            "batch",
            s(&r),
            s(&q),
            "--engine",
            "genome",
            "--aligner",
            s(&fake),
        ])
        .args(["-o", s(&dir.join("x.png"))])
        .env("DOTTIR_CACHE_DIR", dir.join("cache"))
        .env("FAKE_FAIL", "1")
        .output()
        .unwrap();
    assert!(!res.status.success());
    let err = String::from_utf8_lossy(&res.stderr);
    assert!(err.contains("fake failure"), "{err}");
}
