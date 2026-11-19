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
        "--seed-sample",
        "2",
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
    assert!(sidecar.contains("seed_sample = 2"), "{sidecar}");
    let compact: String = sidecar.split_whitespace().collect();
    assert!(compact.contains("depth_sample=[1,1,]"), "{sidecar}");

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

/// Pixels and tEXt chunks of a PNG.
fn read_png(path: &std::path::Path) -> (Vec<u8>, Vec<(String, String)>) {
    let mut reader = png::Decoder::new(std::fs::File::open(path).unwrap())
        .read_info()
        .unwrap();
    let mut buf = vec![0; reader.output_buffer_size()];
    let info = reader.next_frame(&mut buf).unwrap();
    buf.truncate(info.buffer_size());
    let text = reader
        .info()
        .uncompressed_latin1_text
        .iter()
        .map(|c| (c.keyword.clone(), c.text.clone()))
        .collect();
    (buf, text)
}

#[test]
fn coverage_shading_follows_the_ramp_and_keeps_dot_as_score() {
    let dir = tmpdir("cov");
    // A 40 × 500 bp tandem array between unique flanks.
    let unit = dna(500, 11);
    let fa = dir.join("tandem.fa");
    std::fs::write(
        &fa,
        format!(
            ">t\n{}{}{}\n",
            dna(20_000, 12),
            unit.repeat(40),
            dna(20_000, 13)
        ),
    )
    .unwrap();
    let run = |name: &str, extra: &[&str]| {
        let png = dir.join(format!("{name}.png"));
        let dot = dir.join(format!("{name}.dot"));
        let mut args = vec![
            fa.to_str().unwrap().to_string(),
            "-o".into(),
            png.to_str().unwrap().into(),
            "--dot".into(),
            dot.to_str().unwrap().into(),
            "--engine".into(),
            "large".into(),
            "--zoom".into(),
            "1000".into(),
            "--width".into(),
            "0".into(),
            "--depth".into(),
            "none".into(),
        ];
        args.extend(extra.iter().map(|s| s.to_string()));
        let refs: Vec<&str> = args.iter().map(String::as_str).collect();
        let result = batch(&refs);
        assert!(
            result.status.success(),
            "stderr: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let sidecar = std::fs::read_to_string(dir.join(format!("{name}.png.params.toml"))).unwrap();
        let dot = dottir_io::dot_format::DotFile::read(&dot).unwrap();
        let (pixels, text) = read_png(&png);
        (pixels, text, sidecar, dot.pixels)
    };

    let (score_px, score_text, score_sidecar, score_dot) = run("score", &[]);
    assert!(
        score_sidecar.contains("shade = \"score\""),
        "{score_sidecar}"
    );
    assert!(
        !score_sidecar.contains("coverage_ramp_bp"),
        "{score_sidecar}"
    );
    assert!(score_text.iter().all(|(k, _)| k != "shading"));

    let (cov_px, cov_text, cov_sidecar, cov_dot) = run(
        "cov",
        &["--shade", "coverage", "--coverage-ramp", "100,100000"],
    );
    assert!(
        cov_sidecar.contains("shade = \"coverage\""),
        "{cov_sidecar}"
    );
    // The TOML writer may break the array over several lines.
    let flat: String = cov_sidecar.split_whitespace().collect();
    assert!(
        flat.contains("coverage_ramp_bp=[100.0,100000.0,]")
            || flat.contains("coverage_ramp_bp=[100.0,100000.0]"),
        "{cov_sidecar}"
    );
    let chunk = |k: &str| {
        cov_text
            .iter()
            .find(|(key, _)| key == k)
            .map(|(_, v)| v.as_str())
    };
    assert_eq!(chunk("shading"), Some("coverage"));
    assert_eq!(chunk("coverage-white-bp"), Some("100"));
    assert_eq!(chunk("coverage-black-bp"), Some("100000"));
    assert_eq!(cov_dot, score_dot, ".dot keeps the score pixelmap");
    assert_eq!(cov_px.len(), score_px.len());
    assert_ne!(cov_px, score_px);
    // Inside the array a 1 kb pixel crosses about two diagonals of the
    // 500 bp period, so it is darker than the unique main diagonal.
    // Plot at a 50 px axis margin on every side.
    let (w, m) = ((20_000 * 2 + 500 * 40) / 1000, 50);
    let px = |x: usize, y: usize| cov_px[(y + m) * (w + 2 * m) + x + m];
    let (diag, array) = (px(5, 5), px(30, 35));
    assert!(array < diag, "array {array} vs diagonal {diag}");
    assert_eq!(px(30, 5), 255, "unique off-diagonal stays white");

    // `auto` resolves and records the fitted ramp.
    let (_, _, auto_sidecar, _) = run("auto", &["--shade", "coverage"]);
    assert!(
        auto_sidecar.contains("coverage_ramp_bp ="),
        "{auto_sidecar}"
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn coverage_shading_rejects_bad_ramps_and_the_exact_engine() {
    let dir = tmpdir("covbad");
    let fa = dir.join("a.fa");
    std::fs::write(&fa, format!(">a\n{}\n", dna(5000, 21))).unwrap();
    let png = dir.join("a.png");
    let base = [fa.to_str().unwrap(), "-o", png.to_str().unwrap()];
    for ramp in ["5,5", "0,10", "a,b", "100"] {
        let mut args = base.to_vec();
        args.extend([
            "--engine",
            "large",
            "--shade",
            "coverage",
            "--coverage-ramp",
            ramp,
        ]);
        let r = batch(&args);
        assert!(!r.status.success(), "ramp {ramp} accepted");
    }
    let mut args = base.to_vec();
    args.extend(["--engine", "exact", "--shade", "coverage"]);
    let r = batch(&args);
    assert!(!r.status.success());
    assert!(
        String::from_utf8_lossy(&r.stderr).contains("needs the large engine"),
        "{}",
        String::from_utf8_lossy(&r.stderr)
    );
    // `auto` choosing the exact engine falls back to score shading.
    let mut args = base.to_vec();
    args.extend(["--shade", "coverage"]);
    let r = batch(&args);
    assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
    std::fs::remove_dir_all(&dir).ok();
}
