//! Pixelmap goldens from the C dotter (ADR 0010): every case in
//! `tests/golden/c_dotter/cases.tsv` must reproduce C's `.dot` pixelmap,
//! size, window and pixel factor exactly. Zero tolerance (CLAUDE.md).

use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::Command;

use dottir_io::dot_format::DotFile;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

#[test]
fn pixelmaps_match_c_dotter() {
    let root = root();
    let golden = root.join("tests/golden/c_dotter");
    let corpus = root.join("tests/corpora/c_dotter");
    let dir = std::env::temp_dir().join(format!("dottir_cgolden_{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let cases = std::fs::read_to_string(golden.join("cases.tsv")).unwrap();
    let mut failures = Vec::new();
    let mut n = 0;
    for line in cases
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
    {
        let cols: Vec<&str> = line.split('\t').collect();
        let (name, flags, query, subject) = (cols[0], cols[2], cols[3], cols[4]);
        let mut gz = Vec::new();
        flate2::read::GzDecoder::new(
            std::fs::File::open(golden.join(format!("{name}.dot.gz"))).unwrap(),
        )
        .read_to_end(&mut gz)
        .unwrap();
        let want = DotFile::from_reader(&gz[..]).unwrap();

        let dot = dir.join(format!("{name}.dot"));
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_dottir"));
        cmd.arg("batch").arg(corpus.join(query));
        if subject != "-" {
            cmd.arg(corpus.join(subject));
        }
        cmd.args(flags.split_whitespace())
            .arg("--dot")
            .arg(&dot)
            .arg("-o")
            .arg(dir.join(format!("{name}.png")))
            .arg("--no-sidecar");
        let out = cmd.output().unwrap();
        assert!(
            out.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&out.stderr)
        );
        let got = DotFile::read(&dot).unwrap();
        n += 1;

        let header = |d: &DotFile| {
            (
                d.zoom,
                d.width,
                d.height,
                d.pixel_fac,
                d.sliding_window_size,
            )
        };
        if header(&got) != header(&want) {
            failures.push(format!(
                "{name}: (zoom, width, height, pixel factor, window) {:?}, C {:?}",
                header(&got),
                header(&want)
            ));
            continue;
        }
        let diff: Vec<usize> = (0..want.pixels.len())
            .filter(|&i| got.pixels[i] != want.pixels[i])
            .collect();
        if let Some(&i) = diff.first() {
            let w = want.width as usize;
            failures.push(format!(
                "{name}: {} pixels differ; first at x={} y={}: {} (C {})",
                diff.len(),
                i % w,
                i / w,
                got.pixels[i],
                want.pixels[i]
            ));
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
    assert!(n >= 15, "only {n} cases read");
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
