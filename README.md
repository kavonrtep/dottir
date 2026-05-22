# dottir

A modern Rust reimplementation of **[Dotter](https://www.sanger.ac.uk/tool/dotter/)** — the
classic sliding-window dot-matrix plotter by Sonnhammer & Durbin
([Gene 167(2), 1995](https://doi.org/10.1016/0378-1119(95)00714-8)) shipped
in [seqtools-4.28](https://github.com/douglasgscofield/seqtools-4.28).
Same scientific feel — Karlin/Altschul window sizing, score-matrix
scoring, anti-diagonal suppression, live greyramp — in a single
cross-platform binary with no GTK dependency.

## Features

- **CLI and interactive GUI** in one workspace. Same compute engine
  behind both.
- **BLASTN** (forward, reverse, both strands), **BLASTP**, and
  **BLASTX** (three reading frames). `--watson-only` /
  `--crick-only` and `-r` / `-v` reverse-complement options match
  the original Dotter.
- **Self-comparison** with mirror, upper, or lower triangle.
- **Inverted-repeat highlighting** — paint reverse-strand hits in
  magenta against forward-strand greyscale, so palindromes pop.
- **Multi-record FASTA** with **breakline rendering** and
  `record:position` coordinates in the GUI status bar.
- **Live greyramp** — drag the black/white sliders, the image
  updates without re-running compute.
- **Zoom-aware** GUI: smooth pan/zoom, then a debounced recompute
  at a finer tier exposes per-residue detail.
- **Honest memory cap** (`--memory-mib`) — refuses oversized
  pixelmaps rather than silently OOMing.
- **Provenance**: every export writes a `.params.toml` sidecar with
  the SHA-256 of each input and the resolved Karlin parameters.
- **Outputs**: PNG (with `tEXt` provenance), SVG (embedded base64
  PNG + axis labels), and the original C-dotter `.dot` binary
  (read + write).
- **Deterministic**: byte-identical output across runs and across
  rayon thread counts.

## Install

### Conda (recommended)

Pre-built Linux x86_64 packages live on the
[`petrnovak`](https://anaconda.org/petrnovak/dottir) channel:

```sh
mamba install -c petrnovak -c conda-forge dottir
# or:
conda install -c petrnovak -c conda-forge dottir
```

This installs both the `dottir` CLI and the `dottir-gui` binary into
the active environment. The conda package bundles the X11/GL
runtime libraries the GUI needs, so no extra system packages are
required.

### From source

dottir builds with stable Rust. MSRV is **1.85** (the dev
environment uses Rust 1.95 from conda-forge).

```sh
git clone https://github.com/petr/dottir
cd dottir
cargo build --release
# binaries: target/release/dottir (CLI) and target/release/dottir-gui
```

On Linux the GUI needs the usual X11/GL stack:

```sh
sudo apt install \
    libxkbcommon-dev libgl1-mesa-dev libxcb1-dev \
    libxcb-render0-dev libxcb-shape0-dev libxcb-xfixes0-dev \
    libxcursor-dev libxi-dev libxrandr-dev libgtk-3-dev
```

The CLI has no system-library dependencies.

## Quick start

A self-comparison of a single sequence — the classic repeat-finding
plot:

```sh
dottir batch chr4.fa -o chr4.png --auto-zoom 4000
```

(Pass a single FASTA for a self-comparison; pass two for a
pairwise plot.)

A query vs. subject BLASTP at a fixed window:

```sh
dottir batch query.faa target.faa -o p.png \
    --mode blastp --matrix BLOSUM45 -W 12
```

Open the GUI with one sequence pre-loaded (self-comparison) or
both (pairwise):

```sh
dottir-gui chr4.fa -W 25
dottir-gui query.fa subject.fa
```

Full flag reference: `dottir batch --help`, `dottir-gui --help`, or
the [CLI page](./docs/book/src/cli.md) /
[GUI page](./docs/book/src/gui.md) in the user manual.

## Documentation

- [User manual](./docs/book/src/intro.md) (mdBook): installation,
  CLI/GUI reference, algorithm overview, score matrices,
  reproducibility format.
- [Spec](./docs/dottir_specification.md) and
  [implementation plan](./docs/IMPLEMENTATION_PLAN.md): what's
  guaranteed, what's planned.
- [Improvements plan](./docs/IMPROVEMENTS_PLAN.md): progress against
  [REVIEW.md](./docs/REVIEW.md) findings.
- [Changelog](./docs/CHANGELOG.md): per-release notes.
- [Architecture decisions](./docs/adr/): MADR records.

## License

GPL-3.0-or-later (same as upstream Dotter). See [LICENSE](./LICENSE)
and [NOTICE](./NOTICE) for full credits to Sonnhammer, Durbin,
Barson, Scofield, and the NCBI BLAST authors.
