# Installation

dottir is a Rust workspace. The CLI has no system-library dependencies;
the GUI needs the usual X11/GL stack on Linux (no GTK).

## MSRV

Rust **1.85** or newer (the edition2024 stable floor). The dev
environment in this repo uses Rust 1.95 from conda-forge.

## From source — CLI only

```sh
git clone https://github.com/petr/dottir
cd dottir
cargo build --release -p dottir-cli
# binary lands at target/release/dottir
```

## From source — GUI

The GUI is built on `egui` / `eframe`. On Debian / Ubuntu install the
windowing + GL deps first:

```sh
sudo apt install \
    libxkbcommon-dev libgl1-mesa-dev libxcb1-dev \
    libxcb-render0-dev libxcb-shape0-dev libxcb-xfixes0-dev \
    libxcursor-dev libxi-dev libxrandr-dev libgtk-3-dev
```

`libgtk-3-dev` is only needed for the native file picker (`rfd`
crate). Then:

```sh
cargo build --release -p dottir-gui
target/release/dottir-gui            # opens an empty window
target/release/dottir-gui q.fa s.fa  # pre-loads both sequences
```

## With conda-forge Rust

A clean way to pin the toolchain without `rustup`:

```sh
conda create -n dottir -c conda-forge rust
conda activate dottir
cargo build --release
```

## From a release binary

Pending Phase 9 of `docs/IMPLEMENTATION_PLAN.md`. The GitHub Actions
release workflow at `.github/workflows/release.yml` produces static
Linux + Windows artifacts on tag push.
