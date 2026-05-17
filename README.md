# dottir

`dottir` is a Rust reimplementation of **Dotter**, the classic dot-matrix plotter by Sonnhammer & Durbin (1995). It keeps the original idea of score-matrix sliding-window plots, Karlin/Altschul window sizing, and anti-diagonal suppression, while modernizing the codebase and user experience.

## What It Does

- Compute dotplots from FASTA sequences
- Support BLASTN and BLASTP workflows
- Export greyscale PNGs with provenance metadata
- Provide a desktop GUI for interactive inspection
- Keep outputs deterministic across runs and thread counts

## Build

```sh
cargo build --workspace
```

The workspace is pinned to Rust `1.95.0`.

## Run

CLI batch mode:

```sh
cargo run -p dottir-cli -- batch query.fa subject.fa -o plot.png
```

GUI:

```sh
cargo run -p dottir-gui
```

## Project Layout

- `crates/dottir-core/` - algorithmic core
- `crates/dottir-io/` - FASTA, PNG, params sidecar, alignment helpers
- `crates/dottir-cli/` - headless batch binary
- `crates/dottir-gui/` - egui-based interactive frontend
- `docs/book/` - user manual
- `docs/adr/` - architecture decisions

## Inspiration

The project is directly inspired by the original Dotter implementation and its behavior. The goal is to preserve the scientific feel of the original tool while making it easier to maintain, test, and extend in Rust.
