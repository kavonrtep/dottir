# Introduction

`dottir` is a Rust reimplementation of **Dotter** (Sonnhammer & Durbin
1995), a matrix-scored sliding-window sequence dot-matrix plotter with
interactive sensitivity control. The project has three goals, in
priority order:

1. **Faithful preservation** of the original scientific algorithms —
   Karlin/Altschul window-size estimation, score-matrix sliding-window
   scoring, anti-diagonal suppression, and the live greyramp LUT.
2. **Modern cross-platform GUI** (egui/eframe) replacing GTK2. *Not
   yet shipped — see [ADR 0003](./adr.md).*
3. **New features for repeat-oriented genome analysis**: GFF3
   annotation overlays, PAF HSP overlays, region/alignment export,
   inverted-repeat highlighting. *Partially shipped; see the
   architecture decisions for what's pending.*

## What works today

* The `dottir batch` headless CLI computes a dotplot and writes a
  greyscale PNG with full `tEXt` provenance, plus a TOML params
  sidecar with SHA-256 hashes of the inputs and the resolved Karlin
  parameters.
* BLASTN forward, reverse, both-strand; BLASTP forward;
  self-comparison with mirror; inverted-repeat split into separate
  forward/reverse channels.
* Karlin/Altschul λ, K, H, and the window-size derivation produce
  output **bit-identical** to a standalone C reference harness
  extracted from `dotterKarlin.c`.
* Parallel computation via rayon, with output **byte-identical**
  across thread counts.

## What's deferred

See [Architecture decisions](./adr.md) for the full list. The big
ones:

* The egui GUI runtime is deferred until the workspace MSRV bumps
  past Rust 1.75 (the egui ≥ 0.27 transitive-dep chain needs
  edition2024 / Rust 1.85).
* GFF3 / PAF loaders via `noodles-*` have the same blocker.

Until then `dottir batch` is the supported entry point.

## How to read this book

* [Installation](./install.md) — how to build and run.
* [The `dottir batch` CLI](./cli.md) — full flag reference.
* [Algorithm overview](./algorithm.md) — what the kernel does and why.
* [Score matrices](./matrices.md) — which matrices ship, where they
  came from, how to load your own.
* [Reproducibility](./reproducibility.md) — how the params sidecar
  works and why it matters.
* [Crate layout](./crates.md) — the four-crate workspace.
* [Architecture decisions](./adr.md) — index of ADRs.
