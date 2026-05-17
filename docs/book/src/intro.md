# Introduction

`dottir` is a Rust reimplementation of **Dotter** (Sonnhammer & Durbin
1995), a matrix-scored sliding-window sequence dot-matrix plotter with
interactive sensitivity control. The project has three goals, in
priority order:

1. **Faithful preservation** of the original scientific algorithms —
   Karlin/Altschul window-size estimation, score-matrix
   sliding-window scoring, anti-diagonal suppression, and the live
   greyramp LUT.
2. **Modern cross-platform GUI** (egui/eframe) replacing GTK2.
3. **New features for repeat-oriented genome analysis**: GFF3
   annotation overlays, PAF HSP overlays, region/alignment export,
   inverted-repeat highlighting.

## What works today

* `dottir batch` — headless CLI: computes a dotplot from a FASTA pair
  and writes a greyscale PNG (with `tEXt` provenance) plus a TOML
  params sidecar (SHA-256 of inputs, resolved Karlin parameters,
  host info).
* `dottir-gui` — interactive frontend. Light theme by default; pan
  with primary-drag, zoom on cursor with the scroll wheel, click to
  set the crosshair, arrow keys to nudge (Shift = ×10, Ctrl = ×100).
  Live greyramp panel (no recomputation). Settings dialog for
  mode/matrix/W/zoom/strand/self-comparison. Loads via
  File → Open or directly from the command line:
  `dottir-gui query.fa subject.fa -W 25`.
* BLASTN forward / reverse / both-strand; BLASTP forward;
  self-comparison with mirror; optional separate forward/reverse
  channels (for inverted-repeat highlighting).
* Karlin/Altschul λ, K, H, and the window-size derivation produce
  output **bit-identical** to a standalone C reference harness
  extracted from `dotterKarlin.c`.
* Parallel computation via rayon, **byte-identical** across thread
  counts (spec §4.1.11). Memory budget is honest after Phase A1: one
  shared atomic pixelmap per pass regardless of `n_threads`.
* Multi-record FASTA: record boundaries are preserved through the
  load path. The GUI status bar shows `chr4:1234` style coords
  (record + position) when the input has more than one record.
* Eight built-in protein matrices (BLOSUM45/50/62/80/90, PAM30/70/250)
  plus a DNA identity matrix and a `custom_dna(match, mismatch)`
  helper. BLAST-format text matrices load via the strict
  `parse_blast_protein` / `parse_blast_dna` parsers — every cell of
  the alphabet submatrix must be covered, malformed input fails
  loudly.

## What's pending

See [Architecture decisions](./adr.md) for context. The remaining
spec items are tracked in `IMPROVEMENTS_PLAN.md` at the repo root:

* GFF3 / PAF annotation track overlays.
* SVG / PDF / `.dot` exports.
* BLASTX three-frame translation.
* Recompute-on-zoom-settle for sub-pixel detail at high zoom (the
  current GUI is a pure viewport transform over one computed
  pixelmap).
* Sub-dotter spawn from a rubber-band selection.
* Alignment-view dock alongside the canvas.
* Session save/load.

## How to read this book

* [Installation](./install.md) — how to build and run.
* [The `dottir batch` CLI](./cli.md) — full flag reference.
* [Algorithm overview](./algorithm.md) — what the kernel does and
  why.
* [Score matrices](./matrices.md) — which matrices ship, where
  they came from, how to load your own.
* [Reproducibility](./reproducibility.md) — the params sidecar
  format.
* [Crate layout](./crates.md) — the four-crate workspace.
* [Architecture decisions](./adr.md) — index of ADRs.
