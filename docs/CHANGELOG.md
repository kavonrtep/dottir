# Changelog

All notable changes to dottir are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
the project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-05-19

### Added

* New `dottir-render` crate — pure, GUI-free rendering pipeline with
  a tandem-repeat fidelity fixture and a `long_horizontal_runs`
  metric. The `RenderPolicy` enum lets candidate rendering fixes be
  graded numerically against synthetic data instead of GUI
  screenshots.
* GUI: rectangle-select discrete zoom with history stack
  (middle-drag for selection, Esc/Backspace for Back), single-residue
  keyboard crosshair nudge in absolute residue coords (Shift = ×10,
  Ctrl = ×100), double-click 2× zoom, and an explicit Fit button.
* GUI: vector ridge overlay — anti-aliased line segments over
  coherent diagonal runs in the raster, with a single-click toggle
  in the side panel.
* GUI: integer-pixel pan snap (texel-aligned, no sub-pixel sampling
  phase).
* PNG export now includes margins, axis ticks, labels, and
  multi-record sequence names (matches the SVG / CLI output).
* Minor (unlabeled) ticks at step/5 between every labelled tick —
  in the GUI and in PNG/SVG exports.
* Documentation: `docs/reviews/dotter-sizing-model.md` (read of
  C dotter's sizing model) and
  `docs/reviews/gui-nearest-upscale-jaggedness.md` (resampling
  artifact analysis and severity framing).
* CLAUDE.md note on the hermit sandbox; `hermit/` tree ignored.

### Changed

* **GUI render path: the pixmap is now sized to the *physical*
  canvas (one texture texel = one physical framebuffer pixel)**
  instead of the logical canvas. Fixes the diagonal-staircase
  artifact on HiDPI / fractional display scales that injected
  false jogs into pure-diagonal regions — a correctness problem
  since it generated false positives for exactly the indel
  patterns dotplot users scan for. See Option C in
  `docs/reviews/dotter-sizing-model.md`. Cost: pixmap memory and
  compute time scale as `pixels_per_point²` (~2.25× at 1.5×,
  4× at 2×).
* Sequences are now sliced on rect-zoom — compute runs only on
  the selection plus margin instead of the full input, dropping
  recompute cost on deeply-zoomed views.
* Peak memory cut ~3–10× via auto-zoom on load, peak-preserving
  max-pool texture downsample at fit-zoom, and a bounded-memory
  parallel pixelmap allocation that shares one atomic buffer
  across all rayon workers.
* Compute uses auto `pixel_fac`, a bigger pixelmap, and an
  adaptive sampler — intensity output now matches C dotter.
* Axis labels use one decimal in the 1k–10k range so adjacent
  ticks are distinguishable (`9.0k / 9.1k / 9.2k` instead of
  three `9k`s). Both GUI and exports.
* Default greyramp restored to dotter's 40/100; ridge overlay
  default off.
* **API change in `dottir-io`**:
  `png_export::write_grayscale_png_with_axes` and
  `svg_export::write_svg` now take an explicit
  `invert_pixels: bool` parameter. CLI callers pass `true`
  (raw kernel input → invert); GUI callers pass `false`
  (greyramp-mapped input is already in display space).

### Fixed

* SVG export was double-inverting greyramp-mapped pixels — the
  output was flipped relative to the GUI screen. Now matches.
* Scroll-wheel zoom drift on cursor — zoom now stays centred on
  the cursor instead of pulling the view sideways.
* Initial view shows the full data at fit zoom instead of an
  arbitrary zoomed-in portion.
* Crosshair off-by-one at certain compute zoom levels.
* Crosshair coord label background is now transparent so it
  doesn't occlude dotplot diagonals underneath.

## [0.1.2]

> Legacy bundled entry for the 0.1.x line (no per-version history was
> kept). The list below documents the state at the v0.1.2 tag; split
> per-release later if useful.

### Added

* Cargo workspace with `dottir-core`, `dottir-io`, `dottir-cli`, and
  `dottir-gui`.
* Faithful port of Karlin/Altschul λ/K/H statistics
  (`dottir-core::karlin`), bit-identical to a standalone C harness
  extracted from `dotterKarlin.c`.
* Eight built-in protein matrices (BLOSUM45/50/62/80/90 + PAM30/70/250
  vendored from NCBI; BLOSUM62 inline from C-dotter) plus DNA
  identity. Strict `parse_blast_protein` / `parse_blast_dna`
  parsers — every cell must be covered, extra letters via explicit
  allowlist. `ScoreMatrix::validate()` enforces Karlin
  prerequisites. `custom_dna(match, mismatch)` constructor.
* Sliding-window dotplot kernel (`compute_dotplot`):
  - BLASTN forward / reverse / both-strand;
  - BLASTP forward;
  - self-comparison with `Triangle::{Both, Upper, Lower}` and
    optional `disable_mirror`;
  - separate forward/reverse channels for inverted-repeat
    highlighting.
* **Bounded-memory parallel driver**: one shared atomic `PixelMap`
  for the whole pass regardless of `n_threads`. `memory_limit_bytes`
  is now honest. Output byte-identical across thread counts.
* `dottir-io::Sequence` — record-aware FASTA wrapper preserving
  per-record IDs / offsets / source path. `breaks()` for the
  upcoming breakline rendering; `record_at(coord)` powers
  multi-record coord display.
* `dottir-io::fasta::load_fasta_file` — single-pass FASTA read
  returning records + raw bytes (for hashing) in one disk read.
* `dottir-io::png_export` — 8-bit greyscale PNG with `tEXt`
  provenance.
* `dottir-io::params` — TOML sidecar (SHA-256 + resolved parameters +
  Karlin numbers + hostname/OS/UTC timestamp).
* `dottir-io::alignment` — `slice_pair` ±N residues around a coord
  with FASTA-pair / Stockholm / plain output.
* `dottir batch` CLI — positional `<query> <subject> -o <out.png>`
  plus the usual flags. `--auto-zoom` to fit a target max dimension.
* **`dottir-gui` interactive viewer**:
  - Light theme default with View-menu toggle.
  - Command-line pre-load: `dottir-gui q.fa s.fa [options]`
    mirrors the original Dotter.
  - File menu: Open Query / Open Subject (rfd native pickers);
    Save PNG (applies current greyramp LUT).
  - Central canvas: textured pixelmap with primary-drag pan +
    scroll-wheel zoom-on-cursor + click-to-set crosshair + arrow-key
    nudge (Shift = ×10, Ctrl = ×100).
  - Greyramp panel: white/black sliders + Swap/Reset + live LUT
    strip preview.
  - Settings dialog: mode/matrix/W/zoom/pixel-fac/strand/self-comp/
    triangle/memory cap.
  - Status bar: pixelmap dims + W + crosshair coords + pixel value,
    `record_id:position` style when the input has multiple records.
* GitHub Actions: `fmt`, `clippy`, `test` on Linux + Windows;
  `wasm-check` for `dottir-core --no-default-features`; `deny`;
  release workflow producing static Linux+Windows binaries on tag.
* mdBook scaffold at `docs/book/` (intro, install, CLI, GUI,
  algorithm, matrices, reproducibility, crates, ADRs).
* `cargo-deny` configuration.

### Deferred

Tracked in `IMPROVEMENTS_PLAN.md`:

* GFF3 / PAF loaders + annotation track panel in the GUI.
* SVG / PDF / `.dot` exports.
* BLASTX three-frame translation.
* Recompute-on-zoom-settle for sub-pixel detail at high zoom (the
  current GUI is a viewport transform over one computed pixelmap).
* Sub-dotter spawn from rubber-band selection.
* Alignment-view dock.
* Session save/load.

### Goldens

Karlin λ/K/H/window pinned bit-identical to the C reference across
four fixtures (`tests/golden/karlin/values.tsv`). Pixelmap goldens
are not yet pinned against C-dotter output (requires building the
GTK2 binary in `third_party/seqtools/`).

### Test counts

95 tests across the workspace: 41 dottir-core unit, 1 Karlin
integration, 4 parallel determinism, 5 memory budget, 8 phase 1
kernel, 10 phase 2 modes, 22 dottir-io unit, 4 doc-tests.
