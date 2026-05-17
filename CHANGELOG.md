# Changelog

All notable changes to dottir are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
the project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
