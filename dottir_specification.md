# Dottir — Specification and Design Direction

Working name `dottir` (tentative; see §10).
Status: draft v0.1, to be refined with developer input.
Reference implementation: `seqtools-4.28/dotterApp/` (Sonnhammer & Durbin 1995; Barson; Scofield).

## 1. Purpose and scope

A modern Rust reimplementation of Dotter — a matrix-scored, sliding-window sequence
dot-matrix plotter with interactive sensitivity control. The project has three intertwined
goals:

1. **Faithful preservation** of the algorithms and interactive behaviour that make
   the original scientifically distinctive — Karlin/Altschul window-size estimation,
   score-matrix-based sliding-window scoring, anti-diagonal suppression at the pixel
   level, and the live greyramp LUT for sensitivity adjustment.
2. **Modern, cross-platform GUI** replacing the GTK2 frontend, with smoother
   interaction, better performance, and a maintainable codebase.
3. **New features for repeat-oriented genome analysis**, principally GFF3 annotation
   overlays, region/alignment export, PAF-based HSP overlays, and self-comparison
   tooling.

Conventions: requirements use **MUST**, **SHOULD**, **MAY** in the RFC 2119 sense.

## 2. Naming and license

- The working name `dottir` is a placeholder; alternatives include `rotter`, `dot-rs`,
  `seqdot`. Decide before first public release. Crate names must not collide on
  crates.io.
- Original Dotter is GPLv3. The Rust implementation can be a clean-room rewrite under
  any license, **except** for any code translated more-or-less verbatim from
  `dotterKarlin.c` (Karlin/Altschul statistics) and other near-direct ports. Two
  acceptable paths:
  - Adopt GPLv3 for the whole project. Simplest, removes ambiguity, compatible with
    upstream.
  - Adopt MIT/Apache-2.0 dual licensing, but reimplement Karlin/Altschul from the
    original Karlin & Altschul 1990 paper and the public-domain NCBI BLAST sources
    (which have permissive terms). Requires more care; document the lineage in
    `karlin.rs`.
  - Default in this spec: **GPLv3**, with a `LICENSE` and `NOTICE` crediting Sonnhammer,
    Durbin, Barson, Scofield, and NCBI. Revisit if a permissive license becomes important.

## 3. Reference implementation

All file:line references below are against `seqtools-4.28`. Useful entry points:

- `dotterApp/dotplot.c:1476` — `calculateImage()`, the top of the main loop.
- `dotterApp/dotplot.c:1308` — `doCalculateImage()`, the inner sliding-window loop.
- `dotterApp/dotplot.c:1405` — the anti-diagonal suppression rule.
- `dotterApp/dotterKarlin.c:144` — `karlin()`, the λ/K/H computation.
- `dotterApp/dotterKarlin.c:343` — `winsizeFromlambdak()`, window-size derivation.
- `dotterApp/greyramptool.c:265` — `updateGreyMap()`, the LUT generation.
- `dotterApp/dotplot.c:1610` — `loadPlot()` and the `.dot` binary format.
- `dotterApp/dotterMain.c` — CLI option parsing; canonical source for option names.
- `dotter.md` (repo root) — full original CLI documentation.

The original test data in `seqtools-4.28/test/data/` and `examples/` should be reused
as the regression corpus.

## 4. Requirements

### 4.1 Algorithmic — MUST preserve

The following are load-bearing for scientific correctness and the "Dotter feel"; any
deviation must be deliberate and documented.

1. **Karlin/Altschul window-size estimation**. From residue composition of both
   sequences and the chosen score matrix, compute λ, K, H, then the expected MSP
   length for a nominal 100×100 matrix; round to integer. Clamp to `[3, 50]` by default
   (configurable). The user override (`-W <int>`) bypasses estimation.
2. **Score-matrix-based scoring**. Default BLOSUM62 for protein, identity matrix for
   DNA (BLASTN-style). Custom matrices via BLAST format file. The sliding window sums
   *matrix scores*, not match/mismatch counts. This is what gives Dotter sensitivity
   on divergent sequences.
3. **Precomputed score row vector** `scoreVec[residue_type][query_position]`. Required
   for the O(qlen·slen) sum recurrence below. Layout should be flat
   (`row_idx * qlen + col_idx`) for cache locality; this is a permitted internal
   deviation from the C implementation's row-of-pointers.
4. **Sliding-window sum recurrence**. For each subject position `s`, slide over `q`:
   `newsum[q] = oldsum[q-1] + scoreVec[sIndex[s]][q] - scoreVec[sIndex[s-W]][q-W]`
   with a `delrow=zeros` fallback for the warm-up region. Two ping-pong buffers.
   O((qlen·slen)/zoom²) memory output, O(qlen·slen) compute.
5. **Pixel max-merge**. Each pixel covers a `zoom × zoom` block of the score matrix.
   Multiple diagonals fall into the same pixel; keep only the maximum. Score is then
   scaled to a `u8` via `min(255, score * pixelFac / W)`.
6. **Anti-diagonal suppression**. Within each pixel, only keep dots whose local
   sub-pixel position is consistent with the diagonal direction:
   `sPosLocal >= qPosLocal` for forward strand,
   `(zoom - 1 - sPosLocal) >= qPosLocal` for reverse strand.
   This is the rule at `dotplot.c:1405`. Without it the plot acquires diagonal-direction
   noise that obscures true alignments. **Add an explicit unit test for this.**
7. **Strand and frame handling**. BLASTN runs forward and reverse passes and
   max-merges into one pixelmap. BLASTX runs three reading frames similarly. BLASTP
   is one pass.
8. **Self-comparison**. Compute one triangle, mirror across the diagonal. Honour
   `--triangle=u|l` and `--disable-mirror`.
9. **Watson-only / Crick-only**. Skip the corresponding strand pass for DNA.
10. **Reverse-complement options** (`-r`, `-v`). Reverse and complement the horizontal
    or vertical sequence before computation.
11. **Determinism**. Same inputs and parameters MUST produce byte-identical pixelmaps
    across runs and across thread counts.

### 4.2 Visualization — preserve with minor modernization

1. **Greyramp tool**. A 256-byte LUT generated from `(blackPoint, whitePoint)`,
   applied to the pixelmap to produce the displayed image. Swapping black and white
   inverts colors. Reset (undo) returns to defaults. Must update at interactive
   rates (no recomputation of the pixelmap).
2. **Greyramp default points**. White=40, black=100 (matching original defaults
   exposed via `--greyramp-white` / `--greyramp-black`).
3. **HSP overlay modes**: off, solid line, score-colored line, greyscale (replaces
   the dotplot image). Matches original `DOTTER_HSPS_*` enum semantics.
4. **Crosshair**. Click or arrow-key to position. Shows synchronized coordinate in
   both sequences. SHOULD support keyboard nudge by 1 / 10 / 100 (Shift/Ctrl modifiers).
5. **Scale and tick marks**. Horizontal and vertical axes with major/minor ticks and
   coordinate labels. Honour `--suppress-scale`, `--labels-off`, `--labels-size`.
6. **Breaklines** between multiple sequences within a single FASTA (multi-record
   input is concatenated with a coloured break line at boundaries). Honour
   `--breakline-colour`.
7. **Session background colour** (`--session-colour`).
8. **Reversed scale axes** with optional coordinate negation (the original
   `hozScaleRev`, `vertScaleRev`, `negateCoords`).
9. **Tweaks permitted beyond the original**:
   - Alternative colormaps (viridis, magma, inferno) selectable via menu; grey
     remains default. The LUT abstraction stays a `[u8; 256] → Rgba8`.
   - Anti-aliased scale lines and labels (the original is pixel-aligned, which is
     fine but dated).
   - High-DPI awareness.
   - Optional dark theme for the surrounding UI chrome.

### 4.3 GUI — modern reimplementation

Target frontend: **egui / eframe**. Reasons: single binary, no system GTK dependency,
straightforward Linux + Windows + macOS builds, WebAssembly target available for free,
fast immediate-mode rendering well suited to a textured pixelmap plus overlays.

Requirements:

1. **Main window** containing: menu bar, dotplot canvas, greyramp panel (collapsible),
   annotation track panel (collapsible), status bar with coordinate readout.
2. **Dotplot canvas** MUST support:
   - Pan with click-and-drag (middle button or modifier).
   - Zoom with scroll wheel, centered on cursor. Zoom is *display zoom* (pan/zoom over
     the existing pixelmap) — separate from the *computation zoom* (`zoomFactor`,
     which affects the pixelmap resolution and requires recomputation).
   - Crosshair at the active coordinate; updated on click and arrow keys.
   - Selection of a sub-rectangle by drag (with modifier) for export or sub-dotter spawn.
3. **Greyramp panel**: black-point and white-point sliders, swap, reset, numeric
   spinboxes. Live update via callback registration (no recomputation).
4. **Annotation tracks panel**: list of loaded GFF3 tracks with per-track toggle,
   colour, height. See §4.4.
5. **Alignment tool**: a separate dock or floating panel showing a window of ±N
   residues around the crosshair from both sequences, aligned with mismatch
   highlighting. N configurable, default 60. Synchronises with the crosshair.
6. **Settings dialog**: matrix selection, window size override, pixel factor,
   memory limit, colormap, font size.
7. **File operations**: open FASTA (drag-and-drop accepted), load GFF3, load HSPs
   from PAF/GFF, save `.dot`, save session (TOML).
8. **Spawn sub-dotter**: rectangular selection on the main plot → new dotter window
   for that sub-range with a finer computation zoom. Equivalent to the original's
   middle-drag-zoom behaviour.
9. **Cross-platform builds**:
   - **Linux x86_64**: primary target; `cargo build --release` works out of the box.
   - **Windows x86_64**: SHOULD work via standard cross-compilation
     (`cross` / GitHub Actions). Test releases per version.
   - **macOS** (Intel and Apple Silicon): MAY be supported (egui makes it free); not a
     release blocker.
   - **WebAssembly**: MAY be supported as a `--target wasm32-unknown-unknown` build for
     a browser viewer. Sequence I/O via the browser's file picker. Useful for sharing
     interactive dotplots with collaborators. Not a release blocker.

### 4.4 Repeat-focused features (new)

These do not exist in the original and are the practical reason for redoing the tool.

1. **GFF3 annotation overlay**.
   - Load one or more GFF3 files (gzipped accepted) via `noodles-gff`.
   - Each loaded file is a "track" with: name, source attribute filter, color, line
     height, on/off toggle, label rendering toggle.
   - Tracks render alongside the horizontal and vertical axes (one column per axis
     per track). For self-comparison, render on both axes from the same track.
   - Coordinate system: GFF features are mapped to the *full* input sequence range;
     subsetting the dotter window subsets the rendered features.
   - Hover on a feature shows its GFF attributes (a tooltip).
   - Click on a feature jumps the crosshair to the feature's start (or center).
2. **PAF-based HSP overlay**. Load HSPs from a PAF file (e.g. minimap2 output) in
   addition to the original BLAST-style. Render in the existing HSP overlay modes.
3. **Inverted-repeat highlighting**. In self-comparison mode, reverse-strand dots
   are stored in a *separate* channel from forward-strand dots, and rendered in a
   distinct colour (default magenta). Toggle to merge them. Useful for spotting
   inverted repeats, palindromes, satellites with internal symmetry.
4. **Region export**:
   - **Image**: PNG, SVG, PDF of the current view (or of a selected rectangle).
     Embed parameters (window size, pixel factor, zoom, matrix name, input hashes,
     dottir version) in PNG `tEXt` chunks and as an SVG metadata block.
   - **Pixelmap**: dump the raw `Vec<u8>` plus metadata as a `.dot` (original-format,
     read-only compatible) or as a more modern container (`.npz` / `.zarr` / `.h5` —
     pick one; see §10).
   - **Coordinates of selected region**: copy as BED-like text (`seq1<TAB>start<TAB>
     end<TAB>seq2<TAB>start<TAB>end`) for downstream processing.
5. **Alignment export**.
   - At the crosshair (or for a selected diagonal band), export a pairwise alignment
     of the surrounding window as FASTA pair, Stockholm, or simple text. Configurable
     window size (default 100 residues).
   - Optionally invoke an external aligner (system `mafft` / `muscle` / built-in via
     `bio` crate) on the window before export. The simple default just emits the
     ungapped slice.
6. **Multi-sequence input with breaklines** (already in original; called out here
   because it is essential for repeat work on multi-scaffold assemblies). Honour
   the existing breakline rendering and add: hover on a breakline shows the sequence
   name on either side.
7. **Parameter reproducibility sidecar**. Every export action also writes a
   `<output>.params.toml` next to the file with: input file paths, SHA-256 hashes,
   all CLI/GUI parameters at the time of export, dottir version, git SHA, host info.
   Off by default in the GUI (configurable), on by default in CLI batch mode.
8. **CLI batch mode** (`dottir batch ...`):
   - Inputs: query FASTA, subject FASTA, optional GFF3 tracks, optional PAF HSPs.
   - Outputs: PNG/SVG/PDF/.dot, optional alignment dump, params sidecar.
   - All GUI-equivalent parameters exposed as flags or a TOML config file.
   - SHOULD support `--auto-zoom` to target a max output dimension (e.g. 4000px)
     and pick `zoomFactor` accordingly. Avoids surprise OOMs on large inputs.

### 4.5 Non-functional

1. **Performance baseline**: a 1 Mb × 1 Mb BLASTN dotplot at zoom 250 (4000×4000
   pixelmap) MUST complete in ≤ 10 s on a 2024-era 8-core desktop, single-threaded.
   With rayon parallelism, ≤ 2 s. The original takes ~minutes for the same input.
2. **Memory**: pixelmap allocation MUST be checked against `memoryLimit` (default
   0.5 GB, configurable); refuse to allocate larger and suggest a higher zoom factor.
3. **Streaming-friendly input**: FASTA reading via `needletail` or `noodles-fasta`,
   memory-mapped where possible. Gzipped input accepted.
4. **Reproducibility**: documented above (parameter sidecars, deterministic output).
5. **Logging**: `tracing` with `RUST_LOG`-style verbosity. Quiet by default in GUI;
   structured progress in CLI batch.
6. **Error handling**: `anyhow` for application code, `thiserror` for typed errors
   crossing the public API of `dottir-core`. No `panic!` on user input errors.
7. **Documentation**: per-crate `cargo doc` with examples for `dottir-core`'s public
   API. A user manual (mdBook) for the CLI and GUI.
8. **Testing**:
   - Unit tests in each module.
   - Golden tests against the C dotter pixelmaps for a fixed corpus.
   - Property tests (`proptest`) for round-trips: FASTA → pixelmap → save → load →
     equal pixelmap.
   - GUI smoke tests via `egui_kittest` where practical.

## 5. Out of scope (explicitly excluded)

- The exon / intron view from the original. Can be added later if a need emerges; not
  prioritised for repeat analysis.
- Blixem integration. Blixem is a separate tool with separate purpose.
- `libpfetch` (Sanger-internal online sequence fetching).
- The original's bespoke feature-file format; replaced by GFF3.
- Network / cloud functionality.

## 6. Architecture

### 6.1 Crate layout

A Cargo workspace:

```
dottir/
├── Cargo.toml                # workspace
├── crates/
│   ├── dottir-core/          # algorithm, no I/O, no GUI deps
│   ├── dottir-io/            # FASTA, GFF3, PAF, matrix file, .dot, exports
│   ├── dottir-cli/           # batch binary
│   └── dottir-gui/           # interactive binary
├── docs/
│   ├── book/                 # mdBook user manual
│   └── adr/                  # architecture decision records
└── tests/
    ├── golden/               # pinned outputs from C dotter
    └── corpora/              # example inputs
```

Rationale: `dottir-core` having no I/O dependencies means it can be used as a library
from notebooks, other Rust tools, or via Python bindings (`pyo3`) if that becomes
useful. The split also keeps GUI dependencies out of the CLI.

### 6.2 Public API surface (dottir-core)

Stable types and entry points (sketch; final names at developer discretion):

```rust
pub struct ScoreMatrix { /* 24×24 i32, plus name */ }
pub enum BlastMode { Blastn, Blastp, Blastx }
pub enum Strand { Forward, Reverse, Both }

pub struct PlotConfig {
    pub mode: BlastMode,
    pub matrix: ScoreMatrix,
    pub window_size: Option<u32>,   // None => Karlin/Altschul estimate
    pub zoom: u32,
    pub pixel_fac: u32,
    pub strand: Strand,
    pub self_comparison: bool,
    pub memory_limit_bytes: u64,
    // ...
}

pub struct DotPlot {
    pub width: u32,
    pub height: u32,
    pub pixels: Vec<u8>,            // length = width * height
    pub forward_pixels: Option<Vec<u8>>,  // for separate-channel inverted repeats
    pub reverse_pixels: Option<Vec<u8>>,
    pub params: PlotParams,         // resolved window size, λ, K, H, etc.
}

pub fn compute_dotplot(
    query: &[u8],
    subject: &[u8],
    config: &PlotConfig,
) -> Result<DotPlot, DottirError>;

pub fn karlin_window_size(
    matrix: &ScoreMatrix,
    query: &[u8],
    subject: &[u8],
    mode: BlastMode,
) -> Result<KarlinResult, DottirError>;

pub fn greyramp_lut(black: u8, white: u8) -> [u8; 256];
```

### 6.3 Key dependencies

| Concern         | Crate                                |
|-----------------|--------------------------------------|
| FASTA I/O       | `noodles-fasta` (or `needletail`)    |
| GFF3 I/O        | `noodles-gff`                        |
| PAF / SAM I/O   | `noodles-paf`, `noodles-sam`         |
| GUI             | `egui` + `eframe`                    |
| Parallelism     | `rayon`                              |
| CLI parsing     | `clap` (derive)                      |
| Logging         | `tracing`, `tracing-subscriber`      |
| Errors          | `anyhow`, `thiserror`                |
| PNG / colormaps | `image`, `colorgrad`                 |
| SVG             | `tiny-skia` + `resvg`                |
| PDF             | `printpdf` (or render via resvg)     |
| Memory mapping  | `memmap2`                            |
| Config          | `serde`, `toml`                      |
| Hashing         | `sha2`                               |
| Compression     | `flate2` (transparent gzip)          |
| Property tests  | `proptest`                           |

Pin specific versions in `Cargo.toml` once selected; treat new releases as opt-in.

## 7. File formats and interfaces

### 7.1 Inputs

- **FASTA** (`.fa`, `.fasta`, `.fna`, `.faa`), plain or gzipped (`.gz`).
  Multi-record accepted; concatenated with breaklines at record boundaries.
- **GFF3** (`.gff3`, `.gff`), plain or gzipped.
- **PAF** (`.paf`), plain or gzipped — for HSP overlays.
- **BLAST score matrix** files (BLOSUM/PAM format). Built-in: BLOSUM62, BLOSUM50,
  BLOSUM45, BLOSUM80, BLOSUM90, PAM30, PAM70, PAM250. Identity matrix for DNA.
- **Original `.dot` binary** — read for compatibility with archived sessions.

### 7.2 Outputs

- **PNG** (default), **SVG**, **PDF** images.
- **`.params.toml`** sidecar with full provenance.
- **`.dot`** binary (write support optional; see §10).
- **Sequence alignment slices**: FASTA, Stockholm, plain text.
- **BED-style region coordinates** to stdout / clipboard.

### 7.3 Backward compatibility

- The `.dot` binary format from the original (described at `dotplot.c:1610` onward)
  MUST be readable. Format version is encoded in the first byte.
- The original BLAST matrix file format MUST be parsable as-is.

## 8. Validation

### 8.1 Golden tests against C dotter

For a fixed corpus (start with `examples/chr4_ref_seq.fasta` self-comparison,
`examples/Q9H8G1.fasta` × `examples/DA730641.fasta`, plus 3–5 synthetic inputs
including pure tandem repeat, dispersed repeat, inverted repeat), produce
pixelmaps from the C dotter and pin them in `tests/golden/`. The Rust
implementation MUST reproduce these byte-identically when given the same matrix,
window size, zoom, and pixel factor.

Tolerances: zero. If a deliberate algorithmic change is introduced (e.g. fixing an
upstream bug), bump the pixelmap format version and regenerate goldens with a
documented justification.

### 8.2 Property tests

- `pixelmap(reverse(q), reverse(s)) == reverse_xy(pixelmap(q, s))`.
- Self-comparison: `pixelmap(q, q)` is symmetric across the main diagonal (after
  mirror).
- Greyramp LUT is monotone given `black > white`, antitone given `black < white`,
  and saturates at endpoints.

### 8.3 Regression suite

CI runs golden + property tests on Linux x86_64 and Windows x86_64. macOS in CI is
optional. WASM build is gated as a separate CI job that just verifies it compiles
and the test harness runs in headless wasm.

## 9. Phased roadmap

Estimates are rough and assume one developer working part-time; double for full
"production-ready" polish.

| Phase | Deliverable                                                         | Effort |
|-------|---------------------------------------------------------------------|--------|
| 0     | `karlin.rs` ported, with tests pinning λ/K/H against C output.      | 1 wk   |
| 1     | Single-strand BLASTN inner loop, golden-test against C on chr4.     | 1 wk   |
| 2     | BLASTP and BLASTX, dual-strand BLASTN, self-comparison + mirror.    | 1–2 wk |
| 3     | rayon parallelism, SIMD where straightforward, regression unchanged.| 1 wk   |
| 4     | CLI batch mode: PNG/SVG/PDF export, `.dot` read, params sidecar.    | 1–2 wk |
| 5     | GUI MVP: pan/zoom, greyramp, crosshair, alignment view.             | 2–3 wk |
| 6     | GFF3 + PAF overlays, region selection + export, alignment export.   | 2 wk   |
| 7     | Inverted-repeat channel, multi-sequence breaklines, polish, docs.   | 1–2 wk |
| 8     | WASM build (optional).                                              | 1 wk   |
| 9     | Releases: GitHub Actions builds for Linux + Windows.                | 1 wk   |

Total: ~12–18 person-weeks for the full scope, less if Phase 8 is skipped.

The first three phases produce a usable headless tool. The first GUI release is at
end of Phase 5. Phase 6 is where the tool meaningfully exceeds the original.

## 10. Open questions / decisions to make

These should be resolved with the developer early; they affect downstream choices.

1. **Final project name and crate names.** Check crates.io availability.
2. **License.** GPLv3 (simpler, matches upstream) vs. permissive (broader adoption).
3. **FASTA library.** `noodles-fasta` (richer ecosystem, slightly heavier) vs.
   `needletail` (lighter, very fast). Both work; pick after a 30-min benchmark.
4. **Modern pixelmap container.** When dumping the raw matrix for downstream use,
   options are: stay with the original `.dot` format only; add `.npz` (NumPy
   archive — easy interop with Python/R); add HDF5 via `hdf5-rust` (heavyweight);
   add Zarr. **Recommendation: NumPy `.npz` for Python interop, since the dottir-core
   API also gives Rust users `Vec<u8>` directly.**
5. **Alignment-on-the-fly.** Provide a built-in pairwise aligner (e.g. `bio` crate's
   Smith–Waterman) for the alignment export, or shell out to `mafft`/`muscle`?
   **Recommendation: built-in for short windows (<= 1 kb), shell out optional.**
6. **GUI framework.** This spec assumes egui. Alternatives — iced, gtk4-rs, Tauri —
   should only be revisited if a specific need surfaces. **Decision: egui unless
   blocked.**
7. **HSP source.** Built-in lightweight aligner for HSP generation, or only consume
   externally-generated PAF/BLAST output?
   **Recommendation: consume external only.** Building an aligner inside dottir is
   feature creep; minimap2/blast already exist.
8. **WASM.** Worth doing as part of v1 (free with egui), or post-v1?
   **Recommendation: post-v1.**
9. **Python bindings.** A `pyo3` wrapper around `dottir-core` would let you call it
   from notebooks. Defer to post-v1 unless an early need appears.
10. **Style guide and CI.** Adopt `rustfmt` defaults, `clippy --deny warnings` in CI,
    `cargo-deny` for license auditing, `cargo-audit` for advisories. Decide now to
    avoid retrofitting.

## 11. References

- Sonnhammer ELL, Durbin R. *A dot-matrix program with dynamic threshold control
  suited for genomic DNA and protein sequence analysis.* Gene 167(2):GC1-10, 1995.
  (Foundational paper. The window-size-from-Karlin-Altschul reasoning is laid out
  here.)
- Karlin S, Altschul SF. *Methods for assessing the statistical significance of
  molecular sequence features by using general scoring schemes.* PNAS 87:2264-2268,
  1990. (Source of the λ/K/H computation.)
- Altschul SF, Erickson BW. *Significance of nucleotide sequence alignments: a
  method for random sequence permutation that preserves dinucleotide and codon
  usage.* Mol Biol Evol 2:526-538, 1985.
- Seqtools repository: `seqtools-4.28/` — primary source-of-truth for behaviour.
- noodles bioinformatics crates: https://github.com/zaeleus/noodles
- egui: https://www.egui.rs/

---

*This document is a living specification. Treat it as a contract between the project
owner and the developer for v1 scope; deviations from MUST requirements need explicit
agreement and an ADR entry.*
