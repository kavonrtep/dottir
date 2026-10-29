# Changelog

All notable changes to dottir are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
the project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

Performance (findings HIGH 2, HIGH 3, HIGH 4 and MEDIUM 7 of
`docs/efficiency-review-large-genome.md`; output unchanged):

* `--sort-records` no longer prepares the large engine twice: the plot
  reuses the sort pass's index and depth tracks where the sorted input is
  unchanged (`LargeSession::prepare_reusing`). TAIR10 × Col-CC, 8 cores:
  170 → 123 s when the order stays, 180 → 141 s when query records move.
* GUI *View → Sort records* computes its chains on the dotplot's prepared
  session instead of preparing another one, and the compute worker
  re-prepares by reusing what still matches (after sorting, or a window
  size change).
* PAF reading: one pass for HSP overlays instead of one per orientation
  (3 M lines: 3.6 → 2.0 s), no per-line allocation, and deriving
  primaries for `-P` PAFs sorts compact keys and scans only overlapping
  primaries (`read_paf_file`, 3 M lines: 2.4 → 2.0 s, gzipped
  3.0 → 2.1 s). `examples/paf_read_bench.rs` in `dottir-io` times both.
* Genome mode parses the aligner's PAF as it arrives instead of holding
  the whole text first (`paf::PafParser`): TAIR10 × Col-CC `-P` (889 k
  alignments) peaks at 154 MB instead of 225 MB in dottir, and the
  saving grows with `-P` output on large genomes. `--paf-min-len` still
  applies when the cached result is read, after alignments split at
  query chunk cuts are joined.
* The large engine's seed index is built on all cores (chunked windows,
  atomic bucket counters, each bucket sorted afterwards so the index is
  unchanged). Col-CC (142 Mb): 15.5 → 5.2 s per index; TAIR10 × Col-CC
  prepare 61.7 → 42.2 s. `examples/index_build_bench.rs` in
  `dottir-core` times it.

## [0.6.0] - 2026-10-29

Large-sequence and genome-scale comparisons, and pixelmaps identical to
C dotter. (0.5.0 was planned for large mode alone and not released.)

### Added

* **Genome mode** (ADR 0008, spec §4.7): whole-assembly dotplots at
  D-GENIES scale, drawn from alignments instead of matrix cells.
  - `dottir batch REF.fa ASM.fa --engine genome` runs mm2plus or minimap2
    (`--aligner`, `--preset`, `--aligner-args`, `--threads`): query
    records are streamed in `--chunk-size` pieces and re-joined, results
    are cached by input checksums, aligner version and arguments.
    `--paf FILE` draws existing alignments instead. Self-comparisons use
    `-D --dual=no` (`--self-all-copies` for `-P`).
  - Identity-coloured lines (D-GENIES classes), filters (`--min-aln-len`,
    `--min-identity`, `--primary-only`, `--remove-noise`), small-record
    blocks, contig sorting and orientation (`--sort-query`), following the
    D-GENIES source.
  - Exports: `--assoc-out`, `--sorted-query-out`, `--unaligned-out`,
    `--paf-out`; `[genome]` sidecar section.
  - `--engine auto` selects genome mode from `--genome-threshold-cells`
    (10¹⁷) when an aligner is installed.
  - Inputs are indexed (samtools-compatible `.fai`) and never loaded;
    gzipped FASTA is decompressed once into the user cache.
  - `dottir-gui --genome` / `--paf` and *File → Open genome comparison*:
    background alignment with progress and cancel; zoom, pan, hover
    details, display-only filters, sorting, per-record flips, exports;
    drill-down of the view (or a record pair) into the large / exact
    dotplot, with positions reported in the original records; sessions
    store the comparison.
* **Record sorting for dotplots** (spec §4.7.11): `--sort-records
  vertical|horizontal` in `dottir batch` and `dottir-gui`, *View → Sort
  records* (with *File order* to undo) — the D-GENIES contig order and
  orientation from the large engine's synteny chains, for either input
  and either engine. Coordinates and GFF3/BED overlays stay in source
  records; the order is recorded in the sidecar (`[record_sort]`) and in
  sessions. Annotation overlays now also follow genome drill-down
  regions.
* Genome view uses the dotplot layout: controls in a right-hand panel with
  the dotplot's compact fonts, pointer position in the status bar, shared
  menu with *Fit*, *Back* and *Save PNG* acting on the genome plot.
* Genome mode: `--all-copies` (was `--self-all-copies`, kept as an alias)
  now applies to pairwise comparisons too, and the GUI panel has a *Keep
  all copies of repeats (-P)* checkbox (off by default) that re-aligns
  while keeping filters, sort and flips. PAFs without any primary (as
  `-P` writes them) get minimap2-style primaries derived on reading, so
  *Primary only* no longer empties the plot.
* **Large-sequence mode** (ADR 0007, spec §4.6): a seed-based engine for DNA
  comparisons at chromosome/genome scale (tested to 100 Mb × 100 Mb,
  multi-record).
  - Spaced-seed index + two-hit diagonal filter; anchored diagonals are
    scored with the exact sliding-window score, so large-mode pixels never
    exceed the exact engine's and the greyramp keeps its meaning.
  - Repeat classes (`unique`/`low`/`mid`/`high`) per dot from a per-base
    copy-number track — `seed` (diverged repeats), `kmer`, or `none` — plus
    soft-masking and GFF3/BED repeat annotations. Classes are separate
    channels: show, hide or colour them without recomputing.
  - Basic co-linear chaining of unique/low-copy runs into synteny blocks,
    exported as PAF (`--chains-out`, no CIGAR).
  - `dottir batch`: `--engine exact|large|auto`, `--seed-pattern`,
    `--two-hit-dist`, `--max-occ`, `--depth`, `--class-thresholds`,
    `--hide-classes`, `--softmask-class`, `--repeat-annot[-subject]`,
    `--repeat-types`, `--chains-out`, `--depth-out` (bedGraph). Engine and
    parameters are recorded in PNG text chunks and in a new `[large]`
    sidecar section.
  - `dottir-gui`: *Settings → Large sequences*; `Auto` uses the large engine
    for big views and the exact engine once zoomed in. The index is built
    once and reused; progress is shown in the status bar and superseded jobs
    are cancelled. Side panel: class toggles and presets, colour by class,
    synteny-chain overlay. Settings are saved in sessions.
  - GUI sensitivity for large plots: shade by **matching bp per pixel**
    (new optional coverage channel, `LargeConfig::coverage`) with
    white/black sliders in bases and per-view auto-fit; Dotter score
    shading remains selectable.
  - GUI repeat tracks along both axes: the profile of the visible dots
    (aligned copies per base, log scale, split by class) and capped-seed
    marks, with a hover tooltip adding the bases' own copy number
    (`LargePlot::{q,s}_depth_bins`).
  - Synthetic benchmark generator:
    `cargo run --release -p dottir-cli --example simulate_genome`.
  - Seed sampling for overview views (`--seed-sample`, `--depth-sample`,
    GUI settings): one seed in N by seed hash, N from the zoom (plot) and
    sequence length (depth pass); the 100 Mb plot step is about 2× faster
    at 94 % of the lit pixels. Scoring emits once per pixel crossing
    instead of per cell. Large pixelmap format version 2.
  - Chains: minimap2-style masking (a chain mostly inside a better one on
    both axes is dropped) and a 3-anchor minimum. Without them, seed
    sampling fragmented synteny blocks into side chains (100 Mb pair:
    46,064 chains; now 58, against 1,691 before sampling).
  - Coverage no longer saturates: 64-bit counters and an encoding reaching
    ~3.8 Gb per pixel (was ~63 kb), and cells are scored exactly once even
    where runs meet at chunk boundaries (coverage was double-counted there).
* FASTA parsing records soft-masked (lowercase) runs
  (`FastaRecord::soft_masked`, `Sequence::soft_mask`).
* **GUI: repeat annotations for large sequences.** *Settings → Large
  sequences → Repeats (query/subject)* takes a GFF3/BED of known repeats
  (e.g. RepeatMasker) whose bases count as `high`, with a feature-type
  filter; also `dottir-gui --repeat-annot/--repeat-annot-subject/
  --repeat-types`, saved in sessions and passed to *Open in new window*.
* **BLASTX in the GUI alignment view**: the query translated in its
  best-scoring frame around the crosshair (`F1`–`F3`), instead of "not yet
  supported".
* **Genome mode: `--paf-min-len BP`** (GUI: *Skip alignments shorter
  than* in *Open genome comparison*) drops short alignments while the PAF
  is read — PAF files and aligner output — so huge `-P` PAFs of large
  genomes fit in memory (Col-CC self `-P`, 1.0 M lines: 192 → 100 MB at
  5 kb). Recorded in the sidecar.
* **Genome mode: sort the reference along the query**: `--sort-reference`
  (exclusive with `--sort-query`) and *Sort reference along query* in the
  genome view's *Order*; kept in sessions (`sort_reference`) and in the
  sidecar (`sorted_reference`). `arrange` takes an `Option<SortAxis>`.
* **Chain exports** (plan L5): `--chains-out FILE.bedpe` writes BEDPE;
  `--svg` draws the chains as vector lines; GUI *Save chains…* (PAF or
  BEDPE) and chain lines in *Save SVG* when shown.
* **Genome-mode SVG** with vector lines: `--svg` in genome mode and
  *File → Save SVG…* in the genome view (one line per drawn alignment,
  coloured by identity class; separators, mix blocks and record names as
  vectors and text). `dottir_core::genome::plot_lines` gives the PNG
  raster and the SVG the same lines.
* **PDF export** (spec §4.4.2, plan D3): *File → Save PDF…* and
  `dottir batch --pdf PATH` — one page with the plot and its axes, colours
  and HSPs as shown, metadata in the document information. Written
  directly (one Flate image), no new dependency.
* **Selections** (plan D4, G2): a right-button drag selects a rectangle
  with *Zoom to selection*, *Open in new window* and *Copy coordinates
  (BEDPE)*. The new window is `dottir-gui --query-region/--subject-region
  NAME:START-END[:-]`, loading just those regions through the FASTA index
  in their source coordinates.
* **HSP overlays from PAF** (spec §4.2.3, plan E2/E4): *File → Load HSPs
  (PAF)…*, `--hsp-paf` / `--hsp-mode` in `dottir-gui` and `dottir batch`.
  Modes off / line / score (identity-coloured) / greyscale (replaces the
  dots), as the original Dotter. Either PAF orientation; records not
  plotted are skipped; sorted, flipped and region records are followed.
  Kept in sessions and drawn into saved images.
  Alignments with a `cg:Z` CIGAR are drawn along it (gapless stretches,
  ≤ 10 bp diagonal drift per line) instead of start to end, so long
  alignments with indels stay on their dots.
* GUI: while a computation runs, the plot is veiled and a badge in its
  middle shows a spinner and the progress (it was only in the status
  bar, easy to miss after panning a large plot).
* `--export-alignment-at Q,S` (`--alignment-window`, `--alignment-format
  fasta|stockholm|text`) prints the residues around a dotplot position;
  `-o` becomes optional when only exporting (plan D5).

### Changed

* `dottir batch` defaults to `--engine auto`: BLASTN inputs above 10¹¹
  cells (≈ 316 kb × 316 kb, about a minute of exact compute) now use the
  large engine and, without an explicit
  zoom, are fitted to 4000 px. Pass `--engine exact` for the old behaviour.
* The GUI compute worker skips queued requests superseded by a newer one.
* New `DottirError::Cancelled` variant.

### Fixed

* Pixelmaps now match the C dotter exactly, checked against 15 pinned C
  dotter pixelmaps (ADR 0010; pixelmap format version 2). Fixed on the
  way: protein Karlin statistics over the 20 standard amino acids with
  the subject's full length, as C (windows and pixel factors of protein
  and BLASTX plots, and of DNA plots whose subject has `N`s, change); the
  automatic pixel factor truncates; the reverse strand's last row was
  missing; BLASTX plots could be a column narrower; BLOSUM62 `*:*` is −1
  as in C. A Karlin window estimate below 3 now becomes 10, as in C.
* GUI in BLASTX mode treated a pixelmap column (one codon) as one base:
  the crosshair, rectangle zoom, selections, auto zoom, record breaks,
  annotation bands, HSP lines and axis labels were off by 3 on x. They now
  use base positions with 3 bases per column at zoom 1.
* GUI image exports after zooming in (rectangle zoom, wheel past the
  slice) labelled the slice's image as if it started at residue 0 and
  spanned the whole sequence; axes and record names now show the slice
  (`compose_*_with_axes_at`, `write_rgb_svg_with_lines` take an origin).
* SVG exports at zoom > 1 labelled the axes in pixelmap pixels and dropped
  the names of records beyond the pixel width; ticks and names are now in
  residues (`write_svg` / `write_rgb_svg` take the coordinate size, like the
  PNG writers).

* **BLASTX window and pixel factor follow C dotter** (ADR 0009): the
  automatic estimate is the mean over the three reading frames, not frame
  0 alone. `PIXELMAP_FORMAT_VERSION` is now 1; BLASTN/BLASTP output is
  unchanged. (BLASTX itself already worked; the docs said otherwise.)
* Axis tick labels on fine steps no longer repeat (ticks 1050/1100/1150
  were all "1.1k"); they get the decimals they need.

Findings of `docs/reviews/since-v0.4.0.md`:

* Genome mode: joining alignments split at query chunk cuts could drop
  overlapping short hits near a cut (they chose each other as partners).
  Only hits on opposite sides of a cut are joined now; cached runs are
  re-aligned (cache key v2). Col-CC × *A. lyrata* gains 3 alignments.
* `--chains-out` and `--depth-bedgraph` after `--sort-records` (or on
  drill-down regions) reported transformed coordinates under source
  names; they now give source-record positions, lengths and strands.
* Cancelling an aligner run waited for its next output line, and
  indexing could not be cancelled; both now stop within 50 ms.
* GUI PNG/SVG export applied the score greyramp in coverage shading, so
  the coverage sensitivity was lost; exports now use the on-screen
  mapping and record the shading in their metadata. They now also keep
  the screen's colours (repeat classes, inverted repeats in magenta) as
  RGB images.
* FASTA indexing rejected a full-width last line without a newline.
* A PAF made with the axes swapped failed to load when both assemblies
  share record names (`chr1`) with different lengths.
* `LargeSession::compatible` accepted automatic window size / pixel
  factor for a session prepared with manual values.

## [0.4.0] - 2026-09-11

### Added

* **GFF3 with embedded sequences as a sequence input** (ADR 0006). A GFF3
  file that carries its sequences after a `##FASTA` directive can be passed
  directly as a positional argument — one file for a self-comparison, two for
  a pairwise comparison — with no separate FASTA and no `--gff*` flag:

  ```sh
  dottir batch elements.gff3 -o elements.png
  dottir-gui family_a.gff3 family_b.gff3
  ```

  - Applies to `dottir batch`, `dottir periodogram`, the `dottir-gui`
    positional arguments, and **File → Open query/subject sequence…**.
  - In the GUI the features from the same file are bound as that axis's
    annotation overlay automatically (both axes in self-comparison).
    `--gff-query`/`--gff-subject` and a later **File → Load GFF3/BED…** still
    override them. Batch mode has no overlay and logs the feature count
    instead.
  - Plain and gzipped both work; format is chosen by extension, with a
    `##gff-version` content sniff as a fallback.
  - A GFF3 with no `##FASTA` section is rejected as a sequence input with a
    message naming both fixes; it remains a valid annotation input.
  - New corpus `tests/corpora/gff3_with_sequence/`, regenerated by
    `scripts/make_gff3_with_sequence_fixtures.py`.

### Fixed

* `AnnotSet::load` now splits off a `##FASTA` section before parsing, instead
  of passing the embedded sequence lines to the GFF3 record parser.
* An input that is neither FASTA nor GFF3 now names the accepted formats
  instead of reporting `missing '>' header at start`.
* Sequence-load errors no longer repeat the path, which the CLI context and
  the GUI status bar already carry.

## [0.3.0] - 2026-06-30

### Added

* **Annotation overlay (GFF3 / BED)** — load annotation intervals and
  display them on the dotplot as full-span colored bands (ADR 0005). A
  query feature paints as a vertical band across the plot height; a
  subject feature as a horizontal band across the width; overlaps
  self-darken so annotated region-pairs stand out. Works for both
  self-comparison (one file applied to both axes) and pairwise
  comparison (one file per axis).
  - **GFF3** (via `noodles-gff`) is colored by a chosen attribute value
    (default `Name`, plus a synthetic `(type)` key for column 3); a
    feature missing the chosen attribute falls back to its GFF3 type
    rather than a single "none" bucket. **BED** uses a single color per
    file. Plain or gzipped.
  - The right panel carries the **Annotations** controls: a show-bands
    toggle, an **opacity** slider (dims bands so dots show through), and
    a compact **legend** with per-value show/hide. In self-comparison the
    single legend drives both axes. Per-value **color editing** is in
    **Settings → Annotations**.
  - The **alignment view** lists every annotation overlapping the
    crosshair residue on each axis (query and subject) and frames the
    crosshair cell in the feature's color.
  - Load via **File → Load GFF3/BED for query/subject…** (collapsing to
    one entry in self-comparison) or the CLI flags `--gff-query` /
    `--gff-subject` / `--gff` (with `--bed*` aliases; format
    auto-detected by extension).
  - Bundled fixtures under `tests/corpora/annotation_overlay/`
    (`ltr_angela`, `tir_elements`, `tir_simple`) for trying the overlay.

## [0.2.2] - 2026-05-24

### Fixed

* **Conda package on the petrnovak channel** — the 0.2.0 and 0.2.1
  Linux builds were corrupted by conda-build's install-time binary
  prefix replacement. Conda-build pads the build prefix string in
  the binary out to a 255-char placeholder; at install time conda
  does a *variable-length* (not null-padded) byte replacement of
  the placeholder with the actual install prefix. Because
  `x11-dl`'s symbol-name lookup table happens to live in
  `.rodata` immediately after a build-prefix string in our
  binary, every offset in the table ends up ~230 bytes off, and
  at runtime `dlopen()` is called with garbage substrings of the
  symbol buffer ("WMHintsXAll", "tsXAllowEve") instead of
  "libX11.so.6". `dottir-gui` failed at startup with
  `Failed to load one of xlib's shared libraries` even on
  systems with X11 installed. The conda recipe now sets
  `build/detect_binary_files_with_prefix: false` to skip the
  destructive prefix replacement, and bakes RPATH=`$ORIGIN/../lib`
  at link time via RUSTFLAGS so install-time path adjustment
  isn't needed.

  The 0.2.1 release shipped a wrong fix (disabling LTO) that did
  not address the real root cause; superseded by this release.

## [0.2.1] - 2026-05-23

### Fixed

* Attempted conda recipe fix (disable LTO) for the
  `dottir-gui` X11 startup failure on Linux. This was a
  misdiagnosis — see 0.2.2 for the actual fix.

### Added

* README install section: conda install instructions for the
  `petrnovak` channel.

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
