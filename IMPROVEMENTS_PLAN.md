# dottir — Improvements Plan

Concrete follow-up plan addressing the findings in `REVIEW.md`. Estimates
are person-days for one focused developer; double for full polish + tests
+ docs.

Each phase is independently shippable (commit per phase, mirroring how the
original `IMPLEMENTATION_PLAN.md` was executed). The phases are ordered by
priority: A is correctness/scalability we shouldn't ship without; B–C is
parity with what the docs and spec already promise; D–F is feature
completion.

| Phase | Theme | Effort | Blocks |
|-------|-------|--------|--------|
| A | Correctness & scaling | 3-4 d | nothing |
| B | Docs / status sync | 0.5 d | A (so docs reflect post-A state) |
| C | GUI usability gaps | 4-5 d | A3 (record boundaries) |
| D | Export layer | 4-6 d | A3 |
| E | GFF3 / PAF track overlays | 4-5 d | A3, C |
| F | BLASTX + remaining algorithm gaps | 3-4 d | nothing |
| G | Inverted-repeat visualisation & sub-dotter | 2-3 d | C |

---

## Phase A — Correctness and scalability

These are subtle correctness bugs and a memory-cap violation. They block
honest claims about the spec §4.5.2 memory limit and any large-input use.

### A1. Bounded-memory parallel kernel

**Problem.** `plot.rs:329-380` allocates `n_chunks + 1` full pixelmaps
(one per rayon worker, plus the destination) and collects them all before
the final merge. With `target_chunks = 4 × n_threads`, peak memory is
~`(4n + 1) × W × H` bytes. The `memory_limit_bytes` check guards each map
in isolation, so a 256 MB cap with 8 threads can quietly hit 8.5 GB of
peak RSS.

**Fix.**

1. **Replace per-chunk buffers with one shared output and per-chunk
   dot-emit channels.** Each worker emits pixel updates (or a small
   tiled buffer) into the shared map via atomic max on `u8` lanes — or
   into a thread-local *fixed-size* tile that's merged immediately after
   each chunk completes (`reduce_with` instead of `collect`).
2. Reflect the cap honestly: the `memory_limit_bytes` check becomes
   `(W × H) + (n_threads × tile_bytes) ≤ limit`. Document the budget
   in `plot.rs`.
3. Add a `parallel_strategy` knob on `PlotConfig` (`Serial | Parallel { max_workers }`)
   so callers can opt out — the GUI in particular benefits from running
   the compute on a background thread but with fewer rayon workers.

**Files**: `crates/dottir-core/src/{plot.rs,sliding.rs,pixel.rs}`.

**Tests**: the existing `parallel_determinism.rs` already covers byte
identity across thread counts; add a peak-RSS regression that asserts
no pixelmap is allocated more than once per pass at any thread count
(use a global allocator-wrapper for the test, gated behind a `cfg(test)`
counter so it doesn't ship).

**Acceptance**: a 100 Mb × 100 Mb BLASTN run at zoom 250 with
`memory_limit_bytes = 256 MiB` stays under 300 MiB peak RSS regardless of
thread count, and produces byte-identical output to the n=1 baseline.

---

### A2. Strict matrix parsing; separate DNA path

**Problem.** `ScoreMatrix::parse_blast_format` silently drops unknown
header letters (the `J` workaround) and back-fills missing cells with
`-4`. A typo in a user matrix becomes a silent corruption that quietly
changes the dotplot.

**Fix.**

1. Split the parser into two: `parse_blast_protein` and
   `parse_blast_dna`. Protein expects exactly the 24-letter alphabet;
   DNA expects exactly `ACGT` (plus optional `N` row/col that we score
   as zero against everything, explicitly, not via fall-through).
2. Accept extra letters only via an explicit `extra_letters` allow-list
   parameter (default empty). NCBI's `J` then needs an explicit
   `&["J"]` from the caller.
3. Make every cell present. Return `DottirError::InvalidMatrix` with the
   missing (row, col) for any uncovered pair.
4. Add a `Matrix::validate()` method (called from `parse_*` and the
   built-in constructors) that enforces: at least one negative score,
   at least one positive score, symmetric if claimed symmetric.
5. While we're here, expose a `ScoreMatrix::custom_dna(match, mismatch)`
   constructor — there's no built-in for users who want `+1 / -1` or
   `+2 / -3`.

**Files**: `crates/dottir-core/src/matrix.rs`.

**Tests**: parser rejects missing-cell input; parser accepts `J` only
with allow-list; round-trip still works for all built-ins; the existing
NCBI BLOSUM/PAM files load via `parse_blast_protein` with
`extra_letters = &["J"]`.

**Acceptance**: a malformed BLAST matrix file produces a precise error
naming the offending row/column. The "silently fill with -4" path no
longer exists.

---

### A3. Multi-record sequence model: keep the boundaries

**Problem.** `dottir-io::fasta::concatenate` is the only path the rest
of the code uses. It discards offsets, IDs, and descriptions, which
makes spec §4.4.6 (breaklines) and §4.4.5 (alignment export with sane
coords) impossible without re-parsing.

**Fix.**

1. Introduce `dottir_io::Sequence` (struct with `seq: Vec<u8>`,
   `records: Vec<RecordSpan { id, description, range }>`,
   `source_path: PathBuf`). Cheaply produced from `Vec<FastaRecord>`;
   keeps the concatenated buffer that `dottir-core` wants.
2. Add a `breaks(&self) -> &[usize]` helper that returns the
   inter-record offsets in concatenated coords (used by GUI breaklines
   and PNG axis-render).
3. Update `dottir-cli::main::run_batch` and `dottir-gui::app` to load
   into `Sequence` and pass `seq.bytes()` to `compute_dotplot`. Track
   the source path for the params sidecar.
4. Add a helper `Sequence::record_at(coord) -> Option<(&RecordSpan, usize)>`
   returning `(record, position-within-record)`. Used by the GUI status
   bar so it shows `chr4:12345` rather than `12345`.

**Files**: new `crates/dottir-io/src/sequence.rs`; consumers updated.

**Tests**: round-trip a multi-record FASTA → `Sequence` → query by
coordinate; `breaks()` returns the expected offsets; the existing
fasta parser tests stay green.

**Acceptance**: every consumer that needs to know "what record am I in"
asks `Sequence`, not the raw bytes. `concatenate()` is deprecated
(kept for one release, then deleted).

---

### A4. Single-pass FASTA read in CLI

**Problem.** `dottir-cli/src/main.rs:172-180` reads each FASTA file twice
in the common case: once with `std::fs::read` for `String::from_utf8_lossy
+ parse_fasta`, and once via `fasta::read_fasta_file` as a fallback.

**Fix.** Replace with a single streaming load — `read_fasta_file` now
returns the source bytes too (for the sidecar SHA-256) without holding
the whole file twice.

```rust
pub struct LoadedFasta {
    pub records: Vec<FastaRecord>,
    pub bytes: Vec<u8>, // raw on-disk bytes, for hashing
}
```

**Files**: `crates/dottir-io/src/fasta.rs`, `crates/dottir-cli/src/main.rs`.

**Acceptance**: the CLI reads each file exactly once. Peak memory drops
by `2 × file_size` on the load path.

---

## Phase B — Docs / status sync

Reviewer is right that `docs/book/src/intro.md` and `install.md` still
claim "GUI deferred" / "MSRV 1.75". Quick to fix; do this immediately
after Phase A lands so the docs describe the post-A state.

### B1. Rewrite `intro.md`, `install.md`, `cli.md`, `adr.md` to current state

- Remove "GUI deferred" notes; describe `cargo run -p dottir-gui` as the
  default interactive entry point.
- Update MSRV claim to 1.85 (and note that the dev environment uses
  1.95 from conda-forge).
- Update `cli.md` with any new flags from Phases C-D as they land
  (`--memory-limit`, `--svg`, `--dot`, `--gff3`, `--paf`).
- Add a Phase 5 / GUI page describing the panels and the keyboard map.

### B2. CHANGELOG.md: add the "Released MSRV bump + GUI MVP" section

Move the "egui deferred" bullet from the deferred list to the shipped
list. List the new file open / save PNG / greyramp / settings actions.

### B3. ADR 0003 supersession block

`docs/adr/0003-gui-msrv.md` currently says "Superseded by the
MSRV-1.85 bump"; add a one-line `Supersedes-link` to the actual commit
SHA so the trail is unambiguous from the ADR alone.

---

## Phase C — GUI usability gaps

What's in place: file menu, save PNG, greyramp, basic settings, pan/zoom,
crosshair, status bar. What the reviewer flagged + what spec §4.3-§4.4
expect:

### C1. Surface the memory cap as a setting

Currently hardcoded at 1 GiB in `dottir-gui/src/app.rs:197`. Add a
`memory_limit_mib` field to `Settings`, expose as a numeric input in the
Settings window with sensible bounds (8 MiB … available system RAM).
Default to 512 MiB to match the CLI default — surfacing the inconsistency
in the review.

### C2. Zoom quality — recompute on zoom-settle

Implements the review's "Zoom Quality" section.

1. **Debounced recompute**. Wheel-zoom adjusts `display_zoom`
   immediately for a smooth interaction. When wheel events stop for
   `>200 ms`, schedule a recompute with `PlotConfig::zoom = pixel_per_residue`,
   running on a background `std::thread` so the UI stays responsive.
   On completion, swap the texture in.
2. **Multi-resolution cache**. Keep a `BTreeMap<u32, DotPlot>` indexed
   by computation zoom. On zoom settle, look up the closest tier and
   only recompute if more than a 2× factor away. Evict oldest tiers
   beyond a fixed cap (~3 stored plots).
3. **Snap tiers**. Computation-zoom snaps to a power-of-2 ladder
   (1, 2, 4, 8, 16, …) so the cache hits sensibly. Viewport zoom stays
   continuous between snaps so motion is smooth.
4. **Progress indicator**. A spinner + "Recomputing at zoom X…" banner
   while a background recompute is in flight; the current texture stays
   visible (just at lower effective DPI) until the new one swaps in.

**Files**: `crates/dottir-gui/src/app.rs` (the bulk), maybe new
`compute_worker.rs` for the worker thread machinery.

**Acceptance**: zooming into a detail of a 10 Mb plot eventually shows
window-W-sized features cleanly; user never waits on the UI thread.

### C3. Breaklines for multi-record FASTA

Requires A3 (`Sequence::breaks`). Render thin coloured lines on the
canvas at the break offsets in both axes; tooltip on hover shows the
record names. Honour `--breakline-colour` (spec §4.2.6).

### C4. Axes with sequence-coordinate labels

Reviewer's "Recompute overlays and axis labels from sequence
coordinates" point. Render:

- Top axis: query coord, with ticks at every `10^k` boundary and labels
  every major tick.
- Left axis: subject coord, same.
- Tick spacing adapts to `display_zoom` (more ticks when zoomed in).
- Labels are anchored to sequence coords, so they stay readable when
  the underlying pixelmap is replaced at a different computation zoom.

Use egui's painter directly; no extra deps.

### C5. Session save/load

`PlotConfig` + greyramp + crosshair + view transform + input file paths →
TOML at `<name>.dottir-session.toml`. Round-trip via `serde`. File menu
gets "Save Session" / "Open Session".

This is also how the deferred review item "no shared session file format"
gets resolved.

### C6. Background compute & status

Already implied by C2 but worth stating: any compute longer than ~100 ms
moves off the UI thread. Failures (OOM, invalid matrix) surface as a
red toast in the status bar instead of blocking the UI.

### C7. `dottir-gui` command-line: pre-load sequences

Mirror the original `dotter` program's invocation: `dotter q.fa s.fa`
opens with both sequences pre-loaded and the dotplot already computed.
Today dottir-gui takes no CLI args (or three positional args used as a
headless fallback before Phase 5 landed — now obsolete).

Replace the manual `std::env::args` parsing in `dottir-gui/src/main.rs`
with a `clap` parser:

```text
dottir-gui [QUERY] [SUBJECT] [OPTIONS]

Positional:
  QUERY     Optional query FASTA path; loaded at startup.
  SUBJECT   Optional subject FASTA path; loaded at startup.

Options (mirror PlotConfig + the original Dotter flags):
  -W, --window N        Window-size override.
  -z, --zoom N          Computation zoom.
      --mode MODE       blastn | blastp.
      --matrix NAME     Score matrix name (built-in).
      --self            Self-comparison.
      --memory-mib N    Memory cap; surfaces C1 setting.
  -h, --help / -V, --version
```

If `QUERY` and `SUBJECT` are both supplied, the app `recompute()`s
immediately after window creation; otherwise it starts with the
no-plot placeholder (current behaviour).

**Acceptance**: `dottir-gui examples/chr4_ref_seq.fasta examples/chr4_ref_seq.fasta --self -W 25`
opens with a computed self-comparison ready for inspection.

### C8. Light colour scheme as the GUI default

Reviewer's observation: egui's default visuals are dark, but the
plotting area is essentially light (greyscale pixelmap on near-white
background). The mismatch makes axis labels and panel text hard to
read against the plot. Switch the default to
`egui::Visuals::light()` in `DottirApp::new`, and add a "Dark theme"
toggle in the View menu for users who prefer the contrast inversion.

---

## Phase D — Export layer

Reviewer's "Add an explicit export layer" item. Build it once, share
between CLI and GUI.

### D1. SVG export

`dottir-io::svg_export::write_svg(path, plot, params)` using `tiny-skia`
+ `resvg` (or a small hand-rolled SVG writer — pixelmap as
`<image href="data:image/png;base64,...">` plus axes/ticks as SVG
primitives). Embed parameters in an `<svg><metadata>...</metadata>` block.

**CLI**: `--svg <PATH>` (in addition to `-o`).
**GUI**: File → Export → SVG…

### D2. `.dot` binary format

Read C-dotter `.dot` files for backward compatibility (spec §7.3 says
this is a MUST). Write support optional but cheap once read is in.

`dottir-io::dot_format::{read, write}`. Reader maps a `.dot` to a
`DotPlot` so the GUI / CLI can load archived sessions.

### D3. PDF export

Two implementations to compare:

1. Rasterise to a high-DPI PNG and embed in a one-page PDF via
   `printpdf`.
2. Use `resvg` to render the SVG from D1 directly to PDF.

Pick whichever ships smaller and stays cross-platform. Option 1 is
simpler.

### D4. Selection-region export

Rubber-band selection in the GUI (right-click drag) → rectangle in
sequence coords. Actions:

- Copy BED-style coordinates to clipboard.
- Export selection as PNG/SVG/PDF.
- Spawn sub-dotter (see Phase G).

### D5. Alignment-window export

Already partially landed in `dottir_io::alignment::slice_pair`. CLI
flag `--export-alignment-at <q,s> --window N --format fasta|stockholm|text`.

---

## Phase E — GFF3 and PAF track overlays

Now unblocked by the MSRV bump (ADR 0004 was specifically about this).

### E1. GFF3 loader via `noodles-gff`

`dottir-io::gff3::load(path)` returns `Vec<Feature>` with
`{ source, feature_type, range, strand, attributes }`. Handles gzipped
input. Phase tests against the synthetic fixtures called out in the
plan §6.

### E2. PAF loader via `noodles-paf`

Same shape: `dottir-io::paf::load(path)` returns
`Vec<Hsp> { q_start, q_end, s_start, s_end, q_strand, score }`. Used as
HSP overlay input.

### E3. GUI annotation track panel

Per spec §4.4.1: per-track toggle, colour, line height, label rendering.
Tracks render alongside both axes. For self-comparison, the same track
appears on both. Hover tooltips show GFF attributes. Click jumps the
crosshair to feature start.

### E4. PAF HSP overlay modes

Spec §4.2.3: off / solid line / score-coloured line / greyscale-replace.
Existing greyramp panel gains an "HSP mode" segmented control.

---

## Phase F — Algorithmic completion

### F1. BLASTX three-frame translation

Implement the translation table and the three-frame query encoding.
Then `compute_dotplot` for `mode = Blastx` runs the inner kernel three
times (one per frame) and max-merges. Update tests; remove the
`NotImplemented` branch.

**Files**: new `dottir-core::translation`, `plot.rs` extended.

### F2. Reverse-complement query options

Spec §4.1.10's `-r` / `-v` (reverse-complement the horizontal or
vertical sequence before computation). CLI flags + GUI checkboxes.

### F3. Watson-only / Crick-only `--watson-only` / `--crick-only`

Trivially: when set, skip the corresponding strand pass. Already
implementable on top of the existing `Strand` enum, just needs a
top-level flag rather than enum-via-`Both`.

### F4. Pixelmap goldens against real C dotter

Build the GTK2 dotter binary in `third_party/seqtools/` (in a container
with `libgtk2.0-dev`), run it across the corpus, pin the pixelmaps
under `tests/golden/blastn/` etc. Zero-tolerance compare with the Rust
output.

This is mostly tooling; the Karlin golden infrastructure in
`tests/golden_gen/karlin_ref.c` is the template.

---

## Phase G — Inverted-repeat visualisation & sub-dotter

### G1. Inverted-repeat channel rendering

The data is already there (`PlotConfig::separate_strand_channels` →
`DotPlot::reverse_pixels`). Add a `Greyramp::reverse_colour` knob in the
GUI so forward dots render grey and reverse dots render in a distinct
colour (default magenta, spec §4.4.3). Toggle to merge back into a
single channel.

### G2. Spawn sub-dotter from selection

Right-click drag selects a region; "Spawn sub-dotter" → new dottir
window for that sub-range at a finer `PlotConfig::zoom`. Implementation
options:

1. Re-`eframe::run_native` in a separate viewport (egui supports
   multiple viewports out of the box at 0.29).
2. Fork the process with the sub-range as CLI args.

Option 1 keeps everything in-process and gets state-sharing for free.

---

## Cross-cutting

### Shared config layer

The reviewer flagged that CLI and GUI duplicate `PlotConfig` mapping.
Extract a `dottir-core::ConfigBuilder` (or `serde`-derive on
`PlotConfig` directly) so both frontends + session-load go through one
place. Land alongside C5 (sessions).

### Background-job machinery in the GUI

C2 introduces a `std::thread` worker. Generalise so D's exports
(SVG/PDF rasterisation), E's loaders, and any future heavy work all use
the same pool. Avoid pulling in `tokio`; a simple `std::sync::mpsc` +
one worker thread is enough.

### `egui_kittest` smoke tests

A few headless GUI tests: load FASTA → compute → assert texture has
non-zero histogram; greyramp swap inverts; crosshair arrow keys move it.
Worth adding alongside C work.

### `cargo deny` + `cargo audit` on every PR

The workflow at `.github/workflows/release.yml` runs `cargo deny check`
on tag push. Move it into `ci.yml` so every PR is gated.

---

## Suggested execution order

1. **A1 + A2 + A3 + A4** (3-4 days). Honest memory cap, strict matrix
   parsing, record-aware sequences. Land each as its own commit.
2. **B1-B3** (half a day). Docs synced.
3. **C1 + C2** (2-3 days). Memory limit setting + zoom quality. Most
   user-visible improvement after Phase A.
4. **C3 + C4 + C6** (1-2 days). Breaklines, axes, background compute
   plumbing.
5. **D1 + D2** (2-3 days). SVG + `.dot`. PDF (D3) is optional.
6. **E1 + E2 + E3** (3-4 days). GFF3 / PAF / track panel. Last gating
   feature on §4.4.
7. **F1 + F4** (3 days). BLASTX + pixelmap goldens for confidence.
8. **G1 + G2** (1-2 days). Inverted-repeat colour + sub-dotter.
9. **C5 + D4 + D5 + remainder** as time allows.

After step 6 the project meets the spec §11 "Definition of done for
v1.0" except for the optional WASM build (Phase 8 of the original plan,
still mostly free given §C2 already moves compute off the UI thread).

## Out of scope (this plan)

- A Python `pyo3` binding layer.
- The original Dotter's exon/intron view (spec §5 already excludes
  this).
- Built-in HSP generation / aligner (we still consume external
  PAF/BLAST output — spec §10 question #7 was left as "external only").
- A network/cloud sync layer.
- A WASM GUI build (the GUI threading model in C2 makes this harder, not
  easier; revisit only if there's a concrete user ask).
