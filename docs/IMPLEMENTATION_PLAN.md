# dottir — Implementation Plan

A concrete, task-level breakdown expanding §9 of `dottir_specification.md`.
Estimates are rough person-weeks for one part-time developer; double for full
polish. Phases 0–5 produce a usable GUI; Phase 6 is where dottir surpasses the
original; Phases 7–9 are polish and distribution.

Each phase lists deliverables, key tasks, files touched, dependencies, and
acceptance criteria. Spec section references are to `dottir_specification.md`.

---

## Phase −1 — Decisions and scaffolding (prerequisites)

Before any code, resolve enough §10 open questions to unblock work. The rest
can be deferred until first relevant phase.

**Decisions needed now**

1. **Project / crate name.** Verify `dottir-core`, `dottir-io`, `dottir-cli`,
   `dottir-gui` are free on crates.io. Reserve placeholders, or commit to
   final names. Spec §10.1.
2. **License.** Spec default: GPLv3. Confirm with developer. Spec §10.2.
3. **Reference C source location.** Decide whether to vendor a copy of
   `seqtools-4.28` (read-only, used as a side-by-side comparison harness for
   golden-test generation) or to require it as an external tool the developer
   has on `$PATH`. Recommendation: vendor a stripped subset
   (`dotterApp/*`, `examples/`, `test/data/`) under `third_party/seqtools/`
   with the upstream GPLv3 LICENSE preserved.

**Decisions that can be deferred**

- FASTA library (Phase 1). PAF library (Phase 6). Modern pixelmap container
  format (Phase 4). Alignment-on-the-fly approach (Phase 6). WASM (Phase 8).
  pyo3 bindings (post-v1).

**Scaffolding tasks**

- Create Cargo workspace `Cargo.toml` and the four empty crates per §6.1.
- Add `rust-toolchain.toml` pinning a stable Rust (latest stable at start).
- `.gitignore`, `LICENSE` (GPLv3 by default), `NOTICE` crediting Sonnhammer,
  Durbin, Barson, Scofield, NCBI.
- `rustfmt.toml` (defaults), `clippy.toml`, `deny.toml` for `cargo-deny`,
  `.cargo/config.toml` if needed.
- GitHub Actions skeleton: `fmt`, `clippy --deny warnings`, `test` on Linux
  and Windows. Cache `~/.cargo` and `target/`.
- First ADR: `docs/adr/0001-egui-frontend.md` recording the GUI decision.
- First ADR: `docs/adr/0002-license-gplv3.md` if license is confirmed.

**Acceptance**: `cargo check --workspace` is green; CI runs on a no-op PR.

---

## Phase 0 — Karlin/Altschul port (≈ 1 wk)

**Deliverable**: `dottir-core::karlin` module producing λ, K, H and derived
window sizes matching the C output to machine precision for the test corpus.

**Tasks**

- Port `dotterKarlin.c:144` (`karlin()`) to `crates/dottir-core/src/karlin.rs`.
  Preserve numeric algorithm exactly; replace pointer arithmetic with slices.
- Port `dotterKarlin.c:343` (`winsizeFromlambdak()`) including the round-up,
  the 100×100 nominal MSP-length convention, and the `[3, 50]` clamp (with
  the clamp bounds taken from a `KarlinConfig` so they remain configurable).
- Build a residue-composition utility: count residues in a `&[u8]` with the
  matrix's alphabet and produce the probability vector the C code needs.
- Encode the built-in matrices (BLOSUM62, BLOSUM50, BLOSUM45, BLOSUM80,
  BLOSUM90, PAM30, PAM70, PAM250, DNA identity). Source: NCBI matrix files
  shipped with BLAST. Add a `Matrix::parse_blast_format` for user matrices.
- Public API per spec §6.2:
  ```rust
  pub fn karlin_window_size(
      matrix: &ScoreMatrix, query: &[u8], subject: &[u8], mode: BlastMode,
  ) -> Result<KarlinResult, DottirError>;
  ```
- `KarlinResult` exposes `{ lambda, k, h, predicted_msp_length, window_size }`.

**Tests**

- Unit: each built-in matrix loads and round-trips to BLAST format.
- Unit: residue composition on a small fixture (`AAACGT` etc.).
- **Golden**: for 6–8 input pairs (protein and DNA), pin the exact
  `(lambda, k, h, window_size)` produced by C dotter. Compare with
  `f64::to_bits` equality, not `==` with epsilon (spec §4.1.11 determinism).
  These goldens live in `tests/golden/karlin/*.json`.

**Risks**

- The Karlin/Altschul iterative solver in C uses `double` and specific
  termination criteria. To get byte-identical output, mirror the loop
  structure precisely rather than reformulating with a generic root-finder.

**Exit criteria**: golden tests pass; `cargo doc` page exists for
`karlin_window_size` with a worked example.

---

## Phase 1 — Single-strand BLASTN inner loop (≈ 1 wk)

**Deliverable**: `compute_dotplot()` producing pixelmaps byte-identical to C
dotter for forward-strand BLASTN on the chr4 self-comparison corpus.

**Tasks**

- `crates/dottir-core/src/score_vec.rs`: flat
  `scoreVec[residue_type * qlen + col_idx]` builder. Document the layout
  deviation from C in a module-level comment.
- `crates/dottir-core/src/sliding.rs`: the ping-pong sum recurrence with the
  `delrow=zeros` warm-up region. Single-threaded for now.
- `crates/dottir-core/src/pixel.rs`: max-merge into a `Vec<u8>` and the
  `min(255, score * pixelFac / W)` scale step.
- `crates/dottir-core/src/antidiag.rs`: the suppression rule from
  `dotplot.c:1405`. Forward-strand only this phase. **Unit test it
  separately** (spec §4.1.6 mandates this).
- `crates/dottir-core/src/lib.rs`: public `compute_dotplot()` entry point per
  §6.2; for Phase 1, only supports `BlastMode::Blastn` + `Strand::Forward` +
  `self_comparison: false`.
- `crates/dottir-io/src/fasta.rs`: minimal FASTA reader behind a trait so the
  noodles-vs-needletail decision can be made later. Start with whichever you
  picked; benchmark in Phase 3.

**Memory check**

- Honour `memory_limit_bytes` from `PlotConfig`. Refuse allocation and return
  `DottirError::OutOfMemory` with a suggested zoom factor.

**Tests**

- Unit: the anti-diagonal rule on a hand-computed 4×4 sub-pixel grid.
- Unit: sum recurrence on a 12×12 toy input; compare to a naive O(qlen·slen·W)
  reference computed in-test.
- **Golden**: chr4 self-comparison forward strand, zoom 250, W from
  Karlin/Altschul. Pin the `DotPlot.pixels` `Vec<u8>` as a `.bin` under
  `tests/golden/blastn/`. Zero tolerance per spec §8.1.

**Acceptance**: `cargo test --workspace` green; on the corpus, `sha256` of
the pixelmap matches the pinned C output.

---

## Phase 2 — Full BLAST modes, dual strand, self-comparison (≈ 1–2 wk)

**Deliverable**: BLASTN both-strand, BLASTP, BLASTX, self-comparison with
mirror, reverse-complement options.

**Tasks**

- Reverse-strand pass with the second anti-diagonal rule branch.
- Max-merge forward + reverse passes into a single pixelmap (default), with
  an optional separate-channel mode used later by Phase 7's inverted-repeat
  feature. Plumb `DotPlot::{forward_pixels, reverse_pixels}` (spec §6.2).
- BLASTP: single pass over protein alphabet using BLOSUM62 by default.
- BLASTX: three reading frames on the query (or subject) with max-merge.
- Self-comparison: compute upper triangle only, mirror across diagonal.
  Honour `triangle = Upper | Lower` and `disable_mirror`.
- Reverse-complement options `-r`, `-v` from spec §4.1.10.
- Watson-only / Crick-only for DNA (skip the opposite-strand pass).

**Tests**

- Add property tests (`proptest`):
  - `pixelmap(reverse(q), reverse(s)) == reverse_xy(pixelmap(q, s))`.
  - Self-comparison symmetry across the main diagonal after mirror.
- Golden tests for: BLASTN both strands, BLASTP on `Q9H8G1.fasta` ×
  `DA730641.fasta`, BLASTX on a small mRNA × protein pair, self-comparison
  on chr4 (full triangle + mirror).

**Acceptance**: all goldens pass byte-identically; property tests pass with
`PROPTEST_CASES=2048`.

---

## Phase 3 — Parallelism and perf (≈ 1 wk)

**Deliverable**: same outputs, ≤ 2 s on a 1 Mb × 1 Mb BLASTN dotplot at
zoom 250 on a 2024-era 8-core desktop (spec §4.5.1).

**Tasks**

- Partition the subject sequence into chunks along the `s` axis. Each rayon
  worker computes its slice into a thread-local pixelmap; final merge with
  `max` is associative so order-independent → preserves determinism (spec
  §4.1.11).
- Decide on per-thread allocation strategy: thread-local ping-pong buffers
  vs. arena. Document in an ADR.
- Investigate SIMD for the sum recurrence (`std::simd` is stabilising;
  alternatively `wide` crate). Only land if it doesn't break determinism;
  guard behind a feature flag if it changes numeric ordering.
- Benchmark `noodles-fasta` vs `needletail` on a 1 Mb FASTA cold-cache read
  (criterion). Pick one; record in `docs/adr/0003-fasta-library.md`.

**Tests**

- Re-run all Phase 0–2 goldens with `rayon::ThreadPoolBuilder` set to 1, 2,
  4, 8 threads. Output MUST be byte-identical across thread counts.
- `criterion` benches in `crates/dottir-core/benches/dotplot.rs` for the
  target 1 Mb × 1 Mb workload.

**Acceptance**: benches meet §4.5.1 budget; determinism property test runs
on a 100-case proptest matrix of thread counts.

---

## Phase 4 — CLI batch mode (≈ 1–2 wk)

**Deliverable**: `dottir batch …` produces PNG/SVG/PDF + params sidecar and
reads/writes `.dot`.

**Tasks**

- `crates/dottir-cli/src/main.rs`: `clap` derive for the subcommand. Mirror
  the original option names where reasonable (consult `dotter.md`).
- `crates/dottir-io/src/dot.rs`: read the `.dot` binary format
  (`dotplot.c:1610`). Write support optional (decide based on whether anyone
  consumes `.dot` downstream of dottir).
- `crates/dottir-io/src/png.rs`: PNG export. Embed parameters in `tEXt`
  chunks (input paths, SHA-256 hashes, parameters, dottir version, git SHA
  via `build-info` crate or `vergen`).
- `crates/dottir-io/src/svg.rs`: SVG export via `tiny-skia` + `resvg`.
  Metadata block in `<metadata>` element.
- `crates/dottir-io/src/pdf.rs`: PDF via `printpdf`, or rasterise to a
  high-DPI PNG and embed.
- `crates/dottir-io/src/params.rs`: emit `<output>.params.toml` sidecar
  (spec §4.4.7). On by default in CLI batch mode.
- `--auto-zoom <max_dim>` picks `zoomFactor` so the largest output dimension
  fits (spec §4.4.8). Avoids OOMs on large inputs.

**Tests**

- Integration: run `dottir batch` over the corpus; PNG content sha256 pinned
  in `tests/golden/cli/`.
- `.dot` round-trip: read a C-generated `.dot`, render to PNG, compare to a
  golden PNG.
- Params sidecar TOML schema validated by a JSON Schema (or just by
  `toml::from_str` into a typed struct).

**Acceptance**: a headless user can run dottir end-to-end on the corpus
without the GUI crate ever being built. `cargo build -p dottir-cli` is the
only binary needed.

---

## Phase 5 — GUI MVP (≈ 2–3 wk)

**Deliverable**: `dottir-gui` binary with pan/zoom, greyramp, crosshair,
alignment view. Spec §4.3.

**Tasks**

- `crates/dottir-gui/src/main.rs`: `eframe` app skeleton.
- `app.rs`: top-level state (loaded sequences, current `PlotConfig`, current
  `DotPlot`, view transform, crosshair position).
- `canvas.rs`: dotplot canvas as a textured quad with `egui::TextureHandle`.
  Pan with middle-drag, zoom with scroll wheel centred on cursor.
  **Display zoom is separate from computation zoom** (spec §4.3.2) — pan/zoom
  is a viewport transform; recomputing happens only when `zoomFactor` changes
  in settings.
- `greyramp.rs`: black/white sliders + swap/reset/spinboxes. LUT regenerated
  on slider change; pixelmap untouched. Texture rebuilt by applying LUT
  (cheap: `pixels.iter().map(|&p| lut[p as usize])`).
- `crosshair.rs`: keyboard nudge by 1/10/100 (Shift/Ctrl). Status bar shows
  synchronised coordinates in both sequences.
- `axes.rs`: anti-aliased scale lines and ticks. Honour `--suppress-scale`,
  `--labels-off`, `--labels-size`.
- `alignment_view.rs`: dock panel showing ±60 residues around the crosshair
  with mismatch highlighting. Synchronises live.
- `settings.rs`: matrix/window-size/pixel-factor/memory-limit/colormap/font
  dialog. Triggers recompute when params change.
- `breaklines.rs`: render breaklines for multi-record FASTA. Hover shows the
  sequence name on either side (spec §4.4.6 enhancement).
- File ops: drag-and-drop FASTA accepted. `File → Open`, `File → Save .dot`,
  `File → Save session (TOML)`.
- Spawn sub-dotter from a rectangular selection.
- High-DPI awareness. Optional dark theme.

**Tests**

- `egui_kittest` smoke tests where practical: load → render → assert texture
  has non-zero histogram; greyramp swap inverts; crosshair arrow keys move it.
- Manual test plan in `docs/manual_test_plan.md`.

**Acceptance**: developer can run `cargo run -p dottir-gui` and reproduce the
core interactive behaviour of the original Dotter on the corpus.

---

## Phase 6 — GFF3, PAF, region & alignment export (≈ 2 wk)

**Deliverable**: the features that justify rewriting Dotter (spec §4.4).

**Tasks**

- `crates/dottir-io/src/gff.rs`: GFF3 loader via `noodles-gff`. Track =
  `{ name, source_filter, color, height, on, render_labels }`.
  Gzipped accepted (transparent via `flate2`).
- `crates/dottir-gui/src/tracks.rs`: annotation track panel. Each track
  renders alongside both axes (for self-comparison, on both). Tooltip on
  hover; click jumps the crosshair.
- `crates/dottir-io/src/paf.rs`: PAF loader via `noodles-paf` (or hand-rolled
  if noodles-paf lacks something). Render in existing HSP overlay modes
  (off / line / score-coloured / greyscale-replace).
- `crates/dottir-gui/src/selection.rs`: rubber-band selection of a
  sub-rectangle. Hooks: export PNG/SVG of selection; spawn sub-dotter with
  finer zoom; copy BED-like coordinates to clipboard.
- `crates/dottir-io/src/alignment_export.rs`: at the crosshair, slice ±N
  residues from each sequence and write as FASTA pair / Stockholm / text.
  Optionally run a built-in Smith-Waterman from the `bio` crate on windows
  ≤ 1 kb (spec §10.5 recommendation). Otherwise shell out to
  `mafft`/`muscle` if found on `$PATH`.

**Tests**

- GFF3 parser fixture: synthetic gff3 with edge cases (CDS on minus strand,
  multiline features, missing optional cols).
- PAF parser fixture from a small minimap2 run on the corpus.
- Snapshot of GUI overlay: render at a fixed view, hash texture, pin.

**Acceptance**: a user with a fresh assembly + repeat-annotation GFF3 +
minimap2 self-alignment PAF can do meaningful repeat exploration in dottir
that they could not do in the original.

---

## Phase 7 — Inverted-repeat channel, breaklines, polish, docs (≈ 1–2 wk)

**Deliverable**: feature parity + new repeat features fully polished; user
manual exists.

**Tasks**

- Inverted-repeat highlighting: surface `DotPlot::reverse_pixels` in the GUI
  with a separate colour (default magenta); toggle to merge into a single
  channel (spec §4.4.3).
- Breakline hover labels (spec §4.4.6 — possibly already in Phase 5; defer
  here if not).
- Reversed scale axes (`hozScaleRev`, `vertScaleRev`, `negateCoords`) —
  spec §4.2.8.
- `docs/book/`: mdBook user manual covering CLI flags, GUI workflows, GFF3
  track config, params sidecar format.
- `cargo doc` examples for every public item in `dottir-core`.

**Acceptance**: a new user can install dottir and reproduce a published
analysis from the manual.

---

## Phase 8 — WASM build (optional, ≈ 1 wk)

**Deliverable**: `cargo build --target wasm32-unknown-unknown -p dottir-gui`
produces a working browser viewer.

**Tasks**

- Gate I/O behind a `cfg(not(target_arch = "wasm32"))` shim that uses the
  browser file picker on WASM.
- Disable rayon on WASM (no threads in stable wasm32 without coi headers);
  fall back to single-threaded compute.
- Wire up `trunk` build config; add a CI job that just verifies it compiles.

**Acceptance**: a GitHub Pages preview of a static dottir build that accepts
a small FASTA pair via file picker and renders a pixelmap.

---

## Phase 9 — Release engineering (≈ 1 wk)

**Deliverable**: tagged v1.0 with Linux and Windows binaries on GitHub
Releases.

**Tasks**

- GitHub Actions release workflow: build matrix for Linux x86_64 (musl
  static binary) and Windows x86_64 (`cross` or windows-latest runner).
  Upload artifacts on tag push.
- `cargo-deny` and `cargo-audit` gates in CI.
- `CHANGELOG.md` from v0.1 onward.
- Final naming check on crates.io; publish `dottir-core` (the only crate
  that benefits from being on crates.io for library users) if license
  allows.

**Acceptance**: someone with no Rust toolchain can download a Linux or
Windows binary from Releases and run it.

---

## Cross-cutting workstreams

These run alongside the phases, not as separate phases:

**Goldens hygiene**
Goldens live under `tests/golden/{karlin,blastn,blastp,blastx,self,cli}/`.
Each golden has a `<name>.cmd` file recording the exact C dotter invocation
used to generate it. Regeneration requires bumping `PIXELMAP_FORMAT_VERSION`
in `dottir-core` and an ADR.

**ADRs**
Decisions of the form "we deviated from the spec because…" or "we chose X
over Y because…" go in `docs/adr/NNNN-<slug>.md`. Use the
[MADR](https://adr.github.io/madr/) template. Likely first ADRs:
1. egui frontend (Phase −1).
2. License = GPLv3 (Phase −1).
3. FASTA lib choice (Phase 3).
4. Modern pixelmap container format (Phase 4 — `.npz` vs `.zarr` vs HDF5).
5. Alignment-on-the-fly: built-in vs shell-out (Phase 6).

**Documentation**
- `cargo doc` is the API reference. Every public item in `dottir-core` needs
  a doc comment with an example by Phase 7.
- mdBook user manual under `docs/book/` is the end-user reference, built in
  CI and deployed to GitHub Pages on `main`.

**CI**
From Phase −1: `fmt`, `clippy -D warnings`, `test`, `deny`. From Phase 4:
golden integration tests. From Phase 5: WASM compile-only job (if
prioritised earlier than Phase 8).

---

## Open questions tracked here

Mirror of spec §10 with the resolved-at column:

| #  | Question                                  | Resolve by  | Default               |
|----|-------------------------------------------|-------------|-----------------------|
| 1  | Final project/crate names                 | Phase −1    | `dottir-*` (tentative)|
| 2  | License                                   | Phase −1    | GPLv3                 |
| 3  | FASTA library                             | Phase 3     | benchmark + ADR       |
| 4  | Modern pixelmap container                 | Phase 4     | `.npz` (recommended)  |
| 5  | Alignment-on-the-fly                      | Phase 6     | built-in for ≤ 1 kb   |
| 6  | GUI framework                             | Phase −1    | egui (decided)        |
| 7  | HSP source                                | Phase 6     | external only         |
| 8  | WASM in v1?                               | Phase 8     | post-v1               |
| 9  | Python bindings (pyo3)                    | post-v1     | defer                 |
| 10 | Style/CI guide                            | Phase −1    | rustfmt + clippy -D   |

---

## Definition of done for v1.0

- All Phase 4.1 MUST items implemented and golden-tested.
- GUI MVP (Phase 5) usable end-to-end on the corpus.
- Repeat features (Phase 6) functional with at least one documented workflow
  in the manual.
- Performance budget (spec §4.5.1) met.
- Linux + Windows binaries on Releases.
- User manual published.
- All §10 questions either resolved (with ADR) or explicitly deferred to
  post-v1.
