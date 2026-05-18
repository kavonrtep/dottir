# CLAUDE.md — dottir

Guidance for Claude Code sessions working on this repository.

## What this project is

`dottir` is a Rust reimplementation of **Dotter** (Sonnhammer & Durbin 1995), a
matrix-scored, sliding-window sequence dot-matrix plotter with interactive
sensitivity control. Three goals, in priority order:

1. Faithful preservation of the original algorithms (Karlin/Altschul window
   sizing, score-matrix sliding-window, anti-diagonal suppression, live
   greyramp LUT).
2. Modern cross-platform GUI (egui/eframe) replacing GTK2.
3. New features for repeat-oriented genome analysis (GFF3 overlays, PAF HSP
   overlays, region/alignment export, inverted-repeat highlighting).

Authoritative spec: `docs/dottir_specification.md`. RFC-2119
MUST/SHOULD/MAY language applies — treat MUST items in §4.1 as load-bearing
for scientific correctness.

## Reference implementation

The C source-of-truth lives in `seqtools-4.28/` (not currently cloned here).
Upstream: https://github.com/douglasgscofield/seqtools-4.28

Key entry points (file:line in seqtools-4.28):

| What                                         | Where                            |
|----------------------------------------------|----------------------------------|
| Main pixelmap calc loop                      | `dotterApp/dotplot.c:1476`       |
| Inner sliding-window loop                    | `dotterApp/dotplot.c:1308`       |
| Anti-diagonal suppression rule               | `dotterApp/dotplot.c:1405`       |
| Karlin λ/K/H                                 | `dotterApp/dotterKarlin.c:144`   |
| Window size from λ/K                         | `dotterApp/dotterKarlin.c:343`   |
| Greyramp LUT generation                      | `dotterApp/greyramptool.c:265`   |
| `.dot` binary format / `loadPlot()`          | `dotterApp/dotplot.c:1610`       |
| CLI option parsing (canonical option names)  | `dotterApp/dotterMain.c`         |
| Original CLI docs                            | `dotter.md` at seqtools repo root|

Test corpora to reuse: `seqtools-4.28/test/data/` and `examples/`.

## Planned crate layout (Cargo workspace)

```
dottir/
├── crates/
│   ├── dottir-core/   # algorithms; no I/O, no GUI deps
│   ├── dottir-io/     # FASTA, GFF3, PAF, matrix, .dot, exports
│   ├── dottir-cli/    # batch binary
│   └── dottir-gui/    # interactive binary (egui)
├── docs/
│   ├── book/          # mdBook user manual
│   └── adr/           # architecture decision records
└── tests/
    ├── golden/        # pinned pixelmaps from C dotter
    └── corpora/       # inputs
```

The `dottir-core` boundary is load-bearing: it must stay I/O-free so it can be
embedded (notebooks, other Rust tools, future `pyo3` bindings).

## Algorithmic invariants (do not violate without an ADR)

These are the scientific contract. See §4.1 of the spec for the full list.

1. **Karlin/Altschul window-size estimation** from residue composition + score
   matrix. Clamp to `[3, 50]` by default. `-W` overrides.
2. **Score-matrix sliding-window** (not match/mismatch counts). Default
   BLOSUM62 (protein), identity (DNA).
3. **Flat `scoreVec[residue_type * qlen + col_idx]`** layout — permitted
   deviation from the C row-of-pointers for cache locality.
4. **Sum recurrence** with ping-pong buffers and `delrow=zeros` warm-up:
   `newsum[q] = oldsum[q-1] + scoreVec[sIndex[s]][q] - scoreVec[sIndex[s-W]][q-W]`.
5. **Pixel max-merge** over each `zoom × zoom` block; scale to `u8` via
   `min(255, score * pixelFac / W)`.
6. **Anti-diagonal suppression** (the `dotplot.c:1405` rule):
   - Forward: keep iff `sPosLocal >= qPosLocal`.
   - Reverse: keep iff `(zoom - 1 - sPosLocal) >= qPosLocal`.
   Add an explicit unit test — this is the rule most likely to silently drift.
7. **Determinism**: identical pixelmap across runs and across thread counts.
8. Strand/frame handling per spec §4.1.7–10.

Golden tests against the C dotter pixelmaps have **zero tolerance**. Any
deliberate algorithmic deviation requires bumping the pixelmap format version
and writing an ADR.

## Tech stack (pinned in spec §6.3)

- GUI: `egui` + `eframe` (decided; not iced/gtk4-rs/Tauri).
- FASTA: `noodles-fasta` or `needletail` — open, decide after a 30-min bench.
- GFF3 / PAF: `noodles-gff`, `noodles-paf`.
- Parallelism: `rayon`. CLI: `clap` (derive). Logging: `tracing`.
- Errors: `anyhow` (apps) + `thiserror` (across `dottir-core` public API).
- Images: `image` (PNG), `tiny-skia`/`resvg` (SVG), `printpdf` (PDF).
- Colormaps: `colorgrad`. Compression: `flate2` (transparent gzip).
- Hashing: `sha2`. Config: `serde` + `toml`. Property tests: `proptest`.

## Conventions

- **License**: GPLv3 by default (compatible with upstream Dotter). Revisit
  only if permissive licensing becomes important; if so, Karlin/Altschul code
  must be reimplemented from the 1990 paper + permissive NCBI BLAST sources
  rather than ported verbatim.
- **Naming**: `dottir` is tentative. Check crates.io availability before first
  publish.
- **`panic!`**: never on user input. Use `thiserror` for typed errors crossing
  `dottir-core`'s public boundary; `anyhow` inside binaries.
- **Comments**: only when the *why* is non-obvious. Don't restate the code.
- **Determinism**: no `HashMap` iteration order in hot paths producing
  pixelmaps; use `BTreeMap` or sorted vectors where order matters.

## Build / test

The Rust toolchain (1.95.0, matching `rust-toolchain.toml`) lives in a conda
env inside the hermit sandbox — there is no rustup on the host. Activate it
first, then run cargo as normal:

```bash
source /opt/conda/etc/profile.d/conda.sh
conda activate /envs/conda/envs/rust

cargo build --release
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --check
```

`rust-toolchain.toml` is informational only here — without rustup it cannot
auto-install. If the channel pin in that file changes, recreate the env:
`mamba install -p /envs/conda/envs/rust -c conda-forge rust=<version>`.

CI matrix: Linux x86_64 + Windows x86_64. macOS optional. WASM as a separate
"does it compile" job.

## Sandbox environment (`hermit/`)

The `hermit/` subdirectory is **not part of the dottir project**. It is a
self-contained Singularity/Apptainer sandbox for running Claude Code / Codex
CLI agents with read-only data mounts and persistent conda/pip/npm state
under `/envs`. See `hermit/README.md` for layout and `hermit/CLAUDE.md` /
`hermit/AGENTS.md` for the agent-facing context that loads inside the
container.

Practical implications for dottir work performed inside this sandbox:

- The Rust toolchain is provisioned via `mamba create -p /envs/conda/envs/rust`,
  not rustup. Activate it before running cargo (see Build / test above).
- Bare `conda install` / `mamba install` / `pip install` are blocked by
  hermit hooks. For single tools use `htool <name>`; for multi-package envs
  use `mamba create -p /envs/conda/envs/<name> ...` directly.
- The container mounts data paths at their original host locations — paths
  in logs, scripts, and output match between inside and outside the
  container. Do not introduce `/input`→`/output` translation.

## Out of scope (do not propose without checking spec §5)

- Exon/intron view from original Dotter.
- Blixem integration. `libpfetch` (Sanger-internal).
- Built-in HSP generation / aligner (we consume external PAF/BLAST output).
- Original bespoke feature-file format (replaced by GFF3).
- Network / cloud features.

## When making changes

- Algorithm changes in `dottir-core` → bump pixelmap format version + add an
  ADR in `docs/adr/` + regenerate goldens with documented justification.
- New deps → justify in PR description; prefer crates already in §6.3.
- Public API changes in `dottir-core` → update `cargo doc` examples and the
  mdBook user manual under `docs/book/`.
- Anything touching the anti-diagonal rule, pixel max-merge, or sum recurrence
  needs explicit golden-test coverage.

## Documents in this repo

- `docs/dottir_specification.md` — the contract. Read §4 (requirements)
  before starting any non-trivial task.
- `docs/IMPLEMENTATION_PLAN.md` — phased task breakdown.
- `docs/IMPROVEMENTS_PLAN.md` — follow-up plan addressing the review.
- `docs/REVIEW.md` — external code review.
- `docs/CHANGELOG.md` — per-release notes.
- `docs/adr/` — decisions that deviate from spec defaults or resolve §10 open
  questions. Create one when answering: name/license/FASTA-lib/pixelmap
  container/alignment-on-the-fly choices.
