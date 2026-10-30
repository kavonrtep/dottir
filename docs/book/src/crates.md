# Crate layout

```text
dottir/
├── crates/
│   ├── dottir-core/   # algorithms; no I/O, no GUI deps
│   ├── dottir-io/     # FASTA, PNG, TOML sidecar, alignment slice
│   ├── dottir-cli/    # `dottir batch` headless binary
│   └── dottir-gui/    # interactive frontend (egui), incl. the genome view
├── docs/
│   ├── book/          # this manual
│   └── adr/           # architecture decision records
└── tests/
    ├── golden/        # pinned Karlin numbers (+ pixelmap goldens to come)
    └── golden_gen/    # C reference harness extracted from dotterKarlin.c
```

The `dottir-core` boundary is load-bearing per `CLAUDE.md`: it stays
I/O-free so it can be embedded in notebooks, other Rust tools, or
future `pyo3` bindings without dragging in PNG / FASTA / GUI deps.

## Module map

### `dottir-core`

| Module | Role |
|--------|------|
| `alphabet` | ASCII→index tables (protein 24, DNA 4); reverse-complement helpers. |
| `matrix` | `ScoreMatrix` + built-in matrices + BLAST-format parser. |
| `karlin` | Karlin/Altschul λ/K/H and the window-size estimate. |
| `score_vec` | Flat `(n+1) × qlen` precomputed score vector with an explicit zero row for unknown subject residues. |
| `sliding` | The ping-pong sum recurrence kernel, chunkable for rayon parallelism. |
| `antidiag` | Anti-diagonal sub-pixel suppression rule. |
| `pixel` | `PixelMap` with checked allocation, max-merge, and element-wise merge. |
| `plot` | Top-level `compute_dotplot` driver (exact engine). |
| `large` | Large engine ([ADR 0007](./adr.md)): `seed` (spaced seeds over 2-bit DNA), `index` (spaced-seed index of the subject), `hits` (query streaming and the two-hit diagonal filter), `score` (exact sliding-window scoring along the surviving diagonals), `depth` (per-base copy number and repeat classes), `chain` (co-linear synteny chains). |
| `genome` | Genome mode ([ADR 0008](./adr.md)): `alignment` (the alignment model), `filter` (identity, length and D-GENIES noise filters), `order` (sorting query records along the reference), `layout` (record placement on each axis), `arrange` (filters, sorting and layouts for one plot), `raster` (deterministic line rasterisation), `records` (the same sort for ordinary multi-record dotplots). |

### `dottir-io`

| Module | Role |
|--------|------|
| `fasta` | Minimal FASTA reader (plain + gzip). |
| `input` | Positional sequence input: FASTA, or GFF3 with an embedded `##FASTA` section. |
| `png_export` | 8-bit greyscale PNG with `tEXt` provenance. |
| `params` | TOML sidecar struct + SHA-256 helper. |
| `alignment` | ±N residue slice around a crosshair coord. |
| `large` | Large-engine glue: repeat floors from soft-masking and GFF3/BED annotations, depth-track and chain (PAF/BEDPE) exports. |
| `paf` | PAF reader (streaming, either orientation, primaries re-derived for `-P` output) and writer. |
| `aligner` | Runs mm2plus/minimap2 for genome mode, with query chunking. |
| `faidx` | `.fai` index builder and region reader: genome mode reads names and lengths only, and drill-down reads regions. |
| `cache` | Per-user cache directory for aligner results and decompressed inputs. |
| `genome_export` | Genome-mode PNG/SVG, association table, sorted FASTA, unaligned records. |

### `dottir-cli`

The `dottir` binary: `batch` (one-shot plots; `large_batch` and
`genome_batch` hold the large-engine and genome-mode paths),
`periodogram` and `find-peaks`.

### `dottir-gui`

The egui/eframe viewer. `compute_worker` runs the exact and large
engines in the background and picks one per view. `genome_view` is the
alignment-based genome comparison with drill-down to a dotplot, and
`record_sort` sorts records in the background.
