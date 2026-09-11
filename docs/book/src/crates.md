# Crate layout

```text
dottir/
├── crates/
│   ├── dottir-core/   # algorithms; no I/O, no GUI deps
│   ├── dottir-io/     # FASTA, PNG, TOML sidecar, alignment slice
│   ├── dottir-cli/    # `dottir batch` headless binary
│   └── dottir-gui/    # interactive frontend (egui deferred per ADR 0003)
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
| `plot` | Top-level `compute_dotplot` driver. |

### `dottir-io`

| Module | Role |
|--------|------|
| `fasta` | Minimal FASTA reader (plain + gzip). |
| `input` | Positional sequence input: FASTA, or GFF3 with an embedded `##FASTA` section. |
| `png_export` | 8-bit greyscale PNG with `tEXt` provenance. |
| `params` | TOML sidecar struct + SHA-256 helper. |
| `alignment` | ±N residue slice around a crosshair coord. |

### `dottir-cli`

The `dottir batch` binary, wrapping the core + io into a one-shot CLI.

### `dottir-gui`

Currently a headless FASTA-to-PNG fallback. egui runtime pending the
MSRV bump (ADR 0003).
