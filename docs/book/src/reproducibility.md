# Reproducibility

Every `dottir batch` invocation writes two artifacts:

1. **`<output>.png`** — the dotplot itself. The PNG embeds three `tEXt`
   chunks:
   - `dottir-version` — the crate version (e.g. `0.1.0`).
   - `dottir-pixelmap-format-version` — bumps whenever the algorithmic
     contract changes such that goldens need regeneration.
   - `dottir-parameters` — a semicolon-delimited summary like
     `mode=Blastn;matrix=DNA+5/-4;W=25;zoom=1;pixel_fac=50;strand=Both`.
2. **`<output>.png.params.toml`** — the full sidecar (spec §4.4.7).
   Suppress with `--no-sidecar`.

## Sidecar structure

```toml
[dottir]
version = "0.1.0"
pixelmap_format_version = 0

[query]
path = "/path/to/q.fa"
sha256 = "..."           # SHA-256 of the input file as bytes
size_bytes = 12345
n_records = 1
total_residues = 5000

[subject]
path = "/path/to/s.fa"
sha256 = "..."
size_bytes = 12345
n_records = 1
total_residues = 5000

[plot]
mode = "Blastn"
matrix = "DNA+5/-4"
strand = "Both"
window_size = 25         # resolved value (might be the Karlin estimate)
zoom = 1
pixel_fac = 50
self_comparison = false
width = 5000
height = 5000

[plot.karlin]            # present if window was estimated (not -W override)
lambda = 0.191529273986816
k = 0.173345809507352
h = 0.356723794166828
expected_residue_score = 1.862502722123737
expected_msp_score = 38.938557203279061
predicted_msp_length = 21

[host]
hostname = "..."
os = "linux"
timestamp_utc = "2026-05-16T22:33:00Z"
```

## Why this matters

Recreating a published dotplot from scratch requires knowing
*exactly* the window size, zoom, pixel factor, and the SHA-256 of
each input. Without this, a reader can't tell whether they reproduced
the same plot. The sidecar captures all of it in a parseable format.

## Determinism

Given the same inputs and the same sidecar parameters, dottir is
byte-deterministic across runs and across thread counts. This is
spec §4.1.11 and is verified by golden / parallel-determinism tests.
