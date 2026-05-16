# Score matrices

dottir ships eight protein matrices plus a DNA identity matrix as
built-ins:

| Name | Source | Notes |
|------|--------|-------|
| `BLOSUM62` | Inline from `dotter.c:354` (1993-08-09 table) | Used by default for BLASTP. Differs from current NCBI BLOSUM62 in the ambiguous-residue rows (X / B / Z); the C-dotter version is kept verbatim for golden-test reproducibility. |
| `BLOSUM45` / `50` / `80` / `90` | Vendored from NCBI BLAST data files | `crates/dottir-core/data/` |
| `PAM30` / `70` / `250` | Vendored from NCBI BLAST data files | `crates/dottir-core/data/` |
| `DNA+5/-4` | Inline (matches `DNAmatrix()` at `dotter.c:2800`) | `+5` on the diagonal, `-4` off, over the 4-letter ACGT alphabet. |

Look up by name from code:

```rust
use dottir_core::ScoreMatrix;
let m = ScoreMatrix::by_name("BLOSUM80").unwrap();
```

Or load a custom matrix in standard NCBI BLAST text format:

```rust
let text = std::fs::read_to_string("my_matrix.txt")?;
let m = ScoreMatrix::parse_blast_format("MyMatrix", &text)?;
```

The parser tolerates extra columns whose header letters aren't in the
24-residue dotter alphabet (e.g. NCBI's `J` for Leu/Ile), silently
dropping them.
