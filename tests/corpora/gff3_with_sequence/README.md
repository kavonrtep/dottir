# `gff3_with_sequence` corpus

Fixtures for GFF3 input that carries its own sequences after a `##FASTA`
directive, i.e. a self-contained annotated sequence that needs no separate
FASTA.

Regenerate with:

```bash
python3 scripts/make_gff3_with_sequence_fixtures.py
```

## Files

| File                          | Records | Content                                                                 |
|-------------------------------|---------|-------------------------------------------------------------------------|
| `pair_a.gff3`                 | 2       | Synthetic. 600 bp records: 150 bp flank, 300 bp shared element, 150 bp flank. Two features per record. |
| `pair_b.gff3`                 | 2       | Same layout; the element is `pair_a`'s at ~8 % divergence, flanks independent. |
| `no_sequence.gff3`            | 0       | Negative control — valid GFF3, no `##FASTA` section.                     |
| `tir_simple.with_seq.gff3`    | 3       | `annotation_overlay/tir_simple.{gff3,fasta}` merged.                     |
| `tir_elements.with_seq.gff3`  | 10      | `annotation_overlay/tir_elements.{gff3,fasta}` merged.                   |
| `ltr_angela.with_seq.gff3.gz` | 3       | `annotation_overlay/ltr_angela.{gff3,fasta}` merged, gzipped.            |

The `*.with_seq.*` files are byte-preserving merges of their
`annotation_overlay/` sources, so they must load to exactly the same
residues and features as the separate pair. `crates/dottir-io/tests/gff3_with_sequence.rs`
asserts that equivalence.

## Intended uses

- Self-comparison from one file: `dottir batch pair_a.gff3 -o out.png`.
- Pairwise from two files: `dottir batch pair_a.gff3 pair_b.gff3 -o out.png`
  — the shared element gives an off-diagonal block in every record pair.
- GUI overlay: `dottir-gui pair_a.gff3` binds `SharedElement` / `Flank`
  bands to both axes without a `--gff` flag.
- Error path: `no_sequence.gff3` as a positional input must fail with the
  "no embedded sequences" message, not a parse error.
