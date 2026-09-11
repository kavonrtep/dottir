# The `dottir batch` CLI

```text
dottir batch <QUERY> [SUBJECT] -o <OUT.png> [FLAGS]
```

## Required arguments

| Argument | Description |
|----------|-------------|
| `QUERY` | Sequence path; horizontal axis of the dotplot. |
| `SUBJECT` | Sequence path; vertical axis. **Omit for a self-comparison** — subject is set to query and the mirror post-process is enabled. |
| `-o, --output` | Output PNG path. A `<output>.params.toml` sidecar is written next to it. |

Both inputs may be gzipped (`.gz`); the reader auto-detects from the
file's magic bytes.

### Input formats

A positional input is either:

* a **FASTA** file (`.fa`, `.fasta`, `.fna`, `.faa`, …), or
* a **GFF3** file (`.gff`, `.gff3`) that carries its own sequences after
  a `##FASTA` directive.

The second form needs no separate FASTA — one file supplies both the
sequence and its annotation. Multi-record inputs are concatenated in
file order, as always.

```sh
# Self-comparison straight from an annotated GFF3.
dottir batch elements.gff3 -o elements.png

# Pairwise, two annotated GFF3 files, no FASTA anywhere.
dottir batch family_a.gff3 family_b.gff3 -o ab.png
```

Batch mode draws no annotation overlay, so features found in a GFF3
input are counted in the log and otherwise ignored; the pixelmap is
identical to the one you would get from the equivalent FASTA. The
overlay lives in the [GUI](./gui.md).

A GFF3 with no `##FASTA` section cannot serve as a sequence input and is
rejected:

```text
Error: reading query features.gff3

Caused by:
    features.gff3 is a GFF3 file with no embedded sequences: either
    append a `##FASTA` section to it, or pass the FASTA separately and
    give this file with --gff-query/--gff-subject
```

`dottir periodogram` accepts the same two input forms.

## Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--mode {blastn,blastp,blastx}` | `blastn` | BLAST mode. `blastx` is not yet implemented. |
| `--matrix NAME` | `BLOSUM62` for protein, `DNA+5/-4` for BLASTN | Built-in matrix name — see [Score matrices](./matrices.md). |
| `-W, --window N` | Karlin/Altschul estimate | Sliding window size. |
| `-z, --zoom N` | `1` | Pixels per matrix block. Increase to fit larger inputs. |
| `--pixel-fac N` | `50` | Multiplier in `min(255, score * pixel_fac / W)`. |
| `--strand {forward,reverse,both}` | `both` | BLASTP ignores this. |
| `--triangle {both,upper,lower}` | `both` | Self-comparison mirror mode (only meaningful when SUBJECT is omitted). |
| `--disable-mirror` | off | Skip the self-comparison post-process entirely. |
| `--memory-limit-bytes N` | `512 MiB` | Refuse pixelmaps bigger than this. |
| `--auto-zoom MAX_DIM` | off | Auto-pick `--zoom` so the larger output dim is ≤ this. |
| `--no-sidecar` | off | Skip the `.params.toml` sidecar. |

## Examples

Self-comparison of a small repeat region (just one positional):

```sh
dottir batch contig.fa -o contig.png \
    --triangle both --auto-zoom 4000
```

BLASTP between two protein FASTAs at a fixed window:

```sh
dottir batch query.fa target.fa -o p.png \
    --mode blastp --matrix BLOSUM45 -W 12
```

## Output

The PNG is greyscale 8-bit, row-major, with subject on the vertical
axis. `tEXt` chunks embed `dottir-version`,
`dottir-pixelmap-format-version`, and a `dottir-parameters` summary.

The sidecar (`<output>.params.toml`) contains the full provenance —
see [Reproducibility](./reproducibility.md).
