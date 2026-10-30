# Architecture decisions

The full ADRs live in `docs/adr/` (MADR format). Index:

| # | Title | Status |
|---|-------|--------|
| [0001](https://github.com/petr/dottir/blob/main/docs/adr/0001-egui-frontend.md) | egui frontend | Accepted |
| [0002](https://github.com/petr/dottir/blob/main/docs/adr/0002-license-gplv3.md) | License = GPL-3.0-or-later | Accepted |
| [0003](https://github.com/petr/dottir/blob/main/docs/adr/0003-gui-msrv.md) | Defer GUI to MSRV bump | Superseded |
| [0004](https://github.com/petr/dottir/blob/main/docs/adr/0004-defer-gff3-paf.md) | Defer GFF3/PAF to MSRV bump | Accepted |
| [0005](https://github.com/petr/dottir/blob/main/docs/adr/0005-annotation-overlay.md) | Annotation overlay as full-span bands | Accepted |
| [0006](https://github.com/petr/dottir/blob/main/docs/adr/0006-gff3-with-embedded-sequences.md) | GFF3 with embedded sequences as input | Accepted |
| [0007](https://github.com/petr/dottir/blob/main/docs/adr/0007-large-sequence-mode.md) | Seed-based large-sequence mode | Accepted |
| [0008](https://github.com/petr/dottir/blob/main/docs/adr/0008-genome-mode.md) | Genome mode: whole-assembly dotplots from external alignments | Accepted |
| [0009](https://github.com/petr/dottir/blob/main/docs/adr/0009-blastx-karlin-three-frames.md) | BLASTX Karlin estimate over three frames | Accepted |
| [0010](https://github.com/petr/dottir/blob/main/docs/adr/0010-c-dotter-pixelmap-goldens.md) | C dotter pixelmap goldens; matching C where they differed | Accepted |

When making a change that warrants a new ADR (see `docs/adr/README.md`
for criteria), add it under `docs/adr/NNNN-slug.md`, then append a row
here.
