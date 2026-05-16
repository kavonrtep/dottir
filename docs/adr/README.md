# Architecture Decision Records

Decisions of the form "we deviated from the spec because…" or "we chose X
over Y because…" live here. Format: [MADR](https://adr.github.io/madr/).

| #    | Title                       | Status   |
|------|-----------------------------|----------|
| 0001 | egui frontend               | Accepted |
| 0002 | License = GPLv3-or-later    | Accepted |
| 0003 | Defer GUI to MSRV bump      | Accepted |

Add a new ADR whenever the change is one of:

* a deliberate deviation from `dottir_specification.md` requirements,
* a resolution of a §10 open question,
* a tool / dependency choice with non-trivial scope (e.g. FASTA library),
* a change to the pixelmap algorithmic contract (also bump
  `PIXELMAP_FORMAT_VERSION` and regenerate goldens).
