# ADR 0004 — Defer GFF3 and PAF loaders to a future MSRV bump

* Status: Accepted
* Date: 2026-05-16
* Deciders: petr
* Supersedes: nothing; complements ADR 0003.

## Context

Phase 6 of `IMPLEMENTATION_PLAN.md` calls for GFF3 and PAF loaders using
`noodles-gff` and `noodles-paf`. Both pull a similar transitive-dep
chain to `eframe` (winit→toml_edit→indexmap-2.14 family, or the noodles
team's own `bstr`/`csv`-style stack) which triggers the same edition2024
floor: most current `noodles-*` releases target Rust 1.85+.

The Phase 6 GUI features (annotation track panel, rubber-band selection,
PAF HSP overlay) also depend on the deferred Phase 5 egui runtime
(ADR 0003).

## Decision

Defer GFF3 and PAF loaders to a future release that ships against a
post-1.85 MSRV. The Phase 6 piece that does NOT depend on noodles or
egui — the **alignment-export helper** (slice ±N residues around a
crosshair, emit FASTA pair / Stockholm / plain text) — lands now in
`dottir-io::alignment`, since it is useful from the CLI today.

## Consequences

* `dottir-io::alignment` provides `slice_pair` and the three serialisers
  expected by the eventual GUI Phase 6 wiring.
* GFF3/PAF readers are NOT in this release. The dottir-cli has no
  flags for them.
* When ADR 0003's MSRV bump lands, ADR 0004 is re-evaluated: at that
  point `noodles-gff = "0.x"` and `noodles-paf = "0.y"` can almost
  certainly be added directly.

## Alternatives considered

1. Hand-roll a tiny GFF3 parser. GFF3 is line-oriented and a minimal
   subset (attributes, strand, source filter) is ~100 LOC. Tempting but
   accepts technical-debt vs. waiting a release for the well-tested
   noodles parser. Rejected for now.
2. Vendor noodles-gff/paf source under `third_party/` and bypass the
   transitive-dep machinery. Higher maintenance cost than the
   MSRV bump; rejected.
