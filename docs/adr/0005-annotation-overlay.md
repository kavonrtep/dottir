# ADR 0005 — GFF3/BED annotation overlay as full-span bands

* Status: Accepted
* Date: 2026-06-29
* Deciders: petr
* Spec reference: §5 (planned feature), §6.3 (noodles-gff pin)
* Supersedes: revives the GFF3-loader half of ADR 0004 (the MSRV blocker is gone).

## Context

ADR 0004 deferred the GFF3/PAF loaders because the contemporary `noodles-*`
releases needed Rust ≥ 1.85 while the GUI MSRV was pinned lower. The toolchain
now ships 1.95 (conda `rust` env; `rust-version = "1.85"` in the workspace), so
that blocker no longer applies.

We want to overlay annotation intervals on the dotplot for both self-comparison
(one sequence, one annotation file applied to both axes) and pairwise comparison
(one annotation file per axis). A GFF3/BED feature is a 1-D interval on a single
sequence; the dotplot is 2-D (query = X, subject = Y), so the design must define
how a 1-D interval becomes 2-D geometry.

## Decision

1. **Geometry — full-span bands.** A query feature `[a,b)` draws as a vertical
   band spanning the full plot height; a subject feature draws as a horizontal
   band spanning the full width. Bands are drawn *over* the pixelmap texture with
   per-pixel alpha, so dots remain visible. Where a query band and a subject band
   cross, alpha compositing self-darkens the intersection rectangle for free; in
   self-mode each feature's `[a,b]×[a,b]` diagonal square is exactly such a
   crossing.

2. **Formats — GFF3 and BED.**
   - **GFF3** → color by a chosen attribute key (column 9), default `Name`,
     fallback synthetic `(type)` = column 3. Distinct values get an auto-cycled
     palette with a per-value legend (visibility toggle + color override).
   - **BED** → single color per file (BED carries no rich attributes). One color
     picker per loaded BED source.
   - Format chosen by extension (`.gff`/`.gff3`/`.gff3.gz` vs `.bed`/`.bed.gz`).

3. **Loading — per axis.** `File ▸ Load GFF3/BED for query…` / `…for subject…`,
   plus CLI `--gff-query/--bed-query` and `--gff-subject/--bed-subject`. In
   self-mode a single loader/flag feeds both axes from one parse.

4. **Alpha — one global slider.** Color pickers are opaque RGB; the single alpha
   slider is the only source of "darkness", per the requested UX. Overlaps still
   self-darken through compositing.

5. **Parser — `noodles-gff`** (spec §6.3 pin) for GFF3; a tiny hand-rolled BED
   reader in `dottir-io` (BED is a fixed-column TSV; ~40 LOC, avoids a second
   noodles crate). Both emit the same internal `AnnotSet`/`Feature` types.

## Consequences

* New dep: `noodles-gff` (spec-pinned). `flate2` (already present for FASTA)
  covers transparent `.gz`. No `noodles-bed`. Categorical feature colors use a
  small built-in palette (deterministic, dependency-free) rather than the
  spec-pinned `colorgrad`, which targets continuous colormaps.
* `dottir-core` stays I/O-free: parsing lives in `dottir-io`; the band renderer
  lives in `dottir-gui` and consumes plain interval data.
* This is a GUI/IO overlay only — it does **not** touch the pixelmap algorithm,
  so no pixelmap-format version bump and no golden-pixelmap regeneration. The
  interval→pixel coordinate transform gets its own unit test.
* v1 is one annotation source per axis (no layered multi-file). Layering and PAF
  HSP overlays remain future work.

## Alternatives considered

1. **Intersection blocks / gutter tracks** instead of full-span bands. Bands
   chosen because the primary question is "do matches fall inside annotated
   regions?", which bands answer directly while still surfacing region-pairs at
   crossings.
2. **Hand-rolled GFF3 parser** (the ADR 0004 alternative). Rejected: the spec
   pins `noodles-gff` and the MSRV reason to avoid it is gone. BED is kept
   hand-rolled because it is trivial and there is no spec pin for it.
3. **Color GFF3 by feature type (col 3) or single color.** Rejected as the
   default in favor of by-attribute, which is more flexible for repeat families
   (`Name`/`class`); `(type)` remains available as a selectable key.
