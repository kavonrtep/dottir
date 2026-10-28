# ADR 0008 — Genome mode: whole-assembly dotplots from external alignments

* Status: Proposed
* Date: 2026-10-28
* Deciders: petr
* Spec reference: §5 (out of scope: built-in aligner), §10.7 (HSP source:
  consume external output), §4.6 (large mode)
* Plan: [`docs/GENOME_MODE_PLAN.md`](../GENOME_MODE_PLAN.md)
* Related: ADR 0004 (PAF loader deferred; superseded for PAF by this ADR),
  ADR 0007 (large mode)

## Context

Large mode (ADR 0007) scores seed-anchored diagonals with the exact window
score. It is repeat-aware, but its cost grows with the product of the
sequence lengths and with repeat content. Measured on 8 cores with
synthetic, repeat-rich genomes:

| Input | Large engine | minimap2 2.31 `-x asm20` | mm2plus 1.3 `-x asm20` |
|-------|--------------|--------------------------|------------------------|
| 100 Mb × 100 Mb | 82 s, 2.4 GB | 4.4 s, 2.0 GB | 4.9 s, 2.0 GB |
| 600 Mb self | 182 s, 10.3 GB, 6 % positions capped | — | — |
| 600 Mb × 600 Mb | — | 41 s, 11.0 GB | 50 s, 10.3 GB |

The PAF of minimap2 and mm2plus was identical. mm2plus reports up to 7×
on real genomes at 48 threads; the synthetic pairs give its parallel
chaining little to do, so the speedup is unmeasured here.

The target for this mode is D-GENIES scale (Cabanettes & Klopp 2018):
two multi-chromosome assemblies, up to a few Gb and thousands of contigs,
for synteny, rearrangements and assembly comparison. At that scale an
assembly aligner's repeat blindness (high-frequency seeds dropped, only
the best chains kept) is acceptable: repeats are studied one zoom step
down, where large mode takes over.

## Decision

1. **A third engine, "genome mode", above large mode.** It draws
   alignments (PAF records) instead of scoring matrix cells. The three
   engines form one zoom hierarchy:

   | Level  | Scale | Data |
   |--------|-------|------|
   | genome | assemblies, Gb, many records | PAF from an external aligner |
   | large  | ≲ 150 Mb visible | seed engine (ADR 0007) |
   | exact  | ≲ 150 kb visible | Dotter engine (§4.1) |

   Genome mode does not change the exact or large engines, their format
   versions or their goldens.

2. **Alignments come from an external program or a file.** dottir
   never links an aligner. `--paf FILE` loads existing alignments
   (minimap2, unimap, wfmash, MashMap, … — anything writing PAF), and
   `--aligner auto|mm2plus|minimap2|PATH` runs one. `auto` prefers
   `mm2plus` (same CLI and PAF as minimap2, MIT, Bioconda linux-64 only)
   and falls back to `minimap2` (Bioconda linux and macOS, x86_64 and
   arm64). This keeps §5 (no built-in aligner) and §10.7 (consume
   external output); the aligner is invoked, not reimplemented.

3. **Aligner runs are chunked, parallel and cached.**
   * The reference is indexed once (`-d ref.mmi`); query records are cut
     into chunks of `--chunk-size` (10 Mb, as D-GENIES) and aligned
     concurrently within a thread budget. Records split at a cut are
     re-joined when they are co-linear across it (same target and strand,
     gap below a tolerance on both axes).
   * The merged PAF is cached under the user cache directory, keyed by
     SHA-256 of both inputs, the aligner's `--version` output and its full
     argument list. Reopening a comparison does not re-align.
   * The default preset is `asm20` (diverged genomes allowed);
     `--preset asm5|asm10|asm20` and `--aligner-args` pass through.
   * Self-comparison runs with `-D --dual=no`: `-D` drops the trivial
     diagonal so repeat copies become primary alignments, `--dual=no`
     aligns each record pair once, and off-diagonal records are mirrored
     across the diagonal when drawn; the diagonal itself is drawn
     synthetically. On the 10 Mb synthetic genome (750 planted copies)
     plain `asm20` returns only the 4 diagonal records; `-D --dual=no`
     returns 797, including inverted repeats crossing the diagonal.
     `--self-all-copies` adds `-P` (all chains, copy against copy:
     31 k records there) for repeat-family structure.

4. **Display follows D-GENIES, rendered deterministically.**
   * Each alignment is a line from (q_start, t_start) to (q_end, t_end),
     reversed for `-` strand, coloured by identity bin. Identity is
     `1 − de:f` when the record has it (aligner run with `-c`), else
     matches / block length (column 10 / column 11).
   * Identity bins are those of D-GENIES: < 25, 25–50, 50–75, > 75 %,
     on matches / block length unless `de:f` is present. `-c` is not
     passed by default.
   * At overview zoom alignments are rasterised into per-pixel channels
     (max identity, aligned bp per strand); lines are drawn as vectors
     once the visible count falls below a limit. Rasterisation is
     order-independent (max-merge / sums), so output does not depend on
     aligner thread count or PAF order.
   * Filters: minimum length, minimum identity, primary-only, and the
     D-GENIES noise filter (alignments binned by length; bins below 1 %
     of the most populated bin are hidden). Records shorter than 0.2 %
     of the total length are merged into one grey "other" block.

5. **Assembly tools.** Query records can be sorted along the reference:
   each is placed on the reference record holding most of its aligned bp,
   at the length-weighted median of its aligned positions there, and
   oriented by its majority strand; unmatched records go last. Exports:
   association table (TSV), reordered / reoriented query FASTA,
   unaligned records (FASTA), filtered PAF.

6. **Sequences are read by region, not loaded.** Genome mode needs only
   record names and lengths (`.fai`). Drill-down into large or exact mode
   reads just the visible regions through a FASTA index (built on demand;
   bgzip inputs via `.gzi`, plain gzip decompressed once into the cache).
   Two 3 Gb assemblies are never held in memory.

7. **Provenance.** Exports (PNG tEXt, params sidecar) record
   `mode = genome`, the aligner name, version and arguments (or the PAF
   path and its SHA-256), chunk size, filters and sort state.

8. **Engine selection.** `--engine auto` selects genome mode by size:
   when query × subject reaches `--genome-threshold-cells` (10¹⁷, about
   316 Mb × 316 Mb) and an aligner is found. Without an aligner it stays
   on the large engine and logs why. The threshold is tuned on the G0
   corpora.

9. **Packaging.** The conda recipe lists `mm2plus  # [linux and x86_64]`
   and `minimap2  # [not win]` as run requirements. Windows builds ship
   without an aligner; `--paf` and a user-supplied aligner path work
   there.

## Consequences

* New modules: `dottir-core::genome` (I/O-free: alignment model, filters,
  noise filter, record ordering, layout, rasterisation), `dottir-io::paf`
  (reader, plain and gzip), `dottir-io::faidx` (index builder and region
  reader), and a process runner for the aligner outside `dottir-core`.
* Large-mode drill-down from genome mode runs on regions, so its repeat
  classes are estimated within the visible regions, not genome-wide. A
  genome-wide seed-depth track computed during the genome run is adopted
  only if it adds at most 25 % to the genome-mode wall time and fits the
  memory limit on the G0 corpora; otherwise it is skipped and classes
  stay region-local (documented).
* Genome mode is repeat-blind by construction (aligner seed filtering and
  chain selection). The book states this and points to large mode.
* minimap2 cannot index a single sequence ≥ 2³¹ bases; such records are
  rejected with a typed error naming the record.
* New external runtime dependency on Linux and macOS (optional at run
  time: without an aligner only `--paf` is offered).

## Alternatives considered

* **Push large mode further** (more sampling, scoring only at overview
  zoom). Large mode is quadratic in hits and, at 600 Mb, already capped
  6 % of positions at `max_occ = 1000`; matching an assembly aligner on
  Gb assemblies would mean re-deriving minimizer chaining. Rejected; large
  mode stays the repeat-aware middle level.
* **Link minimap2 (C library or the `minimap2` crate).** Adds a C build
  and loses Windows (minimap2 does not build there with MSVC); swapping
  in mm2plus, unimap or wfmash would mean another binding each. Rejected.
* **Reimplement minimizer chaining in Rust.** Violates §5 and duplicates
  a mature tool. Rejected.
* **PAF only, no runner.** Simplest, but caching, chunking and provenance
  are exactly what D-GENIES users get for free; the runner is small.
  Kept as the `--paf` path rather than the only path.
