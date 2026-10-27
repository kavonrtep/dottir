# ADR 0007 — Seed-based large-sequence mode

* Status: Proposed
* Date: 2026-10-27
* Deciders: petr
* Spec reference: §4.1 (algorithmic invariants), §4.5.1 (performance),
  §5 (out of scope: built-in aligner)
* Plan: [`docs/LARGE_MODE_PLAN.md`](../LARGE_MODE_PLAN.md)

## Context

The exact engine scores every cell of the `qlen × slen` matrix. That is the
scientific contract (§4.1) and it reaches the §4.5.1 baseline of 1 Mb × 1 Mb,
but at 100 Mb × 100 Mb (10¹⁶ cells) it is out of reach regardless of SIMD or
thread count. Target use cases at that scale are both **repeat study**
(TE / satellite structure in self-plots) and **synteny** (chromosome vs
chromosome, multi-record assemblies), and the similarity of interest can be
**significantly diverged**, which rules out long exact k-mers or sparse
minimizers on their own.

## Decision

1. **A second, separate engine — "large mode" — for DNA only.** It lives in
   `dottir-core::large` and does not modify the exact engine. The exact
   engine's `PIXELMAP_FORMAT_VERSION` and C-dotter goldens are untouched;
   large mode gets its own `LARGE_PIXELMAP_FORMAT_VERSION` and its own
   deterministic goldens.

2. **Seeding: spaced seeds over a direct-address index.** Default one pattern
   of weight 12 (span ≈ 19, LASTZ-style 12-of-19). Subject positions are
   indexed in a `4^w` count table + `u32` position array (no sort, no
   hashing). Seeds containing `N` or crossing a record boundary are skipped.
   Additional patterns are an opt-in once the sensitivity benchmark justifies
   them.

3. **Two-hit diagonal filter.** A seed hit survives only if another hit lies
   on the same diagonal (± band) within a configurable distance (default
   100 bp). This removes ≈ 99.999 % of random weight-12 hits while keeping
   ≥ 75 %-identity similarity.

4. **Seeded exact scoring.** Surviving anchors define diagonal spans; each
   span is scored with the *exact* sliding-window recurrence (same matrix,
   Karlin window size, pixel factor, anti-diagonal rule, max-merge). Hence
   the invariant: **for every pixel, `large ≤ exact`, with equality on every
   pixel containing an anchor span.** This keeps the greyramp semantics
   identical between modes.

5. **Repeat classes from a separate depth track.** Copy number is taken from
   long k-mers (default k = 21, configurable) counted across all inputs, or
   from GFF3/soft-mask annotations. Each anchor is binned into a class
   (unique / low / mid / high); pixels hold one `u8` channel per
   *(strand × class)*. Showing or hiding repeats is a display operation with
   no recompute. A per-seed occurrence cap bounds runtime; pixels where the
   cap fired are flagged.

6. **Basic co-linear chaining as a display aid.** Anchors are chained per
   *(record pair, strand)* with a minimap2-style simplified DP (bounded
   look-back, gap + diagonal-drift penalty). Chains are drawn as vector lines
   and exportable as PAF-like records **without CIGAR**. No gapped
   extension or base-level alignment is produced. This amends §5: chaining
   is permitted, alignment remains out of scope.

7. **Semantic zoom.** The GUI keeps the index and recomputes only the
   viewport; below a cell-count threshold it hands off to the exact engine.

## Consequences

* No new dependencies are required.
* Large-mode output is not comparable to C dotter except via the `≤ exact`
  invariant, which is property-tested on small inputs.
* Memory budget differs from the exact engine (index ≈ 0.5 GB per pattern
  per 100 Mb); large mode gets its own `memory_limit` default (4 GB).
* Exports (PNG tEXt, params sidecar) MUST record `mode = large` and the
  seed/filter/class parameters so a figure is reproducible.

## Alternatives considered

* **Faster brute force** — insufficient by 3–4 orders of magnitude at 100 Mb.
* **Minimizers alone** — too insensitive for diverged repeats; kept as a
  possible future sampling option on the same index.
* **Suffix array (Gepard)** — allows changing word length without
  re-indexing, but heavier to build and query than a fixed-weight table.
* **Full chaining + alignment** — duplicates minimap2; users can overlay
  minimap2 PAF instead.
