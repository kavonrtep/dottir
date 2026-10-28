# dottir — Large-sequence mode implementation plan

Implements [ADR 0007](./adr/0007-large-sequence-mode.md): seed-based dotplots
for DNA inputs at the 100 Mb scale, multi-record, with repeat classes and
basic chaining. Estimates are person-days for one focused developer.

Target release: **v0.5.0** (phases L0–L3 + CLI), GUI and chaining may land in
v0.5.x.

## Progress

| ID | Item | Status |
|----|------|--------|
| L0 | ADR, spec section, benchmark corpora | ☑ synthetic generator (`examples/simulate_genome.rs`); real plant benchmarks still to run |
| L1 | Seed index + query streaming + two-hit filter | ☑ |
| L2 | Seeded exact scoring → dot layer, CLI `--engine large` | ☑ |
| L3 | Repeat depth track, classes, occurrence cap | ☑ plus `seed` depth method (see ADR 0007) |
| L4 | GUI integration + semantic zoom | ☑ including depth/capped tracks along the axes and GFF3/BED repeat floors |
| L5 | Chaining layer + export | ☑ PAF, BEDPE, SVG chain lines |
| L6 | Docs, book chapter, release | ☑ docs; release pending |

Measured (8 cores, synthetic repeat-rich genomes, seed depth): 100 Mb
self-plot 2 min 14 s / 1.9 GB; 100 Mb × 100 Mb 4 min 7 s / 2.5 GB.

Deviations from the plan below:

* Capped seeds are reported per pixel column/row (`q_capped`, `s_capped`),
  not as a per-pixel flag.
* The CLI flag is `--engine` (`--mode` was already the BLAST mode).
* Large goldens are a pinned FNV hash of all channels + chains in
  `large/tests.rs` rather than PNG files.

## Pipeline overview

```
query, subject (concatenated, RecordSpan ranges — dottir-io::Sequence)
  │
  ├─ L3  depth track: canonical k21 counts over all inputs → u8 per base
  ├─ L1  SeedIndex(subject): spaced seed, direct-address table
  │
  ▼ L1  stream query chunks (rayon), forward + revcomp
seed hits ─► occurrence cap ─► two-hit diagonal filter ─► anchors
                                                            │
             ┌──────────────────────────────┬───────────────┤
             ▼ L2                           ▼ L5            ▼ L3
   diagonal spans → exact          chains per             class per anchor
   sliding-window score →          (record pair, strand)  (from depth track)
   max-merge into channel                                    │
   [strand × class]  ◄───────────────────────────────────────┘
```

Module layout (new, all in `dottir-core`, I/O-free):

```
crates/dottir-core/src/large/
├── mod.rs        # LargeConfig, LargePlot, compute_large_dotplot()
├── seed.rs       # SpacedSeed pattern, 2-bit packing, rolling extraction
├── index.rs      # SeedIndex (count table + positions)
├── hits.rs       # query streaming, two-hit filter, anchors
├── score.rs      # anchor spans → exact window scores → channels
├── depth.rs      # k-mer depth track, RepeatClass
└── chain.rs      # co-linear chaining
```

---

## L0 — Groundwork (≈ 2 d)

**Tasks**

- Accept ADR 0007; add spec §4.6 "Large-sequence mode" (requirements in
  RFC-2119 form, mirroring the ADR decisions) and amend §5 (chaining allowed,
  alignment not).
- **Synthetic corpus generator** (test-only, `dottir-core/tests/support/` or a
  `dev` binary): random background genome with *planted* events at known
  coordinates and divergence — segmental duplications, inversions, tandem
  arrays, TE families (N copies, divergence 0–30 %, inverted copies). Emits
  FASTA + truth BED-pairs. Deterministic from a seed.
- Pick real benchmark pairs (to be supplied by petr): one ~100 Mb
  chromosome self-plot (repeat-rich plant) and one chromosome × chromosome
  synteny pair. Not committed; paths recorded in `tests/golden/large/README.md`.
- Record baseline: exact engine wall-time/RSS at 1, 2, 5 Mb to fix the
  auto-mode threshold.

**Acceptance**: spec + ADR merged; generator reproduces byte-identical FASTA
for a given seed.

---

## L1 — Seed index and hit streaming (≈ 5 d)

**Tasks**

- `seed.rs`: `SpacedSeed { mask: u64, span: u8, weight: u8 }` parsed from a
  `"110110110…"` string; default 12-of-19. Rolling 2-bit encoding of the
  span; extract the masked key (`pext`-free bit gather, table-driven). Any
  non-ACGT in the span invalidates the key. Reverse complement handled by
  encoding the RC stream, not by key transforms (simpler, spaced seeds are
  not RC-symmetric).
- `index.rs`: `SeedIndex::build(subject, records, seed, cfg)`. Two passes:
  count → prefix sum → fill `u32` positions. Positions within a bucket are
  ascending (needed for range queries in L4). Parallel build via per-chunk
  counts, then deterministic placement by chunk order. Memory check against
  `LargeConfig::memory_limit` *before* allocation; typed error otherwise.
- Record boundaries: a seed whose span crosses a `RecordSpan` end is skipped
  (use `Sequence::breaks()`).
- `hits.rs`: split the query into chunks (e.g. 1 Mb + span overlap), rayon
  over chunks, both strands. For each query key, iterate its bucket.
  - **Occurrence cap**: buckets with count > `max_occ` (default 1000,
    configurable) are skipped and counted into a per-pixel *capped* bitmap.
  - **Self-comparison**: emit only `s ≤ q` hits (lower triangle), mirror at
    render time per existing `Triangle`.
  - **Two-hit filter**: per-chunk table keyed by diagonal (`d = s − q`,
    bucketed by `band`) holding the last hit position; accept a hit when a
    prior hit on the same / adjacent bucket is within `two_hit_dist`. Emit
    both hits of the first pair.
- Reverse strand coordinates mapped back to forward subject coordinates
  using the same convention as the exact engine's reverse channel, so
  overlays/crosshair line up.
- Progress + cancellation: `LargeConfig` takes an optional
  `&AtomicBool` cancel flag and a progress callback (`Fn(f32) + Sync`); both
  are plain values, keeping core I/O-free.

**Tests**

- Unit: seed key extraction vs a naïve reference on random strings
  (proptest); N / boundary skipping; bucket ordering.
- Two-hit filter: random sequences of 1 Mb produce < 1e-4 × raw hits;
  planted 75 %-identity 2 kb duplicate is detected with ≥ 95 % of its
  length covered by anchors.
- Determinism: identical anchor set across 1, 4, 16 threads.

**Acceptance**: anchors for a synthetic 10 Mb self-plot in < 10 s on 8 cores.

---

## L2 — Seeded exact scoring and CLI (≈ 5 d) — *first usable deliverable*

**Tasks**

- `score.rs`: group anchors per diagonal; merge overlapping
  `[q − W − ext, q + span + ext]` intervals; score each merged span with the
  exact sliding-window score on that single diagonal (reuse `score_vec` /
  matrix code; a 1-D recurrence along the diagonal). Apply the anti-diagonal
  suppression rule and `zoom × zoom` max-merge exactly as `pixel.rs` does.
- Window size and pixel factor: Karlin on composition (sampled for large
  inputs — deterministic stride sample), same clamps and overrides
  (`-W`, `--pixel-fac`).
- `LargePlot { width, height, channels: Vec<u8> /* [strand][class][px] */,
  capped: BitVec-like Vec<u64>, classes: RepeatClasses, chains: Vec<Chain>,
  params: LargeParams }` with `combined(mask)` → `Vec<u8>` for rendering.
  Until L3 there is a single class.
- Public API:
  ```rust
  pub fn compute_large_dotplot(
      query: &[u8], q_records: &[Range<usize>],
      subject: &[u8], s_records: &[Range<usize>],
      cfg: &LargeConfig,
  ) -> Result<LargePlot, DottirError>;
  ```
  (record ranges passed as plain ranges, so core does not depend on
  `dottir-io::Sequence`.)
- `LARGE_PIXELMAP_FORMAT_VERSION = 1` in `lib.rs`.
- CLI (`dottir batch`): `--mode exact|large|auto` (default `auto`: large
  when `qlen·slen` exceeds the L0-measured threshold; logs the choice).
  Flags: `--seed`, `--two-hit-dist`, `--max-occ`, `--large-memory-limit`.
  PNG tEXt + params sidecar record mode and all large-mode parameters.
- `.dot` export: write the combined visible channel as format 3 with a
  warning that class channels are not preserved (full container deferred).

**Tests**

- **Invariant (proptest, small inputs ≤ 5 kb, zoom 1..8):**
  `large[px] ≤ exact[px]` for all pixels, and equality on every pixel touched
  by an anchor span. This is the key correctness test for the mode.
- Anti-diagonal rule: dedicated test that large mode reproduces the exact
  engine's suppression on a planted anti-diagonal-adjacent pattern.
- Large goldens (`tests/golden/large/`): synthetic 1 Mb corpus PNG hashes,
  pinned against `LARGE_PIXELMAP_FORMAT_VERSION`.
- Sensitivity curve from the synthetic generator: fraction of planted
  events recovered vs divergence (0–30 %), logged to a CSV — decision input
  for multi-pattern seeds.

**Acceptance**

- 100 Mb self-plot (real benchmark) at 4000 px completes in ≤ 5 min on
  16 cores with peak RSS ≤ 3 GB (single pattern).
- Planted events at ≤ 20 % divergence recovered ≥ 90 %.

---

## L3 — Repeat depth track and classes (≈ 5 d)

**Tasks**

- `depth.rs`: count canonical k-mers (default k = 21, configurable 15–31)
  across *all* inputs (query ∪ subject, once if self). Deterministic,
  bounded memory: partition k-mers by top bits into P buckets, sort + count
  each partition, write back saturating counts. Output `Vec<u8>` per input
  position, log-binned copy number (exact up to 15, then log2 bins).
- Alternative depth sources (merged by max):
  - soft-mask: lowercase → at least `mid` class (`--softmask-as repeat|ignore`);
  - GFF3/BED features (existing `AnnotSet`): chosen feature types → `high`
    class (`--repeat-annot <file> --repeat-types LTR,...`). Binding to
    concatenated coordinates is done in `dottir-io` and passed to core as a
    `Vec<u8>` depth vector — core stays I/O-free.
- `RepeatClasses` thresholds: unique (1), low (2–4), mid (5–50), high (> 50);
  configurable. Anchor class = max of its query/subject position classes.
- Channels become `2 strands × 4 classes`. `combined(mask)` merges the
  selected subset by max.
- Optional: soft-mask seeds out entirely (`--softmask-as skip`).
- CLI: `--show-classes unique,low` / `--hide-classes high`, `--depth-k`,
  `--class-thresholds 1,4,50`; `--depth-out <bedgraph>` to export the track.
- Optional colour-by-class rendering (one colour per class) in PNG export.

**Tests**

- Depth counts vs naïve HashMap counter on random + planted-TE corpora.
- Class channels partition the single-class result: max over all class
  channels == L2 output, byte-identical.
- Synthetic: hiding `high` removes planted TE-family dots but keeps planted
  segmental duplication.

**Acceptance**: depth track for 2 × 100 Mb in ≤ 2 min, ≤ 2 GB extra peak.

---

## L4 — GUI integration and semantic zoom (≈ 7 d)

**Tasks**

- Sequences in the GUI are shared as `Arc<[u8]>`; `ComputeRequest` today
  clones `Vec<u8>` per request — unacceptable at 100 Mb.
- `compute_worker.rs`: a `LargeSession` cache holding `SeedIndex` + depth
  track keyed by (input hashes, seed, k); rebuilt only when those change.
- Viewport recompute: stream only the visible query range; restrict bucket
  iteration to the visible subject range via binary search over the sorted
  bucket. Keep the full-view result as the background layer while the
  viewport result streams in.
- **Hand-off**: when visible cells < threshold (default 2·10¹⁰, setting),
  dispatch the exact engine for the viewport instead; display shows which
  engine produced the current view.
- Progress bar + cancel button wired to the core cancel flag.
- Side panel "Large mode": seed, two-hit distance, max occurrences (these
  trigger recompute); class toggles, colour-by-class, capped-pixel overlay
  (display only, instant).
- Depth track drawn as a thin track along both axes (reuses the
  annotation-overlay axis machinery).
- Session TOML stores large-mode parameters.

**Tests**: `dottir-render` fixture for large mode (class toggles change
pixels only where expected); session round-trip.

**Acceptance**: 100 Mb self-plot opens (index build) in ≤ 3 min; zoom to a
10 Mb × 10 Mb viewport refreshes in ≤ 5 s; zoom below threshold shows the
exact engine result.

---

## L5 — Chaining (≈ 5 d)

**Tasks**

- `chain.rs`: per *(q record, s record, strand)*, sort anchors by
  `(q, s)`; DP with bounded look-back `h` (default 50) and max gap `G`
  (default 10 kb): `f(i) = max(w_i, max_j f(j) + α(j,i) − β(j,i))`, with
  α = overlap-adjusted anchor weight, β = gap + diagonal-drift penalty
  (minimap2 2018 form, simplified). Backtrack greedily; anchors used once.
- By default chain only `unique` + `low` class anchors (configurable), so
  repeats do not glue chains.
- Filters: min chain score, min anchors, min span.
- Determinism: stable sort keys, tie-break by anchor index.
- Output `Chain { q_rec, s_rec, strand, q_range, s_range, score,
  n_anchors }`.
- GUI: chains drawn as vector lines over the pixelmap (strand colours from
  the inverted-repeat palette), toggle + min-score slider; reuses the HSP
  overlay draw path where possible.
- Export: `--chains-out file.paf` (PAF columns 1–12, `tp:A:P`, `cm:i`,
  `s1:i`, no CIGAR) and BED-pair; SVG export includes chain lines.

**Tests**: synthetic synteny corpus (blocks with inversions and
translocations) — each planted block recovered as one chain, correct
strand, coordinates within ±500 bp; no chain across records.

**Acceptance**: chaining a 100 Mb × 100 Mb anchor set ≤ 30 s.

---

## L6 — Documentation and release (≈ 2 d)

- mdBook chapter "Large sequences": when to use, parameters, repeat
  classes, interpreting capped pixels, chaining caveats (not an alignment),
  worked synteny + TE examples.
- `cargo doc` example for `compute_large_dotplot`.
- CHANGELOG, spec status lines, ADR 0007 → Accepted, plan progress table.

---

## Risks and mitigations

| Risk | Mitigation |
|------|------------|
| Repeat-rich genomes explode hit counts despite two-hit filter | `max_occ` cap + capped-pixel flag; class-aware skipping (`--skip-classes high` at hit time) as a speed option |
| Single spaced seed too insensitive at > 20 % divergence | L2 sensitivity curve decides; multi-pattern support is designed in (`Vec<SpacedSeed>`), adds ~0.5 GB/pattern |
| Index memory at > 100 Mb subject | Index the shorter sequence; chunked subject indexing (multiple passes) as fallback |
| `u32` positions overflow (> 4.29 Gb concatenated) | Reject with typed error in v0.5; `u64` variant later if needed |
| Divergence between large and exact display semantics | `large ≤ exact` property test is part of CI |
| GUI memory duplication | `Arc<[u8]>` refactor in L4 is a prerequisite, not optional |

## Open decisions (resolve during implementation)

1. Auto-mode threshold for `--mode auto` (measured in L0).
2. Default `max_occ` and class thresholds — tune on the real plant benchmark.
3. Whether a modern multi-channel container (§10 open question) is needed to
   save large plots with classes, or PNG + params is enough for v0.5.
