# dottir — Genome mode implementation plan

Implements [ADR 0008](./adr/0008-genome-mode.md): whole-assembly dotplots
from external alignments (mm2plus / minimap2 or a PAF file), at D-GENIES
scale, with drill-down into the large and exact engines. Estimates are
person-days for one focused developer.

Target release: **v0.6.0** (G0–G4, CLI); GUI drill-down (G5) may land in
v0.6.x.

## Progress

| ID | Item | Status |
|----|------|--------|
| G0 | ADR, spec section, corpora | ☐ ADR proposed |
| G1 | PAF reader, FASTA index, genome layout | ☐ |
| G2 | Filters, ordering, rasterisation, CLI `--engine genome --paf` | ☐ |
| G3 | Aligner runner: detection, chunking, cache, provenance | ☐ |
| G4 | Exports: association table, sorted / unaligned FASTA, filtered PAF | ☐ |
| G5 | GUI genome view and drill-down to large / exact | ☐ |
| G6 | Docs, book chapter, conda recipe, release | ☐ |

## Pipeline overview

```
reference.fa, query.fa ──► .fai (names, lengths; built on demand)
        │
        ├─ G3  aligner runner (or --paf FILE)
        │        ref.mmi once; query cut into 10 Mb chunks, aligned in parallel
        │        ─► merged PAF (split records re-joined) ─► cache keyed by SHA-256
        ▼
  G1  PAF records ─► Alignment {q_rec, q_start, q_end, t_rec, t_start, t_end,
        │                       strand, matches, block_len, identity, primary}
        ▼
  G2  filters (length, identity, primary, noise) ─► record ordering / orientation
        │                                              ─► GenomeLayout (axes)
        ├─► raster: per-pixel max identity + bp per strand ─► PNG / SVG
        ├─► G4 exports
        ▼
  G5  GUI: vector lines when few visible ─► zoom ─► large engine on faidx
                                                     regions ─► exact engine
```

Module layout:

```
crates/dottir-core/src/genome/     # I/O-free
├── mod.rs        # GenomeConfig, GenomePlot, compute_genome_plot()
├── alignment.rs  # Alignment, identity, strand, join across chunk cuts
├── filter.rs     # length / identity / primary filters, noise filter
├── order.rs      # record ordering + orientation along the reference
├── layout.rs     # record order → axis offsets, "other" block
└── raster.rs     # deterministic line rasterisation into channels
crates/dottir-io/src/
├── paf.rs        # PAF reader (plain + gzip), writer for filtered PAF
├── faidx.rs      # .fai builder, region reader (plain, bgzip via .gzi)
└── aligner.rs    # detection, version probe, chunked runs, cache
```

---

## G0 — Groundwork (≈ 2 d)

* ADR 0008 accepted; spec §4.7 "Genome mode" with MUST items (the
  determinism, provenance and filter semantics below) and a §5 note that
  invoking an external aligner is not a built-in aligner.
* Corpora:
  * extend `examples/simulate_genome.rs` to write a **fragmented query
    assembly**: chromosomes cut into contigs, a share of them reversed
    and shuffled, plus unplaced small contigs — gives ground truth for
    ordering and orientation;
  * one real pair (e.g. *A. thaliana* vs *A. lyrata*, the D-GENIES
    reference case: 1 min / 2 GB there) and one plant pair at 1–3 Gb
    for mm2plus vs minimap2 timings.
* Record baseline timings for minimap2 and mm2plus on all corpora.

## G1 — Inputs (≈ 4 d)

* `dottir-io::paf`: tolerant reader (12 mandatory columns, tags `tp:A`,
  `de:f`, `NM:i`, `cg:Z` ignored except for length checks), gzip
  auto-detected as in `fasta.rs`. Decide `noodles-paf` vs a hand-rolled
  parser after checking it handles `de:f` and gzip (CLAUDE.md prefers
  noodles; the format is 12 tab-separated columns).
* Validate PAF against the `.fai`: unknown record names and out-of-range
  coordinates are typed errors listing the first offenders.
* `dottir-io::faidx`: build `.fai` (samtools-compatible) next to the
  input or in the cache when the directory is read-only; region reader
  for plain FASTA and bgzip + `.gzi`; plain gzip is decompressed once into
  the cache. Needed by G4 exports and G5 drill-down.
* `dottir-core::genome::layout`: record order → axis offsets, with
  records below 0.2 % of the total merged into an "other" block.

Tests: PAF round trip, malformed lines, `.fai` byte-identical to samtools
on the test corpora, region reads across line wraps.

## G2 — Core and first CLI deliverable (≈ 5 d) — *first usable deliverable*

* Filters: min length, min identity, primary only (`tp:A:P`), noise
  filter (D-GENIES: length bins = alignments / 10; scan from the shortest
  bin to the most populated one, then on to the first bin with ≤ 1 % of
  its count; shorter alignments are noise).
* Ordering: per query record, the reference record with most aligned bp;
  position = length-weighted median of aligned midpoints there;
  orientation = majority strand by aligned bp; unmatched records last, by
  length. Deterministic tie-breaks by record name.
* Rasterisation: anti-diagonal lines for `-` strand; per pixel the max
  identity and aligned bp per strand, max-merge / sums only, so output
  is independent of PAF order. Identity bins as D-GENIES (< 25, 25–50,
  50–75, > 75 %), configurable.
* CLI: `dottir batch REF.fa QUERY.fa --engine genome --paf aln.paf
  -o plot.png [--min-aln-len] [--min-identity] [--primary-only]
  [--noise-filter] [--sort-query]`; sidecar `mode = genome`, PAF SHA-256,
  filters, sort state. SVG output draws vector lines.
* Goldens: pinned hash of the raster for the synthetic pair, plus a
  property test that shuffling PAF lines leaves the raster unchanged.

## G3 — Aligner runner (≈ 5 d)

* Detection: `--aligner auto` → `mm2plus` on PATH, else `minimap2`, else
  a typed error suggesting `--paf`. `--aligner PATH` accepts any
  minimap2-compatible binary. Version from `--version`; mm2plus prints
  its SIMD launcher lines to stderr, which are logged, not parsed.
* Indexing: `-x PRESET -d ref.mmi REF.fa` once per reference; records
  ≥ 2³¹ bases rejected up front.
* Chunking: query records cut at `--chunk-size` (10 Mb) into a temporary
  FASTA per batch; chunks run concurrently with `-t` split across them so
  total threads ≤ `--threads`. Coordinates are shifted back; records that
  end within a tolerance of a cut are joined with their co-linear partner
  on the other side.
* Self-comparison: `-D --dual=no`, off-diagonal records mirrored at draw
  time, diagonal drawn synthetically; `--self-all-copies` adds `-P`.
  Records with `q_rec == t_rec` and identical coordinates are kept: on
  the `-` strand they are inverted repeats crossing the diagonal.
* `--engine auto`: genome mode when query × subject ≥
  `--genome-threshold-cells` (10¹⁷) and an aligner is found; otherwise
  large engine with a log line saying why.
* Genome-wide depth (decision 4): time the seed-depth pass on both
  inputs alongside the aligner run on the G0 corpora; enable it only if
  it adds ≤ 25 % wall time and fits the memory limit.
* Cache: `$XDG_CACHE_HOME/dottir/genome/<sha256>.paf.gz`, key = SHA-256 of
  both inputs + aligner version + full argv; `--no-cache` and
  `--paf-out FILE` to keep a copy.
* Progress and cancellation (chunks done / total), so G5 can show it and
  stop the child processes.
* Provenance in the sidecar: aligner, version, argv, chunk size.

Tests: a fake aligner script (prints canned PAF per chunk) for chunk
shifting, joining, thread budgeting and cache keys; an integration test
gated on `minimap2` being on PATH.

## G4 — Exports (≈ 2 d)

* `--assoc-out FILE.tsv`: query record, length, assigned reference,
  strand, aligned bp, share of query aligned.
* `--sorted-query-out FILE.fa`: query records in sorted order,
  reverse-complemented where the majority strand is `-` (region reads
  through faidx, streamed).
* `--unaligned-out FILE.fa` (query and reference records with less than
  `--unaligned-threshold` aligned, default 5 %).
* `--paf-out` writes the filtered PAF.

## G5 — GUI (≈ 8 d)

* Open two FASTA files as a genome comparison (or pick the engine in
  *Settings*); the aligner runs in the compute worker with a progress bar
  and cancel. Without an aligner, the dialog asks for a PAF.
* Genome view: record grid lines and names on both axes, "other" block
  greyed; raster at overview, vector lines below a visible-count limit
  (≈ 100 k); hover shows the alignment (records, coordinates, strand,
  length, identity). Filter sliders and the noise toggle are
  display-only; sort / unsort and per-record flip.
* Drill-down: when the visible area falls within the large engine's
  limits (default ≤ 150 Mb visible on both axes), read the visible
  regions through faidx and hand them to the large engine as a
  multi-record input with offsets; the exact engine follows as today.
  Requires a region-input path in `LargeInput` (record spans with
  global offsets) and the `Arc<[u8]>` sequence storage from L4.
* Session file stores the genome view (inputs, PAF or aligner argv,
  filters, order, flips).

## G6 — Documentation and release (≈ 2 d)

* Book chapter "Genome comparisons": workflow, D-GENIES correspondence,
  repeat blindness and when to zoom into large mode, Windows and `--paf`.
* Recipe: `mm2plus  # [linux and x86_64]`, `minimap2  # [not win]` in
  `run`; test command runs a tiny genome comparison when an aligner is
  present.
* CHANGELOG, measured timings table on the G0 corpora.

## Risks and mitigations

| Risk | Mitigation |
|------|------------|
| Aligner memory at Gb scale (600 Mb pair: 11 GB unchunked) | Query chunking (D-GENIES: 81 → 36 GB for human × chimp); reference index reused; memory estimate logged before the run |
| mm2plus linux-64 only, may differ slightly from minimap2 on large genomes | `auto` falls back to minimap2; aligner and version recorded in the sidecar and the cache key |
| Records split at chunk cuts look fragmented | Join co-linear pieces across cuts; test with the fake aligner |
| Millions of alignments slow the GUI | Raster at overview, vector only below a visible-count limit |
| Drill-down classes differ from a whole-genome large run | Classes are region-local; documented and shown in the panel |
| Unindexed gzip inputs | Decompress once into the cache, then faidx |
| PAF from other aligners with missing tags | Identity falls back to columns 10 / 11; tags are optional |

## Decisions

1. Identity bins: D-GENIES settings (< 25, 25–50, 50–75, > 75 % on
   matches / block length; `de:f` used when present; no `-c` by default).
2. Self-comparison: `-D --dual=no` plus mirroring; `-P` opt-in. Measured
   on the 10 Mb synthetic genome with `asm20`:

   | Flags | Records | Aligned | Time |
   |-------|---------|---------|------|
   | none | 4 (diagonal only) | 10 Mb | 1.7 s |
   | `-D` | 756 | 9.8 Mb | 1.1 s |
   | `-D --dual=no` | 797 | 8.9 Mb, mirrored at draw | 1.1 s |
   | `-D -P --dual=no` | 31,123 | 72 Mb | 1.1 s |
   | `-X` (= `-DP --dual=no --no-long-join`) | 31,760 | 62 Mb | 0.7 s |

3. `--engine auto` selects genome mode by size (≥ 10¹⁷ cells, aligner
   present); threshold tuned in G0.
4. Genome-wide depth track only if ≤ 25 % extra wall time and within the
   memory limit; otherwise skipped (region-local classes).
