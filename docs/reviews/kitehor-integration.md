# Combined pipeline: `kitehor kite-periodicity` → `dottir find-peaks`

**Status:** exploratory integration. Documents one way to feed
`kitehor`'s k-mer distance histograms into `dottir find-peaks`
for HOR/founder classification.

`kitehor` ([github.com/kavonrtep/kitehor](https://github.com/kavonrtep/kitehor))
is a separate tool: it does k-mer interval analysis on tandem-repeat
arrays to extract candidate periodicities, then runs its own
rule-based HOR classifier (founder × k = tile). The classifier is
conservative — many real records get tagged `unresolved`.

`dottir find-peaks` operates on a denser periodicity signal
(dottir's windowed-sum periodogram) and produces a richer
classification: fundamentals + harmonics + sub-repeats, with
quantitative scores.

This document explores whether feeding kitehor's per-record `H[d]`
distance histogram into `dottir find-peaks` produces useful
extra information.

## Pipeline

```
FASTA
  ↓ kitehor kite-periodicity --dump-profile DIR
DIR/*.kite.tsv  (one per record; columns: d, H, bg)
  ↓ scripts/kite_profile_to_dottir.py --subtract-bg
synthetic dottir periodogram TSV
  ↓ dottir find-peaks --min-harmonics 0
classified peaks (fundamentals + harmonics + sub-repeats)
```

Concretely:

```bash
# Step 1: kite scan with profile dump (raw H[d] per record).
kitehor kite-periodicity input.fa \
    -o /tmp/kite.tsv \
    --dump-profile /tmp/kite_profiles

# Step 2: convert kite profiles to dottir periodogram-format TSV.
# --subtract-bg removes kite's noise envelope; recommended.
scripts/kite_profile_to_dottir.py /tmp/kite_profiles \
    -o /tmp/kite_as_dottir.tsv \
    --subtract-bg

# Step 3: classify peaks. --min-harmonics 0 is required because
# kite's H[d] for a pure tandem has ONE peak (just the monomer
# period), unlike dottir's periodogram which shows the full
# harmonic ladder. Without this, pure tandems get filtered out.
dottir find-peaks /tmp/kite_as_dottir.tsv \
    --min-score 5 --min-harmonics 0 \
    --subrepeats --subrepeat-min-score 2 \
    -o /tmp/combined_peaks.tsv
```

## Why this works (and why the defaults differ)

`dottir find-peaks` was designed for dottir's own periodogram,
which contains the **autocorrelation-style** signal: a tandem
repeat of period `P` shows peaks at `k = P, 2P, 3P, …` —
fundamental plus a harmonic ladder. The default `--min-harmonics 1`
exploits this to drop lone-peak noise.

Kitehor's `H[d]` is **k-mer-distance histogram**: for a tandem
repeat of period `P`, you typically see ONE peak at `d = P`,
because the most common k-mer distance is one monomer apart. No
harmonic ladder. Pure tandems have `n_harmonics = 0` in the
combined output — this is a feature, not a bug, and the
`--min-harmonics 0` flag is essential.

For **HORs**, kite's `H[d]` does show a partial ladder because
inter-HOR distances reinforce multiples of the founder period.
The hor_k5 smoke test demonstrates this cleanly: founder=150
with `n=5` harmonic at 750 (the HOR tile) plus higher multiples.

## Smoke-test result (`kitehor/test_data/smoke/`)

3-record synthetic fixture: one pure tandem and two HORs.

| Case | Truth | Kite native | Combined |
|---|---|---|---|
| `tandem_pure` | monomer=60 | `tandem`, founder=60 ✓ | fund=60, `n_harmonics=0` ✓ (single peak — correct tandem) |
| `hor_k3` | monomer=100, k=3, tile=300 | `unresolved` ✗ | fund=33 (alias artifact) ✗ |
| `hor_k5` | monomer=150, k=5, tile=750 | `hor`, founder=150, k=5, tile=750 ✓ | fund=150, harmonics `2,3,4,5,6,10,15,20,30` ✓ — n=5 = HOR tile! |

The combined pipeline EXPLICITLY shows the HOR multiplicity as
`harmonic_n` in the ladder (`n=5` for hor_k5 → tile=750), which
kite's `hor_reason` field hides inside a free-text string.

## Real-data result (TRC_7, 32 records)

| Agreement category | Count | Notes |
|---|---|---|
| Same fundamental ±5 bp | 17 | Both pipelines converge on 522/524 (HOR-3 family) or 335/337 (HOR-2 family) — high-confidence call |
| Disagree: kite finds monomer, dottir finds HOR | ~10 | E.g. chr5_4581 (kite 183 / dottir 524); chr7_50208 (kite 186 / dottir 337); chr9_18747 (kite 186 / dottir 524). 174-186 is the **underlying monomer**; 337/524 is the **HOR organization** |
| Disagree: kite finds short artifact | ~5 | E.g. chr11_387 (kite 8, dottir 524). Kite's k-mer distance histogram at short ranges can pick up spurious aliases on records dominated by long HORs |
| Kite native classifier | Most `unresolved`, two `tandem` | Conservative — doesn't surface the rich periodic structure |

For comparison, dottir's own periodogram on the same records
(via `dottir periodogram --find-peaks`) consistently identifies
the HOR-scale fundamentals (524 / 337) across all 32 records.

## When is the combined pipeline useful?

* **Cross-validation of HOR calls** — if kite + dottir agree on
  the fundamental, confidence is high. Worth running both for
  publication-grade calls.
* **Surfacing monomer-vs-HOR structure** — kite tends to surface
  the **monomer period** (~170-186 bp for centromeric satellites),
  dottir tends to surface the **HOR period** (337 / 524 bp). The
  disagreement is informative: it tells you both the building
  block and the higher-order organization.
* **Confirming HOR multiplicity** — for clean HORs, the combined
  pipeline reports `founder` + `harmonic_n = k` (the HOR tile
  multiplicity) in one row. This is the same information kite's
  rule classifier produces, but in a uniform schema.

## Limitations

* **`min_harmonics 0` is required** — kite's H[d] doesn't have a
  classical harmonic ladder, so the default `min_harmonics 1`
  filter from dottir would drop everything. This is documented in
  the help text but easy to forget.
* **Aliasing in short k-mer-distance space** — kite's H[d]
  resolution depends on k-mer size (default 6). Short periods
  (<20 bp) can be aliased; combined-pipeline output at those
  scales is unreliable.
* **The combined pipeline does NOT do what kite's rule classifier
  does** — it doesn't enforce the `d1 = k × p_n` HOR rule. It
  just classifies peaks. For a clean HOR/tandem verdict, run
  kite's own classifier; use the combined pipeline for
  exploratory analysis or cross-validation.

## Files

* `scripts/kite_profile_to_dottir.py` — kite profile dump →
  dottir periodogram-format TSV adapter (stdlib Python).
* `kitehor/` (gitignored) — local clone of the upstream tool.
* `tmp/kite_smoke_*.tsv` and `tmp/kite_TRC_7_*` — exploration
  output (gitignored).
