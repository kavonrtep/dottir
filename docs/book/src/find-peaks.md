# The `dottir find-peaks` CLI

```text
dottir find-peaks <INPUT.tsv> [-o OUT.tsv] [FLAGS]
```

Classifies periodogram or FFT peaks into **fundamentals**,
**harmonics**, and optionally **sub-repeats**. Reads either a
`dottir periodogram` output (`-o`) or its `--fft` companion —
format is auto-detected from the column header.

For the one-pass workflow, you can also invoke this inline from
`dottir periodogram --find-peaks PATH` — same algorithm, same
defaults, no TSV round-trip. The standalone subcommand is for
tuning thresholds without recomputing the (expensive) periodogram.

## Required arguments

| Argument | Description |
|----------|-------------|
| `INPUT` | TSV path. Format auto-detected: periodogram (has `k, signal_mean, z_score` columns) or FFT (has `bin, period_residues, amplitude, peak_rank`). |

## Output flag

| Flag | Default | Description |
|------|---------|-------------|
| `-o, --output PATH` | stdout | Output TSV path. |

## Periodogram-only flags

| Flag | Default | Description |
|------|---------|-------------|
| `--rank-by {signal_mean, z_score}` | `signal_mean` | Periodogram column used for peak detection. `signal_mean` is robust regardless of z-score mode. `z_score` only when source used `--z-score empirical`. |
| `--min-score FLOAT` | 10 (signal_mean) / 5 (z_score) / 5 (auto-threshold floor) | Minimum value in the rank-by column for a peak candidate. Acts as the **floor** when `--auto-threshold` is on. |
| `--boundary-fraction 0..1` | `0.9` | Drop peaks at `period > boundary_fraction × max_period` — kernel edge artifacts saturate at `k ≈ N/2` on large records. Pass `1.0` to disable. |
| `--auto-threshold` | off | Adaptive per-record threshold = `max(floor, median + k × 1.4826 × MAD)` of the rank-by column. Useful for datasets with widely varying per-record signal-to-noise. |
| `--auto-threshold-k FLOAT` | `5` | `k` multiplier for the MAD-based threshold (≈ k-σ equivalent). Higher = stricter. |

## Classification flags (both formats)

| Flag | Default | Description |
|------|---------|-------------|
| `--harmonic-tolerance FLOAT` | `0.02` | Relative tolerance for integer-ratio harmonic detection (2%). |
| `--max-harmonic-n INT` | `30` | Largest integer harmonic considered. Long tandem arrays show clean harmonics out to n>20. |
| `--min-harmonics INT` | `1` | Keep only fundamentals with at least this many detected harmonics. `0` = no filter. A "fundamental" with zero harmonics is almost always noise. |
| `--top-n INT` | unlimited | Limit output rows per record. |
| `--show-harmonics` | off | Include harmonics in the output. By default only fundamentals (and sub-repeats, if enabled) are emitted. |

## Sub-repeat flags (periodogram only)

| Flag | Default | Description |
|------|---------|-------------|
| `--subrepeats` | off | Enable sub-repeat detection: for each fundamental, scan integer divisors at `fundamental / n` and surface any peak found as `kind=subrepeat`. |
| `--subrepeat-min-score FLOAT` | `5` | Minimum score for a sub-repeat (typically lower than `--min-score` — sub-repeats are weaker than their parent HOR). |
| `--max-divisor INT` | `6` | Largest integer divisor scanned (2..N). Most natural HORs are 2- to 5-mers. |
| `--period-tolerance INT` | `2` | `±` bp tolerance when matching sub-repeat candidates around `fundamental / n`. |

---

## Output TSV columns

One row per emitted peak, plus header comments at the top and a
per-record `# diag` line just before each record's data.

| Column | Description |
|--------|-------------|
| `record_id` | FASTA record id. |
| `period_bp` | Peak period in residues. |
| `score` | Value of the rank-by column (signal_mean / z_score / amplitude). |
| `kind` | `fundamental`, `harmonic`, or `subrepeat`. |
| `parent` | Period of the parent fundamental — populated for `harmonic` and `subrepeat`. |
| `harmonic_n` | Integer multiplier `n` such that `period ≈ n × parent`. Only for `kind=harmonic`. |
| `divisor_n` | Integer divisor `n` such that `period ≈ parent / n`. Only for `kind=subrepeat`. |
| `n_harmonics` | Count of distinct integer harmonic positions detected — only for `kind=fundamental`. |
| `harmonics` | Comma-separated list of detected harmonic positions (e.g. `2,3,4,5,6,7`). Gaps are a red flag for partial / local repeats. Only for `kind=fundamental`. |

### Per-record `# diag` line

```
# diag record_id=X length=N threshold=T threshold_mode={fixed|mad} floor=F
```

Visible whether or not the record has any peaks — tells you why an
empty record is empty (threshold too high vs. no candidates above
floor).

---

## How to read it

* **`fundamental`** is the answer for "what's the repeat period
  here?" Look for the strongest fundamental per record.
* **`n_harmonics`** is the signature of a clean tandem array: a
  fundamental with `[2, 3, 4, 5, 6, …]` is a textbook complete
  ladder. A fundamental with `[2, 7, 11, 13, …]` (gaps) is a
  partial or imperfect repeat.
* **`harmonic` rows** (shown only with `--show-harmonics`) explain
  *why* the fundamental was claimed: each lists its `harmonic_n`
  position in the ladder.
* **`subrepeat` rows** (shown only with `--subrepeats`) are the
  monomer underneath an HOR. Example: a 524 bp HOR family with a
  detected sub-repeat at `period=174, parent=524, divisor_n=3`
  indicates a 174 bp monomer organised into 3-mer HORs of 524 bp.
* The output is ordered: each fundamental followed immediately by
  its sub-repeats (when shown), then any harmonics (when shown).

### Interpreting `divisor_n` with HORs

| If you see | It means |
|---|---|
| Fundamental 524, subrepeat 174 div=3 | 524 bp HOR of three ~174 bp monomers |
| Fundamental 338, subrepeat 67 div=5 | 338 bp HOR of five ~67 bp monomers |
| No subrepeat under a fundamental | Either monomers don't self-match (variant HORs), or no integer divisor falls on a strict local maximum |

---

## Algorithm in one paragraph

Strict three-bin local maxima above `min_score` (with the
boundary cap) become candidates. Sort by score descending. Each
candidate is checked against every existing fundamental; if any
match within tolerance, it joins the *best-matching* parent's
ladder (smallest normalised error → smallest distance → smallest
`n` → strongest parent for determinism). Otherwise it's a new
fundamental. After this greedy pass, a consolidation step demotes
any fundamental that turns out to be a harmonic of another (with
its harmonics cascaded), then a reparenting pass re-evaluates
every harmonic against the *final* fundamental set so the result
is independent of discovery order. Sub-repeat detection runs last
on the final fundamentals.

---

## Common workflows

### Default analysis on a saved periodogram

```bash
dottir find-peaks periodogram.tsv -o peaks.tsv
```

Filters to fundamentals only with `min_score=10, min_harmonics=1,
boundary_fraction=0.9`. Same as `dottir periodogram --find-peaks`
inline mode.

### Adaptive threshold for varied-size datasets

```bash
dottir find-peaks periodogram.tsv --auto-threshold --subrepeats -o peaks.tsv
```

Per-record `median + 5×MAD` threshold with `--min-score 5` as
floor. Sub-repeats also surfaced. Works well when records span a
wide range of lengths and per-record signal strengths.

### Strict mode (only clean tandem repeats)

```bash
dottir find-peaks periodogram.tsv --min-score 20 --min-harmonics 5 -o peaks.tsv
```

Only fundamentals with at least five harmonics — filters out
isolated peaks and partial repeats.

### Permissive mode (catch everything, then filter in pandas)

```bash
dottir find-peaks periodogram.tsv \
    --min-score 2 --min-harmonics 0 --show-harmonics --subrepeats \
    -o peaks.tsv
```

Then in pandas:

```python
import pandas as pd
df = pd.read_csv("peaks.tsv", sep="\t", comment="#")
df[df.kind == "fundamental"].sort_values("score", ascending=False).head(20)
```

### FFT input

```bash
dottir find-peaks fft.tsv -o fft_peaks.tsv
```

Auto-detects FFT format. Classification convention flips (FFT:
*larger* period is fundamental, harmonics live at integer
divisors). Sub-repeat detection is periodogram-only — flag is
silently ignored on FFT input.

---

## Limitations

* **Harmonic classification can be wrong when no clean integer
  ladder exists.** A peak with only one detected harmonic (n=2)
  may be coincidence; check the periodogram directly to verify.
* **`--auto-threshold` falls back to floor when MAD=0** (sparse
  periodograms where most positions are zero). The diag line will
  show `threshold = floor` — the adaptive path simply has nothing
  to adapt to.
* **High-`n` near-integer coincidences** are filtered by the
  best-match algorithm but not eliminated entirely. If you see a
  classification that looks wrong, inspect the periodogram rows
  near the candidate's position.
