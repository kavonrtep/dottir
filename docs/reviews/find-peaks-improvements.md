# find_peaks: deferred-item discussion

**Status:** design discussion. Four items deferred from the initial
`dottir find-peaks` port, surfaced during real-data validation
(TRC_7, TRC_2). Each is sketched with options, a recommendation,
and the open question to lock in before implementation.

The items interact (#4 in particular changes the algorithm output
that #1, #2, #3 build on), so suggested implementation order is
**#4 → #2 → #3 → #1**.

---

## 1. Inline `--find-peaks <PATH>` flag on `dottir periodogram`

### Problem

Today the workflow is two commands:

```bash
dottir periodogram input.fa -o p.tsv --fft fft.tsv
dottir find-peaks p.tsv -o peaks.tsv
```

A one-pass workflow would be more convenient for the common case:

```bash
dottir periodogram input.fa -o p.tsv --fft fft.tsv --find-peaks peaks.tsv
```

### Options

* **Minimal** — single `--find-peaks <PATH>` flag, uses defaults
  from the standalone subcommand. For tuning, user runs the
  standalone subcommand on the saved TSV.
* **Full** — mirror every `find-peaks` flag with a `--peaks-`
  prefix (`--peaks-min-score`, `--peaks-min-harmonics`,
  `--peaks-subrepeats`, …). Lots of CLI surface for marginal
  value.

### Decisions to make

* **Input source:** periodogram only, FFT only, or both? Periodogram
  is the primary signal; FFT is diagnostic. Recommend periodogram
  only (operates on the in-memory periodogram before TSV write).
* **Subrepeats default:** off (matches standalone subcommand). User
  opts in via `--peaks-subrepeats` (full mode) or runs standalone
  command (minimal mode).
* **Output requirement:** doesn't require `-o` — can produce peaks
  even when the full periodogram isn't saved.

### Recommendation

**Minimal.** One flag (`--find-peaks <PATH>`), default settings,
periodogram-only. Anyone wanting tuned output runs the standalone
subcommand on the saved TSV. Keeps the inline flag's CLI surface
trivial; power users get the standalone tool.

---

## 2. Auto-detect `--fft-input signal_mean` when source used analytical z-score

### Problem

`dottir periodogram --fft` currently defaults `--fft-input` to
`z_score`. When the user runs with `--z-score analytical`, the
z_score column has uniformly-inflated noise (~120 across all
positions for large records). The FFT then operates on the noise
structure rather than the real periodic signal — spurious peaks
at short periods dominate the output.

Validated on TRC_7: with analytical z-score + default FFT input,
the FFT identified the same 174 bp family in some records, but
the "extra fundamentals" (period 19.76, 22.14, etc.) all came
from FFT-ing the analytical noise.

### Three ways to fix

| Option | Behaviour | Risk |
|---|---|---|
| **(a) Default `--fft-input signal_mean` always** | Robust regardless of z-score mode. signal_mean is always clean. | Slight regression for users who specifically want denoised z-score input — but the empirical-z-score denoising benefit was marginal for FFT peak detection. |
| **(b) Auto-switch based on z-score mode** | `signal_mean` when analytical, `z_score` when empirical. "Smart". | Magic; user might be confused why output changes when toggling z-score mode. Two code paths to maintain. |
| **(c) Keep current default, log a warning** | "z_score with analytical mode is contaminated; consider --fft-input signal_mean". | Footgun stays; warnings are easily ignored. |

### Recommendation

**Option (a).** signal_mean works in both regimes, no surprises,
simplest mental model. The "denoise via empirical z-score" benefit
was theoretical — TRC_7/TRC_2 showed signal_mean works fine for
FFT peak detection. Bonus: removes the need for the
analytical/empirical decision to propagate to FFT-input selection.

### Documentation impact

Behaviour change (default flips). Note in CHANGELOG; mention in
`dottir periodogram --help` that `--fft-input z_score` is now
opt-in.

---

## 3. Adaptive `--min-score`

### Problem

Real datasets show wide variation in per-record signal-to-noise.
A fixed threshold either:

* Under-cuts (lets noise through on small records).
* Over-cuts (misses weak signals on large records).

TRC_2 example: the 335-bp HOR family appears in records ranging
from 13 kbp (signal_mean ~45) to 452 kbp (signal_mean ~10). Same
biological repeat, very different per-bin signal-to-noise.
Default `--min-score 10` catches the small-record peak but the
large-record signal sits right at the threshold.

### Four candidate approaches

| Approach | Formula | Pros | Cons |
|---|---|---|---|
| **(a) Mean + 3σ** | `threshold = mean(scores) + 3 × stdev(scores)` | Standard 3-sigma rule. Easy. | Inflated by real peaks (outliers); over-cuts when many strong peaks present. |
| **(b) MAD-based** | `threshold = median(scores) + N × 1.4826 × MAD(scores)` | Robust to outliers (real peaks don't inflate threshold). | Slightly more code. Median ≈ 0 for sparse periodograms; effectively defaults to N×MAD. |
| **(c) Percentile** | `threshold = scores.percentile(99)` | Conceptually simple: "keep top 1%". | Hard-codes a fraction; long records with many peaks get cut, short records with few peaks under-cut. |
| **(d) Fixed default + adaptive floor** | Keep `--min-score` fixed default, add `--auto-threshold` for adaptive override, with `--min-score` as the absolute floor the adaptive value can't drop below. | Combines safety net + adaptation. | Two flags interact. |

### Recommendation

**(b) MAD-based, opt-in via `--auto-threshold`.** Keep the fixed
default for predictable behaviour; add an opt-in adaptive mode for
varied-size datasets. Default formula: `threshold = median + 5 × 1.4826 × MAD`
(5σ-equivalent, conservative).

Combined with **(d)**: `--min-score` becomes the floor — adaptive
threshold can be above it, never below. Prevents the auto mode
from dropping noise on inputs where the noise floor is genuinely
very low.

### Open question

Should adaptive be the default? Recommend **no** — predictable
thresholds make outputs comparable across runs and across records
of the same dataset. Adaptive should be opt-in via `--auto-threshold`.

### Validation plan

Re-run on TRC_2 with `--auto-threshold` before locking in. Confirm:

* All 23 records get a fundamental (today some need explicit
  low `--min-score 5`).
* No false-positive noise peaks in small records (which previously
  passed at `--min-score 10` but might now fail at adaptive
  threshold).

---

## 4. Stronger harmonic-ratio matching at high `n`

### Problem

At high integer `n` (n=20+), the relative tolerance (2%) allows
accidental matches between unrelated periods. The test case from
the port:

* Family A: period 7 (peaks at 7, 14, 21, 28, …)
* Family B: period 13 (peaks at 13, 26, 39, 52, 78, …)

Peak at 78 (= 6 × 13):
* `78 / 7 = 11.143`, n=11, err=0.013 < 2% → **matches 7** as 11th harmonic
* `78 / 13 = 6.000`, n=6, err=0 → **matches 13** as 6th harmonic (correct)

Current algorithm picks the **first** match in iteration order,
so 78 ends up classified as 11th harmonic of 7 — wrong.

### Two-part fix

#### Part A — best-match instead of first-match

Among all fundamentals a candidate matches, pick the **lowest n**.
Break ties by smallest absolute distance to expected position.

Biological rationale: a peak is most likely a direct multiple of
the closest-matching fundamental, not a distant 11th harmonic of
an unrelated one.

For period 78:
* vs 7: n=11, distance 1
* vs 13: n=6, distance 0
* Lowest n → match 13 ✓

For period 1572 in real HOR data (524 bp fundamental, 174 bp monomer):
* vs 524: n=3, distance 0
* vs 174: n=9, distance 6
* Lowest n → match 524 ✓ (the HOR is the closer parent)

#### Part B — tighter tolerance

Combine relative + absolute:

```
tolerance_bp = max(absolute_floor_bp, relative_pct × n × fundamental)
```

* `absolute_floor_bp = 1` (default)
* `relative_pct = 0.02` (current default)
* For 7-bp fundamental at n=11: relative = 0.02 × 11 × 7 = 1.54 bp, absolute floor = 1 → use 1.54 bp
* For 524-bp at n=3: relative = 0.02 × 3 × 524 = 31 bp → use 31 bp

Distance threshold scales sensibly with `n × fundamental` (the
expected absolute period) instead of just `fundamental`.

#### Combined effect

Part A fixes the family-assignment problem. Part B tightens
marginal matches that come from rounding noise. Both together is
the robust answer; Part A alone fixes the TRC test case but
leaves some borderline matches in real data.

### Recommendation

**Both A and B.** Default `absolute_floor_bp = 1`, keep
`relative_pct = 0.02` for backwards compatibility.

### Regression test

Add a test fixture: two coprime fundamentals with overlapping
high-harmonic positions. Today's test (`two_independent_periods_kept_separate`)
uses periods 13 and 23 specifically to avoid this — should add a
companion test with 7 and 13 (the bug case) and confirm it now
passes.

---

## Suggested implementation order

1. **#4 (matching fix)** — affects algorithm output, should land
   first so #1/#2/#3 build on the corrected behaviour.
2. **#2 (FFT default)** — small (default flip + CHANGELOG note).
3. **#3 (adaptive threshold)** — opt-in `--auto-threshold`, validated
   on TRC_2.
4. **#1 (inline flag)** — wires the now-improved find-peaks into
   the periodogram subcommand.

Roughly **3 small commits** (#1, #2, #4) and **1 medium commit**
(#3).

---

## Open questions to lock in before implementation

1. **#1:** minimal flag (`--find-peaks <PATH>` only) or full flag set
   (mirror every find-peaks tunable with `--peaks-*` prefix)?
2. **#2:** OK to change default `--fft-input` to `signal_mean`?
   Behaviour change; note in CHANGELOG.
3. **#3:** Adaptive opt-in (`--auto-threshold`) with MAD-based
   default? Keep `--min-score` as the floor the auto value can
   never drop below?
4. **#4:** Best-match (lowest n) + tightened tolerance both?
   Default `absolute_floor_bp = 1`, keep `relative_pct = 0.02`?
