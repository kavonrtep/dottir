# Review: find_peaks improvement plan

**Status:** response to `docs/reviews/find-peaks-improvements.md`
**Date:** 2026-05-19

This is a review of the four deferred `find-peaks` items, with concrete
opinions on the decision points and a few extra suggestions before
implementation.

## Summary position

I agree with the overall direction and the proposed implementation order:
fix harmonic assignment first, then change the FFT default, then add adaptive
thresholding, then wire in inline peak finding.

Two details should be corrected before implementation:

1. The proposed high-`n` tolerance formula in item #4 is effectively already
   what the current code does. `harmonic_ratio()` computes
   `abs(ratio - n) / n`, which is the same as
   `abs(period - n * fundamental) / (n * fundamental)`. In bp terms, the
   current acceptance window is already
   `relative_pct * n * fundamental`. Adding
   `max(absolute_floor_bp, relative_pct * n * fundamental)` only adds a low-end
   floor; it does not tighten high-`n` matching.
2. `--auto-threshold` plus "`--min-score` is the floor" will not recover weak
   records if the default floor remains 10. If TRC_2 records need
   `--min-score 5`, then `--auto-threshold` with an implicit floor of 10 cannot
   fix them. Either auto mode needs a lower default floor, or users must pass
   `--min-score 5 --auto-threshold`.

## Decision Points

| Item | Recommendation | Rationale |
|---|---|---|
| #1 inline `--find-peaks <PATH>` | Yes, minimal mode only | Useful for routine workflows, but the standalone command should remain the tuning surface. |
| #1 input source | Periodogram only | The periodogram is the primary signal; FFT peak finding is diagnostic and more parameter-sensitive. |
| #1 subrepeats default | Off | Keep parity with standalone defaults. Users who need subrepeats should run `dottir find-peaks` explicitly. |
| #2 FFT default | Change default to `signal_mean` | It is robust across z-score modes and avoids analytical z-score noise. |
| #3 adaptive threshold | Add opt-in MAD threshold | Good idea, but make floor semantics explicit and report the resolved threshold per record. |
| #4 harmonic matching | Implement first, but revise the matching rule | Best-parent assignment is needed. The tolerance change as written is not a high-`n` tightening. |

## #1 Inline `--find-peaks <PATH>`

I would implement the minimal version:

```bash
dottir periodogram input.fa -o p.tsv --find-peaks peaks.tsv
```

The inline mode should:

- operate on the in-memory periodogram signal, not by parsing the TSV it just
  wrote;
- use the same default ranking as standalone `find-peaks`, currently
  `signal_mean`;
- write only periodogram-derived peaks;
- not require `-o`, because users may want only `peaks.tsv`;
- keep subrepeats and custom thresholds out of the inline CLI surface for now.

One practical CLI detail: avoid allowing two outputs to compete for stdout. If
`--find-peaks -` is ever supported, it should be rejected when the periodogram
output is also stdout. The simpler first implementation is to require an actual
path for `--find-peaks`.

I would not mirror all `find-peaks` options behind `--peaks-*`. That creates a
second API surface for the same classifier and makes future changes harder.

## #2 Default `--fft-input signal_mean`

I agree with changing the default from `z_score` to `signal_mean`.

The current default is surprising because `periodogram --fft` defaults to
FFT-ing `z_score`, while default greyramp sensitivity makes `--z-score auto`
resolve to empirical mode. That is both expensive and, under analytical mode,
can produce the short-period artifacts described in the review.

Recommended behavior:

- `--fft-input signal_mean` by default;
- `--fft-input z_score` remains available for users who explicitly want it;
- help text should say that `z_score` is opt-in and is most meaningful with
  empirical z-scores;
- examples for quick period discovery should probably include `--z-score off`
  when the z-score column is not needed.

This is a behavior change, but it is the right one. Put it in the changelog.

## #3 Adaptive `--min-score`

I agree with adding adaptive thresholding, but I would make it opt-in and
record the resolved threshold in output comments.

Recommended first version:

```text
threshold = max(floor, median(scores) + k * 1.4826 * MAD(scores))
```

with:

- `--auto-threshold` to enable it;
- `k = 5` as the initial conservative default;
- `--min-score` acting as the floor;
- a per-record comment such as `# record_id=X threshold=12.34 threshold_mode=mad floor=5`.

The floor needs a clear decision. My preference:

- if `--auto-threshold` is set and the user did not pass `--min-score`, use an
  auto-mode floor of 5 for `signal_mean`;
- if the user passes `--min-score`, honor it as the floor exactly.

That preserves predictable behavior for explicit users while making
`--auto-threshold` useful for the motivating TRC_2 case. Keeping the default
floor at 10 is safer but undercuts the main reason to add auto mode.

Additional implementation notes:

- Compute the robust baseline after dropping non-finite scores.
- Consider computing the threshold over all eligible offsets after the boundary
  cap, not over offsets that will never be candidates.
- If `MAD == 0`, fall back to the floor instead of inventing a threshold.
- Keep the existing strict local-maximum extraction for the first version.
  Prominence and smoothing can come later.

## #4 Harmonic Matching

This is the most important item, and I agree it should land first. I would
adjust the proposed solution.

### What is correct

The current first-match behavior is fragile. If a candidate peak matches more
than one fundamental, parent selection should be based on match quality, not
iteration order.

There should also be a regression test for overlapping harmonic families, for
example periods 7 and 13 with a shared/ambiguous high harmonic.

### What I would change

Do not make "lowest n" the primary rule. It fixes the `78 = 6 * 13` example,
but it can misassign an exact high-order harmonic to a worse low-order near
match.

Use a best-match score instead:

1. accept only candidates within tolerance;
2. choose the parent with the smallest normalized error, for example
   `distance_bp / allowed_bp`;
3. break ties by smaller absolute distance;
4. then prefer lower `n`;
5. then prefer stronger parent score for determinism.

This still assigns 78 to 13 because the distance is zero for 13 and non-zero
for 7. It also keeps exact high-order matches from being stolen by approximate
low-order matches.

### Important implementation detail

Best-match should run against the final set of fundamentals, not only the
fundamentals that happened to be discovered earlier in the greedy pass.

The current flow can classify a peak as a harmonic before a better parent is
later discovered. `consolidate()` demotes fundamentals, but it does not fully
reparent every already-classified harmonic to the best final parent. I would add
a reparenting pass after consolidation:

1. run greedy classification;
2. consolidate fundamentals;
3. collect final fundamentals;
4. for every harmonic candidate, recompute its best parent among all final
   fundamentals;
5. rebuild each fundamental's harmonic list from those final assignments.

That makes the outcome independent of candidate ordering except where score
ordering is intentionally used to decide which peaks become candidate
fundamentals.

### Tolerance

The proposed `max(absolute_floor_bp, relative_pct * n * fundamental)` formula is
not a high-`n` tightening relative to the current code. It is the current rule
expressed in bp, plus a floor.

I would not add a new tolerance knob until the best-parent/reparenting fix is
validated on TRC_7 and TRC_2. If tightening is still needed after that, consider
one of these explicit policies:

- an absolute cap, such as `allowed_bp = min(cap_bp, relative_pct * expected)`;
- a lower default relative tolerance for high `n`;
- a parent-selection score that penalizes high `n` without rejecting exact high
  harmonics.

The first implementation can keep the existing `harmonic_tolerance = 0.02`.

## Additional Suggestions

### 1. Add candidate diagnostics

Peak finding is parameter-sensitive. It will be much easier to validate real
records if the output header includes per-record diagnostics:

- resolved threshold;
- number of strict local maxima;
- number passing threshold;
- number passing boundary filter;
- number of final fundamentals;
- rank column used.

This can be comments in the TSV, not extra data columns.

### 2. Improve boundary semantics

Standalone `find-peaks` currently infers the boundary cap from the maximum
period present in the TSV. That is fine for default periodograms, but it is not
the same as using true sequence length when `--max-offset` was explicit.

The periodogram record header already contains `length=...`. The standalone
loader could parse it and prefer true record length for the boundary cap, with
the current max-period inference as a fallback.

### 3. Add real-ish regression fixtures

The current tests are useful but synthetic. Before changing defaults, add small
fixtures that mimic the observed TRC_2/TRC_7 failure modes:

- analytical z-score FFT artifact case;
- weak large-record peak near the fixed threshold;
- overlapping harmonic families, including the 7/13 class of bug;
- HOR plus monomer/subrepeat case.

These do not need to be large FASTA files. Synthetic score vectors are enough
for core tests, plus one CLI integration fixture for TSV parsing/output.

### 4. Keep strict local maxima for now, but plan prominence

Strict three-bin maxima are deterministic and simple, but broad or plateaued
peaks can be missed. I would not change this in the same batch as harmonic
matching. A later improvement could add optional prominence-based detection,
possibly after a tiny smoothing window, but only after the current behavior is
covered by regression tests.

### 5. Be explicit about expensive defaults

Default `periodogram` with GUI greyramp currently resolves `--z-score auto` to
empirical z-score. That can be expensive. If routine peak detection moves toward
`signal_mean`, the docs/examples should show the faster path when z-scores are
not needed:

```bash
dottir periodogram input.fa -o p.tsv --z-score off
dottir find-peaks p.tsv -o peaks.tsv
```

or, after item #1:

```bash
dottir periodogram input.fa --z-score off --find-peaks peaks.tsv
```

## Revised Implementation Order

1. Implement #4 as best-parent matching plus post-consolidation reparenting.
   Add overlapping-family regression tests.
2. Change #2 default FFT input to `signal_mean`, update help text and changelog.
3. Add #3 `--auto-threshold`, including per-record threshold diagnostics.
   Validate on TRC_2 before considering any default change.
4. Add #1 minimal inline `--find-peaks <PATH>` using the improved classifier.

## Final Recommendation

Lock in these decisions:

- minimal inline peak output, periodogram-only;
- FFT default becomes `signal_mean`;
- adaptive threshold is opt-in, MAD-based, and reports resolved thresholds;
- harmonic matching gets best-parent selection and a final reparenting pass;
- do not implement the proposed tolerance formula as a "tightening" until the
  current-code equivalence is resolved.

