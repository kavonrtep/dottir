# The `dottir periodogram` CLI

```text
dottir periodogram <INPUT.fa> [-o OUT.tsv] [--fft FFT.tsv] [--find-peaks PEAKS.tsv] [FLAGS]
```

Computes a per-record self-comparison periodogram from a DNA FASTA,
optionally with the FFT magnitude spectrum and inline peak
classification. One pass over the input; up to three outputs.

A **periodogram** sums the windowed, greyramp-sensitized dotplot
signal along every anti-diagonal at offset `k = s - q`. Each row
of the output corresponds to one offset; the row's value tells you
how much aligned signal exists between residues `i` and `i+k` for
all valid `i`. Periodic structure shows as peaks at the periods of
the repeats and their integer multiples (harmonics).

The FFT and find-peaks outputs are derived from the periodogram and
are designed to be useful in their own right — see below.

## Required arguments

| Argument | Description |
|----------|-------------|
| `INPUT` | DNA FASTA path. Multi-record inputs produce one periodogram block per record in the same output TSV. Gzip auto-detected from magic bytes. |

## Output flags

| Flag | Default | Description |
|------|---------|-------------|
| `-o, --output PATH` | stdout | Periodogram TSV path. |
| `--fft PATH` | off | Opt-in FFT magnitude spectrum (one row per frequency bin per record). |
| `--find-peaks PATH` | off | Opt-in inline peak classification. Uses standalone `dottir find-peaks` defaults; for tuning save the periodogram with `-o` and run `dottir find-peaks` on it. |

## Scoring flags

| Flag | Default | Description |
|------|---------|-------------|
| `--match-score INT` | `5` | DNA match score. |
| `--mismatch-score INT` | `-4` | DNA mismatch score. |
| `-W, --window N` | Karlin/Altschul auto | Sliding window size in residues. Auto-derived per record from Karlin/Altschul statistics (clamped to `[3, 50]`). |
| `--pixel-fac N` | Karlin auto | Multiplier in `min(255, score * pixel_fac / W)`. `0` = auto: `0.2 * 256 / E[M]`. |
| `--greyramp-white 0..255` | `40` | Noise floor — pixels at or below this contribute 0. |
| `--greyramp-black 0..255` | `100` | Saturation — pixels at or above this contribute 255. |
| `--min-offset N` | `3` | Smallest period `k` reported. Skips trivial `k=1, 2` homopolymer spikes. |
| `--max-offset N` | `floor(N/2)` | Largest period reported. |

The greyramp acts exactly like the GUI's noise filter: values below
`white` map to 0 in the periodogram signal, values above `black` map
to 255, between is a linear ramp. Lowering `white` lets through more
faint signal; raising it suppresses noise.

## z-score flags

| Flag | Default | Description |
|------|---------|-------------|
| `--z-score MODE` | `auto` | `auto`, `analytical`, `empirical`, or `off`. |
| `--z-shuffles N` | `200` | Number of shuffles for empirical mode. |
| `--seed N` | `0` | RNG seed for empirical shuffles (`0` → deterministic default). |

* `auto` → `analytical` if greyramp is identity (white=0, black=255), `empirical` otherwise.
* `analytical` → closed-form, fast, but biased (noise floor uniformly inflated).
* `empirical` → per-record shuffle-based null; clean z-scores but `200 × N²` cost.
* `off` → emit `nan` in the z_score column. Recommended for fast survey runs that don't need z-scores (the `--fft` and `--find-peaks` defaults don't need them).

**Recommended for routine analysis:** `--z-score off` plus
`--find-peaks` is the fast path.

## FFT flags

| Flag | Default | Description |
|------|---------|-------------|
| `--fft-input {signal_mean, z_score}` | `signal_mean` | Periodogram column to FFT. `signal_mean` is robust regardless of z-score mode. `z_score` is opt-in and only meaningful with `--z-score empirical`. |
| `--fft-top-peaks N` | `10` | Number of local-maximum bins to mark with a rank in the FFT output's `peak_rank` column. `0` disables annotation. |

## Compute flags

| Flag | Default | Description |
|------|---------|-------------|
| `--threads N` | rayon default (one per CPU) | `1` forces single-threaded. |
| `--memory-limit BYTES` | `1 GiB` | Per-record cap. Streaming algorithm is O(N) memory; this guards against pathological inputs. |

---

## Output: periodogram TSV (`-o`)

One row per `(record, offset)` pair, plus header comments and a
per-record `# record_id=…` block above each record's rows.

| Column | Description |
|--------|-------------|
| `record_id` | FASTA record id (text before first whitespace in the header). |
| `k` | Period in residues (offset between query and subject position). |
| `raw_sum` | `sum_i min(255, score(i, i+k) * pixel_fac / W)` — sum of pixel-fac-scaled raw kernel values along the diagonal at offset `k`. Sensitivity-independent. |
| `signal_sum` | `sum_i sensitivity(raw(i, i+k))` — sum after the greyramp sensitivity ramp. The primary signal. |
| `signal_mean` | `signal_sum / n_pairs(k)` where `n_pairs(k) = N - W - k`. Per-pair average, comparable across `k`. |
| `z_score` | z-score relative to the chosen null (analytical or empirical). `nan` if `--z-score off`. |

### How to read it

* **Peaks at `k = P` and `k = 2P, 3P, …`** → tandem repeat of period
  `P`. The harmonic ladder is the signature.
* **One peak at `k = P` alone (no `2P`, `3P`)** → likely noise, a
  short partial repeat, or a near-record-length artifact.
* **`signal_mean` is the most comparable column** — bins with the
  same biological strength but different `k` show similar
  `signal_mean` values. Use this for ranking.
* **`raw_sum` is biased toward small `k`** because more pairs fit at
  small offsets — `n_pairs(k) = N - W - k` decreases linearly with
  `k`. Use `signal_mean` for ranking unless you specifically want the
  raw sum.
* **`z_score` under `analytical` mode is uniformly inflated** on
  large records (noise floor ~110-140 across all positions). Don't
  use it as a hard threshold unless you ran `--z-score empirical`.

### Per-record metadata header

Just before each record's rows you'll see:

```
# record_id=chr1 length=12345 window=27 pixel_fac=34 freq_A=… freq_C=… freq_G=… freq_T=… null_mean=… null_var=… z=…
```

`window` and `pixel_fac` are the resolved values (after Karlin
auto-derivation if `auto` was requested). The residue frequencies
are what Karlin/Altschul saw; useful when reproducing or
re-analyzing.

---

## Output: FFT TSV (`--fft`)

One row per `(record, bin)`. The spectrum is the FFT of the
`signal_mean` (default) or `z_score` column of the periodogram,
detrended (mean-subtracted), Hann-windowed, zero-padded to the next
power of two.

| Column | Description |
|--------|-------------|
| `record_id` | FASTA record id. |
| `bin` | Frequency bin index, `0..=padded_len/2`. Bin 0 is DC; bin `padded_len/2` is Nyquist. |
| `frequency` | Cycles per residue (`bin / padded_len`). |
| `period_residues` | Period in residues (`padded_len / bin`), `inf` for bin 0. |
| `amplitude` | `\|FFT[bin]\|` magnitude. |
| `peak_rank` | `1` for the brightest non-DC local maximum, `2` for the next, … up to `--fft-top-peaks`. Empty otherwise. |

### How to read it

* **Peak at frequency `f`** ↔ **periodic structure at period
  `1/f` residues** (`= padded_len / bin`).
* The FFT of a periodogram is a *spectrum of a denoised
  autocorrelation* — it inherits the periodogram's window-aware
  scoring and greyramp denoising. Peaks correspond to dominant
  periodicities in the dotplot view.
* For an ideal impulse train at period `P`, the FFT shows peaks at
  `1/P, 2/P, 3/P, …` — the fundamental **plus** its harmonics. The
  top-ranked peak isn't always the fundamental; check the
  ranked-peaks set for the longest period in the family.
* **Frequency resolution = `1 / padded_len`**. For a periodogram
  with `max_offset = 60` the padded length is 64, so adjacent
  bins differ by ~0.5 residues in period near the long end.
  Raise `--max-offset` for finer resolution.

---

## Output: inline find-peaks TSV (`--find-peaks`)

See [`dottir find-peaks`](./find-peaks.md) for the full schema and
options. Inline mode uses defaults (`signal_mean` rank,
`min_score=10`, `min_harmonics=1`, no subrepeats). For tuning, save
the periodogram with `-o` and run `dottir find-peaks` on it.

---

## Common workflows

### Survey a FASTA for tandem repeats (fast path)

```bash
dottir periodogram input.fa --z-score off \
    -o periodogram.tsv \
    --find-peaks peaks.tsv
```

`peaks.tsv` is usually all you need to read. `periodogram.tsv`
is the raw signal in case you want to inspect a specific record
or re-run with different `find-peaks` thresholds.

### With FFT for period-detection cross-check

```bash
dottir periodogram input.fa --z-score off \
    -o periodogram.tsv \
    --fft fft.tsv \
    --find-peaks peaks.tsv
```

The FFT view is most useful as a complement to the periodogram
peaks — strong agreement between the two is high-confidence
periodic structure.

### Full statistical analysis (slow)

```bash
dottir periodogram input.fa --z-score empirical --z-shuffles 200 \
    -o periodogram.tsv \
    --find-peaks peaks.tsv
```

Empirical z-scores are reliable but cost `200 × N²` per record.
For records over ~30 kbp this gets expensive fast; consider
running on a subset of records or with `--threads 1` for
predictable resource use.

### Tuned peak classification on a saved periodogram

```bash
# One-time compute (expensive part):
dottir periodogram input.fa -o periodogram.tsv

# Iterate cheaply on the saved TSV:
dottir find-peaks periodogram.tsv --min-score 5 --min-harmonics 2 --subrepeats
dottir find-peaks periodogram.tsv --auto-threshold --subrepeats
dottir find-peaks periodogram.tsv --rank-by z_score --min-score 50
```

---

## Memory and time notes

* **Periodogram is O(N²) per record, O(N) memory** — the streaming
  algorithm doesn't materialize a pixmap. Records up to a few Mbp
  are tractable; very-large records (>100 kbp) dominate wall time.
* **Records are processed sequentially** within one `dottir
  periodogram` invocation; the kernel parallelizes within each
  record via rayon. Use `--threads 1` for single-threaded
  predictability.
* **`--z-score empirical` is the most expensive option** — `200 × N²`
  per record at the default shuffle count. Skip it for routine
  surveys.
