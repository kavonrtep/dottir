#!/usr/bin/env python3
"""find-peaks: extract ranked period peaks from a `dottir periodogram`
TSV or its FFT TSV (`--fft <PATH>` output).

Auto-detects which file format you handed it (by column header). For
periodogram input, finds local-max bins above a z-score threshold.
For FFT input, uses the already-annotated `peak_rank` column.

In both cases the output is then de-duplicated: peaks that look like
integer harmonics of a stronger fundamental are marked as such. By
default only fundamentals are emitted; pass `--show-harmonics` to see
the full harmonic ladder.

Output TSV columns:

    record_id    sequence id (one per FASTA record in the original input)
    period_bp    period of this peak in residues / base-pairs
    score        z_score (periodogram) or amplitude (FFT)
    kind         'fundamental' | 'harmonic'
    parent       period (bp) of the fundamental, when kind == 'harmonic'
    harmonic_n   integer multiple, when kind == 'harmonic'
    n_harmonics  number of harmonics detected for this fundamental
                 (helpful: real tandem repeats have long ladders;
                 isolated peaks have 0)
    harmonics    comma-separated list of (harmonic_n) values seen
                 for this fundamental (e.g. "2,3,4,5,6" — gaps in
                 the ladder are a red flag for partial / local repeats)

Algorithm: greedy strongest-first. The strongest peak is always a
fundamental; subsequent peaks at near-integer multiples of an already-
claimed fundamental are marked as harmonics.

This is a stdlib-only prototype — no numpy / pandas. Once the
heuristics are settled, the same logic will land as a `dottir
find-peaks` Rust subcommand.

Examples:
    scripts/find_peaks.py tmp/test7.tsv
    scripts/find_peaks.py tmp/fft7.tsf --top-n 20
    scripts/find_peaks.py tmp/test7.tsv --z-threshold 50 --show-harmonics
"""

import argparse
import csv
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Format detection + loading
# ---------------------------------------------------------------------------

PERIODOGRAM_COLS = {"record_id", "k", "z_score", "signal_mean"}
FFT_COLS = {"record_id", "bin", "period_residues", "amplitude", "peak_rank"}


def detect_format(path: Path) -> str:
    """Return 'periodogram' or 'fft' by inspecting the column header."""
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = set(line.rstrip("\n").split("\t"))
            if PERIODOGRAM_COLS.issubset(cols):
                return "periodogram"
            if FFT_COLS.issubset(cols):
                return "fft"
            raise ValueError(
                f"unrecognised columns in {path}: {sorted(cols)}\n"
                f"expected periodogram ({sorted(PERIODOGRAM_COLS)}) "
                f"or fft ({sorted(FFT_COLS)})"
            )
    raise ValueError(f"no header row found in {path}")


def load(path: Path) -> dict[str, list[dict]]:
    """Load a TSV, skipping `#` comment lines. Returns
    {record_id: [row_dict, ...]} preserving file order within each
    record."""
    by_record: dict[str, list[dict]] = {}
    with open(path) as f:
        rows = csv.DictReader(
            (line for line in f if not line.startswith("#")), delimiter="\t"
        )
        for row in rows:
            by_record.setdefault(row["record_id"], []).append(row)
    return by_record


# ---------------------------------------------------------------------------
# Peak finding
# ---------------------------------------------------------------------------


def periodogram_peaks(
    rows: list[dict], rank_by: str, min_score: float
) -> list[dict]:
    """Strict three-bin local maxima in the chosen `rank_by` column
    above `min_score`.

    `rank_by` must be a column name in the periodogram TSV — typically
    ``signal_mean`` (recommended default; robust against the analytical
    z-score bias) or ``z_score`` (only reliable when the source used
    ``--z-score empirical``).

    Returns: list of {'period_bp': int, 'score': float, 'signal_mean': float,
    'z_score': float}. `score` is whatever column was used for ranking;
    `signal_mean` / `z_score` are always populated for downstream context.
    """
    peaks = []
    n = len(rows)
    if n < 3:
        return peaks
    scores = [float(r[rank_by]) for r in rows]
    for i in range(1, n - 1):
        if scores[i] < min_score:
            continue
        if scores[i] > scores[i - 1] and scores[i] > scores[i + 1]:
            peaks.append(
                {
                    "period_bp": int(rows[i]["k"]),
                    "score": scores[i],
                    "signal_mean": float(rows[i]["signal_mean"]),
                    "z_score": float(rows[i]["z_score"]),
                }
            )
    return peaks


def find_subrepeats(
    fundamentals: list[dict],
    rows: list[dict],
    rank_by: str,
    subrepeat_min_score: float,
    max_divisor: int,
    period_tolerance: int = 2,
) -> list[dict]:
    """For each fundamental, scan for subrepeats at integer divisors of
    its period.

    Biologically: a tandem array of monomer length `m` organized into
    higher-order repeats (HORs) of `n` monomers shows the strongest
    periodogram signal at the HOR length `n*m` (because adjacent HORs
    are more conserved than adjacent monomers). The monomer itself
    appears at the divisor position `n*m / n = m` — usually with a
    weaker signal than the HOR.

    We scan each fundamental's period `P` at positions `P/n` for
    `n in 2..max_divisor`. Tolerates `period_tolerance` bp of slack
    (the actual peak may not sit at the exact integer divisor when
    monomers have irregular size). Requires the candidate to be a
    strict three-bin local maximum AND above `subrepeat_min_score`.

    Returns: list of {'period_bp', 'score', 'kind'='subrepeat',
    'parent', 'divisor_n'} dicts.
    """
    # Lookup: integer offset -> score in the rank-by column.
    period_score: dict[int, float] = {}
    for r in rows:
        period_score[int(r["k"])] = float(r[rank_by])

    subrepeats = []
    for fund in fundamentals:
        period_f = fund["period_bp"]
        # Only chase divisors of integer periodogram peaks (FFT
        # fundamentals are floats and have different semantics).
        if not isinstance(period_f, int):
            continue
        for n in range(2, max_divisor + 1):
            target = round(period_f / n)
            if target < 3:
                break  # subsequent divisors only smaller; out of range
            # Find the strongest score within ±tolerance of the target;
            # the actual monomer peak can drift a couple of bp.
            best_score = -1.0
            best_period: int | None = None
            for p in range(target - period_tolerance, target + period_tolerance + 1):
                if p < 3 or p not in period_score:
                    continue
                if period_score[p] > best_score:
                    best_score = period_score[p]
                    best_period = p
            if best_period is None or best_score < subrepeat_min_score:
                continue
            # Local-max check at best_period.
            prev_s = period_score.get(best_period - 1, -1e9)
            next_s = period_score.get(best_period + 1, -1e9)
            if best_score <= prev_s or best_score <= next_s:
                continue
            subrepeats.append(
                {
                    "period_bp": best_period,
                    "score": best_score,
                    "kind": "subrepeat",
                    "parent": period_f,
                    "divisor_n": n,
                }
            )
    return subrepeats


def fft_peaks(rows: list[dict], top_n_limit: int | None) -> list[dict]:
    """Pre-annotated peaks from the FFT TSV's `peak_rank` column.

    Returns: list of {'period_bp': float, 'score': float, 'rank': int}.
    """
    peaks = []
    for row in rows:
        rank_s = row["peak_rank"].strip()
        if not rank_s:
            continue
        rank = int(rank_s)
        if top_n_limit is not None and rank > top_n_limit:
            continue
        peaks.append(
            {
                "period_bp": float(row["period_residues"]),
                "score": float(row["amplitude"]),
                "rank": rank,
            }
        )
    return peaks


# ---------------------------------------------------------------------------
# Harmonic dedup
# ---------------------------------------------------------------------------


def classify_harmonics(
    peaks: list[dict],
    fmt: str,
    tolerance: float,
    max_n: int,
) -> list[dict]:
    """Greedy strongest-first classification.

    Each peak that's a near-integer multiple/divisor of an earlier
    fundamental is marked as a harmonic of it. Returns the input list
    (mutated): each dict gains `kind`, `parent`, `harmonic_n` keys.

    `fmt` controls the harmonic direction:
      * periodogram: peak at LARGER period (k) is harmonic of one at
        SMALLER period — `peak.period == n * fundamental.period`
      * fft:         peak at SMALLER period is harmonic of one at
        LARGER period — `fundamental.period == n * peak.period`
    """
    # Sort by score descending — the strongest peak is always a
    # fundamental, and later weaker peaks are checked against the
    # already-claimed set.
    peaks.sort(key=lambda p: -p["score"])
    fundamentals: list[dict] = []
    for p in peaks:
        matched = None
        for f in fundamentals:
            n = _harmonic_n(p["period_bp"], f["period_bp"], fmt, tolerance, max_n)
            if n is not None:
                matched = (f, n)
                break
        if matched is None:
            p["kind"] = "fundamental"
            p["parent"] = None
            p["harmonic_n"] = None
            p["_harmonic_ns"] = []
            fundamentals.append(p)
        else:
            f, n = matched
            p["kind"] = "harmonic"
            p["parent"] = f["period_bp"]
            p["harmonic_n"] = n
            f["_harmonic_ns"].append(n)
    # Annotate fundamentals with harmonic-ladder stats so the caller
    # can filter out spike-only peaks. `n_harmonics` is the count of
    # distinct integer ladder positions seen; `_harmonics_sorted` is
    # the sorted unique ladder for display.
    for f in fundamentals:
        unique_ns = sorted(set(f["_harmonic_ns"]))
        f["n_harmonics"] = len(unique_ns)
        f["_harmonics_sorted"] = unique_ns
    # Consolidation pass: greedy strongest-first can flip when a
    # harmonic's amplitude exceeds its fundamental's (common in FFT
    # where the 6th–9th harmonic of a long-period repeat can beat
    # the fundamental). Fix by checking each fundamental against the
    # others and demoting the one on the wrong side of the
    # biological-convention direction.
    _consolidate_fundamentals(peaks, fmt, tolerance, max_n)
    return peaks


def _consolidate_fundamentals(
    peaks: list[dict],
    fmt: str,
    tolerance: float,
    max_n: int,
) -> None:
    """Demote 'fundamentals' that are actually harmonics of another
    'fundamental' on the wrong side of the convention direction.

    Convention:
      * Periodogram: shorter period is the fundamental.
      * FFT:         longer period is the fundamental.

    Cascades the demoted peak's harmonics to the surviving
    fundamental, multiplying their harmonic_n by the demotion ratio
    (h(n_old) of demoted = h(n_old * demotion_n) of survivor).
    """
    # Iterate until stable — a single pass can leave secondary cases
    # unresolved (e.g., F1 → F2 then F2 → F3).
    while True:
        fundamentals = [p for p in peaks if p["kind"] == "fundamental"]
        demoted = False
        for f in fundamentals:
            survivor = _find_promoting_fundamental(f, fundamentals, fmt, tolerance, max_n)
            if survivor is None:
                continue
            # Compute demotion ratio. For periodogram, demoted (f) has
            # larger period: f.period = n * survivor.period.
            # For FFT, demoted (f) has smaller period:
            # survivor.period = n * f.period.
            if fmt == "periodogram":
                n_demote = round(f["period_bp"] / survivor["period_bp"])
            else:
                n_demote = round(survivor["period_bp"] / f["period_bp"])
            # Demote f → harmonic of survivor.
            f["kind"] = "harmonic"
            f["parent"] = survivor["period_bp"]
            f["harmonic_n"] = n_demote
            f.pop("_harmonics_sorted", None)
            f.pop("n_harmonics", None)
            # Re-attribute peaks that were harmonics of the demoted
            # fundamental: they're now harmonics of `survivor` with
            # multiplied n.
            new_ns_for_survivor = [n_demote]
            for p in peaks:
                if p["kind"] == "harmonic" and p.get("parent") == f["period_bp"]:
                    old_n = p["harmonic_n"]
                    new_n = old_n * n_demote
                    p["parent"] = survivor["period_bp"]
                    p["harmonic_n"] = new_n
                    new_ns_for_survivor.append(new_n)
            survivor["_harmonics_sorted"] = sorted(
                set(survivor["_harmonics_sorted"]) | set(new_ns_for_survivor)
            )
            survivor["n_harmonics"] = len(survivor["_harmonics_sorted"])
            demoted = True
            break  # restart the outer loop with the updated fundamental set
        if not demoted:
            return


def _find_promoting_fundamental(
    f: dict,
    fundamentals: list[dict],
    fmt: str,
    tolerance: float,
    max_n: int,
) -> dict | None:
    """Return another fundamental `g` for which `f` is a harmonic
    (in the biological-convention direction), or None.

    Picks the survivor that gives the cleanest integer ratio if
    multiple candidates match — that's a stand-in for "the real
    fundamental"; concretely: the smallest absolute ratio error.
    """
    best: tuple[float, dict] | None = None
    for g in fundamentals:
        if g is f:
            continue
        # Direction matters: f is demoted if g is the fundamental
        # under our convention.
        if fmt == "periodogram":
            # Periodogram: smaller-period is fundamental. f is demoted
            # only if g.period < f.period.
            if g["period_bp"] >= f["period_bp"]:
                continue
            ratio = f["period_bp"] / g["period_bp"]
        else:  # fft
            # FFT: larger-period is fundamental. f is demoted only if
            # g.period > f.period.
            if g["period_bp"] <= f["period_bp"]:
                continue
            ratio = g["period_bp"] / f["period_bp"]
        n = round(ratio)
        if n < 2 or n > max_n:
            continue
        err = abs(ratio - n) / n
        if err > tolerance:
            continue
        if best is None or err < best[0]:
            best = (err, g)
    return best[1] if best is not None else None


def _harmonic_n(
    peak_period: float,
    fundamental_period: float,
    fmt: str,
    tolerance: float,
    max_n: int,
) -> int | None:
    """If `peak_period` is the nth harmonic of `fundamental_period`,
    return n. Otherwise None.

    Periodogram convention: harmonics live at *larger* k values.
        peak_period / fundamental_period ≈ n  (n in 2..=max_n)
    FFT convention: harmonics live at *smaller* periods.
        fundamental_period / peak_period ≈ n  (n in 2..=max_n)
    """
    if fmt == "periodogram":
        ratio = peak_period / fundamental_period
    else:  # fft
        ratio = fundamental_period / peak_period
    n = round(ratio)
    if n < 2 or n > max_n:
        return None
    if abs(ratio - n) / n > tolerance:
        return None
    return n


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "input",
        type=Path,
        help="periodogram or FFT TSV (auto-detected by column header)",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="output TSV path (default: stdout)",
    )
    parser.add_argument(
        "--rank-by",
        choices=["signal_mean", "z_score"],
        default="signal_mean",
        help=(
            "periodogram: column to use as the peak-detection score. "
            "Default 'signal_mean' — robust against the analytical "
            "z-score bias that uniformly inflates noise across all "
            "offsets. Pass 'z_score' only when the source used "
            "--z-score empirical (clean z) or you want the old "
            "behaviour."
        ),
    )
    parser.add_argument(
        "--min-score",
        type=float,
        default=10.0,
        help=(
            "periodogram: minimum value of the --rank-by column for a "
            "peak (default 10). Suggested values: 10 for signal_mean "
            "(default), 5 for z_score."
        ),
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=None,
        help=(
            "FFT: limit input to peaks with rank <= TOP_N "
            "(default: all annotated). For periodogram, limits the "
            "output (per record) after classification."
        ),
    )
    parser.add_argument(
        "--harmonic-tolerance",
        type=float,
        default=0.02,
        help=(
            "relative tolerance for harmonic-ratio matching (default "
            "0.02 = 2%%). Larger = stricter dedup; smaller = catches "
            "more harmonics."
        ),
    )
    parser.add_argument(
        "--max-harmonic-n",
        type=int,
        default=30,
        help=(
            "largest integer harmonic to consider (default 30). "
            "Long tandem arrays show harmonics out to n>20; the "
            "default covers most natural-sequence cases."
        ),
    )
    parser.add_argument(
        "--show-harmonics",
        action="store_true",
        help="include harmonics in output (default: fundamentals only)",
    )
    parser.add_argument(
        "--subrepeats",
        action="store_true",
        help=(
            "periodogram only: also scan each fundamental's integer "
            "divisors and report any peaks found there as 'subrepeat'. "
            "Useful for spotting the monomer of an HOR family (e.g. "
            "the 174 bp monomer underneath a 524 bp = 3*174 HOR)."
        ),
    )
    parser.add_argument(
        "--subrepeat-min-score",
        type=float,
        default=5.0,
        help=(
            "minimum --rank-by value for a subrepeat peak (default 5). "
            "Subrepeats are typically weaker than their parent HOR, "
            "so this defaults below --min-score."
        ),
    )
    parser.add_argument(
        "--max-divisor",
        type=int,
        default=6,
        help=(
            "largest integer divisor to scan for subrepeats (default 6). "
            "Most natural HORs are 2-5-mers; 6 covers all common cases."
        ),
    )
    parser.add_argument(
        "--min-harmonics",
        type=int,
        default=0,
        help=(
            "fundamentals filter: require at least this many detected "
            "harmonics (default 0). e.g. 3 keeps fundamentals with a "
            "clear ladder (period + 2× + 3× + ...) and drops "
            "single-peak artifacts."
        ),
    )
    args = parser.parse_args()

    fmt = detect_format(args.input)
    by_record = load(args.input)
    out = open(args.output, "w") if args.output else sys.stdout
    try:
        _write_header(out, args, fmt)
        out.write(
            "\t".join(
                [
                    "record_id",
                    "period_bp",
                    "score",
                    "kind",
                    "parent",
                    "harmonic_n",
                    "divisor_n",
                    "n_harmonics",
                    "harmonics",
                ]
            )
            + "\n"
        )
        for rid, rows in by_record.items():
            if fmt == "periodogram":
                peaks = periodogram_peaks(rows, args.rank_by, args.min_score)
            else:
                peaks = fft_peaks(rows, args.top_n)
            classify_harmonics(
                peaks, fmt, args.harmonic_tolerance, args.max_harmonic_n
            )
            # Subrepeat scan — only meaningful for periodogram (integer
            # offsets; FFT peaks are float-period). Runs AFTER the
            # harmonic consolidation so it operates on the final
            # fundamental set.
            subrepeats: list[dict] = []
            if args.subrepeats and fmt == "periodogram":
                fundamentals_for_subrepeats = [
                    p for p in peaks if p["kind"] == "fundamental"
                ]
                subrepeats = find_subrepeats(
                    fundamentals_for_subrepeats,
                    rows,
                    args.rank_by,
                    args.subrepeat_min_score,
                    args.max_divisor,
                )
                # De-dup: drop subrepeat candidates that are already
                # in the peaks list (e.g., as a fundamental in their
                # own right).
                existing_periods = {p["period_bp"] for p in peaks}
                subrepeats = [
                    s for s in subrepeats if s["period_bp"] not in existing_periods
                ]
            # Final filter / ordering.
            if not args.show_harmonics:
                peaks = [p for p in peaks if p["kind"] == "fundamental"]
            if args.min_harmonics > 0:
                peaks = [
                    p
                    for p in peaks
                    if p["kind"] != "fundamental"
                    or p.get("n_harmonics", 0) >= args.min_harmonics
                ]
            peaks.sort(key=lambda p: -p["score"])
            if fmt == "periodogram" and args.top_n is not None:
                peaks = peaks[: args.top_n]
            for p in peaks:
                _emit(out, rid, p, fmt)
                # Subrepeats are emitted directly after their parent
                # so they're easy to scan visually.
                for s in subrepeats:
                    if s["parent"] == p["period_bp"]:
                        _emit(out, rid, s, fmt)
    finally:
        if args.output:
            out.close()
    return 0


def _write_header(out, args, fmt: str) -> None:
    out.write(f"# find_peaks.py (prototype)\n")
    out.write(f"# input: {args.input}\n")
    out.write(f"# input_format: {fmt}\n")
    if fmt == "periodogram":
        out.write(f"# rank_by: {args.rank_by}\n")
        out.write(f"# min_score: {args.min_score}\n")
    else:
        out.write(f"# top_n_input: {args.top_n if args.top_n else 'all-annotated'}\n")
    out.write(f"# harmonic_tolerance: {args.harmonic_tolerance}\n")
    out.write(f"# max_harmonic_n: {args.max_harmonic_n}\n")
    out.write(f"# show_harmonics: {args.show_harmonics}\n")
    if fmt == "periodogram":
        out.write(f"# subrepeats: {args.subrepeats}\n")
        if args.subrepeats:
            out.write(f"# subrepeat_min_score: {args.subrepeat_min_score}\n")
            out.write(f"# max_divisor: {args.max_divisor}\n")


def _emit(out, rid: str, p: dict, fmt: str) -> None:
    period = p["period_bp"]
    period_s = str(period) if isinstance(period, int) else f"{period:.2f}"
    score_s = f"{p['score']:.2f}"
    kind = p["kind"]
    parent_s = ""
    harm_n_s = ""
    div_n_s = ""
    n_harm_s = ""
    harm_list_s = ""
    if kind == "harmonic":
        parent = p["parent"]
        parent_s = str(parent) if isinstance(parent, int) else f"{parent:.2f}"
        harm_n_s = str(p["harmonic_n"])
    elif kind == "subrepeat":
        parent = p["parent"]
        parent_s = str(parent) if isinstance(parent, int) else f"{parent:.2f}"
        div_n_s = str(p["divisor_n"])
    else:  # fundamental
        n_harm_s = str(p.get("n_harmonics", 0))
        harm_list_s = ",".join(str(n) for n in p.get("_harmonics_sorted", []))
    out.write(
        f"{rid}\t{period_s}\t{score_s}\t{kind}\t{parent_s}\t{harm_n_s}\t"
        f"{div_n_s}\t{n_harm_s}\t{harm_list_s}\n"
    )


if __name__ == "__main__":
    sys.exit(main())
