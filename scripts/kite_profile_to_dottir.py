#!/usr/bin/env python3
"""Convert one or more kitehor profile TSVs (`<case_id>.kite.tsv` with
columns `d, H, bg`) into a dottir-periodogram-format TSV that
`dottir find-peaks` can consume directly.

Usage:
    kite_profile_to_dottir.py KITE_DIR -o out.tsv [--min-offset N] [--subtract-bg]

Each kite profile is treated as a separate record. The `H` column maps
to dottir's `signal_mean` (and `raw_sum` / `signal_sum`). The `bg`
column can optionally be subtracted to give a noise-floor-corrected
signal via `--subtract-bg`. dottir's `z_score` column is emitted as
`nan` — find-peaks defaults to ranking by `signal_mean` anyway.

Output schema matches `dottir periodogram -o`:

    record_id  k  raw_sum  signal_sum  signal_mean  z_score
"""

import argparse
import sys
from pathlib import Path


def convert(profile_path: Path, min_offset: int, subtract_bg: bool):
    """Yield (record_id, rows) where rows are tuples ready for TSV.

    Pads with zero-rows for any missing `d` values so the output is
    dense (find-peaks requires `k = min_offset + i`).
    """
    record_id = profile_path.stem.replace(".kite", "")
    raw: dict[int, float] = {}
    max_d = -1
    with open(profile_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        if header != ["d", "H", "bg"]:
            raise ValueError(f"unexpected header in {profile_path}: {header}")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            d = int(parts[0])
            h = float(parts[1])
            bg = float(parts[2])
            value = max(0.0, h - bg) if subtract_bg else h
            raw[d] = value
            if d > max_d:
                max_d = d
    rows = []
    for d in range(min_offset, max_d + 1):
        rows.append((record_id, d, raw.get(d, 0.0)))
    return record_id, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("kite_dir", type=Path,
                        help="directory of <case_id>.kite.tsv files (kitehor --dump-profile output)")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="output TSV (default: stdout)")
    parser.add_argument("--min-offset", type=int, default=3,
                        help="smallest d to include (matches dottir periodogram default; default 3)")
    parser.add_argument("--subtract-bg", action="store_true",
                        help="emit max(0, H - bg) instead of raw H — subtracts kite's noise envelope")
    args = parser.parse_args()

    profile_files = sorted(args.kite_dir.glob("*.kite.tsv"))
    if not profile_files:
        sys.exit(f"no *.kite.tsv files in {args.kite_dir}")

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        out.write(f"# converted from kite profiles in {args.kite_dir}\n")
        out.write(f"# min_offset: {args.min_offset}\n")
        out.write(f"# subtract_bg: {args.subtract_bg}\n")
        out.write("record_id\tk\traw_sum\tsignal_sum\tsignal_mean\tz_score\n")
        for path in profile_files:
            record_id, rows = convert(path, args.min_offset, args.subtract_bg)
            for record_id, k, v in rows:
                # raw_sum and signal_sum get same value (integer-ish) as
                # placeholders; signal_mean is what find-peaks ranks by
                # under default --rank-by signal_mean.
                out.write(f"{record_id}\t{k}\t{int(round(v))}\t{int(round(v))}\t{v:.4f}\tnan\n")
    finally:
        if args.output:
            out.close()


if __name__ == "__main__":
    main()
