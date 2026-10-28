#!/bin/sh
# Regenerate the C dotter pixelmap goldens (ADR 0010). Needs `dotter`
# (SeqTools 4.44.1 was used) on PATH; batch mode needs no display.
# Run from the repository root.
set -eu
corpus=tests/corpora/c_dotter
out=tests/golden/c_dotter
grep -v '^#' "$out/cases.tsv" | while IFS="$(printf '\t')" read -r name cflags _dflags query subject; do
    [ "$subject" = "-" ] && subject=$query
    # shellcheck disable=SC2086
    dotter $cflags -b "$out/$name.dot" "$corpus/$query" "$corpus/$subject" >/dev/null 2>&1
    gzip -n -9 -f "$out/$name.dot"
done
dotter --version 2>&1 | head -n 1 > "$out/dotter-version.txt"
