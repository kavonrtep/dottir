#!/usr/bin/env python3
"""Build small dottir test fixtures: extract DANTE_TIR hAT + Mutator elements
plus a flank on each side, and re-coordinate all overlapping GFF3 features into
the extracted window.

Output (under tests/corpora/annotation_overlay/):
  tir_elements.fasta  - one record per window, id "<chrom>_<wstart>-<wend>" (1-based incl.)
  tir_elements.gff3   - every feature from the source GFF3 overlapping a window,
                        coordinates shifted (and clipped) to that window
  tir_elements.manifest.tsv - window <-> source-element mapping

Coordinates stay in genome (+) orientation; feature strands are preserved as-is.
Deterministic: no randomness, stable selection order.
"""
import argparse
import subprocess
import sys
from pathlib import Path

REPO = Path("/home/petr/PycharmProjects/dottir")
GENOME = REPO / "tmp" / "genome_cleaned.fasta"
GFF = REPO / "tmp" / "Repeat_Annotation_Unified.gff3"
OUTDIR = REPO / "tests" / "corpora" / "annotation_overlay"

FLANK = 5000
# pick mid-sized elements so the central element is clearly visible in a dotplot
MIN_LEN, MAX_LEN = 1500, 8000


def contig_lengths(fai: Path) -> dict[str, int]:
    lengths = {}
    for line in fai.read_text().splitlines():
        name, ln = line.split("\t")[:2]
        lengths[name] = int(ln)
    return lengths


def load_dante_tir(gff: Path, kind: str) -> list[dict]:
    """All DANTE_TIR rows whose classification contains `kind`."""
    out = []
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[1] != "DANTE_TIR" or kind not in f[8]:
                continue
            out.append(
                dict(chrom=f[0], start=int(f[3]), end=int(f[4]),
                     strand=f[6], attrs=f[8])
            )
    return out


def select(cands: list[dict], clen: dict[str, int], n: int) -> list[dict]:
    """Mid-sized elements with >=FLANK clearance from both ends, spread across
    contigs, deterministic."""
    ok = [
        c for c in cands
        if MIN_LEN <= (c["end"] - c["start"] + 1) <= MAX_LEN
        and c["start"] - 1 >= FLANK
        and clen[c["chrom"]] - c["end"] >= FLANK
    ]
    ok.sort(key=lambda c: (c["chrom"], c["start"]))
    picked, seen = [], set()
    # first pass: one per contig for diversity
    for c in ok:
        if c["chrom"] not in seen:
            picked.append(c)
            seen.add(c["chrom"])
        if len(picked) == n:
            return picked
    # second pass: fill remainder in order
    for c in ok:
        if c not in picked:
            picked.append(c)
        if len(picked) == n:
            break
    return picked


def index_features(gff: Path) -> dict[str, list[tuple]]:
    """chrom -> list of (start, end, raw_fields) for fast overlap lookup."""
    by_chrom: dict[str, list[tuple]] = {}
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            by_chrom.setdefault(f[0], []).append((int(f[3]), int(f[4]), f))
    return by_chrom


def faidx(region: str) -> str:
    return subprocess.run(
        ["samtools", "faidx", str(GENOME), region],
        check=True, capture_output=True, text=True,
    ).stdout


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prefix", default="tir_elements",
                    help="output basename under tests/corpora/annotation_overlay/")
    ap.add_argument("--hat", type=int, default=5, help="number of hAT windows")
    ap.add_argument("--mutator", type=int, default=5, help="number of Mutator windows")
    ap.add_argument("--dante-tir-only", action="store_true",
                    help="keep only DANTE_TIR features (source_tool=DANTE_TIR) "
                         "in the overlay, instead of every overlapping feature")
    args = ap.parse_args()

    fai = GENOME.with_suffix(GENOME.suffix + ".fai")
    if not fai.exists():
        subprocess.run(["samtools", "faidx", str(GENOME)], check=True)
    clen = contig_lengths(fai)

    picks = (
        [("hAT", c) for c in select(load_dante_tir(GFF, "TIR/hAT"), clen, args.hat)]
        + [("Mutator", c) for c in select(load_dante_tir(GFF, "MuDR_Mutator"), clen, args.mutator)]
    )

    feat_idx = index_features(GFF)
    OUTDIR.mkdir(parents=True, exist_ok=True)

    fasta_lines, gff_lines, manifest = [], ["##gff-version 3"], []
    manifest.append("\t".join(
        ["window_id", "source_chrom", "elem_start", "elem_end", "elem_strand",
         "class", "window_start", "window_end", "n_features"]))

    for kind, c in picks:
        chrom = c["chrom"]
        wstart = c["start"] - FLANK          # 1-based inclusive
        wend = min(c["end"] + FLANK, clen[chrom])
        wid = f"{chrom}_{wstart}-{wend}"

        # sequence
        rec = faidx(f"{chrom}:{wstart}-{wend}").splitlines()
        seq = "".join(rec[1:])
        fasta_lines.append(f">{wid}")
        fasta_lines.extend(seq[i:i + 80] for i in range(0, len(seq), 80))

        # overlapping features, shifted to window (1-based), clipped to window
        nfeat = 0
        for fstart, fend, f in feat_idx.get(chrom, []):
            if fend < wstart or fstart > wend:
                continue
            if args.dante_tir_only and not (
                f[1] == "DANTE_TIR" and "source_tool=DANTE_TIR" in f[8]
            ):
                continue
            ns = max(fstart, wstart) - wstart + 1
            ne = min(fend, wend) - wstart + 1
            row = list(f)
            row[0] = wid
            row[3] = str(ns)
            row[4] = str(ne)
            gff_lines.append("\t".join(row))
            nfeat += 1

        manifest.append("\t".join(map(str, [
            wid, chrom, c["start"], c["end"], c["strand"],
            kind, wstart, wend, nfeat])))
        print(f"{wid}\t{kind}\telem {c['start']}-{c['end']} ({c['strand']})\t"
              f"{wend - wstart + 1} bp\t{nfeat} features")

    (OUTDIR / f"{args.prefix}.fasta").write_text("\n".join(fasta_lines) + "\n")
    (OUTDIR / f"{args.prefix}.gff3").write_text("\n".join(gff_lines) + "\n")
    (OUTDIR / f"{args.prefix}.manifest.tsv").write_text("\n".join(manifest) + "\n")
    print(f"\nWrote fixtures to {OUTDIR} (prefix '{args.prefix}')")
    return 0


if __name__ == "__main__":
    sys.exit(main())
