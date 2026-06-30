#!/usr/bin/env python3
"""Build a dottir LTR test fixture: extract complete Angela (Ty1/copia) elements
from DANTE_LTR output, *including their child features* (LTRs, PBS, TSD, protein
domains), plus a flank on each side, and re-coordinate everything into the
extracted window.

DANTE_LTR GFF3 is hierarchical: a `transposable_element` parent (ID=TE_...) with
`long_terminal_repeat`, `primer_binding_site`, `target_site_duplication`
(source=dante_ltr) and `protein_domain` (source=dante) children that point back
via `Parent=TE_...`. We group by element ID, keep the whole group, and only
rewrite column 1 (seqid) + the coordinates — ID/Parent links are left intact so
the overlay can still reconstruct the hierarchy.

Output (under tests/corpora/annotation_overlay/):
  ltr_angela.fasta         - one record per window, id "<chrom>_<wstart>-<wend>"
  ltr_angela.gff3          - parent + children per window, coords shifted to window
  ltr_angela.manifest.tsv  - window <-> source-element mapping

Coordinates stay in genome (+) orientation; feature strands are preserved as-is.
Deterministic: stable selection order, no randomness.
"""
import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path("/home/petr/PycharmProjects/dottir")
GENOME = REPO / "tmp" / "genome_cleaned.fasta"
GFF = REPO / "tmp" / "DANTE_LTR.gff3"
OUTDIR = REPO / "tests" / "corpora" / "annotation_overlay"

FLANK = 5000
LINEAGE = "Ty1/copia|Angela"
RANK = "DLTP"  # most complete: Domains + LTRs + TSD + PBS

ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")


def contig_lengths(fai: Path) -> dict[str, int]:
    out = {}
    for line in fai.read_text().splitlines():
        name, ln = line.split("\t")[:2]
        out[name] = int(ln)
    return out


def group_elements(gff: Path):
    """Return (parents, children): parents[elem_id] = parent row fields;
    children[elem_id] = list of child rows."""
    parents, children = {}, {}
    with gff.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            attrs = f[8]
            if f[2] == "transposable_element":
                m = ID_RE.search(attrs)
                if m:
                    parents[m.group(1)] = f
            else:
                m = PARENT_RE.search(attrs)
                if m:
                    children.setdefault(m.group(1), []).append(f)
    return parents, children


def select(parents: dict, clen: dict[str, int], n: int) -> list[tuple[str, list]]:
    cands = [
        (eid, f) for eid, f in parents.items()
        if LINEAGE in f[8] and f"Rank={RANK}" in f[8]
        and int(f[3]) - 1 >= FLANK
        and clen[f[0]] - int(f[4]) >= FLANK
    ]
    cands.sort(key=lambda kv: (kv[1][0], int(kv[1][3])))
    picked, seen = [], set()
    for eid, f in cands:                       # one per contig for diversity
        if f[0] not in seen:
            picked.append((eid, f))
            seen.add(f[0])
        if len(picked) == n:
            return picked
    for kv in cands:
        if kv not in picked:
            picked.append(kv)
        if len(picked) == n:
            break
    return picked


def faidx(region: str) -> str:
    return subprocess.run(
        ["samtools", "faidx", str(GENOME), region],
        check=True, capture_output=True, text=True,
    ).stdout


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--prefix", default="ltr_angela")
    ap.add_argument("-n", "--count", type=int, default=3)
    ap.add_argument("--flank", type=int, default=FLANK)
    args = ap.parse_args()
    flank = args.flank

    fai = GENOME.with_suffix(GENOME.suffix + ".fai")
    if not fai.exists():
        subprocess.run(["samtools", "faidx", str(GENOME)], check=True)
    clen = contig_lengths(fai)

    parents, children = group_elements(GFF)
    picks = select(parents, clen, args.count)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fasta_lines, gff_lines = [], ["##gff-version 3"]
    manifest = ["\t".join(
        ["window_id", "source_chrom", "elem_id", "elem_start", "elem_end",
         "elem_strand", "rank", "window_start", "window_end", "n_features"])]

    for eid, pf in picks:
        chrom = pf[0]
        estart, eend, strand = int(pf[3]), int(pf[4]), pf[6]
        wstart = estart - flank                      # 1-based inclusive
        wend = min(eend + flank, clen[chrom])
        wid = f"{chrom}_{wstart}-{wend}"

        rec = faidx(f"{chrom}:{wstart}-{wend}").splitlines()
        seq = "".join(rec[1:])
        fasta_lines.append(f">{wid}")
        fasta_lines.extend(seq[i:i + 80] for i in range(0, len(seq), 80))

        # parent first, then children in genomic order
        group = [pf] + sorted(children.get(eid, []), key=lambda r: int(r[3]))
        for f in group:
            row = list(f)
            row[0] = wid
            row[3] = str(int(f[3]) - wstart + 1)
            row[4] = str(int(f[4]) - wstart + 1)
            gff_lines.append("\t".join(row))

        rank = (re.search(r"Rank=([A-Z]+)", pf[8]) or [None, "?"])[1]
        manifest.append("\t".join(map(str, [
            wid, chrom, eid, estart, eend, strand, rank,
            wstart, wend, len(group)])))
        print(f"{wid}\t{eid}\telem {estart}-{eend} ({strand})\t"
              f"{wend - wstart + 1} bp\t{len(group)} features "
              f"({len(group) - 1} children)")

    (OUTDIR / f"{args.prefix}.fasta").write_text("\n".join(fasta_lines) + "\n")
    (OUTDIR / f"{args.prefix}.gff3").write_text("\n".join(gff_lines) + "\n")
    (OUTDIR / f"{args.prefix}.manifest.tsv").write_text("\n".join(manifest) + "\n")
    print(f"\nWrote fixtures to {OUTDIR} (prefix '{args.prefix}')")
    return 0


if __name__ == "__main__":
    sys.exit(main())
