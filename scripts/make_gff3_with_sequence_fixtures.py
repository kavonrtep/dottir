#!/usr/bin/env python3
"""Generate the `tests/corpora/gff3_with_sequence/` fixtures.

Two kinds of fixture are produced:

* **Derived** — the existing `tests/corpora/annotation_overlay/` pairs
  (`*.gff3` + `*.fasta`) merged into one self-contained GFF3 by appending
  a `##FASTA` section. Byte-identical features and residues, so a test can
  assert that the merged file and the separate pair load to the same thing.
* **Synthetic** — a small deterministic pairwise pair (`pair_a`, `pair_b`)
  with a shared element, cheap enough for an end-to-end `dottir batch` run,
  plus a negative control with no `##FASTA` section.

Deterministic: the synthetic sequences come from a seeded `random.Random`,
so re-running reproduces the checked-in bytes.

Usage: python3 scripts/make_gff3_with_sequence_fixtures.py
"""

from __future__ import annotations

import gzip
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "tests" / "corpora" / "annotation_overlay"
DST = ROOT / "tests" / "corpora" / "gff3_with_sequence"

LINE_WIDTH = 80


def merge(gff: Path, fasta: Path) -> bytes:
    """Concatenate features and sequences into one GFF3."""
    features = gff.read_bytes()
    if not features.endswith(b"\n"):
        features += b"\n"
    return features + b"##FASTA\n" + fasta.read_bytes()


def wrap(seq: str) -> str:
    return "\n".join(seq[i : i + LINE_WIDTH] for i in range(0, len(seq), LINE_WIDTH))


def random_dna(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(rng: random.Random, seq: str, rate: float) -> str:
    out = []
    for base in seq:
        if rng.random() < rate:
            out.append(rng.choice([b for b in "ACGT" if b != base]))
        else:
            out.append(base)
    return "".join(out)


def synthetic_pair() -> tuple[str, str]:
    """Two GFF3-with-FASTA files sharing a diverged element.

    Layout (identical in both files, so the dotplot has an obvious
    off-diagonal block):

        record  flank(150) element(300) flank(150)

    `pair_b`'s element is the `pair_a` element at ~8 % divergence; the
    flanks are independent random sequence.
    """
    rng = random.Random(20260911)
    element = random_dna(rng, 300)

    def build(name: str, element_seq: str) -> str:
        records = []
        features = []
        for i in (1, 2):
            seqid = f"{name}_ctg{i}"
            left = random_dna(rng, 150)
            right = random_dna(rng, 150)
            seq = left + element_seq + right
            records.append((seqid, seq))
            # GFF3 is 1-based inclusive: the element occupies 151..450.
            features.append(
                f"{seqid}\tdottir_fixture\trepeat_region\t151\t450\t.\t+\t."
                f"\tID={name}_elem{i};Name=SharedElement;family=synthetic"
            )
            features.append(
                f"{seqid}\tdottir_fixture\tregion\t1\t150\t.\t+\t."
                f"\tID={name}_flank{i};Name=Flank;family=synthetic"
            )
        header = [
            "##gff-version 3",
            *(f"##sequence-region {seqid} 1 {len(seq)}" for seqid, seq in records),
        ]
        fasta = [f">{seqid}\n{wrap(seq)}" for seqid, seq in records]
        return "\n".join([*header, *features, "##FASTA", *fasta]) + "\n"

    a = build("pair_a", element)
    b = build("pair_b", mutate(rng, element, 0.08))
    return a, b


def no_sequence_fixture() -> str:
    """Negative control: valid GFF3, no `##FASTA` section."""
    return (
        "##gff-version 3\n"
        "##sequence-region pair_a_ctg1 1 600\n"
        "pair_a_ctg1\tdottir_fixture\trepeat_region\t151\t450\t.\t+\t."
        "\tID=pair_a_elem1;Name=SharedElement;family=synthetic\n"
    )


def main() -> None:
    DST.mkdir(parents=True, exist_ok=True)

    for stem in ("tir_simple", "ltr_angela", "tir_elements"):
        gff = SRC / f"{stem}.gff3"
        fasta = SRC / f"{stem}.fasta"
        if not (gff.exists() and fasta.exists()):
            print(f"skip {stem}: source pair missing")
            continue
        merged = merge(gff, fasta)
        if stem == "ltr_angela":
            # One gzipped fixture to exercise the transparent-.gz path.
            out = DST / f"{stem}.with_seq.gff3.gz"
            out.write_bytes(gzip.compress(merged, mtime=0))
        else:
            out = DST / f"{stem}.with_seq.gff3"
            out.write_bytes(merged)
        print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size} bytes)")

    a, b = synthetic_pair()
    (DST / "pair_a.gff3").write_text(a)
    (DST / "pair_b.gff3").write_text(b)
    (DST / "no_sequence.gff3").write_text(no_sequence_fixture())
    for name in ("pair_a.gff3", "pair_b.gff3", "no_sequence.gff3"):
        p = DST / name
        print(f"wrote {p.relative_to(ROOT)} ({p.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
