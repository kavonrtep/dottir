# Genome comparisons

For whole assemblies — many chromosomes or thousands of contigs, up to a
few gigabases — dottir has a third engine, **genome mode**
([ADR 0008](./adr.md)). Instead of scoring matrix cells it draws the
alignments of an assembly aligner, much like
[D-GENIES](https://dgenies.toulouse.inra.fr/): one line per alignment,
coloured by identity, with contig sorting and assembly exports.

Genome mode is **repeat-blind by construction**: assembly aligners skip
very frequent seeds and keep only the best chains, so transposable
elements, satellites and old duplications mostly do not appear. Use it
for synteny, rearrangements and assembly comparison; use the
[large engine](./large.md) on a chromosome or region to see the repeats.

## Alignments

Genome mode either reads a PAF file or runs an aligner:

```sh
# Run an aligner (mm2plus if installed, else minimap2).
dottir batch reference.fa assembly.fa --engine genome -o genome.png

# Use alignments you already have (any PAF: minimap2, unimap, wfmash…).
dottir batch reference.fa assembly.fa --paf aln.paf -o genome.png
```

The first input is the **horizontal axis (reference)**, the second the
vertical one. A PAF made the other way round is recognised from its record
names and swapped. Inputs may be gzipped; they are decompressed once into
the cache (see below) and indexed (`.fai`), and only names and lengths are
read — the sequences are never loaded.

| Flag | Default | Meaning |
|------|---------|---------|
| `--aligner` | `auto` | `auto` (mm2plus, else minimap2), a name on PATH, or a path |
| `--preset` | `asm20` | minimap2 `-x` preset (`asm5`, `asm10`, `asm20`) |
| `--aligner-args` | — | extra arguments, e.g. `"-c"` for exact identities |
| `--threads` | all cores | aligner threads |
| `--chunk-size` | 10 Mb | query records go to the aligner in pieces this long |
| `--self-all-copies` | off | self-comparison: keep every chain (`-P`) |
| `--no-cache` | off | re-run the aligner even if a cached result exists |

[mm2plus](https://github.com/at-cg/mm2-plus) is an accelerated minimap2
with the same options and output (Bioconda, Linux x86-64 only). On an
8-thread laptop it was not faster than minimap2 on *Arabidopsis*
assemblies; its gains are reported for large genomes on many cores.

**Chunking.** Long query records are cut into `--chunk-size` pieces before
alignment and the pieces of an alignment that crosses a cut are joined
again. Like D-GENIES, this keeps the aligner's memory bounded (Col-CC ×
*A. lyrata*: 1.9 GB chunked, 3.0 GB unchunked, same time).

**Self-comparison** (one input) runs with `-D --dual=no`: the trivial
diagonal is not aligned, so duplicated regions show up as alignments in
their own right; the plot mirrors them and draws the diagonal itself.
`--self-all-copies` adds `-P`, which aligns every copy of every repeat
against every other — the Col-CC centromeres alone give a million
alignments.

**Cache.** Results are kept in `~/.cache/dottir/genome/` (or
`$XDG_CACHE_HOME/dottir`, `%LOCALAPPDATA%\dottir\cache` on Windows,
`$DOTTIR_CACHE_DIR` to override), keyed by both inputs' SHA-256, the
aligner name and version, its arguments and the chunk size. Re-plotting
with other filters or sorting takes well under a second.

`--engine auto` uses genome mode when the first × second input length
reaches `--genome-threshold-cells` (10¹⁷, about 316 Mb × 316 Mb) and an
aligner is installed; otherwise it stays with the large engine and says
why in the log.

## Display

Each alignment is a line, coloured by identity in the D-GENIES classes:
below 25 %, 25–50 %, 50–75 % and above 75 %. Identity is matches / block
length (PAF columns 10 / 11), or 1 − divergence when the PAF has `de:f`
tags (aligner run with `-c`, several times slower). Without `-c` the
identities are approximate and low for diverged genomes; compare colours
within one plot rather than across aligner settings.

| Flag | Default | Meaning |
|------|---------|---------|
| `--min-aln-len` | 0 | hide alignments shorter than this (bp) |
| `--min-identity` | 0 | hide alignments below this identity (0–1) |
| `--primary-only` | off | hide secondary alignments |
| `--remove-noise` | off | the D-GENIES noise filter (below) |
| `--line-width` | 2 | line thickness, pixels |
| `--auto-zoom` | 4000 | largest plot dimension, pixels |

**Noise filter.** With more than 1000 alignments, their lengths are put
into equal-width bins (one per ten alignments); from the most populated
bin onward, the first bin holding fewer than 1/50 of its count sets the
limit, and shorter alignments are hidden.

Runs of five or more consecutive records shorter than 0.2 % of the
assembly are drawn as one shaded block labelled "*n* small".

## Sorting and exports

`--sort-query` sorts the vertical assembly along the reference, as
D-GENIES does: every alignment weighs (1 + length)²; a record goes to the
reference record with the largest total weight, at the weighted mean
position of its alignments there. A record is reverse-complemented when its
long alignments run against the reference. Records without alignments go
last. Reversed records are labelled with a trailing `-`.

The same sort is available for ordinary dotplots of multi-record FASTA
files, with either input sorted: see
[Sorting records](./large.md#sorting-records).

| Flag | Writes |
|------|--------|
| `--assoc-out FILE.tsv` | per record in plot order: best reference record, strand, aligned extents (D-GENIES columns) |
| `--sorted-query-out FILE.fa` | the vertical assembly in plot order, reverse-complemented where shown reversed (header ` dottir:revcomp`) |
| `--unaligned-out PREFIX` | records without any alignment: `PREFIX.reference.fa`, `PREFIX.query.fa` |
| `--paf-out FILE.paf` | the alignments drawn (after filters), query = vertical axis |

```sh
# Order and orient a contig assembly along a reference.
dottir batch reference.fa contigs.fa --engine genome --sort-query \
    --assoc-out contigs.assoc.tsv --sorted-query-out contigs.sorted.fa \
    -o contigs.png
```

The params sidecar has a `[genome]` section with the aligner, its version
and arguments (or the PAF path), the PAF's SHA-256, filters and sort state.

## In the GUI

```sh
dottir-gui reference.fa assembly.fa --genome          # run an aligner
dottir-gui reference.fa assembly.fa --paf aln.paf     # or read a PAF
```

or *File → Open genome comparison…*. Indexing and alignment run in the
background with a progress bar and *Cancel*; cached results open at once.

The layout is the dotplot's: menu bar with the two inputs, controls on
the right, pointer position in the status bar. *Fit* and *Back* (also in
*View*, and *File → Save PNG*) act on the genome view. The right panel
has the same filters as the CLI (applied instantly, no re-alignment),
*Sort query along reference*, line width, the identity legend and exports
(plot PNG, association table, query FASTA in plot order, drawn alignments
as PAF). Dotplot-only menu items are disabled. The plot is redrawn at
screen resolution for every view, so lines stay sharp at any zoom.

| Mouse | |
|-------|-|
| Wheel | zoom at the pointer |
| Drag | pan |
| Shift-drag or middle-drag | zoom into a rectangle |
| Double-click | fit |
| Hover | alignment coordinates, strand, length, identity |
| Right-click | dot-plot the view or the record pair under the pointer; flip a query record |
| Esc / Backspace | back |

**Drill-down.** *Dot-plot this view* (enabled once both axes show at most
150 Mb, adjustable) reads the visible regions from the FASTA indexes and
opens them in the usual dotplot, where the large and exact engines take
over — with the repeats the aligner skipped. Regions keep the orientation
shown in the genome view (flipped records are reverse-complemented), and
the status bar reports positions in the original records. Repeat classes
there are computed from the regions alone. *View → Genome view* returns,
and *View → Dot-plot view* goes back the other way.

Sessions store the comparison (inputs, PAF or aligner settings), filters,
sorting and flips; re-opening a session re-runs the aligner from its
cache.

## Performance

8-thread laptop, *Arabidopsis* assemblies (Col-CC 142 Mb, 5 chromosomes;
*A. lyrata* 207 Mb, 696 scaffolds):

| Comparison | Aligner run | Plot from cache |
|------------|-------------|-----------------|
| Col-CC × *A. lyrata* (`asm20`) | 8.6 s, 1.9 GB | 0.1 s |
| Col-CC self (`asm20`) | 52 s, 4.2 GB | 0.1 s |

Run time and memory depend on which assembly is the reference: the
reference is indexed whole, the other is streamed in chunks.

## Limitations

* DNA only. Records of 2³¹ bp or more cannot be aligned (minimap2).
* No aligner is packaged for Windows; use `--paf` with alignments made
  elsewhere, or `--aligner PATH`.
* Drill-down repeat classes are region-local (a genome-wide depth track
  costs more than the alignment itself and is not computed).
