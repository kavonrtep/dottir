# Large sequences

The exact engine scores every cell of the query × subject matrix. That is
what makes Dotter's pictures trustworthy, but it stops scaling at a few
megabases. For chromosome- and genome-scale DNA comparisons dottir has a
second engine, the **large engine** ([ADR 0007](./adr.md)).

## How it works

1. The subject is indexed with a **spaced seed** (default 12-of-19,
   `1110100110010101111`). Spaced seeds tolerate scattered mismatches much
   better than contiguous k-mers, which matters for diverged repeats.
2. The query is streamed on both strands. A seed hit survives only if a
   second hit lies on the same diagonal within `--two-hit-dist` (100 bp).
   This removes essentially all random hits while keeping similarity down
   to roughly 75–80 % identity.
3. The surviving diagonal runs are scored with **the same sliding-window
   score as the exact engine** (matrix, window size, pixel factor,
   anti-diagonal rule). A large-mode pixel is therefore never darker than
   the exact pixel, and equal on the drawn diagonals; the greyramp means the
   same in both engines.
4. Every dot gets a **repeat class** — `unique`, `low`, `mid`, `high` — from
   the copy number of the bases it connects, so repeats can be shown,
   hidden or coloured without recomputing.
5. Optionally, unique/low-copy runs are **chained** into co-linear synteny
   blocks (coordinates and a score; not alignments).

## When it is used

`dottir batch --engine auto` (the default) switches to the large engine
when query × subject exceeds `--large-threshold-cells` (10¹¹, about
316 kb × 316 kb — roughly a minute of exact compute; a 1 Mb self-plot
takes ~10 min exact versus 0.5 s large). Use `--engine large` or `--engine exact` to force either.
Without `--zoom`/`--auto-zoom`, a large-engine plot is fitted to 4000 px.

In the GUI (*Settings → Large sequences*), `Auto` uses the large engine for
any view with more than 2·10¹⁰ visible cells (about 140 kb × 140 kb) and the
exact engine once you zoom in further, so detail views stay fully
Dotter-faithful. The first large view builds the index (status bar:
*Indexing for large mode…*); later views reuse it.

## Repeat classes

Copy number is estimated per base, within its own sequence:

| `--depth` | What it counts | Use |
|-----------|----------------|-----|
| `seed` (default) | 1 + diagonal runs covering the base in a self-comparison | Diverged repeats (old TE families). One extra seeding pass per sequence. |
| `kmer` | exact canonical k-mer copies (`--depth-k`, default 21) | Faster; sees only young repeats (< ~5 % divergence). |
| `none` | nothing | Classes come only from soft-masking / annotations. |

Thresholds `--class-thresholds 2,5,51` mean: fewer than 2 copies is
`unique`, 2–4 `low`, 5–50 `mid`, 51 and more `high`. A dot takes the higher
class of its two bases. In addition:

* soft-masked (lowercase) bases are at least `mid`
  (`--softmask-class none|low|mid|high`);
* features from `--repeat-annot` / `--repeat-annot-subject` (GFF3 or BED,
  filtered by `--repeat-types LTR,...`) are `high`. A GFF3-with-FASTA input
  supplies its own features when `--repeat-types` is given.

Seeds occurring more than `--max-occ` times (1000) are skipped to bound run
time — typically satellites. They count as `high`, and the number of capped
positions is reported in the log and the params sidecar.

## Recipes

```sh
# Repeat landscape of a chromosome: everything, classes coloured in the GUI.
dottir-gui chr1.fa

# Synteny between two assemblies, repeats hidden, chains exported.
dottir batch asmA.fa asmB.fa -o synteny.png \
    --hide-classes mid,high --chains-out synteny.paf

# Only the repeats, using RepeatMasker output as an extra mask.
dottir batch chr1.fa -o repeats.png \
    --hide-classes unique,low --repeat-annot chr1.rm.gff3

# Export the copy-number track (binned to the zoom) for a genome browser.
dottir batch chr1.fa -o chr1.png --depth-out chr1.depth
```

In the GUI the side panel shows the four classes with checkboxes, the
**Synteny** (unique + low) / **Repeats** (mid + high) / **All** presets,
**Colour by class** (black / blue / orange / red) and **Show synteny
chains** with a minimum-score slider. These are display-only changes.

### Sensitivity

At overview zoom a pixel covers thousands of bases on each axis, so the
best window score inside it is almost always the maximum and the Dotter
greyramp has nothing to work with. The GUI therefore shades large-engine
plots by **matching bp**: the number of diagonal bases with a positive
window score inside the pixel (Gepard-like density). Two sliders set the
ramp, in bases per pixel and on a log scale:

* **white below** — less matching sequence than this is hidden; raise it
  to drop short, old repeat fragments;
* **black at** — this much or more is drawn black.

**Auto-fit each view** (on by default) sets white at the 2nd and black at
the 95th percentile of the lit pixels whenever a new view arrives; moving a
slider turns it off. **Window score** switches back to Dotter shading with
the usual greyramp, which is the right choice once zoomed in to a few
bases per pixel.

### Axis tracks

With **Repeat tracks on axes**, a strip along the top (query) and left
(subject) edge is the profile of the dots you can see: for each pixel
column/row, the matching bp of its visible dots divided by the bases the
column spans — **aligned copies per base**, on a log scale (1 ≈ a tenth
of the strip, 1000 = full). In a self-comparison unique sequence reads 1
(the diagonal) and a 50-copy element about 50. The bar is split by the
dots' class: unique (grey), low (blue), mid (orange), high (red). Hidden
classes and dots below the white threshold do not count, so the strip
changes with the class toggles and the sensitivity sliders.

A red mark on the outer edge flags seeds skipped by `--max-occ`; those
copies are not drawn. Hovering the strip shows the coordinates, the
matching bp and copies per base, the class split and the copy number of
the bases themselves (median and maximum).

## Seed sampling

At overview zoom a pixel spans tens or hundreds of kilobases, and looking
up every seed is mostly wasted work. The large engine therefore looks up
only one seed in *N*, chosen by a hash of the seed itself: every copy of a
sequence carries the same seeds, so a match keeps 1/*N* of its hits instead
of disappearing, and the two-hit distance is scaled by *N* to compensate.
What sampling loses is short or diverged fragments, which may be left with
fewer than two sampled hits; bigger pixels hold many fragments, so few lit
pixels go dark.

| Flag | Default | Meaning |
|------|---------|---------|
| `--seed-sample` | `auto` | plot rate: zoom / 6000, at most 1/32 |
| `--depth-sample` | `auto` | copy-number pass rate: length / 25 Mb, at most 1/16 |

Below 6 kb per pixel nothing is sampled, so zoomed-in views are unaffected.
On the 100 Mb benchmark at 25 kb per pixel, rate 1/4 kept 94 % of the lit
pixels (1/8: 83 %, 1/25: 48 %) and halved the plot time. Pass
`--seed-sample 1 --depth-sample 1` for a full-sensitivity figure. The rates
used are logged and recorded in the params sidecar; in the GUI they are
under *Settings → Large sequences* and shown in the side panel.

## Chains

`--chains-out FILE.paf` writes one PAF line per chain: record-local
coordinates, strand, `tp:A:P`, `cm:i:` (anchors) and `s1:i:` (score). There
is no CIGAR. By default only `unique` and `low` runs are chained
(`--chain-max-class`), so repeats do not glue unrelated blocks together;
`--chain-min-score` (300) drops short chains, and chains with fewer than 3
anchors are dropped too. As in minimap2 (`--mask-level 0.5`), a chain is
not reported when half or more of both its query and subject span lies
inside a better chain of the same record pair and strand; this removes
side chains along small indels of a block while keeping duplicated blocks
elsewhere. For base-level alignments,
run minimap2 and load its PAF.

A name ending in `.bedpe` writes BEDPE instead: query interval, subject
interval (0-based, end-exclusive, source records), `chainN`, score, `+`
and the chain's relative strand. With `--svg` the chains are also drawn
over the image as vector lines (forward green, reverse magenta, mirrored
for a self-comparison). In the GUI, *Save chains…* under *Show synteny
chains* writes the chains above the min-score slider in either format,
and *Save SVG* includes the chains when they are shown.

## Sorting records

Contigs, scaffolds or BAC clones in file order make a scattered plot.
`--sort-records vertical` reorders the second input's records along the
first, and reverse-complements those that run against it, so the plot
becomes one diagonal per reference record; `horizontal` sorts the first
input along the second.

```sh
dottir batch chromosome.fa contigs.fa --sort-records vertical -o sorted.png
```

The rule is the D-GENIES contig sort used in
[genome mode](./genome.md#sorting-and-exports), with the large engine's
synteny chains of the whole plot standing in for alignments (one extra
large-engine pass over the inputs, also when the exact engine draws the
plot). Records without chains go last. It works for both engines, DNA
only, with two inputs.

Only the input order changes, not the scores. Coordinates in the status
bar, the alignment view and annotation overlays stay in the original
records (`ctg12:1,234 (−)` for a reverse-complemented record). The new
order, with reversed records marked `-`, is in the sidecar's
`[record_sort]` section.

In the GUI: *View → Sort records*, or `dottir-gui A.fa B.fa --sort-records
vertical`; *File order* undoes it, and sessions redo it on opening.

On the synthetic contig corpus (42 contigs from a 10 Mb genome, 13
reversed, plus 4 unrelated) every contig is placed on the right
chromosome, in the right order and orientation; sorting took 9 s.

## Parameters

| Flag | Default | Effect |
|------|---------|--------|
| `--engine` | `auto` | `exact`, `large` or `auto` |
| `--large-threshold-cells` | 1e11 | `auto` switch point |
| `--seed-pattern` | `1110100110010101111` | weight 6–13, span ≤ 32 |
| `--two-hit-dist` | 100 | max distance between two hits on one diagonal |
| `--max-occ` | 1000 | skip seeds more frequent than this |
| `--seed-sample` / `--depth-sample` | `auto` | seed sampling (see above) |
| `--extend` | W | cells scored beyond each run |
| `--depth` / `--depth-k` | `seed` / 21 | copy-number method |
| `--class-thresholds` | `2,5,51` | class boundaries |
| `--hide-classes` | — | classes left out of the image |
| `--large-memory-limit-bytes` | 6 GiB | index + depth + channels |

## Performance

On an 8-core desktop with synthetic, repeat-rich genomes:

| Input | Wall time | Peak RSS |
|-------|-----------|----------|
| 10 Mb self | 4–7 s | 0.5 GB |
| 100 Mb self | 57 s | 1.8 GB |
| 100 Mb × 100 Mb | 1 min 22 s | 2.4 GB |
| 600 Mb self (`--large-memory-limit-bytes` 16 GiB) | 3 min 2 s | 10.3 GB |

These use the default automatic seed sampling.

Run time grows with the number of seed hits, i.e. with repeat content. If
a repeat-rich genome is slow, lower `--max-occ` or use a heavier seed.

## Limitations

* DNA only; `-r`/`-v` and `--repeat-uniform-rendering` are not supported.
* The `.dot` export stores the shown classes merged; class channels are not
  preserved. PNG/SVG plus the params sidecar are the reproducible outputs.
* Similarity below roughly 75 % identity is not seeded; zoom in (exact
  engine) to see it.
