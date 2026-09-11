# The `dottir-gui` interactive viewer

The GUI is built on `egui` / `eframe`. Cross-platform; same binary
on Linux, macOS, and Windows.

## Command line

```text
dottir-gui [OPTIONS] [QUERY] [SUBJECT]
```

Mirrors the original `dotter [options] <horizontal> [<vertical>]`
invocation. Both positional arguments are optional — run with none to
open an empty window and load sequences via **File → Open**. Pass only
one for a self-comparison.

Each positional is a FASTA file, or a GFF3 file with an embedded
`##FASTA` section (see [Sequences from GFF3](#sequences-from-gff3)).

| Flag | Default | Description |
|------|---------|-------------|
| `-W, --window N` | Karlin/Altschul estimate | Sliding window size. |
| `-z, --zoom N` | `1` | Computation zoom (pixels per matrix block). |
| `-p, --pixel-fac N` | `50` | Multiplier in `min(255, score * pixel_fac / W)`. |
| `--mode {blastn,blastp}` | `blastn` | BLAST mode. |
| `--matrix NAME` | DNA+5/-4 or BLOSUM62 | Built-in score matrix. |
| `--strand {forward,reverse,both}` | `both` | BLASTN strand selection. |
| `-m, --memory-mib N` | `512` | Pixelmap memory cap, MiB. |
| `--gff PATH` | — | Annotation overlay (GFF3/BED) applied to **both** axes. `--bed` alias. |
| `--gff-query PATH` | — | Annotation overlay for the query (horizontal) axis. `--bed-query` alias. |
| `--gff-subject PATH` | — | Annotation overlay for the subject (vertical) axis. `--bed-subject` alias. |

Per-axis flags override `--gff`, which in turn overrides features
embedded in a GFF3 positional input. Format is auto-detected; both plain
and gzipped files are accepted.

### Sequences from GFF3

A GFF3 file may carry the sequences it annotates after a `##FASTA`
directive. Such a file is self-contained and can be given directly as a
positional argument — no FASTA, no `--gff*` flag:

```sh
# Self-comparison: sequence and overlay from one file.
dottir-gui elements.gff3

# Pairwise: one annotated file per axis.
dottir-gui family_a.gff3 family_b.gff3
```

The features from the same file are bound as that axis's overlay
automatically. **File → Open query/subject sequence…** accepts the same
form. A GFF3 with no `##FASTA` section is not a sequence input — load it
through the `--gff*` flags or **File → Load GFF3/BED…** instead.

## Panels

* **Top menu**: File → Open query / Open subject / Save PNG; View →
  Reset pan/zoom / Reset greyramp / Switch theme / Settings.
* **Central canvas**: textured pixelmap with the crosshair overlay.
* **Right panel — Greyramp**: white/black sliders + Swap/Reset + a
  live LUT strip preview. The LUT is applied on every redraw; the
  underlying pixelmap is not recomputed (spec §4.2.1).
* **Right panel — Annotations** (when an overlay is loaded): a
  show-bands toggle, an opacity slider, and the color legend. See
  [Annotations](#annotations-gff3--bed) below.
* **Bottom alignment view**: the residue-level alignment around the
  crosshair (forward and, for BLASTN, reverse strand), plus the list of
  annotations overlapping the crosshair on each axis.
* **Bottom status bar**: pixelmap dimensions + window size +
  crosshair coordinates + pixel value. For multi-record FASTAs, the
  coordinates are rendered as `record_id:position` rather than
  bare concatenated offsets.

## Mouse / keyboard

| Input | Action |
|-------|--------|
| Primary-button drag on the canvas | Pan |
| Scroll wheel on the canvas | Zoom on cursor |
| Click on the canvas | Set crosshair |
| Arrow keys | Nudge crosshair by 1 |
| Shift + arrow keys | Nudge by 10 |
| Ctrl + arrow keys | Nudge by 100 |
| `,` / `.` | Step along the main diagonal (both axes ±1). Matches original Dotter. |
| `[` / `]` | Step along the anti-diagonal (q±1, s∓1). Matches original Dotter. |
| `Space` | Snap crosshair to the brightest pixel within ±64 px |

Shift / Ctrl multiply the diagonal and anti-diagonal steps the same
way as the arrow keys (×10 / ×100). The full keymap is also visible
in the GUI under **View → Keyboard shortcuts…**.

## Theme

Light theme is the default — the plotting area renders greyscale on
a near-white background, and a dark surround would muddle axis
labels. **View → Switch to dark theme** flips to egui's dark visuals.

## Settings dialog

Behind **View → Settings…**. Changes recompute the dotplot.

* Mode: Blastn or Blastp.
* Matrix: one of the eight built-ins (protein) or DNA+5/-4
  (nucleotide).
* Window size: explicit value or "auto (Karlin)".
* Zoom (computation): 1..64.
* Pixel factor: 1..255.
* Strand: Forward / Reverse / Both (BLASTN only).
* Self-comparison: triangle Both / Upper / Lower.

## Annotations (GFF3 / BED)

Load interval annotations and overlay them on the dotplot (ADR 0005).
A **query** feature paints as a vertical band spanning the plot height;
a **subject** feature as a horizontal band spanning the width; their
intersection (a region-pair of interest) self-darkens. Dots stay visible
through the bands.

![dottir-gui with an Angela LTR element overlaid: colored annotation
bands, the side-panel legend, and the alignment view listing every
feature under the crosshair](../../dottir_screenshot.png)

**Loading.** Use **File → Load GFF3/BED for query / subject…** (a single
entry in self-comparison), or the `--gff*` flags above — or skip both
and open a GFF3 that carries its own sequences, which binds its features
on load. GFF3 is parsed
via `noodles-gff`, BED as plain intervals; both plain and gzipped.

**Coloring.** GFF3 features are colored by a chosen attribute, selected
under **View → Settings… → Annotations → Color by** (default `Name`,
plus a synthetic `(type)` key for column 3). A feature missing the
chosen attribute **falls back to its GFF3 type** rather than a single
"none" bucket, so e.g. unnamed `long_terminal_repeat` / `primer_binding_site`
children still get distinct colors. BED files use one color per file.

**Right-panel controls.**

* **Show bands** — master on/off.
* **Opacity** — dims every band so the underlying dots show through
  (`0` = invisible).
* **Legend** — one row per value: a show/hide checkbox, a color swatch,
  and the feature count. In self-comparison the single legend drives
  both axes. Per-value **color editing** lives in **Settings →
  Annotations**.

**Alignment view.** When the crosshair is set, the bottom dock lists the
annotations overlapping the crosshair residue for the query (`q:`) and
subject (`s:`) — all of them when features nest — and frames the
crosshair cell in each alignment row in the feature's color.

## Out of scope (this release)

These were called out in the spec but aren't shipped yet — see the
[ADR index](./adr.md) and `docs/IMPROVEMENTS_PLAN.md`:

* PAF / HSP alignment overlays.
* Sub-dotter spawn from a rubber-band selection.
