# The `dottir-gui` interactive viewer

The GUI is built on `egui` / `eframe`. Cross-platform; same binary
on Linux, macOS, and Windows.

## Command line

```text
dottir-gui [OPTIONS] [QUERY] [SUBJECT]
```

Mirrors the original `dotter [options] <horizontal> <vertical>`
invocation. Both positional arguments are optional — run with none to
open an empty window and load FASTAs via **File → Open**.

| Flag | Default | Description |
|------|---------|-------------|
| `-W, --window N` | Karlin/Altschul estimate | Sliding window size. |
| `-z, --zoom N` | `1` | Computation zoom (pixels per matrix block). |
| `-p, --pixel-fac N` | `50` | Multiplier in `min(255, score * pixel_fac / W)`. |
| `--mode {blastn,blastp}` | `blastn` | BLAST mode. |
| `--matrix NAME` | DNA+5/-4 or BLOSUM62 | Built-in score matrix. |
| `--strand {forward,reverse,both}` | `both` | BLASTN strand selection. |
| `--self-comparison` | off | Treat the inputs as a self pair. |
| `-m, --memory-mib N` | `512` | Pixelmap memory cap, MiB. |

## Panels

* **Top menu**: File → Open query / Open subject / Save PNG; View →
  Reset pan/zoom / Reset greyramp / Switch theme / Settings.
* **Central canvas**: textured pixelmap with the crosshair overlay.
* **Right panel — Greyramp**: white/black sliders + Swap/Reset + a
  live LUT strip preview. The LUT is applied on every redraw; the
  underlying pixelmap is not recomputed (spec §4.2.1).
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

## Out of scope (this release)

These were called out in the spec but aren't shipped yet — see the
[ADR index](./adr.md) and `IMPROVEMENTS_PLAN.md`:

* GFF3 / PAF annotation track overlays.
* Recompute-on-zoom-settle for sub-pixel detail.
* Sub-dotter spawn from a rubber-band selection.
* Alignment-view dock.
* Session save/load.
* SVG / PDF export.
