# Dotter's Plot-Sizing Model — Architectural Note

**Status:** finding from reading the C dotter source. Confirms the user's
observation that dotter's plot area changes size with zoom rather than
fitting a fixed canvas. Relevant context for choosing the GUI rendering
fix described in
[`gui-nearest-upscale-jaggedness.md`](./gui-nearest-upscale-jaggedness.md):
dotter avoids the GPU-resampling problem entirely by never resampling,
because its plot is always its natural size.

## The model in one sentence

**Dotter sizes the pixmap from sequence length and zoom; the window
adapts via scrollbars or whitespace. The image is always drawn 1:1.**

There is no "fit to canvas" code path. There is no resampling. The
notion of "display zoom separate from compute zoom" doesn't exist —
there is one zoom, and it determines both the pixmap size and the
on-screen size simultaneously.

## Source references

1. **Image dimension is fixed by sequence length and zoom**
   (`dotterApp/dotplot.c:724` `getImageDimension`):
   ```c
   int imageLen = (int)ceil((double)seqLen / getScaleFactor(properties, horizontal));
   ```
   Per-axis: `imageWidth = ceil(qlen / zoom)`, `imageHeight = ceil(slen / zoom)`.
   The window dimensions never enter this calculation.

2. **`getScaleFactor` is just `zoomFactor`** for non-BLASTX
   (`dotterApp/dotplot.c:3262`):
   ```c
   result = zoomFactor * resFactor;   // resFactor == 1 for BLASTN/P
   ```
   Same semantics as `PlotConfig::zoom` in dottir-core: residues per
   pixel.

3. **Initial zoom comes from the memory budget, not the window**
   (`dotterApp/dotter.c:474` `getInitZoomFactor`):
   ```c
   result = (int)sqrt((qLen / numFrames / 1e6 * sLen - 1e-6) / memoryLimit) + 1;
   ```
   Picks the smallest zoom that keeps the pixmap under
   `dc->memoryLimit` (MB). If you pass `-z N` explicitly via the CLI,
   that's used instead. **The screen / window size is never an input.**

4. **Hard cap at MAX_IMAGE_DIMENSION = 16000 px** per axis
   (`dotterApp/dotplot.c:53, 779-802`). If the natural dimension
   exceeds 16000 the zoom is bumped up to fit:
   ```c
   if (properties->imageWidth > MAX_IMAGE_DIMENSION) {
       int origLen = properties->imageWidth * dwc->zoomFactor;
       properties->imageWidth = MAX_IMAGE_DIMENSION;
       dwc->zoomFactor = ceil((double)origLen / (double)properties->imageWidth);
       properties->imageWidth = getImageDimension(properties, TRUE);
       properties->imageHeight = getImageDimension(properties, FALSE);
       g_warning("Image too wide - setting zoom to %d\n", (int)dwc->zoomFactor);
   }
   ```
   This is the *only* place dotter ever revises the user's choice of
   zoom — and it's a safety cap on the X11 image, not a fit-to-canvas.

5. **Plot widget lives inside a GtkScrolledWindow**
   (`dotterApp/dotplot.c:932-934`):
   ```c
   GtkWidget *scrollWin = gtk_scrolled_window_new(NULL, NULL);
   gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrollWin), GTK_WIDGET(table));
   gtk_scrolled_window_set_policy(GTK_SCROLLED_WINDOW(scrollWin),
                                  GTK_POLICY_AUTOMATIC, GTK_POLICY_AUTOMATIC);
   ```
   `GTK_POLICY_AUTOMATIC` means: scrollbars when needed, gone when not.
   If the image is bigger than the viewport you scroll; smaller, there's
   whitespace around it.

6. **Interactive zoom change just recomputes**
   (`dotterApp/dotter.c:2275-2292` `onZoomFactorChanged`,
   `dotterApp/dotter.c:2571` zoom entry field):
   ```c
   properties->dotterWinCtx->zoomFactor = newValue;
   /* triggers recompute + relayout */
   ```
   The user types a new zoom into a text entry. The recompute produces
   a new pixmap at the new `ceil(seqLen / zoom)` size, the widget
   resizes to match, and the scrollbars adjust automatically. There
   are no zoom buttons, no scroll-wheel zoom, no rectangle-select
   zoom. The only ways to change the view are:
   - Type a number into the zoom entry.
   - Scroll the scrollbar.

7. **The final blit is pure raster, 1:1** (`dotterApp/dotplot.c:2761`
   `drawImage`, already documented in
   [`tandem-repeat-square-variation.md`](./tandem-repeat-square-variation.md)):
   ```c
   gdk_draw_image(drawable, gc, properties->image,
                  0, 0, properties->plotRect.x, properties->plotRect.y,
                  properties->image->width, properties->image->height);
   ```
   Source and destination dimensions are identical. No sampling.

## Comparison with dottir today

| Dimension | dotter | dottir (current) |
|---|---|---|
| Pixmap size | `ceil(seqLen / zoom)` — fixed by zoom | targets canvas size (display-matched zoom from `854aca9`) |
| What "zoom" means | residues/pixel; sets pixmap size and on-screen size simultaneously | two zooms: `PlotConfig::zoom` (compute) and `display_zoom` (viewport, currently pinned to 1.0) |
| Drawn at | native pixmap size, 1:1 blit | logical canvas size, GPU samples to physical |
| Window vs image | scroll/whitespace if mismatch | always fills the canvas |
| Resample ratio in pixels | always exactly 1.0 | almost never exactly 1.0 (HiDPI, fractional ppp) |
| Interactive navigation | type zoom number, scroll | wheel zoom, rect zoom, pan, Fit, Back, etc. |
| Per-pixel sampling artifacts | none (no sampler) | staircase or blur, depending on regime |
| Memory tied to | sequence length and zoom only | canvas size × `ppp²` |

## Why dotter looks better at every zoom

This is the underlying reason for the user's observation that "dotter
performs better at any zoom level":

- Dotter doesn't have an "any zoom" continuous resample regime.
- Dotter has a finite set of discrete zooms (whatever the user typed),
  and each one is rendered at exactly the pixmap's natural size with a
  1:1 blit.
- Pixel-perfect rendering means every diagonal that's 1 texel wide in
  the pixmap is 1 screen pixel wide on display, end of story.
- The pixmap-level discretization (the "non-divisor zoom aliasing"
  effect documented in
  [`tandem-repeat-square-variation.md`](./tandem-repeat-square-variation.md))
  is still present in dotter's pixmaps, but it stays small because the
  gradient greyramp smooths it visually and there's no GPU resampler
  amplifying it.

In contrast, dottir's "display zoom matched to canvas + GPU sample to
fill" approach guarantees a fractional resample ratio on basically any
input + window combination, and the sampler artifacts swamp the
underlying pixmap quality.

## Implications for the rendering fix

This finding adds a third option to the fix space documented in
[`gui-nearest-upscale-jaggedness.md`](./gui-nearest-upscale-jaggedness.md):

| Option | UX | Implementation cost | Visual quality |
|---|---|---|---|
| **A. CPU-resample pixmap to physical canvas size** (previous review's recommendation) | unchanged — fit-to-canvas with wheel-zoom, rect-zoom | new CPU resampler (nearest-expand for upscale, max-pool for downscale, 1:1 for matched); `ppp²` memory overhead | excellent (1:1 with physical pixels) |
| **B. Compute pixmap at physical canvas size** (revert `080d72f`) | unchanged | one-line revert + `Vec2(pw/ppp, ph/ppp)` draw rect | fixes upscale only; downscale still degrades |
| **C. Dotter-style scrollbox** (new) | different — pixmap is its natural size, scroll to navigate, fewer zooms | drop wheel-zoom-as-display-zoom; rect-zoom recomputes at new zoom (already does); embed plot in egui `ScrollArea`; remove display-zoom entirely | excellent (1:1 always, identical to dotter) |

Option C is structurally the cleanest match to dotter and would
produce visually identical output. It would also remove a lot of
code (display-zoom math, wheel-zoom handler, GPU sampler selection,
max-pool fallback, integer-pixel render path). The cost is a UX
change: continuous-feeling pan/zoom interactions get replaced with
discrete zoom + scrollbar.

Option A keeps the modern UX but requires writing and maintaining
a correct CPU resampler. Option B is the cheapest partial fix and
acceptable as a stepping-stone.

## Recommendation

Make this a deliberate UX decision before writing code. Options A
and C produce visually equivalent output (pixel-perfect 1:1) but
imply different products. C is dottir-faithful; A is what a modern
user expects from a GUI app. Either is defensible. B is a
half-measure and shouldn't be the long-term answer.

Worth prototyping option C in a branch — wrapping the plot canvas in
an `egui::ScrollArea` and seeing how the pan/zoom interactions feel —
before committing to the more invasive option A. If the dotter-style
UX feels acceptable, the whole rendering subsystem gets dramatically
simpler.
