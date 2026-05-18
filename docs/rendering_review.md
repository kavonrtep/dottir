# Dotplot Rendering Review

This note summarizes the current rendering problem, observations from the
original Dotter implementation, and recommended changes for `dottir`.

## Problem

The core dotplot calculation is faithful to Dotter, but the GUI rendering path
can introduce display artifacts when a large sparse pixelmap is shown in a much
smaller canvas. Typical failure modes are:

- dashed or broken diagonals when sparse one-pixel ridges are averaged away
- stair-step or "scaly" diagonals when peak-preserving pooling is later sampled
  with a visible lattice
- moire-like patterns from regular source diagonals interacting with regular
  downsample grids

These are display artifacts. They should not be fixed by changing the
scientific pixelmap contract in `dottir-core`.

## Current `dottir` Path

The GUI currently treats the computed pixelmap as a texture and applies a
viewport transform over it:

- compute one pixelmap at `PlotConfig::zoom`
- apply the greyramp LUT into an RGBA texture
- optionally max-pool at extreme zoom-out to preserve sparse bright/dark cells
- upload the texture to egui
- draw it into the plot rectangle with GPU texture sampling

Recent experiments improved the situation:

- pool target was made HiDPI-aware
- pooled textures are sampled with `LINEAR` during downsample
- `NEAREST` is kept for high zoom inspection

That is directionally correct, but the remaining artifacts are likely caused by
the non-overlapping integer max-pool grid itself. Max-pooling preserves peaks,
but because its block boundaries are fixed in source pixelmap space, it can
impose a regular lattice on diagonal patterns.

## Original Dotter Observations

Dotter does not appear to solve overview rendering with a higher-quality image
filter. It mostly avoids the problem by making the computed raster match the
display raster.

Relevant code paths:

- `third_party/seqtools-4.28/dotterApp/dotplot.c:724`
  `getImageDimension()` computes the image size as
  `ceil(seq_len / getScaleFactor(...))`.
- `third_party/seqtools-4.28/dotterApp/dotter_.h:218`
  `zoomFactor` is documented as the factor used to scale the dotplot:
  values greater than 1 zoom out, values less than 1 zoom in.
- `third_party/seqtools-4.28/dotterApp/dotplot.c:775`
  `createImage()` allocates a `GdkImage` at the computed image dimensions.
- `third_party/seqtools-4.28/dotterApp/dotplot.c:3002`
  `transformGreyRampImage()` applies the greyramp on CPU directly into the
  `GdkImage`.
- `third_party/seqtools-4.28/dotterApp/dotplot.c:2761`
  `drawImage()` calls `gdk_draw_image()` using the image's native width and
  height.
- `third_party/seqtools-4.28/dotterApp/dotplot.c:912`
  the dotplot widget is put in a scrolled window. The image is not continuously
  scaled to fit the viewport.

The important invariant is:

> Dotter computes the image at the intended displayed scale, then draws it 1:1.

This is different from the current `dottir` approach, where the GUI often
computes a large image and then asks the GPU to downsample it heavily.

Dotter does use vector line drawing, but primarily for overlays such as HSPs,
breaklines, scale marks, and crosshair. The computed dotplot itself remains a
raster image.

## Assessment Of Candidate Fixes

### GPU Sampling Tweaks

Using `LINEAR` for downsampled textures is better than `NEAREST`, and using
mipmaps may help further on the glow backend:

```rust
TextureOptions::LINEAR.with_mipmap_mode(Some(egui::TextureFilter::Linear))
```

This is worth testing because it gives a cheap second-stage smoothing pass after
max-pooling. It should be considered a fallback improvement, not the main
architecture.

### Max-Pool Then Antialias

This is the best texture-only compromise:

1. Max-pool to preserve sparse diagonal evidence.
2. Antialias after that, preferably in display/ink space.

It avoids pure averaging, which can erase one-pixel diagonals, and avoids pure
max-pooling, which can produce blocky lattices.

However, non-overlapping max-pool still has a fixed source-space phase. For
regular diagonal fields, that phase can remain visible.

### Screen-Space Resampling

A better raster renderer would build a texture for the current viewport where
each output pixel samples the exact source footprint it covers. This aligns the
filter with the screen, pan, and zoom instead of with arbitrary source block
boundaries.

This is more work than the current texture cache, but it directly targets the
remaining artifact. It also lets the renderer combine strategies per output
pixel:

- preserve peak ink within the source footprint
- blend local coverage for antialiasing
- avoid depending on a global pool grid

### Vector Ridge Overlay

Vector rendering is promising as an optional enhancement, but it should not
replace the raster dotplot.

The dotplot is a raster score field. Every pixel value can matter, including
weak dots, local clusters, gaps, and noise. Turning it entirely into lines would
require thresholding and joining, which becomes interpretation.

A safer long-term design is a hybrid:

- keep the raster pixelmap as the authoritative base layer
- optionally extract coherent diagonal runs
- render those runs as anti-aliased vector strokes at screen resolution
- only vectorize runs with enough support, for example minimum length, score
  threshold, and limited gap bridging

This would help long biological diagonals remain visually continuous while
leaving ambiguous or noisy detail in the raster layer.

## Recommendation

The primary fix should follow Dotter's invariant:

> Choose the computation zoom so the produced pixelmap is close to the physical
> display size, then render it close to 1:1.

For initial fit-to-window, do not compute a `16000 x 16000` pixelmap and then
downsample it to a `1500 x 1500` canvas. Instead, choose a computation zoom that
directly produces roughly the canvas-sized image.

For a full-plot view, a practical target is:

```text
target_zoom_x = ceil(query_len   / physical_plot_width)
target_zoom_y = ceil(subject_len / physical_plot_height)
target_zoom   = max(target_zoom_x, target_zoom_y, 1)
```

where `physical_plot_width = logical_plot_width * ctx.pixels_per_point()`, and
similarly for height.

When the user zooms or pans:

- keep smooth viewport interaction for immediate feedback
- debounce expensive recomputation
- after the view settles, recompute at a scale where `display_zoom` returns near
  1.0 in physical pixels
- cache recent computed zoom levels to avoid repeating expensive work

This makes texture filtering a preview/fallback path, not the normal way to
create overview quality.

## Suggested Implementation Plan

1. **Initial fit uses display-matched computation zoom**
   - On startup or file load, estimate the physical plot area.
   - Pick `PlotConfig::zoom` from sequence lengths and physical canvas size.
   - Compute that pixelmap first.
   - Render near 1:1.

2. **Revise zoom-settle recompute**
   - Current tiering by powers of two is useful, but it should aim for
     `display_zoom * pixels_per_point()` near 1.0.
   - Allow recompute targets based on actual viewport scale, not only the
     `[0.5, 2.0]` logical threshold.

3. **Keep texture fallback quality improvements**
   - Use `LINEAR` for downsample.
   - Test mipmaps on glow.
   - Snap image rects to physical pixels before painting to reduce subpixel
     phase artifacts.

4. **If artifacts remain, add screen-space raster cache**
   - Build a viewport-sized CPU texture from the visible source footprint.
   - Use peak-preserving plus coverage-aware filtering in display/ink space.
   - This should replace global max-pool for overview rendering.

5. **Optional later enhancement: vector ridge overlay**
   - Extract supported diagonal runs from the pixelmap.
   - Render them as anti-aliased overlay strokes.
   - Keep it user-toggleable and clearly separate from the base raster view.

## Non-Recommendations

- Do not change `dottir-core` pixel emission to make the GUI prettier.
- Do not make overlapping/dilating max-pool the default; it thickens noise and
  can merge nearby features.
- Do not make Lanczos/sinc filtering the first solution. It is more code, may
  ring on sparse high-contrast data, and does not preserve peaks semantically.
- Do not replace the dotplot with vector lines. Use vectorization only as an
  optional ridge enhancement.

## Bottom Line

Original Dotter looks stable because it computes the raster at the displayed
scale and draws it 1:1 in a scrolled widget. `dottir` should adopt that
principle first. GPU filtering, max-pooling, mipmaps, and vector overlays are
secondary tools for preview quality and feature enhancement, not substitutes for
scale-matched computation.
