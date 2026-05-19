# GUI Screenshot Jaggedness vs GUI PNG Export

**Status:** root cause identified from direct comparison of the GUI
screen capture and the GUI PNG export at max zoom.

This revises the earlier emphasis in
[`tandem-repeat-square-variation.md`](./tandem-repeat-square-variation.md).
The non-divisor compute-zoom effect is real and can create small
tile-to-tile differences in the underlying dotplot raster. It does not
explain the remaining strong jaggedness seen on screen after the recent
period-aware zoom changes. The dominant remaining issue is the GUI
display path.

## Evidence

The user saved two images from the same GUI view at max zoom:

| File | Source | Dimensions | Mode | Visual result |
|---|---:|---:|---:|---|
| `tmp/max_zoom.png` | GUI PNG export | 171 x 171 | 8-bit grayscale | clean 1-pixel diagonals |
| `tmp/max_zoom_screenshot.png` | OS screenshot of GUI | 254 x 259 | RGB | jagged/staircased diagonals |

Both are derived from the same `DotPlot` data. The PNG export path maps
`plot.pixels` through the current greyramp and writes the native raster
directly:

- `crates/dottir-gui/src/app.rs:3760` builds the greyramp LUT.
- `crates/dottir-gui/src/app.rs:3768` maps `plot.pixels`.
- `crates/dottir-gui/src/app.rs:3769` writes the PNG at
  `plot.width x plot.height`.

The on-screen path uploads the same mapped raster as an egui texture:

- `crates/dottir-gui/src/app.rs:1131` builds the same greyramp LUT.
- `crates/dottir-gui/src/app.rs:1163` maps `plot.combined()`.
- `crates/dottir-gui/src/app.rs:1169` creates a `ColorImage` at
  `[plot.width, plot.height]`.
- `crates/dottir-gui/src/app.rs:1174` uploads it with
  `TextureOptions::NEAREST`.
- `crates/dottir-gui/src/app.rs:2436` draws it with logical size
  `Vec2::new(pw, ph)`.

The important difference is not the dotplot data. It is the final
logical-point to physical-pixel conversion done by egui/wgpu/windowing.

## Measurements

Image metadata:

```text
tmp/max_zoom.png            171 x 171, L,   min/max 0..255
tmp/max_zoom_screenshot.png 254 x 259, RGB, min/max 0..255
```

Dark-pixel statistics using threshold `< 128`:

```text
export:
  dark pixels: 174
  dark bbox:   x=14..157, y=13..156
  row dark-run lengths: 168 runs of width 1, 3 runs of width 2

screenshot:
  dark pixels: 358
  dark bbox:   x=20..235, y=20..235
  row dark-run lengths: 157 runs of width 1, 93 runs of width 2, 5 runs of width 3
```

The screenshot roughly doubles the number of dark pixels and introduces
many 2-pixel and 3-pixel horizontal runs. That is the visible staircase.

When `tmp/max_zoom.png` is resized to the screenshot dimensions with
nearest-neighbour sampling, it matches the screenshot very closely:

```text
nearest-resized export vs screenshot:
  mean absolute error: 2.41 grey levels
  exact pixel match:   97.24 %
  pixels differing by >32 grey levels: 1.07 %
```

This is decisive: the screen image is essentially a nearest-neighbour
resample of the clean exported raster.

## Mechanism

The GUI code currently treats "1:1" as:

```rust
plot pixel count == egui logical point count
```

That is not the same as:

```text
plot pixel count == physical screen pixel count
```

On a fractional or HiDPI display, egui logical points are scaled by
`ctx.pixels_per_point()`. The max-zoom export is 171 x 171 raster
pixels, but the screenshot is about 254 x 259 physical pixels. The
ratios are approximately:

```text
x: 254 / 171 = 1.485
y: 259 / 171 = 1.515
```

Those are close to a 1.5x display scale, but not exact integer
multiples after window/crop/clipping effects. With nearest sampling,
the physical-pixel mapping becomes position-dependent:

```text
physical pixel x -> source texel floor(x / 1.5)
```

That means a one-texel-wide diagonal is displayed as an alternating
sequence of one-pixel and two-pixel runs. The source raster is clean;
the screen representation is not.

The current comments around `app.rs:1170` and `app.rs:2427` say the
render is 1:1 and therefore avoids GPU resampling. That is only true in
egui logical coordinates. It is false in physical framebuffer
coordinates whenever `pixels_per_point != 1.0`.

## Relationship To The Previous Finding

There are now two distinct effects:

1. **Dotplot raster phase aliasing:** non-divisor compute zoom can put
   identical record tiles on different residue-to-pixel phases. This is
   a real data-raster issue and is what the period-aware zoom snapping
   and boundary-aware downsampling address.
2. **GUI display resampling:** even a clean raster becomes jagged when
   egui draws a 171-pixel texture into about 254-259 physical pixels
   with nearest sampling. This is now the dominant visible issue in the
   max-zoom screenshot/export comparison.

The second effect explains why the GUI export looks like the desired
image while the screenshot still looks jagged. The export bypasses the
logical-to-physical display scale. The screen path does not.

## Best Fix Direction

The robust GUI fix should make the displayed texture land on physical
pixels, not just logical egui points.

Recommended approach:

1. Continue computing the authoritative dotplot raster in core as today.
2. For GUI display, create a screen raster whose dimensions are the
   intended physical framebuffer dimensions for the visible plot area.
3. Resample on CPU using an explicit policy before texture upload.
4. Upload that already-sized raster and draw it so one texture texel
   maps to one physical framebuffer pixel.

For this dotplot use case, the CPU policy should be deterministic and
peak-preserving:

- At or above 1:1: use integer-aligned nearest expansion only when the
  scale is exactly integer; otherwise prefer recomputing/choosing a
  compute zoom that matches the physical viewport.
- When reducing: use max-pool or boundary-aware max-pool so thin
  diagonals are not averaged away.
- Avoid relying on GPU nearest sampling for fractional scales.

Alternative lower-effort fixes:

- Use linear sampling for fractional on-screen scales. This reduces
  staircasing but blurs one-pixel diagonals and is less dotter-faithful.
- Snap the draw rectangle to an integer physical multiple of the source
  raster and accept unused padding/clipping. This avoids fractional
  nearest sampling but may waste plot area.
- Target compute zoom against physical canvas size instead of logical
  canvas size. This is closer to dotter's display model, but still must
  handle cases where the resulting pixmap and physical canvas differ by
  a non-integer ratio.

## Conclusion

The remaining jaggedness is not primarily caused by alignment scoring,
greyramp, window size, or repeat-boundary downsampling. The GUI PNG
export proves that the underlying max-zoom raster can look clean.

The screenshot differs because the GUI display path scales that raster
from 171 x 171 logical texture pixels to roughly 254 x 259 physical
screen pixels. Nearest-neighbour sampling at that fractional scale
creates the staircase pattern.
