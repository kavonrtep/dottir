# Per-Square Variation on Tandem-Repeat Self-Comparison

**Status:** root cause confirmed for the reported fixture; proposed fix
is to make auto-zoom divisor-aware for repeat arrays.

## Symptom

When dottir is run on `test-data/tandem_repeat2.fasta` (12 identical
records of `KF293390.1`, 451 bp each) in self-comparison mode, the
resulting plot shows a 12×12 grid of "monomer-vs-monomer" squares. Each
square *should* render identically — every record is byte-equal, so the
local alignment structure in any (i, j) square is the same alignment
structure as any other (i', j') square.

Observed: adjacent squares in the rendered grid show visibly different
diagonal patterns. Some squares look clean; others look "staircase";
some look thicker than others.

Reference: dotter on the same input + same compute zoom shows
uniform-looking squares "at any zoom level" (user observation).

See `tmp/Selection_023.png` (full plot, dottir left vs dotter right)
and `tmp/Selection_024.png` (3×3 zoomed view of the corner).

## What was ruled out

### Greyramp default

Initial hypothesis: switching the default greyramp from `40/100`
(gradient) to `100/100` (hard threshold) collapses the LUT to binary
white/black and exposes per-pixel intensity discretization that the
gradient was previously smoothing. Reverting to `40/100` (dotter's
default) made the diagonals look smoother but did **not** eliminate
the per-square variation — user reports dotter still looks
significantly more uniform even after the revert.

So the greyramp choice contributes to perceived "jaggedness" of
individual lines but is not the cause of the per-square asymmetry.

### Ridge overlay

The vector ridge overlay (`a006c34`, `dottir-core::ridges`) is now
default-off. Disabling it does not change the per-square pattern,
confirming that the asymmetry is in the underlying raster, not the
overlay.

### Kernel math vs. dotter

The dotter inner loop at `third_party/seqtools/dotterApp/dotplot.c:1308`
(`doCalculateImage`) is mathematically identical to dottir's
`emit_pixel` in `crates/dottir-core/src/sliding.rs:234`. Both:

- score every (q_idx, s_idx) residue pair via the sliding window
- compute `dotposq = (qIdx - win2) / zoom`, `dotposs` similarly
- compute per-block-local positions
  `qPosLocal = qIdx - win2 - dotposq * zoom`, `sPosLocal` likewise
- apply the anti-diagonal suppression rule `sPosLocal >= qPosLocal`
  (forward; mirror for reverse)
- max-merge `min(255, score * pixelFac / W)` into `pixelmap[dotposs *
  width + dotposq]`

There is no asymmetry in the rule itself — it depends only on per-block
local positions, which tile perfectly. So **if the kernel produces a
different pixmap than dotter for the same input, the difference is
not in the suppression rule.**

## Two candidate causes

### A. Pixmap-level: non-integer record length in pixmap pixels

For input with R records of N residues each at compute zoom Z,
pixel boundaries are at residue multiples of Z while record boundaries
are at multiples of N. When N is not a multiple of Z, record
boundaries fall *inside* pixmap pixels, and the pixmap cannot tile.

Concretely for this input (N = 451):

- Pixel (qp, sp) is the max-merge over the 22×22 residue block
  `[qp*Z, qp*Z+Z) × [sp*Z, sp*Z+Z)`.
- The "monomer alignment" diagonal in square (i, j) sits at residue
  offset `Δ = (i - j) × 451`.
- In pixmap coords this offset is `Δ / Z`. With Z=22, `451 / 22 = 20.5`
  — non-integer. The diagonal of square (i, j) for |i-j| odd lies on
  *half-integer* pixel rows and must zigzag between two adjacent rows;
  for |i-j| even it lies on integer pixel rows and renders cleanly.

So squares alternate by parity of `(i - j) × 451 mod Z`:

| Z   | 451 mod Z | clean (i-j) values        | aliased (i-j) values |
|-----|-----------|---------------------------|----------------------|
| 11  | 0         | all (Z divides N)         | none                 |
| 22  | 11        | even                      | odd                  |
| 6   | 1         | multiples of 6            | rest                 |
| 7   | 3         | multiples of 7            | rest                 |

**This is mathematically inherent.** Both dotter and dottir have it,
because both kernels are identical. But it predicts *positionally
correlated* variation (clean rows and aliased rows alternate in a
predictable pattern), not random per-square noise.

Verification: at compute_zoom = 1 (one residue per pixel), the
pixmap *must* tile perfectly because no aggregation happens. If
dottir at zoom=1 still shows per-square variation, this hypothesis
is wrong (or incomplete).

### B. Render-pipeline: HiDPI NEAREST upscale phase

Dotter's rendering is a pure `gdk_draw_image` blit at 1:1
(`dotplot.c:2761`, `dotplot.c:3002`). It never resamples.

Dottir uploads the pixmap as an egui texture, then draws it at
*logical* canvas size, and the GPU upscales to physical pixels via
`TextureOptions::NEAREST` (`display_zoom >= 1.0`) or `LINEAR`
(`< 1.0`).

On a HiDPI display with `pixels_per_point = 2`:

- Pixmap size: `pixmap_dim = ceil(seq_len / compute_zoom)`
- Logical canvas: `logical_dim = some value derived from screen size`
- Physical canvas: `physical_dim = logical_dim × 2`
- Texel-to-physical ratio: depends on whether
  `pixmap_dim × ppp` equals `physical_dim` exactly

If `pixmap_dim × ppp ≠ physical_dim`, NEAREST sampling produces a
non-uniform upscale: most texels become a 2×2 block of physical
pixels, but some become 1×2 or 2×1 or 3×2 depending on the
fractional ratio. **The pattern of "which texels get which physical
size" is position-dependent.** A 1-pixel-wide diagonal can render as
two physical pixels in one region and one physical pixel in another,
even though the underlying pixmap pixels are identical.

This is exactly the kind of artifact that disappears in dotter
(which never resamples) and persists in dottir at any zoom level
(because the HiDPI upscale always happens, regardless of compute zoom).

Verification: set `WAYLAND_DISPLAY=` / force `WINIT_X11_SCALE_FACTOR=1`
(or run on a non-HiDPI display), or temporarily hardcode
`pixels_per_point = 1.0` in the egui context, then take the same
screenshot. If the per-square variation disappears at ppp=1 but
returns at ppp=2, the HiDPI upscale is the cause.

## Recommended verification sequence

In order from cheapest / most informative to most invasive:

1. **CLI PNG export.** Run dottir CLI on `tandem_repeat2.fasta`
   self-comparison and write a PNG. CLI has no GPU, no HiDPI upscale —
   it's a CPU raster of the pixmap with the greyramp LUT. Examine the
   PNG at 100 % in an image viewer.
   - If the CLI PNG shows per-square variation: cause is **A** (pixmap-level).
   - If the CLI PNG looks uniform but the GUI shows variation: cause is
     **B** (render-pipeline) or some other GUI-only effect.

2. **GUI at `--zoom 1`.** Force compute zoom = 1 via the CLI flag,
   disable Auto-fit in Settings, then re-test. At zoom=1 the pixmap
   tiles perfectly and cause **A** vanishes. Any remaining variation
   is render-pipeline.

3. **GUI at ppp=1.** Run the GUI on a non-HiDPI display (or with
   pixels-per-point coerced to 1.0). If variation present at ppp=2
   disappears at ppp=1, confirms cause **B**.

4. **Diff against dotter's raw pixmap.** Dotter can save its `.dot`
   binary (format: `dotplot.c:1610` `loadPlot()`). Add a temporary
   reader in dottir-io to parse it, then compute our pixmap at the same
   compute zoom and `diff` byte-by-byte. Byte-identical = our pixmap
   matches dotter's exactly; the remaining difference is purely render.


## Confirmed finding

The additional checks on `test-data/tandem_repeat2.fasta` support the
pixmap-level explanation:

- `zoom=22` produced visibly different monomer squares, and the raw
  `.dot` export showed square-level variation in dark-pixel counts.
- `zoom=11`, `zoom=41`, and `zoom=1` all produced uniform square
  statistics on the same fixture.

That is exactly the signature expected from non-divisor zoom aliasing:
the repeat length is 451 bp, so when the compute zoom does not divide
451, record boundaries and pixel boundaries drift against each other
and identical monomers land on different sub-pixel phases.

This means the reported asymmetry is not a greyramp bug, not a ridge
overlay bug, and not a self-comparison mirror bug. It is a consequence
of how the repeat length interacts with the compute zoom.

## Possible fixes (when root cause is confirmed)

For cause **A** (non-divisor zoom aliasing):

- **Snap auto-fit zoom to a sequence divisor.** When picking the
  display-matched zoom for multi-record input, prefer the nearest
  zoom that divides each record length (or at least the shortest
  record). For 451-bp records this would prefer Z ∈ {1, 11, 41, 451}
  near the canvas-derived target. Cheap, no kernel change. Falls back
  to "any zoom" when no divisor is close enough.
- **CPU downsample with sub-pixel-precise tile alignment.** Compute the
  pixmap at zoom = 1, then CPU-downsample with a kernel aware of
  record boundaries (one downsample window per record, padded
  internally). Tile-perfect by construction. Costs ~Z² more memory at
  compute time.
- **Live with it + document.** It's mathematically inherent to the
  pixmap representation; the user can avoid it by zooming such that
  Z divides the record length.

## Best robust solution

The robust implementation should have two layers:

1. **Divisor-aware auto-zoom** as the default path.
2. **Boundary-aware downsampling** as an exact visual mode for repeat
   arrays where identical units must tile identically regardless of the
   chosen display zoom.

The default should still preserve the scientific pixelmap contract:
dottir-core computes the Dotter-style raster at `PlotConfig::zoom` with
the existing sliding-window, anti-diagonal suppression, and max-merge
rules. The divisor-aware zoom path fixes the common case cheaply by
choosing a zoom that aligns record and pixel boundaries.

Boundary-aware downsampling should be implemented as a rendering/export
stage, not as a replacement for the authoritative pixelmap. It should:

- compute or reuse a higher-resolution source raster, ideally `zoom=1`
  when memory permits;
- split the source into record-aligned tiles using FASTA record
  boundaries;
- downsample each tile in tile-local coordinates so every identical
  monomer-vs-monomer square uses the same sampling phase;
- use a peak-preserving reducer such as max or dark-ink max so sparse
  diagonal evidence is not averaged away;
- keep the result clearly marked as a display raster, not the scientific
  core pixelmap.

Recommended behavior:

- Keep the normal display-matched auto-zoom heuristic.
- When the input is a multi-record repeat array or self-comparison of
  concatenated records, search nearby zooms for divisors of the repeat
  period or, at minimum, the shortest record length.
- Prefer the closest divisor that stays near the requested visual size.
- If no reasonable divisor exists, fall back to boundary-aware
  downsampling when the user enables "repeat-uniform rendering" or when
  the export path explicitly asks for visual normalization.
- Fall back to the existing zoom and raster path when neither strategy
  is applicable.

This keeps the normal path fast and faithful while providing an exact
visual normalization path for tandem-repeat interpretation.

For cause **B** (HiDPI NEAREST upscale phase):

- **Compute the pixmap at *physical* canvas size, not logical.** This
  was the pre-`080d72f` behaviour; that commit deliberately switched
  to logical-sized pixmaps to cut HiDPI memory by ppp². The trade-off
  was "slight NEAREST upscale". If the upscale turns out to be the
  dominant visual artifact, the memory saving is the wrong call.
- **CPU-upscale before texture upload.** Resample the pixmap to
  physical canvas size on the CPU using a uniform NEAREST policy,
  then upload at physical size and disable GPU scaling. Removes the
  position-dependent texel-to-pixel mapping at the cost of a larger
  texture and a CPU resample step.
- **Mipmaps + LINEAR everywhere.** Drop the NEAREST/LINEAR crossover
  and let the GPU do trilinear filtering. Smooths everything but
  blurs sharp 1-pixel diagonals (the dot-plot's primary content).

## Notes for whoever picks this up

- The greyramp default is back to `40/100` (dotter's gradient) as of
  the latest commit — that's the better baseline for any visual A/B.
- The ridge overlay is default-off; flip it on in Settings if you want
  to see what "vector lines from residue evidence" looks like for
  comparison (it sidesteps both **A** and **B** by drawing in
  screen-space rather than sampling the raster).
- `test-data/tandem_repeat2.fasta` and `tmp/Selection_023.png` /
  `tmp/Selection_024.png` are the reference fixtures for this issue.
- The "non-divisor → aliasing" story is testable on synthetic inputs:
  generate N records of any length L, vary L (or vary the input
  canvas size to pick different auto-zooms), and check which
  configurations produce clean grids vs. checkerboard parity.
