//! Misleading-artifact regression test.
//!
//! Tandem-repeat dotplots should render as a set of parallel,
//! one-pixel-wide diagonals at every scale. Any rendering policy
//! that introduces wider horizontal runs (the staircase) is
//! producing visual signal that is not in the data — the failure
//! mode documented in `docs/reviews/gui-nearest-upscale-jaggedness.md`.

use dottir_render::{
    fixtures::tandem_repeat, metrics::long_horizontal_runs, render_to_pixels, Greyramp,
    RenderPolicy,
};

/// Sanity check on the fixture itself: the native 1:1 raster
/// already has the expected geometry (no horizontal jogs). If this
/// regresses, every other test in this file becomes meaningless.
#[test]
fn native_1to1_has_no_horizontal_jogs() {
    let plot = tandem_repeat(b"ACGT", 250);
    let g = Greyramp::default();
    let (w, h, img) = render_to_pixels(&plot, &g, RenderPolicy::Native1To1, 0, 0);
    let extra = long_horizontal_runs(&img, w, h, 128, 1);
    assert_eq!(
        extra, 0,
        "native 1:1 raster of a tandem-repeat dotplot should have no row \
         runs >1 px; got {extra}"
    );
}

/// Documents the current bug. The GUI's two-stage-nearest pipeline
/// at a typical HiDPI scale (1.5×) introduces horizontal jogs that
/// the source raster does not contain. This test should be flipped
/// to `assert_eq!(extra, 0, …)` once the GUI is migrated to a
/// non-resampling render policy.
#[test]
fn current_nearest_at_1_5x_introduces_horizontal_jogs() {
    let plot = tandem_repeat(b"ACGT", 250);
    let g = Greyramp::default();
    let pw = (plot.width as f32 * 1.5) as u32;
    let ph = (plot.height as f32 * 1.5) as u32;
    let (w, h, img) = render_to_pixels(&plot, &g, RenderPolicy::CurrentNearest, pw, ph);
    let extra = long_horizontal_runs(&img, w, h, 128, 1);
    assert!(
        extra > 0,
        "expected the current nearest-sample pipeline at 1.5× to reproduce \
         the staircase artifact (run length >1 px); got {extra} long runs. \
         If this assertion fires, the resampling bug may have been fixed — \
         flip the test to assert_eq!(extra, 0)."
    );
}
