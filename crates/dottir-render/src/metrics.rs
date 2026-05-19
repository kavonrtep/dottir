//! Geometric fidelity metrics for rendered dotplot images.
//!
//! Each metric is a pure function over the 8-bit grayscale output of
//! [`crate::render_to_pixels`]. The intent is that the same metric
//! can be applied to any candidate render policy, so we can grade
//! candidate fixes numerically.

/// For each row, count dark-pixel runs (value `< threshold`) whose
/// length exceeds `max_expected`.
///
/// On a clean diagonal of slope ±1 in the source raster, each row
/// contributes at most one dark pixel, so `max_expected = 1` should
/// yield zero. The previous review's manual measurement
/// (`docs/reviews/gui-nearest-upscale-jaggedness.md` "Measurements")
/// found 93 width-2 and 5 width-3 runs in the GUI screenshot vs. 3
/// width-2 runs in the export — this metric is the automated form
/// of exactly that comparison.
pub fn long_horizontal_runs(
    image: &[u8],
    width: u32,
    height: u32,
    threshold: u8,
    max_expected: u32,
) -> usize {
    assert_eq!(image.len(), (width as usize) * (height as usize));
    let w = width as usize;
    let h = height as usize;
    let mut count = 0;
    for y in 0..h {
        let row = &image[y * w..(y + 1) * w];
        let mut x = 0;
        while x < w {
            if row[x] < threshold {
                let mut run = 1;
                while x + run < w && row[x + run] < threshold {
                    run += 1;
                }
                if run as u32 > max_expected {
                    count += 1;
                }
                x += run;
            } else {
                x += 1;
            }
        }
    }
    count
}
