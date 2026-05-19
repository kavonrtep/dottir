//! Pure rendering pipeline for dottir `DotPlot`s.
//!
//! The GUI's render is conceptually:
//!
//!   raw u8 raster → greyramp LUT → resample to physical canvas size
//!
//! Today that last step happens implicitly via egui (texture nearest-
//! sample) + the display compositor (logical→physical scale). On any
//! non-integer `pixels_per_point` this introduces the staircase /
//! false-jog artifact documented in
//! `docs/reviews/gui-nearest-upscale-jaggedness.md`.
//!
//! This crate isolates the pipeline as a pure function so rendering
//! correctness can be measured by unit tests, against synthetic
//! fixtures with known geometric structure, without standing up a GUI
//! or depending on a particular display scale.
//!
//! Each [`RenderPolicy`] variant corresponds to one fix option from
//! `docs/reviews/dotter-sizing-model.md`. The intent is that the same
//! metrics on the same fixtures grade every candidate policy, so the
//! choice between Options A / B / C can be made on numbers rather
//! than on screenshots.

pub mod fixtures;
pub mod metrics;

use dottir_core::DotPlot;

/// Greyramp LUT — raw pixel value (0..=255) → displayed grey
/// (0..=255). Mirrors the GUI's private `Greyramp` struct
/// (`crates/dottir-gui/src/app.rs:163`). Kept inline here for now;
/// once the GUI is migrated to consume this crate the canonical
/// definition should move here and the GUI copy deleted.
#[derive(Clone, Copy, Debug)]
pub struct Greyramp {
    pub white: u8,
    pub black: u8,
    pub swap: bool,
}

impl Default for Greyramp {
    fn default() -> Self {
        Self {
            white: 40,
            black: 100,
            swap: false,
        }
    }
}

impl Greyramp {
    pub fn lut(&self) -> [u8; 256] {
        let mut lut = [0u8; 256];
        let (lo, hi) = if self.white <= self.black {
            (self.white as i32, self.black as i32)
        } else {
            (self.black as i32, self.white as i32)
        };
        for (i, slot) in lut.iter_mut().enumerate() {
            let i = i as i32;
            let v: u8 = if i <= lo {
                if self.swap {
                    0
                } else {
                    255
                }
            } else if i >= hi {
                if self.swap {
                    255
                } else {
                    0
                }
            } else {
                let t = (i - lo) as f32 / (hi - lo).max(1) as f32;
                let g = 255.0 * (1.0 - t);
                if self.swap {
                    (255.0 - g) as u8
                } else {
                    g as u8
                }
            };
            *slot = v;
        }
        lut
    }
}

/// Which rendering policy to apply. Each variant maps to one option
/// in the fix-space discussion in the review docs.
#[derive(Clone, Copy, Debug)]
pub enum RenderPolicy {
    /// What the GUI does today: build a texture at the plot's natural
    /// size, let the display compositor nearest-sample it to the
    /// physical canvas dimensions. This is the policy under test —
    /// it should *fail* the fidelity metrics on any fixture at a
    /// non-integer scale.
    CurrentNearest,
    /// Option C (dotter-style): the rendered image *is* the natural
    /// plot raster; no scaling. `physical_w/h` are ignored and the
    /// returned image has the plot's native dimensions. Should pass
    /// every fidelity metric trivially.
    Native1To1,
    // Future:
    // /// Option B: compute the pixmap at the physical canvas size,
    // /// draw 1:1.
    // ComputeAtPhysical,
    // /// Option A: CPU resample with max-pool on downscale and
    // /// integer-aligned nearest on upscale.
    // CpuMaxPool,
}

/// Render a `DotPlot` to an 8-bit grayscale image.
///
/// Returns `(width, height, pixels)` where `pixels.len() == width *
/// height`. For [`RenderPolicy::Native1To1`] the returned dimensions
/// are `plot.width × plot.height` and `physical_w/h` are ignored.
pub fn render_to_pixels(
    plot: &DotPlot,
    greyramp: &Greyramp,
    policy: RenderPolicy,
    physical_w: u32,
    physical_h: u32,
) -> (u32, u32, Vec<u8>) {
    let lut = greyramp.lut();
    let combined = plot.combined();
    let raster: Vec<u8> = combined.iter().map(|&v| lut[v as usize]).collect();
    let (sw, sh) = (plot.width, plot.height);

    match policy {
        RenderPolicy::Native1To1 => (sw, sh, raster),

        RenderPolicy::CurrentNearest => {
            assert!(physical_w > 0 && physical_h > 0);
            let mut out = vec![255u8; (physical_w as usize) * (physical_h as usize)];
            // Mirrors egui's NEAREST texture sampling under a
            // non-integer display scale: each physical pixel reads
            // the source texel at `floor(x * sw / pw)`.
            for y in 0..physical_h {
                let src_y = ((y as u64 * sh as u64) / physical_h as u64).min(sh as u64 - 1) as u32;
                let src_row = (src_y * sw) as usize;
                let dst_row = (y * physical_w) as usize;
                for x in 0..physical_w {
                    let src_x =
                        ((x as u64 * sw as u64) / physical_w as u64).min(sw as u64 - 1) as u32;
                    out[dst_row + x as usize] = raster[src_row + src_x as usize];
                }
            }
            (physical_w, physical_h, out)
        }
    }
}
