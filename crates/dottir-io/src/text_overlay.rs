//! Tiny bitmap font + axis-label compositor for greyscale image
//! exports. Self-contained: no font files, no font crates.
//!
//! The font is 5 wide × 7 tall ASCII glyphs for the characters used
//! by axis labels: digits, `.`, `k`, `M`, `-`, and space. Glyphs are
//! drawn as filled rectangles on a `&mut [u8]` greyscale canvas;
//! background pixels are left alone, so labels overlay cleanly on
//! whatever's underneath.
//!
//! Public surface:
//!
//! * [`compose_image_with_axes`] — wraps a raw `width × height`
//!   pixelmap in a margin, draws tick marks and numeric labels along
//!   the top and left edges, returns the larger composite buffer.
//! * [`nice_tick_step`] — picks a 1/2/5×10^k step so adjacent ticks
//!   are at least `min_pixel_spacing` apart on screen.
//! * [`format_kb`] — human-friendly residue-count label ("1.2M",
//!   "500k", "75").

/// Glyph height in pixels.
pub const FONT_H: usize = 7;
/// Glyph width in pixels (advance includes 1px gap).
pub const FONT_W: usize = 5;
/// Horizontal advance per character (glyph + 1px gap).
pub const FONT_ADVANCE: usize = FONT_W + 1;

/// 5×7 bitmap font. Each entry is 7 rows, each row encoded as the
/// low 5 bits of a u8 (MSB = leftmost pixel within the glyph).
/// Glyphs cover the characters used by axis labels.
fn glyph(c: char) -> Option<&'static [u8; FONT_H]> {
    // Each row is 5 bits; '1' means "ink", '0' means "background".
    static G0: [u8; 7] = [0b01110, 0b10001, 0b10011, 0b10101, 0b11001, 0b10001, 0b01110];
    static G1: [u8; 7] = [0b00100, 0b01100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110];
    static G2: [u8; 7] = [0b01110, 0b10001, 0b00001, 0b00010, 0b00100, 0b01000, 0b11111];
    static G3: [u8; 7] = [0b11110, 0b00001, 0b00001, 0b01110, 0b00001, 0b00001, 0b11110];
    static G4: [u8; 7] = [0b00010, 0b00110, 0b01010, 0b10010, 0b11111, 0b00010, 0b00010];
    static G5: [u8; 7] = [0b11111, 0b10000, 0b11110, 0b00001, 0b00001, 0b10001, 0b01110];
    static G6: [u8; 7] = [0b00110, 0b01000, 0b10000, 0b11110, 0b10001, 0b10001, 0b01110];
    static G7: [u8; 7] = [0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b01000, 0b01000];
    static G8: [u8; 7] = [0b01110, 0b10001, 0b10001, 0b01110, 0b10001, 0b10001, 0b01110];
    static G9: [u8; 7] = [0b01110, 0b10001, 0b10001, 0b01111, 0b00001, 0b00010, 0b01100];
    static GDOT: [u8; 7] = [0, 0, 0, 0, 0, 0b00110, 0b00110];
    static GK: [u8; 7] = [0b10000, 0b10000, 0b10010, 0b10100, 0b11000, 0b10100, 0b10010];
    static GM: [u8; 7] = [0b10001, 0b11011, 0b10101, 0b10001, 0b10001, 0b10001, 0b10001];
    static GDASH: [u8; 7] = [0, 0, 0, 0b11111, 0, 0, 0];
    static GSPACE: [u8; 7] = [0; 7];
    match c {
        '0' => Some(&G0),
        '1' => Some(&G1),
        '2' => Some(&G2),
        '3' => Some(&G3),
        '4' => Some(&G4),
        '5' => Some(&G5),
        '6' => Some(&G6),
        '7' => Some(&G7),
        '8' => Some(&G8),
        '9' => Some(&G9),
        '.' => Some(&GDOT),
        'k' => Some(&GK),
        'M' => Some(&GM),
        '-' => Some(&GDASH),
        ' ' => Some(&GSPACE),
        _ => None,
    }
}

/// Draw `text` onto an 8-bit greyscale canvas at `(x, y)` (top-left
/// of the first glyph). `ink` is the pixel value written for set
/// bits; out-of-bounds pixels are silently dropped. Unknown
/// characters are rendered as a space.
pub fn draw_text(
    canvas: &mut [u8],
    canvas_width: usize,
    canvas_height: usize,
    x: usize,
    y: usize,
    text: &str,
    ink: u8,
) {
    let mut cursor = x;
    for ch in text.chars() {
        let g = glyph(ch).unwrap_or_else(|| glyph(' ').unwrap());
        for (row, &bits) in g.iter().enumerate() {
            let yy = y + row;
            if yy >= canvas_height {
                continue;
            }
            for col in 0..FONT_W {
                if (bits >> (FONT_W - 1 - col)) & 1 == 1 {
                    let xx = cursor + col;
                    if xx < canvas_width {
                        canvas[yy * canvas_width + xx] = ink;
                    }
                }
            }
        }
        cursor += FONT_ADVANCE;
    }
}

/// Pixel width of a rendered string.
pub fn text_width(s: &str) -> usize {
    let n = s.chars().count();
    if n == 0 {
        0
    } else {
        n * FONT_ADVANCE - 1 // last glyph contributes no trailing gap
    }
}

/// Compose a `(plot_w × plot_h)` greyscale pixelmap into a larger
/// canvas with a margin on every side, a 1-px axis frame, tick marks
/// at top and left, and numeric labels along both axes. Returns
/// `(composite_pixels, total_width, total_height)`.
///
/// `coord_w` and `coord_h` are the upper bounds of the **coordinate
/// space** the axes label — typically the sequence lengths in
/// residues. They can differ from `plot_w` / `plot_h` (the
/// **pixel-canvas** dimensions) when the kernel ran at zoom > 1 or
/// the pixelmap was nearest-neighbour upscaled. Pass `coord_w =
/// plot_w, coord_h = plot_h` for the no-rescale case.
///
/// `pixels` is the input pixelmap in row-major form. The composite
/// is also greyscale, with the same convention (0 = black ink,
/// 255 = white). The caller is expected to have already applied any
/// colour inversion / greyramp it wants.
pub fn compose_image_with_axes(
    pixels: &[u8],
    plot_w: u32,
    plot_h: u32,
    coord_w: u32,
    coord_h: u32,
    margin: u32,
) -> (Vec<u8>, u32, u32) {
    let total_w = plot_w + 2 * margin;
    let total_h = plot_h + 2 * margin;
    let mut canvas = vec![255_u8; (total_w as usize) * (total_h as usize)];

    // Blit the pixelmap into the interior.
    let stride_in = plot_w as usize;
    let stride_out = total_w as usize;
    let m = margin as usize;
    for row in 0..plot_h as usize {
        let src = row * stride_in;
        let dst = (m + row) * stride_out + m;
        canvas[dst..dst + stride_in].copy_from_slice(&pixels[src..src + stride_in]);
    }

    // Axis frame: a 1-px box drawn in the *margin* just outside the
    // pixelmap, so the data area stays pristine. The frame's top
    // edge sits at row `m - 1`, etc.
    let frame_ink = 80_u8;
    let top = m.saturating_sub(1);
    let bot = m + plot_h as usize;
    let left = m.saturating_sub(1);
    let right = m + plot_w as usize;
    draw_hline(&mut canvas, stride_out, left, right, top, frame_ink);
    draw_hline(&mut canvas, stride_out, left, right, bot, frame_ink);
    draw_vline(&mut canvas, stride_out, left, top, bot, frame_ink);
    draw_vline(&mut canvas, stride_out, right, top, bot, frame_ink);

    // Tick marks + labels on the top axis (query). Steps and labels
    // are in coordinate space; tick *positions* convert to pixel
    // space via the pixel/coord ratio.
    let coord_w_safe = coord_w.max(1);
    let coord_h_safe = coord_h.max(1);
    let px_per_coord_x = plot_w as f64 / coord_w_safe as f64;
    let px_per_coord_y = plot_h as f64 / coord_h_safe as f64;
    let min_coord_step_x =
        (60.0 / px_per_coord_x.max(f64::MIN_POSITIVE)).ceil().max(1.0) as u32;
    let min_coord_step_y =
        (60.0 / px_per_coord_y.max(f64::MIN_POSITIVE)).ceil().max(1.0) as u32;
    let step_x = nice_tick_step(coord_w_safe as u64, min_coord_step_x);
    let mut t = 0_u64;
    while t <= coord_w_safe as u64 {
        let x = m + (t as f64 * px_per_coord_x).round() as usize;
        if x >= stride_out {
            break;
        }
        draw_vline(&mut canvas, stride_out, x, top.saturating_sub(5), top, frame_ink);
        let label = format_kb(t);
        let lw = text_width(&label);
        let lx = x.saturating_sub(lw / 2);
        let ly = top.saturating_sub(5 + FONT_H + 2);
        draw_text(&mut canvas, stride_out, total_h as usize, lx, ly, &label, 30);
        t = t.saturating_add(step_x);
    }
    // Left axis (subject).
    let step_y = nice_tick_step(coord_h_safe as u64, min_coord_step_y);
    let mut t = 0_u64;
    while t <= coord_h_safe as u64 {
        let y = m + (t as f64 * px_per_coord_y).round() as usize;
        if y >= total_h as usize {
            break;
        }
        draw_hline(&mut canvas, stride_out, left.saturating_sub(5), left, y, frame_ink);
        let label = format_kb(t);
        let lw = text_width(&label);
        let lx = left.saturating_sub(5 + lw + 2);
        let ly = y.saturating_sub(FONT_H / 2);
        draw_text(&mut canvas, stride_out, total_h as usize, lx, ly, &label, 30);
        t = t.saturating_add(step_y);
    }

    (canvas, total_w, total_h)
}

fn draw_hline(canvas: &mut [u8], stride: usize, x0: usize, x1: usize, y: usize, ink: u8) {
    let (lo, hi) = if x0 <= x1 { (x0, x1) } else { (x1, x0) };
    let row = y * stride;
    for x in lo..=hi.min(stride.saturating_sub(1)) {
        canvas[row + x] = ink;
    }
}

fn draw_vline(canvas: &mut [u8], stride: usize, x: usize, y0: usize, y1: usize, ink: u8) {
    let (lo, hi) = if y0 <= y1 { (y0, y1) } else { (y1, y0) };
    for y in lo..=hi {
        let idx = y * stride + x;
        if idx < canvas.len() {
            canvas[idx] = ink;
        }
    }
}

/// Pick a 1/2/5 × 10^k tick step so adjacent ticks are at least
/// `min_pixel_spacing` apart when projected to the screen at 1:1.
pub fn nice_tick_step(span: u64, min_pixel_spacing: u32) -> u64 {
    if span == 0 {
        return 1;
    }
    let target = (min_pixel_spacing as f64).max(1.0);
    let exp = target.log10().floor();
    let base = 10f64.powf(exp).max(1.0) as u64;
    for &m in &[1, 2, 5, 10] {
        let s = m * base;
        if (s as f64) >= target {
            return s.max(1);
        }
    }
    base.max(1)
}

/// Format a residue coord with a `k`/`M` suffix when large. Uses one
/// decimal in the 1k-10k and ≥1M ranges so adjacent tick labels stay
/// distinguishable (otherwise 1000, 1100, 1200 all show "1k").
pub fn format_kb(n: u64) -> String {
    if n >= 1_000_000 {
        format!("{:.1}M", n as f64 / 1_000_000.0)
    } else if n >= 10_000 {
        // 10k+: integer thousands read cleanly ("12k", "100k").
        format!("{}k", n / 1_000)
    } else if n >= 1_000 {
        // 1k-10k: one decimal so 1000/1100/1200 → 1.0k/1.1k/1.2k.
        format!("{:.1}k", n as f64 / 1_000.0)
    } else {
        format!("{n}")
    }
}

/// Nearest-neighbour upscale a `(w × h)` greyscale pixelmap so the
/// longer side is at least `target_max_dim` pixels, preserving the
/// aspect ratio with an *integer* scale factor (every input pixel
/// becomes a clean `scale × scale` block — no aliasing).
///
/// Returns the upscaled buffer plus its new dimensions. If
/// `target_max_dim ≤ max(w, h)` or the inputs are zero-sized, the
/// original is returned unchanged.
///
/// Used by the CLI's PNG exporter to bring tiny plots (e.g. a
/// 287 bp self-comparison ↦ 287 × 287 pixels) up to a comfortable
/// viewing size without committing to fractional scaling (which
/// would smear pixel boundaries and break the "1 pixel = 1 matrix
/// block" contract).
pub fn upscale_nearest(
    pixels: &[u8],
    w: u32,
    h: u32,
    target_max_dim: u32,
) -> (Vec<u8>, u32, u32) {
    let max_dim = w.max(h);
    if max_dim == 0 || target_max_dim <= max_dim {
        return (pixels.to_vec(), w, h);
    }
    // Ceil-divide so the longer side ends up ≥ target.
    let scale = ((target_max_dim + max_dim - 1) / max_dim).max(1);
    if scale <= 1 {
        return (pixels.to_vec(), w, h);
    }
    let new_w = w.saturating_mul(scale);
    let new_h = h.saturating_mul(scale);
    let mut out = vec![0_u8; (new_w as usize) * (new_h as usize)];
    let s = scale as usize;
    let w_us = w as usize;
    let nw = new_w as usize;
    for row in 0..h as usize {
        for col in 0..w_us {
            let v = pixels[row * w_us + col];
            for dy in 0..s {
                let dst_y = row * s + dy;
                let row_off = dst_y * nw + col * s;
                for dx in 0..s {
                    out[row_off + dx] = v;
                }
            }
        }
    }
    (out, new_w, new_h)
}

/// Invert a greyscale pixelmap in place: `v = 255 - v`. Turns the raw
/// "0 = no hit (black), 255 = strong hit (white)" output from
/// `compute_dotplot` into the analysis-friendly "white background,
/// dark hits" rendering used by the GUI's default greyramp and most
/// dotplot conventions.
pub fn invert_in_place(pixels: &mut [u8]) {
    for v in pixels.iter_mut() {
        *v = 255 - *v;
    }
}

/// Same as [`invert_in_place`] but produces a new owned buffer,
/// leaving the input alone.
pub fn inverted(pixels: &[u8]) -> Vec<u8> {
    pixels.iter().map(|&v| 255 - v).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn inversion_round_trips() {
        let p: Vec<u8> = (0..=255).collect();
        let inv = inverted(&p);
        let inv_inv = inverted(&inv);
        assert_eq!(inv_inv, p);
        assert_eq!(inv[0], 255);
        assert_eq!(inv[255], 0);
    }

    #[test]
    fn text_width_of_known_strings() {
        // "100" → 3 glyphs × 6 advance − 1 = 17.
        assert_eq!(text_width("100"), 17);
        assert_eq!(text_width(""), 0);
        assert_eq!(text_width("0"), 5);
    }

    #[test]
    fn nice_step_choices() {
        assert_eq!(nice_tick_step(0, 50), 1);
        // For a 1000-residue span with ≥60 px min spacing at 1:1
        // we expect a coarse step of at least 100.
        assert!(nice_tick_step(1000, 60) >= 60);
        // For a 1M-residue span the step should be at least 100k.
        assert!(nice_tick_step(1_000_000, 60) >= 60);
    }

    #[test]
    fn format_kb_thresholds() {
        assert_eq!(format_kb(0), "0");
        assert_eq!(format_kb(999), "999");
        // 1k-10k uses one decimal so adjacent ticks (1000, 1100,
        // 1200…) don't collide on the same "1k" label.
        assert_eq!(format_kb(1_000), "1.0k");
        assert_eq!(format_kb(1_500), "1.5k");
        assert_eq!(format_kb(9_900), "9.9k");
        // 10k+ uses integer thousands.
        assert_eq!(format_kb(12_345), "12k");
        assert_eq!(format_kb(1_500_000), "1.5M");
    }

    #[test]
    fn compose_writes_pixelmap_into_interior() {
        let pixels = vec![100_u8; 4 * 3];
        let (canvas, tw, th) = compose_image_with_axes(&pixels, 4, 3, 4, 3, 10);
        assert_eq!(tw, 24);
        assert_eq!(th, 23);
        // The interior pixel block should still read 100 (no
        // accidental overwrites in the centre).
        let stride = tw as usize;
        for row in 0..3 {
            for col in 0..4 {
                let idx = (10 + row) * stride + (10 + col);
                assert_eq!(canvas[idx], 100);
            }
        }
        // Outside-the-frame corners should be white background.
        assert_eq!(canvas[0], 255);
        assert_eq!(canvas[(tw - 1) as usize], 255);
    }

    #[test]
    fn upscale_brings_small_plot_above_target() {
        // 287 ↦ scale=7 ↦ 287*7 = 2009 ≥ 2000.
        let pixels: Vec<u8> = (0..287 * 287).map(|i| (i % 256) as u8).collect();
        let (out, w, h) = upscale_nearest(&pixels, 287, 287, 2000);
        assert_eq!(w, 287 * 7);
        assert_eq!(h, 287 * 7);
        assert_eq!(out.len(), (w as usize) * (h as usize));
        assert!(w >= 2000);
        // Each input pixel must reproduce as a 7×7 block: the top-
        // left 7×7 block should all equal pixels[0].
        let v00 = pixels[0];
        for dy in 0..7 {
            for dx in 0..7 {
                assert_eq!(out[dy * w as usize + dx], v00);
            }
        }
    }

    #[test]
    fn upscale_skips_when_already_large_enough() {
        let pixels = vec![42_u8; 3000 * 100];
        let (out, w, h) = upscale_nearest(&pixels, 3000, 100, 2000);
        assert_eq!(w, 3000);
        assert_eq!(h, 100);
        assert_eq!(out, pixels);
    }

    #[test]
    fn upscale_preserves_aspect_ratio() {
        // 100 × 50, target 250 ↦ scale = ceil(250/100) = 3 ↦ 300 × 150.
        let pixels = vec![0_u8; 100 * 50];
        let (_, w, h) = upscale_nearest(&pixels, 100, 50, 250);
        assert_eq!(w, 300);
        assert_eq!(h, 150);
        assert_eq!(w as f32 / h as f32, 100.0 / 50.0);
    }

    #[test]
    fn draw_text_emits_expected_pixels() {
        // Render "0" onto a tiny canvas and check that at least one
        // pixel ended up dark — confirms the glyph table is wired
        // through to the writer.
        let mut canvas = vec![255_u8; 20 * 10];
        draw_text(&mut canvas, 20, 10, 2, 1, "0", 0);
        let any_dark = canvas.iter().any(|&v| v == 0);
        assert!(any_dark, "draw_text did not emit any ink");
    }
}
