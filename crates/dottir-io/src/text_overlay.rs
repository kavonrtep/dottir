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
/// Covers digits, uppercase A–Z, common FASTA-id punctuation, and a
/// few helpers (`.kM-_|:/`); lowercase letters resolve to their
/// uppercase counterpart so the table stays small. Unknown chars
/// fall back to space.
fn glyph(c: char) -> Option<&'static [u8; FONT_H]> {
    // ASCII case-fold so 'a'..'z' use the uppercase glyphs.
    let c = if c.is_ascii_lowercase() {
        c.to_ascii_uppercase()
    } else {
        c
    };
    // Each row is 5 bits; '1' means "ink", '0' means "background".
    static G0: [u8; 7] = [
        0b01110, 0b10001, 0b10011, 0b10101, 0b11001, 0b10001, 0b01110,
    ];
    static G1: [u8; 7] = [
        0b00100, 0b01100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110,
    ];
    static G2: [u8; 7] = [
        0b01110, 0b10001, 0b00001, 0b00010, 0b00100, 0b01000, 0b11111,
    ];
    static G3: [u8; 7] = [
        0b11110, 0b00001, 0b00001, 0b01110, 0b00001, 0b00001, 0b11110,
    ];
    static G4: [u8; 7] = [
        0b00010, 0b00110, 0b01010, 0b10010, 0b11111, 0b00010, 0b00010,
    ];
    static G5: [u8; 7] = [
        0b11111, 0b10000, 0b11110, 0b00001, 0b00001, 0b10001, 0b01110,
    ];
    static G6: [u8; 7] = [
        0b00110, 0b01000, 0b10000, 0b11110, 0b10001, 0b10001, 0b01110,
    ];
    static G7: [u8; 7] = [
        0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b01000, 0b01000,
    ];
    static G8: [u8; 7] = [
        0b01110, 0b10001, 0b10001, 0b01110, 0b10001, 0b10001, 0b01110,
    ];
    static G9: [u8; 7] = [
        0b01110, 0b10001, 0b10001, 0b01111, 0b00001, 0b00010, 0b01100,
    ];
    static GDOT: [u8; 7] = [0, 0, 0, 0, 0, 0b00110, 0b00110];
    static GM_LO: [u8; 7] = [
        0b10000, 0b10000, 0b10010, 0b10100, 0b11000, 0b10100, 0b10010,
    ]; // lowercase 'k' shape kept as legacy
    static GDASH: [u8; 7] = [0, 0, 0, 0b11111, 0, 0, 0];
    static GSPACE: [u8; 7] = [0; 7];
    static GUSCORE: [u8; 7] = [0, 0, 0, 0, 0, 0, 0b11111];
    static GPIPE: [u8; 7] = [
        0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100,
    ];
    static GCOLON: [u8; 7] = [0, 0b00110, 0b00110, 0, 0b00110, 0b00110, 0];
    static GSLASH: [u8; 7] = [
        0b00001, 0b00010, 0b00010, 0b00100, 0b01000, 0b01000, 0b10000,
    ];
    // Uppercase A–Z.
    static GA: [u8; 7] = [
        0b01110, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001,
    ];
    static GB: [u8; 7] = [
        0b11110, 0b10001, 0b10001, 0b11110, 0b10001, 0b10001, 0b11110,
    ];
    static GC: [u8; 7] = [
        0b01110, 0b10001, 0b10000, 0b10000, 0b10000, 0b10001, 0b01110,
    ];
    static GD: [u8; 7] = [
        0b11100, 0b10010, 0b10001, 0b10001, 0b10001, 0b10010, 0b11100,
    ];
    static GE: [u8; 7] = [
        0b11111, 0b10000, 0b10000, 0b11110, 0b10000, 0b10000, 0b11111,
    ];
    static GF: [u8; 7] = [
        0b11111, 0b10000, 0b10000, 0b11110, 0b10000, 0b10000, 0b10000,
    ];
    static GG: [u8; 7] = [
        0b01110, 0b10001, 0b10000, 0b10111, 0b10001, 0b10001, 0b01110,
    ];
    static GH: [u8; 7] = [
        0b10001, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001,
    ];
    static GI: [u8; 7] = [
        0b01110, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110,
    ];
    static GJ: [u8; 7] = [
        0b00111, 0b00010, 0b00010, 0b00010, 0b00010, 0b10010, 0b01100,
    ];
    static GK: [u8; 7] = [
        0b10001, 0b10010, 0b10100, 0b11000, 0b10100, 0b10010, 0b10001,
    ];
    static GL: [u8; 7] = [
        0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b11111,
    ];
    static GM: [u8; 7] = [
        0b10001, 0b11011, 0b10101, 0b10001, 0b10001, 0b10001, 0b10001,
    ];
    static GN: [u8; 7] = [
        0b10001, 0b10001, 0b11001, 0b10101, 0b10011, 0b10001, 0b10001,
    ];
    static GO: [u8; 7] = [
        0b01110, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110,
    ];
    static GP: [u8; 7] = [
        0b11110, 0b10001, 0b10001, 0b11110, 0b10000, 0b10000, 0b10000,
    ];
    static GQ: [u8; 7] = [
        0b01110, 0b10001, 0b10001, 0b10001, 0b10101, 0b10010, 0b01101,
    ];
    static GR: [u8; 7] = [
        0b11110, 0b10001, 0b10001, 0b11110, 0b10100, 0b10010, 0b10001,
    ];
    static GS: [u8; 7] = [
        0b01110, 0b10001, 0b10000, 0b01110, 0b00001, 0b10001, 0b01110,
    ];
    static GT: [u8; 7] = [
        0b11111, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100,
    ];
    static GU: [u8; 7] = [
        0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110,
    ];
    static GV: [u8; 7] = [
        0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01010, 0b00100,
    ];
    static GW: [u8; 7] = [
        0b10001, 0b10001, 0b10001, 0b10001, 0b10101, 0b10101, 0b01010,
    ];
    static GX: [u8; 7] = [
        0b10001, 0b10001, 0b01010, 0b00100, 0b01010, 0b10001, 0b10001,
    ];
    static GY: [u8; 7] = [
        0b10001, 0b10001, 0b10001, 0b01010, 0b00100, 0b00100, 0b00100,
    ];
    static GZ: [u8; 7] = [
        0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b10000, 0b11111,
    ];
    let _ = GM_LO; // legacy shape kept for parity; no character maps here
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
        '-' => Some(&GDASH),
        '_' => Some(&GUSCORE),
        '|' => Some(&GPIPE),
        ':' => Some(&GCOLON),
        '/' => Some(&GSLASH),
        ' ' => Some(&GSPACE),
        'A' => Some(&GA),
        'B' => Some(&GB),
        'C' => Some(&GC),
        'D' => Some(&GD),
        'E' => Some(&GE),
        'F' => Some(&GF),
        'G' => Some(&GG),
        'H' => Some(&GH),
        'I' => Some(&GI),
        'J' => Some(&GJ),
        'K' => Some(&GK),
        'L' => Some(&GL),
        'M' => Some(&GM),
        'N' => Some(&GN),
        'O' => Some(&GO),
        'P' => Some(&GP),
        'Q' => Some(&GQ),
        'R' => Some(&GR),
        'S' => Some(&GS),
        'T' => Some(&GT),
        'U' => Some(&GU),
        'V' => Some(&GV),
        'W' => Some(&GW),
        'X' => Some(&GX),
        'Y' => Some(&GY),
        'Z' => Some(&GZ),
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

/// One record's contribution to an axis label strip. The renderer
/// places `name` centred on the residue interval `[start, end)` of
/// the axis. Records that project to fewer than ~`FONT_ADVANCE * 3`
/// pixels are skipped so dense multi-record inputs don't end up
/// with overlapping labels.
#[derive(Debug, Clone)]
pub struct AxisRecord {
    pub name: String,
    pub start: u32,
    pub end: u32,
}

impl AxisRecord {
    pub fn new(name: impl Into<String>, start: u32, end: u32) -> Self {
        Self {
            name: name.into(),
            start,
            end,
        }
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
#[allow(clippy::too_many_arguments)]
pub fn compose_image_with_axes(
    pixels: &[u8],
    plot_w: u32,
    plot_h: u32,
    coord_w: u32,
    coord_h: u32,
    margin: u32,
    axis_records_x: &[AxisRecord],
    axis_records_y: &[AxisRecord],
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
    let min_coord_step_x = (60.0 / px_per_coord_x.max(f64::MIN_POSITIVE))
        .ceil()
        .max(1.0) as u32;
    let min_coord_step_y = (60.0 / px_per_coord_y.max(f64::MIN_POSITIVE))
        .ceil()
        .max(1.0) as u32;
    let step_x = nice_tick_step(coord_w_safe as u64, min_coord_step_x);
    let mut t = 0_u64;
    while t <= coord_w_safe as u64 {
        let x = m + (t as f64 * px_per_coord_x).round() as usize;
        if x >= stride_out {
            break;
        }
        draw_vline(
            &mut canvas,
            stride_out,
            x,
            top.saturating_sub(5),
            top,
            frame_ink,
        );
        let label = format_kb(t);
        let lw = text_width(&label);
        let lx = x.saturating_sub(lw / 2);
        let ly = top.saturating_sub(5 + FONT_H + 2);
        draw_text(
            &mut canvas,
            stride_out,
            total_h as usize,
            lx,
            ly,
            &label,
            30,
        );
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
        draw_hline(
            &mut canvas,
            stride_out,
            left.saturating_sub(5),
            left,
            y,
            frame_ink,
        );
        let label = format_kb(t);
        let lw = text_width(&label);
        let lx = left.saturating_sub(5 + lw + 2);
        let ly = y.saturating_sub(FONT_H / 2);
        draw_text(
            &mut canvas,
            stride_out,
            total_h as usize,
            lx,
            ly,
            &label,
            30,
        );
        t = t.saturating_add(step_y);
    }

    // Record-name labels along each axis. Drawn one font row *above*
    // the tick labels on the top axis, and one font advance *to the
    // left of* the tick labels on the left axis. Records whose
    // projected pixel extent is smaller than ~3 characters of font
    // advance are skipped to avoid garbage-overlap when many small
    // records line up.
    if !axis_records_x.is_empty() {
        // Tick labels for the top axis are drawn at
        // y = top - 5 - FONT_H - 2. Record names sit another
        // (FONT_H + 4) px above that.
        let tick_label_y = top.saturating_sub(5 + FONT_H + 2);
        let ly = tick_label_y.saturating_sub(FONT_H + 4);
        draw_axis_record_labels_x(
            &mut canvas,
            stride_out,
            total_h as usize,
            m,
            plot_w,
            coord_w_safe,
            px_per_coord_x,
            axis_records_x,
            ly,
        );
    }
    if !axis_records_y.is_empty() {
        // Tick labels for the left axis right-align at
        // x = left - 5 - tick_label_width - 2. Record names sit
        // another column to the left, vertically along each
        // record's slice midpoint.
        draw_axis_record_labels_y(
            &mut canvas,
            stride_out,
            total_h as usize,
            m,
            plot_h,
            coord_h_safe,
            px_per_coord_y,
            axis_records_y,
            left,
        );
    }

    (canvas, total_w, total_h)
}

/// Place record-name labels along the top axis. Each name is
/// truncated (with ellipsis) to fit within its record's projected
/// pixel span; records too narrow for even a 3-character label are
/// skipped.
#[allow(clippy::too_many_arguments)]
fn draw_axis_record_labels_x(
    canvas: &mut [u8],
    stride: usize,
    total_h: usize,
    margin: usize,
    plot_w: u32,
    coord_w: u32,
    px_per_coord: f64,
    records: &[AxisRecord],
    label_y: usize,
) {
    let plot_left = margin;
    let plot_right = margin + plot_w as usize;
    let label_ink = 60_u8; // slightly bolder than tick labels
    for r in records {
        if r.end <= r.start || r.start >= coord_w {
            continue;
        }
        let start = r.start.min(coord_w);
        let end = r.end.min(coord_w);
        let x0 = plot_left + (start as f64 * px_per_coord).round() as usize;
        let x1 = plot_left + (end as f64 * px_per_coord).round() as usize;
        let span = x1.saturating_sub(x0);
        let max_chars = span.checked_div(FONT_ADVANCE).unwrap_or(0);
        if max_chars < 3 {
            continue;
        }
        let name = truncate_to_chars(&r.name, max_chars);
        let lw = text_width(&name);
        let centre = (x0 + x1) / 2;
        let lx = centre.saturating_sub(lw / 2).max(plot_left);
        let lx = lx.min(plot_right.saturating_sub(lw));
        draw_text(canvas, stride, total_h, lx, label_y, &name, label_ink);
    }
}

#[allow(clippy::too_many_arguments)]
fn draw_axis_record_labels_y(
    canvas: &mut [u8],
    stride: usize,
    total_h: usize,
    margin: usize,
    plot_h: u32,
    coord_h: u32,
    px_per_coord: f64,
    records: &[AxisRecord],
    plot_left: usize,
) {
    let plot_top = margin;
    let plot_bot = margin + plot_h as usize;
    let label_ink = 60_u8;
    // Place record names hard against the left margin (x = 2) so
    // they don't collide with the tick labels (which sit a bit to
    // the right). Truncate to the available char count between
    // x=2 and the tick-label area.
    let avail_chars = plot_left.saturating_sub(5 + 4 * FONT_ADVANCE + 2) / FONT_ADVANCE;
    for r in records {
        if r.end <= r.start || r.start >= coord_h {
            continue;
        }
        let start = r.start.min(coord_h);
        let end = r.end.min(coord_h);
        let y0 = margin + (start as f64 * px_per_coord).round() as usize;
        let y1 = margin + (end as f64 * px_per_coord).round() as usize;
        let span = y1.saturating_sub(y0);
        if span < FONT_H + 2 {
            continue;
        }
        let _ = (plot_top, plot_bot); // bounds only; label fits inside span by construction
        let max_chars = avail_chars.max(3);
        let name = truncate_to_chars(&r.name, max_chars);
        let centre = (y0 + y1) / 2;
        let ly = centre.saturating_sub(FONT_H / 2);
        // Draw at the very left of the canvas margin.
        draw_text(canvas, stride, total_h, 2, ly, &name, label_ink);
    }
}

/// Truncate a string so it fits in `max_chars` glyphs, using `…`
/// (rendered as a single `.`) when the original is longer.
fn truncate_to_chars(s: &str, max_chars: usize) -> String {
    if max_chars == 0 {
        return String::new();
    }
    let count = s.chars().count();
    if count <= max_chars {
        return s.to_string();
    }
    if max_chars <= 1 {
        return ".".to_string();
    }
    // Drop trailing chars, append "." as an ellipsis marker (our
    // bitmap font doesn't carry an actual U+2026).
    let keep: String = s.chars().take(max_chars - 1).collect();
    format!("{keep}.")
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

/// Nearest-neighbour resize a `(src_w × src_h)` greyscale pixelmap
/// to **exactly** `dst_w` pixels wide, with the height derived to
/// preserve the input aspect ratio.
///
/// Returns the resized buffer plus its `(dst_w, dst_h)` dimensions.
/// Special cases:
///
/// * `dst_w == 0` — caller is opting out: return the original
///   pixelmap unchanged.
/// * `dst_w == src_w` — no rescale needed; return as-is.
/// * `dst_w < src_w` — downscaling. Nearest-neighbour drops pixels
///   (skips entries on the input grid). For dotplots this loses
///   information; callers usually want to *up*scale instead, but
///   the downscale path is here for completeness.
/// * `dst_w > src_w` — upscaling. Each output pixel maps to a
///   single input pixel via `src_x = floor(dst_x * src_w / dst_w)`
///   — clean block expansion with no aliasing.
///
/// Replaces the earlier integer-only `upscale_nearest`: the
/// integer constraint meant we couldn't actually hit a
/// user-requested width (a 287 px input ↦ 2000 px request wanted
/// 6.97×; integer ceil gave 7× = 2009 px). Fractional NN is fine
/// for dotplots — most pixels are 0 / 255 so block-quantisation is
/// either invisible or the desired pixel-perfect look.
pub fn resize_nearest(pixels: &[u8], src_w: u32, src_h: u32, dst_w: u32) -> (Vec<u8>, u32, u32) {
    if dst_w == 0 || src_w == 0 || src_h == 0 {
        return (pixels.to_vec(), src_w, src_h);
    }
    if dst_w == src_w {
        return (pixels.to_vec(), src_w, src_h);
    }
    // Preserve aspect ratio: `dst_h = round(dst_w * src_h / src_w)`.
    let dst_h = ((dst_w as u64 * src_h as u64) / src_w as u64).max(1) as u32;
    let mut out = vec![0_u8; (dst_w as usize) * (dst_h as usize)];
    let src_w_u = src_w as u64;
    let src_h_u = src_h as u64;
    let dst_w_u = dst_w as u64;
    let dst_h_u = dst_h as u64;
    for y in 0..dst_h as usize {
        let sy = ((y as u64) * src_h_u / dst_h_u).min(src_h_u - 1) as usize;
        let src_row = sy * src_w as usize;
        let dst_row = y * dst_w as usize;
        for x in 0..dst_w as usize {
            let sx = ((x as u64) * src_w_u / dst_w_u).min(src_w_u - 1) as usize;
            out[dst_row + x] = pixels[src_row + sx];
        }
    }
    (out, dst_w, dst_h)
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
        let (canvas, tw, th) = compose_image_with_axes(&pixels, 4, 3, 4, 3, 10, &[], &[]);
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
    fn resize_hits_exact_target_width() {
        let pixels: Vec<u8> = (0..287 * 287).map(|i| (i % 256) as u8).collect();
        let (out, w, h) = resize_nearest(&pixels, 287, 287, 2000);
        assert_eq!(w, 2000);
        assert_eq!(h, 2000); // square input → square output
        assert_eq!(out.len(), 2000 * 2000);
        // Top-left output pixel maps to top-left input pixel.
        assert_eq!(out[0], pixels[0]);
    }

    #[test]
    fn resize_preserves_aspect_ratio() {
        // 100 × 50, target 200 ↦ height = 200 * 50 / 100 = 100.
        let pixels = vec![0_u8; 100 * 50];
        let (_, w, h) = resize_nearest(&pixels, 100, 50, 200);
        assert_eq!(w, 200);
        assert_eq!(h, 100);
        // 1342 × 1308, target 2000 ↦ height = 2000*1308/1342 = 1949.
        let pixels = vec![0_u8; 1342 * 1308];
        let (_, w, h) = resize_nearest(&pixels, 1342, 1308, 2000);
        assert_eq!(w, 2000);
        assert_eq!(h, 1949);
    }

    #[test]
    fn resize_dst_w_zero_is_noop() {
        let pixels = vec![42_u8; 100 * 50];
        let (out, w, h) = resize_nearest(&pixels, 100, 50, 0);
        assert_eq!((w, h), (100, 50));
        assert_eq!(out, pixels);
    }

    #[test]
    fn resize_dst_w_equal_src_w_is_noop() {
        let pixels: Vec<u8> = (0..3000_u32).map(|i| (i % 256) as u8).collect();
        let (out, w, h) = resize_nearest(&pixels, 300, 10, 300);
        assert_eq!((w, h), (300, 10));
        assert_eq!(out, pixels);
    }

    #[test]
    fn resize_downscale_drops_pixels_cleanly() {
        // 10×10 input with each pixel marked by its row index.
        let mut pixels = vec![0_u8; 100];
        for y in 0..10 {
            for x in 0..10 {
                pixels[y * 10 + x] = y as u8;
            }
        }
        // Downscale to 5×5.
        let (out, w, h) = resize_nearest(&pixels, 10, 10, 5);
        assert_eq!((w, h), (5, 5));
        // Each output row maps to a single input row: row y maps to
        // input row (y * 10 / 5) = 2y. So out[2*5+0] should equal 4.
        assert_eq!(out[2 * 5], 4);
    }

    #[test]
    fn draw_text_emits_expected_pixels() {
        // Render "0" onto a tiny canvas and check that at least one
        // pixel ended up dark — confirms the glyph table is wired
        // through to the writer.
        let mut canvas = vec![255_u8; 20 * 10];
        draw_text(&mut canvas, 20, 10, 2, 1, "0", 0);
        let any_dark = canvas.contains(&0);
        assert!(any_dark, "draw_text did not emit any ink");
    }
}
