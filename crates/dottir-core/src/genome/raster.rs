//! Deterministic rasterisation of alignments into a genome plot
//! (spec §4.7.7). Each pixel keeps the highest identity drawn into it, so
//! the result does not depend on alignment order.

use super::alignment::Alignment;
use super::layout::AxisLayout;

/// D-GENIES identity class limits (`Paf.limit_idy`): classes are
/// `< 0.25`, `< 0.5`, `< 0.75` and the rest.
pub const IDENTITY_LIMITS: [f64; 3] = [0.25, 0.5, 0.75];

/// Pixel code for an identity: 0 = nothing drawn, else `1 + 254·identity`.
pub fn identity_code(identity: f64) -> u8 {
    1 + (identity.clamp(0.0, 1.0) * 254.0).round() as u8
}

/// Identity class (0..=3) of a non-zero pixel code.
pub fn identity_class(code: u8) -> Option<usize> {
    if code == 0 {
        return None;
    }
    let id = f64::from(code - 1) / 254.0;
    Some(IDENTITY_LIMITS.iter().take_while(|&&l| id >= l).count())
}

#[derive(Debug, Clone, PartialEq)]
pub struct RasterConfig {
    /// Bases per pixel on both axes.
    pub bp_per_px: u64,
    /// Line thickness in pixels.
    pub line_width: u32,
    /// Self-comparison: also draw every alignment mirrored across the
    /// diagonal, and the diagonal of each record (spec §4.7.3).
    pub self_comparison: bool,
}

/// A rendered genome plot: x = target (reference), y = query, origin top
/// left.
#[derive(Debug, Clone, PartialEq)]
pub struct GenomeRaster {
    pub width: u32,
    pub height: u32,
    pub bp_per_px: u64,
    /// Row-major [`identity_code`]s.
    pub pixels: Vec<u8>,
}

impl GenomeRaster {
    fn stamp(&mut self, x: i64, y: i64, w: i64, code: u8) {
        let off = (w - 1) / 2;
        for yy in (y - off).max(0)..(y - off + w).min(i64::from(self.height)) {
            for xx in (x - off).max(0)..(x - off + w).min(i64::from(self.width)) {
                let p = &mut self.pixels[yy as usize * self.width as usize + xx as usize];
                *p = (*p).max(code);
            }
        }
    }

    /// Draw a line between two axis points (bp).
    fn line(&mut self, from: (u64, u64), to: (u64, u64), w: i64, code: u8) {
        let px = |v: u64, max: u32| ((v / self.bp_per_px) as i64).min(i64::from(max) - 1);
        let (x0, y0) = (px(from.0, self.width), px(from.1, self.height));
        let (x1, y1) = (px(to.0, self.width), px(to.1, self.height));
        let (dx, dy) = (x1 - x0, y1 - y0);
        let n = dx.abs().max(dy.abs());
        if n == 0 {
            self.stamp(x0, y0, w, code);
            return;
        }
        // Round to nearest with integer arithmetic only.
        let step = |d: i64, i: i64| d.signum() * ((2 * d.abs() * i + n) / (2 * n));
        for i in 0..=n {
            self.stamp(x0 + step(dx, i), y0 + step(dy, i), w, code);
        }
    }
}

/// Rasterise `alns` onto the target (x) and query (y) layouts.
pub fn rasterize<'a>(
    alns: impl IntoIterator<Item = &'a Alignment>,
    x: &AxisLayout,
    y: &AxisLayout,
    cfg: &RasterConfig,
) -> GenomeRaster {
    let bp = cfg.bp_per_px.max(1);
    let mut r = GenomeRaster {
        width: x.total.div_ceil(bp).max(1) as u32,
        height: y.total.div_ceil(bp).max(1) as u32,
        bp_per_px: bp,
        pixels: Vec::new(),
    };
    r.pixels = vec![0; r.width as usize * r.height as usize];
    let w = i64::from(cfg.line_width.max(1));
    for a in alns {
        let code = identity_code(a.identity());
        let start = (
            x.to_axis(a.t_rec, a.t_at_q_start()),
            y.to_axis(a.q_rec, a.q_start),
        );
        let end = (
            x.to_axis(a.t_rec, a.t_at_q_end()),
            y.to_axis(a.q_rec, a.q_end),
        );
        r.line(start, end, w, code);
        if cfg.self_comparison {
            let start = (
                x.to_axis(a.q_rec, a.q_start),
                y.to_axis(a.t_rec, a.t_at_q_start()),
            );
            let end = (
                x.to_axis(a.q_rec, a.q_end),
                y.to_axis(a.t_rec, a.t_at_q_end()),
            );
            r.line(start, end, w, code);
        }
    }
    if cfg.self_comparison {
        let code = identity_code(1.0);
        for s in &x.segments {
            let (a, b) = (x.to_axis(s.rec, 0), x.to_axis(s.rec, s.len));
            let (c, d) = (y.to_axis(s.rec, 0), y.to_axis(s.rec, s.len));
            r.line((a, c), (b, d), w, code);
        }
    }
    r
}

#[cfg(test)]
mod tests {
    use super::super::alignment::GenomeRecord;
    use super::*;

    fn layout(lens: &[u64]) -> AxisLayout {
        let recs: Vec<GenomeRecord> = lens
            .iter()
            .enumerate()
            .map(|(i, &len)| GenomeRecord {
                name: format!("r{i}"),
                len,
            })
            .collect();
        AxisLayout::new(&recs, None, &[], 0.0)
    }

    fn aln(q: (u64, u64), t: (u64, u64), reverse: bool, id: f64) -> Alignment {
        Alignment {
            q_rec: 0,
            q_start: q.0,
            q_end: q.1,
            t_rec: 0,
            t_start: t.0,
            t_end: t.1,
            reverse,
            matches: ((q.1 - q.0) as f64 * id) as u64,
            block_len: q.1 - q.0,
            mapq: 60,
            primary: true,
            divergence: None,
        }
    }

    fn cfg() -> RasterConfig {
        RasterConfig {
            bp_per_px: 10,
            line_width: 1,
            self_comparison: false,
        }
    }

    #[test]
    fn identity_classes_follow_dgenies_limits() {
        assert_eq!(identity_class(0), None);
        assert_eq!(identity_class(identity_code(0.2)), Some(0));
        assert_eq!(identity_class(identity_code(0.25)), Some(1));
        assert_eq!(identity_class(identity_code(0.74)), Some(2));
        assert_eq!(identity_class(identity_code(0.99)), Some(3));
    }

    #[test]
    fn forward_and_reverse_lines() {
        let (x, y) = (layout(&[100]), layout(&[100]));
        let fwd = aln((0, 100), (0, 100), false, 0.9);
        let r = rasterize([&fwd], &x, &y, &cfg());
        assert_eq!((r.width, r.height), (10, 10));
        for i in 0..10 {
            assert_ne!(r.pixels[i * 10 + i], 0, "diagonal pixel {i}");
        }
        let rev = aln((0, 100), (0, 100), true, 0.9);
        let r = rasterize([&rev], &x, &y, &cfg());
        // Query start pairs with target end: top-right to bottom-left.
        assert_ne!(r.pixels[9], 0);
        assert_ne!(r.pixels[9 * 10], 0);
        assert_eq!(r.pixels[0], 0);
    }

    #[test]
    fn raster_does_not_depend_on_order_and_keeps_max_identity() {
        let (x, y) = (layout(&[1000]), layout(&[1000]));
        let alns = [
            aln((0, 500), (0, 500), false, 0.3),
            aln((0, 500), (0, 500), false, 0.9),
            aln((200, 900), (900, 1000), true, 0.6),
        ];
        let a = rasterize(alns.iter(), &x, &y, &cfg());
        let b = rasterize(alns.iter().rev(), &x, &y, &cfg());
        assert_eq!(a, b);
        assert_eq!(a.pixels[0], identity_code(0.9));
    }

    /// Pinned output for a fixed pseudo-random alignment set; tied to
    /// `GENOME_FORMAT_VERSION`.
    #[test]
    fn golden_raster_hash() {
        let recs: Vec<GenomeRecord> = [3_000_000u64, 1_500_000, 700_000, 40_000]
            .iter()
            .enumerate()
            .map(|(i, &len)| GenomeRecord {
                name: format!("r{i}"),
                len,
            })
            .collect();
        let x = AxisLayout::new(&recs, None, &[], 0.002);
        let y = AxisLayout::new(
            &recs,
            Some(&[2, 0, 3, 1]),
            &[false, true, false, false],
            0.002,
        );
        let mut state = 7u64;
        let mut next = |n: u64| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) % n
        };
        let alns: Vec<Alignment> = (0..500)
            .map(|_| {
                let (q, t) = (next(4) as usize, next(4) as usize);
                let len = 1 + next(recs[q].len.min(recs[t].len) / 4);
                let qs = next(recs[q].len - len);
                let ts = next(recs[t].len - len);
                Alignment {
                    q_rec: q as u32,
                    q_start: qs,
                    q_end: qs + len,
                    t_rec: t as u32,
                    t_start: ts,
                    t_end: ts + len,
                    reverse: next(2) == 1,
                    matches: next(len + 1),
                    block_len: len,
                    mapq: 0,
                    primary: true,
                    divergence: None,
                }
            })
            .collect();
        let c = RasterConfig {
            bp_per_px: 2000,
            line_width: 2,
            self_comparison: false,
        };
        let r = rasterize(alns.iter(), &x, &y, &c);
        let hash = r.pixels.iter().fold(0xcbf2_9ce4_8422_2325u64, |h, &b| {
            (h ^ u64::from(b)).wrapping_mul(0x100_0000_01b3)
        });
        assert_eq!(super::super::GENOME_FORMAT_VERSION, 1);
        assert_eq!((r.width, r.height), (2620, 2620));
        assert_eq!(
            hash, 1885063123674135242,
            "update the pin only with a GENOME_FORMAT_VERSION bump"
        );
    }

    #[test]
    fn self_comparison_mirrors_and_draws_the_diagonal() {
        let l = layout(&[100]);
        let c = RasterConfig {
            self_comparison: true,
            ..cfg()
        };
        let off = aln((0, 20), (60, 80), false, 0.9);
        let r = rasterize([&off], &l, &l, &c);
        assert_ne!(r.pixels[6], 0); // (x = 6, y = 0)
        assert_ne!(r.pixels[6 * 10], 0); // mirrored (x = 0, y = 6)
        assert_ne!(r.pixels[5 * 10 + 5], 0); // diagonal
    }
}
