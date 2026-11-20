//! Cost estimate of an exact-engine plot (docs/GUI_NAVIGATION_PLAN.md N2).

use dottir_core::{BlastMode, Strand};

/// Cells per second assumed before an exact run has been measured.
pub const DEFAULT_RATE: f64 = 1e9;
/// An exact plot asked for explicitly starts without asking below this.
pub const ASK_ABOVE_EXPLICIT_S: f64 = 5.0;
/// Any other exact plot (global Exact engine, or zooming an exact view)
/// asks first above this.
pub const ASK_ABOVE_S: f64 = 30.0;
/// Only runs at least this long update the measured rate.
pub const MEASURE_ABOVE_S: f64 = 1.0;

/// What an exact computation is expected to cost.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExactCost {
    pub cells: f64,
    pub seconds: f64,
    pub bytes: u64,
}

/// The work and memory of one exact run.
#[derive(Clone, Copy, Debug)]
pub struct ExactJob {
    pub mode: BlastMode,
    pub strand: Strand,
    /// Residues on each axis (bases for a BLASTX query).
    pub qlen: usize,
    pub slen: usize,
    pub self_comparison: bool,
    pub zoom: u32,
    pub separate_strands: bool,
    pub threads: usize,
}

impl ExactJob {
    /// Query cells visited: one pass per strand (BLASTX: three frames of
    /// `qlen / 3`); a self-comparison pass stops at the diagonal.
    pub fn cells(&self) -> f64 {
        let (q, s) = (self.qlen as f64, self.slen as f64);
        let per_pass = if self.self_comparison {
            let n = q.min(s);
            n * (n + 1.0) / 2.0 + (s - n).max(0.0) * q
        } else {
            q * s
        };
        match self.mode {
            BlastMode::Blastx => 3.0 * (q / 3.0).floor() * s,
            BlastMode::Blastp => per_pass,
            BlastMode::Blastn => match self.strand {
                Strand::Both => 2.0 * per_pass,
                _ => per_pass,
            },
        }
    }

    /// Pixelmap(s), the score vector and the per-thread sum buffers.
    pub fn bytes(&self) -> u64 {
        let z = self.zoom.max(1) as u64;
        let qcols = match self.mode {
            BlastMode::Blastx => self.qlen as u64 / 3,
            _ => self.qlen as u64,
        };
        let (w, h) = (qcols.div_ceil(z), (self.slen as u64).div_ceil(z));
        let maps = 1 + (self.separate_strands && self.strand == Strand::Both) as u64;
        let alphabet = match self.mode {
            BlastMode::Blastn => 5,
            _ => 25,
        };
        w * h * maps + 4 * alphabet * qcols + 2 * 4 * qcols * self.threads.max(1) as u64
    }

    pub fn cost(&self, rate: f64) -> ExactCost {
        let cells = self.cells();
        ExactCost {
            cells,
            seconds: cells / rate.max(1.0),
            bytes: self.bytes(),
        }
    }
}

/// Running estimate of the exact engine's speed on this machine.
#[derive(Clone, Copy, Debug, Default)]
pub struct Rate {
    measured: Option<f64>,
}

impl Rate {
    pub fn get(&self) -> f64 {
        self.measured.unwrap_or(DEFAULT_RATE)
    }

    pub fn is_measured(&self) -> bool {
        self.measured.is_some()
    }

    /// Fold in a finished run; short runs are dominated by overhead and
    /// are ignored.
    pub fn record(&mut self, cells: f64, seconds: f64) {
        if seconds < MEASURE_ABOVE_S || cells <= 0.0 {
            return;
        }
        let r = cells / seconds;
        self.measured = Some(match self.measured {
            None => r,
            Some(old) => 0.5 * old + 0.5 * r,
        });
    }
}

/// `1.2 s`, `4 min`, `2.5 h`, `3.1 days`.
pub fn fmt_duration(s: f64) -> String {
    if s < 60.0 {
        format!("{s:.1} s")
    } else if s < 3600.0 {
        format!("{:.0} min", s / 60.0)
    } else if s < 86_400.0 {
        format!("{:.1} h", s / 3600.0)
    } else {
        format!("{:.1} days", s / 86_400.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn job() -> ExactJob {
        ExactJob {
            mode: BlastMode::Blastn,
            strand: Strand::Forward,
            qlen: 1000,
            slen: 2000,
            self_comparison: false,
            zoom: 10,
            separate_strands: false,
            threads: 4,
        }
    }

    #[test]
    fn cells_count_strands_and_the_self_triangle() {
        assert_eq!(job().cells(), 2e6);
        let both = ExactJob {
            strand: Strand::Both,
            ..job()
        };
        assert_eq!(both.cells(), 4e6);
        let self_cmp = ExactJob {
            slen: 1000,
            self_comparison: true,
            ..job()
        };
        assert_eq!(self_cmp.cells(), 1000.0 * 1001.0 / 2.0);
        let blastx = ExactJob {
            mode: BlastMode::Blastx,
            qlen: 3000,
            ..job()
        };
        assert_eq!(blastx.cells(), 3.0 * 1000.0 * 2000.0);
    }

    #[test]
    fn bytes_count_pixelmaps_and_buffers() {
        // 100 × 200 pixels, score vector 5 × 4 B × 1000, buffers 2 × 4 B × 1000 × 4.
        assert_eq!(job().bytes(), 20_000 + 20_000 + 32_000);
        let sep = ExactJob {
            strand: Strand::Both,
            separate_strands: true,
            ..job()
        };
        assert_eq!(sep.bytes(), 2 * 20_000 + 20_000 + 32_000);
    }

    #[test]
    fn rate_ignores_short_runs_and_averages() {
        let mut r = Rate::default();
        assert_eq!(r.get(), DEFAULT_RATE);
        r.record(1e9, 0.5);
        assert!(!r.is_measured());
        r.record(4e9, 2.0);
        assert_eq!(r.get(), 2e9);
        r.record(8e9, 2.0);
        assert_eq!(r.get(), 3e9);
        assert_eq!(job().cost(2e6).seconds, 1.0);
    }

    #[test]
    fn durations_read_naturally() {
        assert_eq!(fmt_duration(1.24), "1.2 s");
        assert_eq!(fmt_duration(240.0), "4 min");
        assert_eq!(fmt_duration(9000.0), "2.5 h");
        assert_eq!(fmt_duration(3.0 * 86_400.0), "3.0 days");
    }
}
