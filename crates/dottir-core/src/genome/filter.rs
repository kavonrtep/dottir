//! Alignment filters (spec §4.7.5), including the D-GENIES noise filter.

use super::alignment::Alignment;

/// D-GENIES only filters noise in plots with more than this many
/// alignments.
pub const NOISE_MIN_ALIGNMENTS: usize = 1000;

#[derive(Debug, Clone, PartialEq)]
pub struct FilterConfig {
    /// Minimum alignment length (longer of the two spans).
    pub min_len: u64,
    /// Minimum identity in `[0, 1]`.
    pub min_identity: f64,
    pub primary_only: bool,
    /// Drop alignments shorter than [`noise_threshold`].
    pub remove_noise: bool,
}

impl Default for FilterConfig {
    fn default() -> Self {
        Self {
            min_len: 0,
            min_identity: 0.0,
            primary_only: false,
            remove_noise: false,
        }
    }
}

/// Length of an alignment's line in plot coordinates (bp), as D-GENIES
/// measures it for the noise filter and contig sorting.
pub fn line_len(a: &Alignment) -> f64 {
    (a.q_len() as f64).hypot(a.t_len() as f64)
}

/// D-GENIES `Paf.parse_paf` noise limit: histogram the line lengths into
/// `n / 10` equal-width bins, find the most populated bin, then the first
/// bin from there on holding fewer than 1/50 of its count; lines shorter
/// than that bin's lower edge are noise. `None` when there are at most
/// [`NOISE_MIN_ALIGNMENTS`] lines or no bin qualifies.
pub fn noise_threshold(lens: &[f64]) -> Option<f64> {
    if lens.len() <= NOISE_MIN_ALIGNMENTS {
        return None;
    }
    let n_bins = lens.len() / 10;
    let lo = lens.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = lens.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if hi <= lo {
        return None;
    }
    let width = (hi - lo) / n_bins as f64;
    let mut counts = vec![0u64; n_bins];
    for &l in lens {
        let b = (((l - lo) / width) as usize).min(n_bins - 1);
        counts[b] += 1;
    }
    let max = *counts.iter().max()?;
    let peak = counts.iter().position(|&c| c == max)?;
    let limit = (peak..n_bins).find(|&i| (counts[i] as f64) < max as f64 / 50.0)?;
    Some(lo + limit as f64 * width)
}

/// Which alignments pass `cfg`. The noise threshold is computed over the
/// alignments passing the other filters; it is returned when applied.
pub fn keep_mask(alns: &[Alignment], cfg: &FilterConfig) -> (Vec<bool>, Option<f64>) {
    let mut keep: Vec<bool> = alns
        .iter()
        .map(|a| {
            a.len() >= cfg.min_len
                && a.identity() >= cfg.min_identity
                && (!cfg.primary_only || a.primary)
        })
        .collect();
    if !cfg.remove_noise {
        return (keep, None);
    }
    let lens: Vec<f64> = alns
        .iter()
        .zip(&keep)
        .filter(|(_, &k)| k)
        .map(|(a, _)| line_len(a))
        .collect();
    let threshold = noise_threshold(&lens);
    if let Some(t) = threshold {
        for (a, k) in alns.iter().zip(&mut keep) {
            *k &= line_len(a) >= t;
        }
    }
    (keep, threshold)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn aln(len: u64, id: f64, primary: bool) -> Alignment {
        Alignment {
            q_rec: 0,
            q_start: 0,
            q_end: len,
            t_rec: 0,
            t_start: 0,
            t_end: len,
            reverse: false,
            matches: (len as f64 * id) as u64,
            block_len: len,
            mapq: 0,
            primary,
            divergence: None,
        }
    }

    #[test]
    fn noise_threshold_cuts_after_the_peak() {
        // 2000 short lines around 100 bp, a tail of 30 long ones.
        let mut lens: Vec<f64> = (0..2000).map(|i| 100.0 + (i % 50) as f64).collect();
        lens.extend((0..30).map(|i| 10_000.0 + 1000.0 * i as f64));
        let t = noise_threshold(&lens).unwrap();
        assert!(t > 150.0 && t < 10_000.0, "{t}");
        assert_eq!(noise_threshold(&lens[..1000]), None);
    }

    #[test]
    fn filters_combine() {
        let alns = vec![
            aln(500, 0.9, true),
            aln(5000, 0.5, true),
            aln(5000, 0.9, false),
        ];
        let cfg = FilterConfig {
            min_len: 1000,
            min_identity: 0.6,
            ..FilterConfig::default()
        };
        assert_eq!(keep_mask(&alns, &cfg).0, vec![false, false, true]);
        let cfg = FilterConfig {
            primary_only: true,
            ..FilterConfig::default()
        };
        assert_eq!(keep_mask(&alns, &cfg).0, vec![true, true, false]);
    }
}
