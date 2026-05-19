//! Power-spectrum analysis of 1-D real signals, intended to find
//! periodicities in [`crate::periodogram::Periodogram`] output.
//!
//! Why FFT a periodogram? A tandem repeat of period P makes the
//! periodogram itself periodic: P(k) has peaks at k = P, 2P, 3P, …
//! The FFT of P(k) (as a function of k) therefore has a peak at
//! frequency 1/P, which directly identifies the underlying period.
//! Going through the periodogram first means the FFT inherits the
//! periodogram's denoising (windowed scoring + sensitivity ramp)
//! instead of seeing every base-level base-pair frequency.
//!
//! Pipeline:
//! 1. Detrend (default: subtract the input mean — otherwise the
//!    `f = 0` bin dominates and hides everything else).
//! 2. Apply a window function (default: Hann) to suppress
//!    spectral leakage from the finite signal edges.
//! 3. Zero-pad to the next power of two (default: on) for FFT
//!    efficiency and finer frequency resolution.
//! 4. Real FFT via [`realfft`].
//! 5. Find local-maximum bins, rank by amplitude, mark the top N.
//!
//! Each spectrum bin `i` corresponds to a period of
//! `padded_length / i` units of the input signal's index. For a
//! periodogram indexed in residues, that's a period in residues.
//!
//! # Example
//!
//! ```
//! use dottir_core::spectrum::{spectrum_of_signal, SpectrumConfig};
//!
//! // Synthetic periodogram: peaks every 7 indices.
//! let mut signal = vec![0.0; 256];
//! for k in (7..signal.len()).step_by(7) {
//!     signal[k] = 1.0;
//! }
//! let s = spectrum_of_signal(&signal, &SpectrumConfig::default()).unwrap();
//! // For an ideal impulse train, the fundamental (period 7) and its
//! // harmonics (3.5, 2.33, …) have equal amplitude — the top peak
//! // may land on any of them. The fundamental will be in the ranked
//! // peaks though.
//! let periods: Vec<f64> = s.ranked_peaks().iter().map(|&b| s.period(b)).collect();
//! assert!(periods.iter().any(|p| (p - 7.0).abs() < 0.5), "periods: {periods:?}");
//! ```

use realfft::num_complex::Complex;
use realfft::RealFftPlanner;

use crate::error::DottirError;

// ---------------------------------------------------------------------------
// Knobs
// ---------------------------------------------------------------------------

/// Window function applied to the detrended signal before FFT.
///
/// Windowing tapers the signal toward zero at both endpoints,
/// suppressing the high-frequency artefacts caused by treating a
/// finite-length sample as one period of an infinite signal
/// (spectral leakage).
#[derive(Clone, Copy, Debug, Default)]
pub enum WindowFn {
    /// No window — rectangular. Sharpest frequency resolution but
    /// worst leakage; mostly useful as a control.
    None,
    /// Hann (raised-cosine). Good general-purpose default: ~32 dB
    /// side-lobe rejection, moderate main-lobe width.
    #[default]
    Hann,
}

impl WindowFn {
    fn coefficient(self, n: usize, len: usize) -> f64 {
        if len <= 1 {
            return 1.0;
        }
        match self {
            WindowFn::None => 1.0,
            WindowFn::Hann => {
                let denom = (len - 1) as f64;
                0.5 * (1.0 - (2.0 * std::f64::consts::PI * n as f64 / denom).cos())
            }
        }
    }
}

/// What to subtract from the signal before windowing.
#[derive(Clone, Copy, Debug, Default)]
pub enum DetrendMode {
    /// Leave the signal alone. The `f = 0` (DC) bin will contain
    /// the squared mean and dominate the spectrum.
    None,
    /// Subtract the arithmetic mean. Cheap and almost always what
    /// you want — moves the signal's DC component out of the way.
    #[default]
    SubtractMean,
}

/// Configuration for [`spectrum_of_signal`].
#[derive(Clone, Copy, Debug)]
pub struct SpectrumConfig {
    pub window: WindowFn,
    pub detrend: DetrendMode,
    /// Zero-pad the windowed signal up to the next power of two
    /// before FFT. Defaults to `true` — gives finer frequency
    /// resolution and the FFT runs at peak speed.
    pub pad_to_pow2: bool,
    /// How many local-maximum bins to mark with a rank (1 = brightest
    /// non-DC peak). `0` disables peak detection.
    pub top_peaks: usize,
}

impl Default for SpectrumConfig {
    fn default() -> Self {
        Self {
            window: WindowFn::default(),
            detrend: DetrendMode::default(),
            pad_to_pow2: true,
            top_peaks: 10,
        }
    }
}

// ---------------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------------

/// Result of [`spectrum_of_signal`].
#[derive(Clone, Debug)]
pub struct Spectrum {
    /// Length of the user-supplied signal, before padding.
    pub input_len: usize,
    /// Length of the buffer the FFT actually ran on (≥ `input_len`,
    /// possibly padded up to a power of two).
    pub padded_len: usize,
    /// `|X[i]|` for `i in 0..=padded_len/2`. Bin 0 is DC; bin
    /// `padded_len/2` is the Nyquist frequency.
    pub amplitude: Vec<f64>,
    /// `Some(rank)` (1-indexed) for the top-N local maxima by
    /// amplitude; `None` everywhere else. Same length as
    /// [`Self::amplitude`].
    pub peak_ranks: Vec<Option<u32>>,
}

impl Spectrum {
    /// Frequency at bin `i`, in cycles per input-index unit.
    /// `frequency(0) = 0` (DC); `frequency(padded_len / 2) = 0.5`
    /// (Nyquist).
    #[inline]
    pub fn frequency(&self, bin: usize) -> f64 {
        bin as f64 / self.padded_len as f64
    }

    /// Period at bin `i`, in input-index units (= residues if the
    /// input was a per-residue periodogram). `period(0) = +inf`.
    #[inline]
    pub fn period(&self, bin: usize) -> f64 {
        if bin == 0 {
            f64::INFINITY
        } else {
            self.padded_len as f64 / bin as f64
        }
    }

    /// The bin index of the brightest peak (rank 1). `None` if peak
    /// detection was disabled or the signal had no local maxima.
    pub fn top_peak(&self) -> Option<usize> {
        self.peak_ranks
            .iter()
            .enumerate()
            .find_map(|(i, r)| if *r == Some(1) { Some(i) } else { None })
    }

    /// All peak bins in rank order (rank 1 first).
    pub fn ranked_peaks(&self) -> Vec<usize> {
        let mut peaks: Vec<(u32, usize)> = self
            .peak_ranks
            .iter()
            .enumerate()
            .filter_map(|(i, r)| r.map(|rk| (rk, i)))
            .collect();
        peaks.sort_unstable_by_key(|&(rk, _)| rk);
        peaks.into_iter().map(|(_, i)| i).collect()
    }
}

// ---------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------

/// Compute the magnitude spectrum of a real-valued signal.
///
/// See module docs for the pipeline. Returns
/// [`DottirError::EmptySequence`] for an empty input.
pub fn spectrum_of_signal(input: &[f64], config: &SpectrumConfig) -> Result<Spectrum, DottirError> {
    if input.is_empty() {
        return Err(DottirError::EmptySequence);
    }
    let input_len = input.len();

    // 1. Detrend.
    let mut buf: Vec<f64> = match config.detrend {
        DetrendMode::None => input.to_vec(),
        DetrendMode::SubtractMean => {
            let mean = input.iter().sum::<f64>() / input_len as f64;
            input.iter().map(|x| x - mean).collect()
        }
    };

    // 2. Window.
    for (n, x) in buf.iter_mut().enumerate() {
        *x *= config.window.coefficient(n, input_len);
    }

    // 3. Pad.
    let padded_len = if config.pad_to_pow2 {
        next_pow2(input_len.max(2))
    } else {
        input_len
    };
    if padded_len > input_len {
        buf.resize(padded_len, 0.0);
    }

    // 4. Real FFT.
    let mut planner = RealFftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(padded_len);
    let mut spectrum_buf: Vec<Complex<f64>> = fft.make_output_vec();
    fft.process(&mut buf, &mut spectrum_buf)
        .map_err(|e| DottirError::InvalidConfig(format!("realfft failed: {e}")))?;

    let amplitude: Vec<f64> = spectrum_buf.iter().map(|c| c.norm()).collect();

    // 5. Peaks.
    let peak_ranks = if config.top_peaks > 0 {
        rank_local_maxima(&amplitude, config.top_peaks)
    } else {
        vec![None; amplitude.len()]
    };

    Ok(Spectrum {
        input_len,
        padded_len,
        amplitude,
        peak_ranks,
    })
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn next_pow2(n: usize) -> usize {
    if n <= 1 {
        return 1;
    }
    // For usize values where the high bit isn't set, this is exact.
    let mut x = 1usize;
    while x < n {
        x <<= 1;
    }
    x
}

/// Find local-maximum bins in `amp` (strict three-bin maxima,
/// excluding bin 0), sort by amplitude descending, take the first
/// `top_n`, and return a `Vec<Option<u32>>` with rank labels at
/// peak positions.
fn rank_local_maxima(amp: &[f64], top_n: usize) -> Vec<Option<u32>> {
    let mut ranks = vec![None; amp.len()];
    if amp.len() < 3 {
        return ranks;
    }
    // Skip bin 0 (DC) — even after detrending, numeric noise leaves
    // a small DC residue that occasionally looks like a peak.
    let mut peaks: Vec<(usize, f64)> = (1..amp.len() - 1)
        .filter(|&i| amp[i] > amp[i - 1] && amp[i] > amp[i + 1])
        .map(|i| (i, amp[i]))
        .collect();
    // Stable-sort by amplitude descending; ties broken by lower bin
    // index (lower frequency / longer period) so the output is
    // deterministic across runs.
    peaks.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.0.cmp(&b.0))
    });
    for (rank, &(bin, _)) in peaks.iter().take(top_n).enumerate() {
        ranks[bin] = Some(rank as u32 + 1);
    }
    ranks
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn impulse_train(len: usize, period: usize) -> Vec<f64> {
        let mut s = vec![0.0; len];
        for k in (period..len).step_by(period) {
            s[k] = 1.0;
        }
        s
    }

    #[test]
    fn next_pow2_basics() {
        assert_eq!(next_pow2(0), 1);
        assert_eq!(next_pow2(1), 1);
        assert_eq!(next_pow2(2), 2);
        assert_eq!(next_pow2(3), 4);
        assert_eq!(next_pow2(255), 256);
        assert_eq!(next_pow2(256), 256);
        assert_eq!(next_pow2(257), 512);
    }

    #[test]
    fn empty_input_errors() {
        let err = spectrum_of_signal(&[], &SpectrumConfig::default()).unwrap_err();
        assert!(matches!(err, DottirError::EmptySequence));
    }

    #[test]
    fn dc_only_input_has_no_peak_after_detrend() {
        let signal = vec![42.0; 64];
        let s = spectrum_of_signal(&signal, &SpectrumConfig::default()).unwrap();
        // After mean subtraction the signal is identically zero, so
        // all bins are zero and there are no local maxima.
        assert!(s.amplitude.iter().all(|&v| v.abs() < 1e-9));
        assert!(s.peak_ranks.iter().all(|r| r.is_none()));
    }

    #[test]
    fn pure_sinusoid_peaks_at_its_frequency() {
        // 128-sample sinusoid completing 8 cycles → bin 8 should
        // dominate (period = 128/8 = 16 samples per cycle).
        let n = 128;
        let cycles = 8.0;
        let signal: Vec<f64> = (0..n)
            .map(|k| (2.0 * std::f64::consts::PI * cycles * k as f64 / n as f64).sin())
            .collect();
        let cfg = SpectrumConfig {
            pad_to_pow2: false, // exact bin alignment
            ..SpectrumConfig::default()
        };
        let s = spectrum_of_signal(&signal, &cfg).unwrap();
        let top = s.top_peak().expect("expected at least one peak");
        // Window broadening + leakage can shift the apparent peak
        // by 1 bin; allow a small tolerance.
        assert!(
            (top as i64 - cycles as i64).abs() <= 1,
            "expected peak near bin {cycles}, got {top}"
        );
    }

    #[test]
    fn impulse_train_recovers_period() {
        // Period-7 impulse train. An ideal impulse train has equal
        // amplitude at the fundamental f=1/7 and all its harmonics
        // 2/7, 3/7, …, so the top peak may land on a harmonic (e.g.
        // period 3.5). We don't yet dedup harmonics, so the realistic
        // check is that the fundamental period 7 appears among the
        // top peaks.
        let signal = impulse_train(256, 7);
        let s = spectrum_of_signal(&signal, &SpectrumConfig::default()).unwrap();
        let peaks = s.ranked_peaks();
        let near_seven = peaks.iter().any(|&bin| (s.period(bin) - 7.0).abs() < 0.5);
        assert!(
            near_seven,
            "expected period ~7 among ranked peaks, got periods {:?}",
            peaks.iter().map(|&b| s.period(b)).collect::<Vec<_>>()
        );
        // And every top peak should be a (sub-)harmonic of 7 —
        // its period divides 7 (within a small tolerance) or vice
        // versa. Equivalent statement: frequency is a rational
        // multiple of 1/7.
        for &bin in peaks.iter().take(3) {
            let freq = s.frequency(bin);
            let cycles_per_7 = freq * 7.0;
            assert!(
                (cycles_per_7 - cycles_per_7.round()).abs() < 0.05,
                "peak at frequency {freq} (period {}) is not a harmonic of 1/7",
                s.period(bin)
            );
        }
    }

    #[test]
    fn ranked_peaks_are_in_amplitude_order() {
        let signal = impulse_train(512, 11);
        let s = spectrum_of_signal(&signal, &SpectrumConfig::default()).unwrap();
        let peaks = s.ranked_peaks();
        assert!(!peaks.is_empty());
        // Rank 1 has the largest amplitude of any ranked peak.
        let amps: Vec<f64> = peaks.iter().map(|&i| s.amplitude[i]).collect();
        for w in amps.windows(2) {
            assert!(
                w[0] >= w[1] - 1e-12,
                "ranks should be descending in amplitude: {amps:?}"
            );
        }
    }

    #[test]
    fn period_and_frequency_are_inverses() {
        let s = Spectrum {
            input_len: 100,
            padded_len: 128,
            amplitude: vec![0.0; 65],
            peak_ranks: vec![None; 65],
        };
        for bin in [1, 5, 32, 64] {
            let f = s.frequency(bin);
            let p = s.period(bin);
            assert!(
                (f * p - 1.0).abs() < 1e-12,
                "f * p should be 1 at bin {bin}"
            );
        }
        assert!(s.period(0).is_infinite());
    }

    #[test]
    fn top_peaks_zero_disables_peak_detection() {
        let signal = impulse_train(256, 7);
        let cfg = SpectrumConfig {
            top_peaks: 0,
            ..SpectrumConfig::default()
        };
        let s = spectrum_of_signal(&signal, &cfg).unwrap();
        assert!(s.peak_ranks.iter().all(|r| r.is_none()));
        assert!(s.top_peak().is_none());
    }

    #[test]
    fn ranks_are_unique_and_dense() {
        // Construct a spectrum-like signal where we know how many
        // local maxima there are, then check the ranks 1..=k cover
        // exactly k bins.
        let signal: Vec<f64> = (0..200)
            .map(|i| (i as f64 * 0.13).sin() + (i as f64 * 0.41).sin())
            .collect();
        let cfg = SpectrumConfig {
            top_peaks: 5,
            ..SpectrumConfig::default()
        };
        let s = spectrum_of_signal(&signal, &cfg).unwrap();
        let mut ranks: Vec<u32> = s.peak_ranks.iter().filter_map(|r| *r).collect();
        ranks.sort_unstable();
        assert!(!ranks.is_empty(), "expected at least one peak");
        assert!(ranks.len() <= 5);
        // Dense 1..=k.
        for (i, &r) in ranks.iter().enumerate() {
            assert_eq!(r as usize, i + 1);
        }
    }
}
