//! Large mode against the exact engine.

use super::*;
use crate::alphabet::reverse_complement_dna;
use crate::plot::{compute_dotplot, PlotConfig};
use proptest::prelude::*;

fn lcg_dna(n: usize, mut state: u64) -> Vec<u8> {
    (0..n)
        .map(|_| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(state >> 62) as usize]
        })
        .collect()
}

/// Substitute roughly `rate_pct` % of bases, deterministically.
fn mutate(seq: &[u8], rate_pct: u64, mut state: u64) -> Vec<u8> {
    seq.iter()
        .map(|&b| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            if (state >> 32) % 100 < rate_pct {
                b"ACGT"[((state >> 20) & 3) as usize]
            } else {
                b
            }
        })
        .collect()
}

/// Query/subject pair with a forward and an inverted diverged copy.
fn planted_pair(len: usize, unit_len: usize, rate: u64, seed: u64) -> (Vec<u8>, Vec<u8>) {
    let unit = lcg_dna(unit_len, seed ^ 0xabc);
    let mut q = lcg_dna(len, seed);
    let mut s = lcg_dna(len, seed ^ 0x55);
    let at = len / 4;
    q.splice(at..at + unit_len, unit.iter().copied());
    let m = mutate(&unit, rate, seed ^ 7);
    s.splice(at / 2..at / 2 + unit_len, m.iter().copied());
    let rc = reverse_complement_dna(&mutate(&unit, rate, seed ^ 9));
    let at2 = len / 2;
    s.splice(at2..at2 + unit_len, rc.iter().copied());
    (q, s)
}

fn exact(q: &[u8], s: &[u8], zoom: u32, window: Option<u32>, self_comp: bool) -> DotPlot {
    let mut cfg = PlotConfig::default_blastn(ScoreMatrix::dna_identity());
    cfg.zoom = zoom;
    cfg.window_size = window;
    cfg.separate_strand_channels = true;
    cfg.self_comparison = self_comp;
    compute_dotplot(q, s, &cfg).unwrap()
}

fn large_cfg(zoom: u32, window: Option<u32>) -> LargeConfig {
    LargeConfig {
        zoom,
        window_size: window,
        depth: DepthMethod::None,
        ..LargeConfig::default()
    }
}

fn large(q: &[u8], s: &[u8], cfg: &LargeConfig) -> LargePlot {
    compute_large_dotplot(&LargeInput::pair(q, &[], s, &[]), cfg).unwrap()
}

fn assert_le(large: &[u8], exact: &[u8], what: &str) {
    assert_eq!(large.len(), exact.len());
    for (i, (&l, &e)) in large.iter().zip(exact).enumerate() {
        assert!(l <= e, "{what}: pixel {i}: large {l} > exact {e}");
    }
}

#[test]
fn never_exceeds_exact_and_matches_on_planted_copies() {
    let (q, s) = planted_pair(6000, 800, 10, 3);
    for zoom in [1, 3, 8] {
        let e = exact(&q, &s, zoom, None, false);
        let l = large(&q, &s, &large_cfg(zoom, None));
        assert_eq!((l.width, l.height), (e.width, e.height));
        assert_eq!(l.params.window_size, e.params.window_size);
        assert_eq!(l.params.pixel_fac, e.params.pixel_fac);
        let lf = l.merged(Some(false), ClassMask::ALL);
        let lr = l.merged(Some(true), ClassMask::ALL);
        let er = e.reverse_pixels.as_ref().unwrap();
        assert_le(&lf, &e.pixels, "forward");
        assert_le(&lr, er, "reverse");
        // Inside the planted copy's box, the pixels large mode lit carry
        // the exact engine's value (cells off the anchored diagonal could
        // only raise the exact pixel where the copy itself is weak).
        let at = 6000 / 4;
        let w = l.width as usize;
        let z = zoom as usize;
        let (mut lit, mut equal) = (0, 0);
        for y in (at / 2 + 50) / z..(at / 2 + 750) / z {
            for x in (at + 50) / z..(at + 750) / z {
                let i = y * w + x;
                if lf[i] > 0 {
                    lit += 1;
                    equal += (lf[i] == e.pixels[i]) as usize;
                }
            }
        }
        assert!(lit * z >= 600, "zoom {zoom}: only {lit} lit");
        assert!(equal * 10 >= lit * 9, "zoom {zoom}: {equal}/{lit} equal");
        assert!(lr.iter().any(|&v| v > 0));
    }
}

#[test]
fn self_comparison_matches_exact_on_diagonal() {
    let unit = lcg_dna(500, 77);
    let mut seq = lcg_dna(5000, 78);
    seq.splice(3000..3500, mutate(&unit, 5, 1));
    seq.splice(1000..1500, unit.iter().copied());
    let e = exact(&seq, &seq, 4, None, true);
    let mut cfg = large_cfg(4, None);
    cfg.self_comparison = true;
    let l = compute_large_dotplot(&LargeInput::self_comparison(&seq, &[]), &cfg).unwrap();
    let lf = l.merged(Some(false), ClassMask::ALL);
    assert_le(&lf, &e.pixels, "self forward");
    assert_le(
        &l.merged(Some(true), ClassMask::ALL),
        e.reverse_pixels.as_ref().unwrap(),
        "self reverse",
    );
    let w = l.width as usize;
    for x in 20..(l.width as usize - 20) {
        assert_eq!(lf[x * w + x], e.pixels[x * w + x], "diagonal pixel {x}");
    }
    // Mirrored: symmetric like the exact engine.
    for y in 0..w {
        for x in 0..y {
            assert_eq!(lf[y * w + x], lf[x * w + y]);
        }
    }
}

#[test]
fn random_sequences_stay_nearly_blank() {
    let q = lcg_dna(100_000, 1);
    let s = lcg_dna(100_000, 2);
    let l = large(&q, &s, &large_cfg(50, None));
    let lit = l
        .merged(None, ClassMask::ALL)
        .iter()
        .filter(|&&v| v > 0)
        .count();
    assert!(lit < 20, "{lit} pixels lit on random input; {:?}", l.stats);
}

#[test]
fn viewport_equals_crop_of_full_plot() {
    let (q, s) = planted_pair(20_000, 3000, 15, 11);
    let zoom = 5;
    let full = large(&q, &s, &large_cfg(zoom, Some(20)));
    let mut cfg = large_cfg(zoom, Some(20));
    cfg.viewport = Some(Viewport {
        q: 4000..12_003,
        s: 1500..13_000,
    });
    let v = large(&q, &s, &cfg);
    assert_eq!((v.origin_q, v.origin_s), (4000, 1500));
    let (vw, fw) = (v.width as usize, full.width as usize);
    let (x0, y0) = (v.origin_q / 5, v.origin_s / 5);
    for strand in [false, true] {
        let fm = full.merged(Some(strand), ClassMask::ALL);
        let vm = v.merged(Some(strand), ClassMask::ALL);
        for y in 0..v.height as usize {
            for x in 0..vw {
                assert_eq!(
                    vm[y * vw + x],
                    fm[(y + y0) * fw + x + x0],
                    "strand {strand} pixel ({x},{y})"
                );
            }
        }
        assert!(vm.iter().any(|&p| p > 0));
    }
}

#[test]
fn deterministic_across_thread_counts() {
    let (q, s) = planted_pair(200_000, 5000, 15, 5);
    let mut cfg = large_cfg(40, None);
    cfg.chain = Some(ChainConfig::default());
    let run = |threads| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| large(&q, &s, &cfg))
    };
    let a = run(1);
    let b = run(4);
    assert_eq!(a.channels, b.channels);
    assert_eq!(a.segments, b.segments);
    assert_eq!(a.chains, b.chains);
    assert_eq!(a.q_capped, b.q_capped);
}

#[test]
fn record_breaks_and_n_runs_are_respected() {
    let unit = lcg_dna(1000, 9);
    let mut q = lcg_dna(500, 1);
    q.extend_from_slice(&unit);
    let mut s = unit[..500].to_vec();
    s.extend_from_slice(b"NNNNNNNNNNNNNNNNNNNNNNNNN");
    s.extend_from_slice(&unit[525..]);
    // Second subject record starts mid-copy.
    let s_breaks = [700usize];
    let mut cfg = large_cfg(1, Some(15));
    cfg.chain = Some(ChainConfig {
        min_score: 50.0,
        ..ChainConfig::default()
    });
    let l = compute_large_dotplot(&LargeInput::pair(&q, &[], &s, &s_breaks), &cfg).unwrap();
    assert!(!l.segments.is_empty());
    for seg in &l.segments {
        let rec_lo = (seg.s_start as usize) >= 700;
        let rec_hi = (seg.s_end as usize - 1) >= 700;
        assert_eq!(rec_lo, rec_hi, "segment crosses record break: {seg:?}");
        assert_eq!(seg.s_rec, rec_lo as u32);
    }
    assert!(l
        .chains
        .iter()
        .all(|c| (c.s_start >= 700) == (c.s_end > 700)));
}

#[test]
fn repeat_classes_partition_the_unclassified_plot() {
    let te = lcg_dna(400, 4);
    let mut seq = lcg_dna(30_000, 3);
    for i in 0..12 {
        let at = 1000 + i * 2300;
        seq.splice(at..at + 400, mutate(&te, 3, i as u64));
    }
    let unit = lcg_dna(1500, 8);
    seq.splice(5000..6500, unit.iter().copied());
    seq.splice(20_000..21_500, mutate(&unit, 5, 2));

    let mut cfg = large_cfg(20, None);
    cfg.self_comparison = true;
    let flat = compute_large_dotplot(&LargeInput::self_comparison(&seq, &[]), &cfg).unwrap();
    cfg.depth = DepthMethod::Kmer { k: 15 };
    cfg.classes = ClassThresholds {
        low: 2,
        mid: 5,
        high: 8,
    };
    let classed = compute_large_dotplot(&LargeInput::self_comparison(&seq, &[]), &cfg).unwrap();
    assert_eq!(
        flat.merged(None, ClassMask::ALL),
        classed.merged(None, ClassMask::ALL)
    );
    // The 12-copy TE lands in `high`; the 2-copy duplication does not.
    let high = classed.merged(None, ClassMask::from_classes(&[RepeatClass::High]));
    let rest = classed.merged(None, ClassMask::ALL.with(RepeatClass::High, false));
    let w = classed.width as usize;
    let (dx, dy) = (5500 / 20, 20_500 / 20);
    assert!(rest[dy * w + dx] > 0 && high[dy * w + dx] == 0);
    let (tx, ty) = ((1000 + 200) / 20, (1000 + 2300 * 5 + 200) / 20);
    assert!(high[ty * w + tx] > 0);
}

#[test]
fn seed_depth_sees_diverged_repeats_that_kmers_miss() {
    // A 20-copy family at ~18 % divergence: almost no shared 21-mers.
    let te = lcg_dna(1500, 40);
    let mut seq = lcg_dna(80_000, 41);
    for i in 0..20u64 {
        let at = 2000 + i as usize * 3900;
        let copy = mutate(&te, 18, 100 + i);
        seq.splice(at..at + 1500, copy);
    }
    let breaks: [usize; 0] = [];
    let input = LargeInput::self_comparison(&seq, &breaks);
    let mut cfg = large_cfg(1, None);
    cfg.self_comparison = true;
    let fraction_high = |cfg: &LargeConfig| {
        let s = LargeSession::prepare(&input, cfg, &Control::default()).unwrap();
        let d = &s.query_depth().depth;
        let in_te = (0..20).flat_map(|i| 2000 + i * 3900 + 100..2000 + i * 3900 + 1400);
        let n = in_te.clone().count();
        let hi = in_te
            .filter(|&p| cfg.classes.classify(d[p]) >= RepeatClass::Mid)
            .count();
        let bg = (1..20)
            .flat_map(|i| 2000 + i * 3900 - 2000..2000 + i * 3900 - 600)
            .filter(|&p| cfg.classes.classify(d[p]) >= RepeatClass::Low)
            .count();
        (hi as f64 / n as f64, bg)
    };
    cfg.depth = DepthMethod::Kmer { k: 21 };
    let (kmer, _) = fraction_high(&cfg);
    cfg.depth = DepthMethod::Seed;
    let (seed, bg) = fraction_high(&cfg);
    assert!(kmer < 0.3, "kmer depth unexpectedly sensitive: {kmer}");
    assert!(
        seed > 0.8,
        "seed depth classified only {seed} of TE bases as repeats"
    );
    assert!(bg < 500, "{bg} background bases called repetitive");
}

#[test]
fn capped_seeds_are_reported() {
    let mut seq = lcg_dna(20_000, 1);
    let sat = b"ACGTTGCAAGGCTTACGATCGA".repeat(200);
    seq.splice(5000..5000 + sat.len(), sat.iter().copied());
    let mut cfg = large_cfg(100, None);
    cfg.self_comparison = true;
    cfg.max_occ = 50;
    let l = compute_large_dotplot(&LargeInput::self_comparison(&seq, &[]), &cfg).unwrap();
    assert!(l.stats.capped_query_positions > 3000);
    assert!(l.q_capped[52] > 0 && l.q_capped[10] == 0);
    assert!(l.s_capped[52] > 0 && l.s_capped[10] == 0);
}

#[test]
fn session_reuse_and_cancel() {
    let (q, s) = planted_pair(20_000, 2000, 10, 1);
    let cfg = large_cfg(10, None);
    let input = LargeInput::pair(&q, &[], &s, &[]);
    let session = LargeSession::prepare(&input, &cfg, &Control::default()).unwrap();
    let a = session.compute(&cfg, &Control::default()).unwrap();
    let b = compute_large_dotplot(&input, &cfg).unwrap();
    assert_eq!(a.channels, b.channels);

    let mut other = cfg.clone();
    other.seed = SpacedSeed::contiguous(11).unwrap();
    assert!(!session.compatible(&other));
    assert!(session.compute(&other, &Control::default()).is_err());

    let flag = AtomicBool::new(true);
    let ctl = Control {
        cancel: Some(&flag),
        progress: None,
    };
    assert!(matches!(
        session.compute(&cfg, &ctl),
        Err(DottirError::Cancelled)
    ));
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(24))]

    #[test]
    fn prop_large_le_exact(
        seed in 0u64..10_000,
        rate in 0u64..25,
        zoom in 1u32..7,
        window in prop::option::of(8u32..40),
        unit_len in 100usize..900,
    ) {
        let (q, s) = planted_pair(3000, unit_len, rate, seed);
        let e = exact(&q, &s, zoom, window, false);
        let mut cfg = large_cfg(zoom, window);
        cfg.seed = SpacedSeed::parse("11011011011").unwrap();
        let l = large(&q, &s, &cfg);
        assert_le(&l.merged(Some(false), ClassMask::ALL), &e.pixels, "forward");
        assert_le(&l.merged(Some(true), ClassMask::ALL), e.reverse_pixels.as_ref().unwrap(), "reverse");
    }
}

/// Pinned output for a fixed input. A change here is an output change:
/// bump `LARGE_PIXELMAP_FORMAT_VERSION`, note it in the changelog, and
/// re-pin (the assertion message prints the new value).
#[test]
fn golden_channels_hash() {
    let (q, s) = planted_pair(40_000, 4000, 12, 2024);
    let mut cfg = large_cfg(16, None);
    cfg.depth = DepthMethod::Kmer { k: 15 };
    cfg.chain = Some(ChainConfig::default());
    let l = large(&q, &s, &cfg);
    // FNV-1a over every channel, then chain coordinates.
    let mut h: u64 = 0xcbf29ce484222325;
    let mut eat = |b: u8| {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    };
    for ch in &l.channels {
        ch.iter().for_each(|&b| eat(b));
        eat(0xff);
    }
    for c in &l.chains {
        for v in [c.q_start, c.q_end, c.s_start, c.s_end] {
            v.to_le_bytes().iter().for_each(|&b| eat(b));
        }
    }
    assert_eq!(LARGE_PIXELMAP_FORMAT_VERSION, 1);
    assert_eq!(h, GOLDEN, "large-mode output changed: new hash {h:#018x}");
}

const GOLDEN: u64 = 0x7f156c95cbc4c727;
