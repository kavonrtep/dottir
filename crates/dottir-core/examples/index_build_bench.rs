//! Times `SeedIndex::build` on a FASTA file (records concatenated, breaks
//! at record starts). Usage: `index_build_bench SEQ.fa [REPEATS]`.

use std::time::Instant;

use dottir_core::alphabet::{encode, AlphabetKind};
use dottir_core::large::{SeedIndex, SpacedSeed};

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args
        .next()
        .expect("usage: index_build_bench SEQ.fa [REPEATS]");
    let repeats: usize = args.next().map_or(3, |r| r.parse().expect("REPEATS"));
    let text = std::fs::read_to_string(&path).expect("read FASTA");
    let mut seq = Vec::new();
    let mut breaks = Vec::new();
    for line in text.lines() {
        if line.starts_with('>') {
            if !seq.is_empty() {
                breaks.push(seq.len());
            }
        } else {
            seq.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    let enc = encode(&seq, AlphabetKind::Dna);
    let seed = SpacedSeed::default();
    println!("{} bases, {} records", enc.len(), breaks.len() + 1);
    for _ in 0..repeats {
        let t = Instant::now();
        let idx = SeedIndex::build(&seed, &enc, &breaks).unwrap();
        println!(
            "build {:.2} s, {} positions",
            t.elapsed().as_secs_f64(),
            idx.n_positions()
        );
    }
}
