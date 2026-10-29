//! Time the PAF readers on a real file:
//! `cargo run --release -p dottir-io --example paf_read_bench -- ALN.paf[.gz] QUERY.fai TARGET.fai`
//! (record names and lengths come from the samtools `.fai` indexes).

use std::path::Path;
use std::time::Instant;

use dottir_core::genome::GenomeRecord;
use dottir_io::paf::{read_paf_file, read_paf_subset_oriented};

fn records(fai: &str) -> Vec<GenomeRecord> {
    std::fs::read_to_string(fai)
        .expect("reading .fai")
        .lines()
        .map(|l| {
            let mut f = l.split('\t');
            GenomeRecord {
                name: f.next().unwrap().to_string(),
                len: f.next().unwrap().parse().unwrap(),
            }
        })
        .collect()
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [paf, q, t] = args.as_slice() else {
        eprintln!("usage: paf_read_bench ALN.paf QUERY.fai TARGET.fai");
        std::process::exit(2);
    };
    let (q, t) = (records(q), records(t));
    let path = Path::new(paf);

    let start = Instant::now();
    let alns = read_paf_file(path, &q, &t, 0).expect("read_paf_file");
    println!(
        "read_paf_file: {} alignments in {:.2}s",
        alns.len(),
        start.elapsed().as_secs_f64()
    );
    let mut alns = alns;
    for a in alns.iter_mut() {
        a.primary = false;
    }
    let start = Instant::now();
    dottir_core::genome::alignment::assign_primary(
        &mut alns,
        dottir_core::genome::alignment::MASK_LEVEL,
    );
    println!(
        "assign_primary: {} primaries in {:.2}s",
        alns.iter().filter(|a| a.primary).count(),
        start.elapsed().as_secs_f64()
    );
    drop(alns);

    let start = Instant::now();
    let (alns, skipped, swapped, _) =
        read_paf_subset_oriented(path, &t, &q).expect("read_paf_subset_oriented");
    println!(
        "read_paf_subset_oriented: {} alignments ({skipped} skipped, swapped {swapped}) in {:.2}s",
        alns.len(),
        start.elapsed().as_secs_f64()
    );
}
