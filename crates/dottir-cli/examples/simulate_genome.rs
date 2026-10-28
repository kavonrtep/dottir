//! Deterministic synthetic genomes for large-mode benchmarks (plan L0).
//!
//! ```text
//! cargo run --release -p dottir-cli --example simulate_genome -- \
//!     OUT_PREFIX LEN_MB [RECORDS] [SEED] [DIVERGENCE_PCT]
//! ```
//!
//! Writes `OUT_PREFIX.a.fa` (genome A: random background plus TE families
//! of 5–500 copies at 0–25 % divergence, a tandem satellite array and
//! segmental duplications), `OUT_PREFIX.b.fa` (genome B: A rearranged —
//! inversions, a translocation — and diverged by `DIVERGENCE_PCT`), and
//! `OUT_PREFIX.truth.bed` with the planted events in A coordinates.
//!
//! For genome mode it also writes `OUT_PREFIX.contigs.fa`: B cut into
//! contigs of 20–200 % of `len / (10 · RECORDS)`, 1 % diverged, shuffled,
//! a third reverse-complemented, plus unplaced random contigs; and
//! `OUT_PREFIX.contigs.tsv` (contig, B record, start, end, strand; `*` for
//! unplaced). Sorting the contigs along B must recover the B order.

use std::io::Write;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let mut x = self.0;
        x ^= x >> 33;
        x = x.wrapping_mul(0xff51afd7ed558ccd);
        x ^ (x >> 33)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
    fn base(&mut self) -> u8 {
        b"ACGT"[self.below(4) as usize]
    }
    fn dna(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| self.base()).collect()
    }
    /// Substitutions plus ~10 % as many 1-bp indels.
    fn mutate(&mut self, seq: &[u8], pct: f64) -> Vec<u8> {
        let p = (pct / 100.0 * 1_000_000.0) as u64;
        let mut out = Vec::with_capacity(seq.len());
        for &b in seq {
            if self.below(1_000_000) < p {
                match self.below(20) {
                    0 => {}
                    1 => {
                        out.push(b);
                        out.push(self.base());
                    }
                    _ => out.push(self.base()),
                }
            } else {
                out.push(b);
            }
        }
        out
    }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

fn write_fasta(path: &str, records: &[(String, Vec<u8>)]) -> std::io::Result<()> {
    let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);
    for (id, seq) in records {
        writeln!(w, ">{id}")?;
        for line in seq.chunks(80) {
            w.write_all(line)?;
            w.write_all(b"\n")?;
        }
    }
    w.flush()
}

fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: simulate_genome OUT_PREFIX LEN_MB [RECORDS] [SEED] [DIVERGENCE_PCT]");
        std::process::exit(2);
    }
    let prefix = &args[1];
    let len = (args[2].parse::<f64>().expect("LEN_MB") * 1e6) as usize;
    let n_rec: usize = args.get(3).map_or(4, |s| s.parse().expect("RECORDS"));
    let mut rng = Rng(args.get(4).map_or(1, |s| s.parse().expect("SEED")));
    let div: f64 = args
        .get(5)
        .map_or(10.0, |s| s.parse().expect("DIVERGENCE_PCT"));

    let mut a = rng.dna(len);
    let mut truth = Vec::new();

    // TE families: (unit length, copies, divergence %). Copies overwrite
    // background at random places, a third of them inverted.
    let scale = (len as f64 / 1e7).max(0.01);
    let families = [
        (6000usize, (40.0 * scale) as usize, 2.0),
        (3000, (200.0 * scale) as usize, 8.0),
        (1500, (500.0 * scale) as usize, 18.0),
        (8000, (5.0 * scale).max(2.0) as usize, 25.0),
    ];
    for (fi, &(ulen, copies, d)) in families.iter().enumerate() {
        let unit = rng.dna(ulen);
        for _ in 0..copies.max(2) {
            let mut c = rng.mutate(&unit, d);
            let inv = rng.below(3) == 0;
            if inv {
                c = revcomp(&c);
            }
            let at = rng.below((len - c.len()) as u64) as usize;
            a[at..at + c.len()].copy_from_slice(&c);
            truth.push(format!(
                "{at}\t{}\tTE{fi}\t{d}\t{}",
                at + c.len(),
                if inv { '-' } else { '+' }
            ));
        }
    }
    // Satellite array: 180-bp monomers, 2 % divergence, 100 kb.
    let mono = rng.dna(180);
    let sat_len = (100_000).min(len / 20);
    let at = len / 3;
    let mut sat = Vec::new();
    while sat.len() < sat_len {
        sat.extend(rng.mutate(&mono, 2.0));
    }
    sat.truncate(sat_len);
    a[at..at + sat_len].copy_from_slice(&sat);
    truth.push(format!("{at}\t{}\tsatellite\t2\t+", at + sat_len));
    // Segmental duplications, 50 kb at 5 % divergence.
    for k in 0..4 {
        let l = 50_000.min(len / 50);
        let from = rng.below((len - l) as u64) as usize;
        let to = rng.below((len - l) as u64) as usize;
        let mut c = rng.mutate(&a[from..from + l], 5.0);
        c.truncate(l);
        let inv = k % 2 == 1;
        if inv {
            c = revcomp(&c);
        }
        let n = c.len();
        a[to..to + n].copy_from_slice(&c);
        truth.push(format!(
            "{to}\t{}\tsegdup_of_{from}\t5\t{}",
            to + n,
            if inv { '-' } else { '+' }
        ));
    }

    // Genome B: A cut into blocks, some inverted, two swapped, diverged.
    let n_blocks = 20;
    let bl = len / n_blocks;
    let mut blocks: Vec<Vec<u8>> = (0..n_blocks)
        .map(|i| a[i * bl..if i + 1 == n_blocks { len } else { (i + 1) * bl }].to_vec())
        .collect();
    for i in [3usize, 11, 16] {
        blocks[i] = revcomp(&blocks[i]);
    }
    blocks.swap(6, 13);
    let b: Vec<u8> = blocks.concat();
    let b = rng.mutate(&b, div);

    let split = |seq: &[u8], name: &str| -> Vec<(String, Vec<u8>)> {
        let step = seq.len().div_ceil(n_rec);
        seq.chunks(step)
            .enumerate()
            .map(|(i, c)| (format!("{name}{}", i + 1), c.to_vec()))
            .collect()
    };
    let b_recs = split(&b, "chrB");
    write_fasta(&format!("{prefix}.a.fa"), &split(&a, "chrA"))?;
    write_fasta(&format!("{prefix}.b.fa"), &b_recs)?;

    let mean = (len / (10 * n_rec)).max(1000) as u64;
    let mut contigs = Vec::new();
    for (rec, seq) in &b_recs {
        let mut at = 0usize;
        while at < seq.len() {
            let l = (mean / 5 + rng.below(mean * 9 / 5)) as usize;
            let end = (at + l).min(seq.len());
            let rev = rng.below(3) == 0;
            let mut c = rng.mutate(&seq[at..end], 1.0);
            if rev {
                c = revcomp(&c);
            }
            let strand = if rev { '-' } else { '+' };
            contigs.push((format!("{rec}\t{at}\t{end}\t{strand}"), c));
            at = end;
        }
    }
    for _ in 0..(contigs.len() / 10).max(1) {
        let l = (mean / 10 + rng.below(mean / 5)) as usize;
        contigs.push(("*\t*\t*\t*".to_string(), rng.dna(l)));
    }
    for i in (1..contigs.len()).rev() {
        contigs.swap(i, rng.below(i as u64 + 1) as usize);
    }
    let mut t = std::fs::File::create(format!("{prefix}.contigs.tsv"))?;
    let mut named = Vec::with_capacity(contigs.len());
    for (i, (truth, seq)) in contigs.into_iter().enumerate() {
        let id = format!("ctg{:05}", i + 1);
        writeln!(t, "{id}\t{truth}")?;
        named.push((id, seq));
    }
    write_fasta(&format!("{prefix}.contigs.fa"), &named)?;
    let mut t = std::fs::File::create(format!("{prefix}.truth.bed"))?;
    for line in truth {
        writeln!(t, "genomeA\t{line}")?;
    }
    eprintln!(
        "wrote {prefix}.a.fa ({} bp), {prefix}.b.fa ({} bp)",
        a.len(),
        b.len()
    );
    Ok(())
}
