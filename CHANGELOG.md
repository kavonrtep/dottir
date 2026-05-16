# Changelog

All notable changes to dottir are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/);
the project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
* Cargo workspace with `dottir-core`, `dottir-io`, `dottir-cli`, and
  `dottir-gui` crates.
* Faithful port of Karlin/Altschul λ/K/H statistics
  (`dottir-core::karlin`), bit-identical to a standalone C harness
  extracted from `dotterKarlin.c`.
* Built-in score matrices: BLOSUM62 (inline from C-dotter), BLOSUM45,
  BLOSUM50, BLOSUM80, BLOSUM90, PAM30, PAM70, PAM250 (vendored from
  NCBI), and DNA identity (`+5 / -4`).
* BLAST-format text matrix parser and round-trip serialiser.
* Sliding-window dotplot kernel
  (`compute_dotplot`) supporting:
  - BLASTN forward, reverse, both-strand;
  - BLASTP forward;
  - self-comparison with `Triangle::{Both, Upper, Lower}` and an
    optional `disable_mirror`;
  - separate forward/reverse channels for inverted-repeat highlighting.
* rayon-parallel driver with byte-identical determinism across thread
  counts.
* `dottir-io::fasta` minimal FASTA reader (plain + gzip).
* `dottir-io::png_export` 8-bit greyscale PNG writer with `tEXt`
  provenance chunks.
* `dottir-io::params` TOML sidecar (SHA-256 hashes + resolved
  parameters + Karlin numbers + host info).
* `dottir-io::alignment` ±N-residue slice helper around a crosshair
  coord, with FASTA-pair / Stockholm / plain output.
* `dottir batch` CLI with flags mirroring upstream Dotter where
  reasonable; `--auto-zoom` to fit a target max dimension.
* GitHub Actions workflow: fmt, clippy, test on Linux + Windows.
* mdBook scaffold at `docs/book/`.
* `cargo-deny` configuration.

### Deferred
* egui/eframe GUI runtime — blocked on Rust 1.85 MSRV bump
  (see `docs/adr/0003-gui-msrv.md`).
* GFF3 / PAF loaders via `noodles-*` — same MSRV blocker
  (see `docs/adr/0004-defer-gff3-paf.md`).
* SVG / PDF / `.dot` exports — Phase 4-extra.
* BLASTX three-frame translation — Phase 2-extra.

### Goldens
Karlin λ/K/H/window pinned bit-identical to the C reference across four
fixtures (`tests/golden/karlin/values.tsv`). Pixelmap goldens are not
yet pinned against C-dotter output (requires building the GTK2 binary
in `third_party/seqtools/`).
