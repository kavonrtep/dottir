# ADR 0002 — License under GPLv3-or-later

* Status: Accepted
* Date: 2026-05-16
* Deciders: petr
* Spec reference: §2, §10 question #2

## Context

The original Dotter is GPLv3. Spec §2 lists two paths for the Rust port:

1. Adopt **GPLv3** for the whole project. Simplest and removes any doubt
   about derivative status, particularly for `dotterApp/dotterKarlin.c`
   which we port more-or-less verbatim.
2. Adopt a permissive license (MIT / Apache-2.0) and reimplement Karlin/
   Altschul cleanly from Karlin & Altschul 1990 and the public-domain NCBI
   BLAST sources.

## Decision

We license dottir under **GPLv3-or-later** for v1.0.

## Rationale

* The Karlin/Altschul implementation in `crates/dottir-core/src/karlin.rs`
  is a structural port of the GPL'd reference code (loop shape, MAXIT,
  SUMLIMIT, fudge-to-0.1 behaviour all carry over). Trying to argue this
  is "clean-room" would be specious.
* Sequence-analysis users overwhelmingly distribute and consume GPL tools.
  Compatibility with `seqtools`, `samtools`, `htslib`, `BLAST+` matters
  more than the marginal adoption boost from a permissive license.
* If permissive licensing becomes load-bearing later (e.g. for inclusion
  in a non-GPL pipeline), Karlin can be reimplemented from primary sources
  in a separate crate at that point and the rest of the codebase can be
  relicensed with author consent — it is clean-room except for that one
  module.

## Consequences

* `LICENSE` at repo root is GPLv3.
* `Cargo.toml` workspace declares `license = "GPL-3.0-or-later"`.
* `NOTICE` credits Sonnhammer, Durbin, Barson, Scofield, and NCBI.
* Dependencies must be GPL-compatible. `cargo-deny` (Phase 9) will enforce
  this.
