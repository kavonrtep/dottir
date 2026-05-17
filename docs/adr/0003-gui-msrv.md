# ADR 0003 — Defer egui/eframe runtime to a future MSRV bump

* Status: **Superseded** by commit
  [cd3924e](https://github.com/petr/dottir/commit/cd3924e) ("Bump
  MSRV to 1.85 and unpin transitive deps", 2026-05-17). Phase 5
  GUI MVP shipped in [cc1fbf1](https://github.com/petr/dottir/commit/cc1fbf1)
  (egui/eframe at 0.29), light theme + CLI pre-load follow-ups in
  [0de5f72](https://github.com/petr/dottir/commit/0de5f72) and
  [436a809](https://github.com/petr/dottir/commit/436a809).
* Date: 2026-05-16
* Deciders: petr

## Context

ADR 0001 fixed `egui/eframe` as the GUI framework for Phase 5+. Phase 5
of `IMPLEMENTATION_PLAN.md` budgets 2-3 weeks to land the interactive
MVP.

When we tried to add `eframe = "0.27.2"` to `crates/dottir-gui/Cargo.toml`
the dependency resolution failed:

```
eframe 0.27.2
  → winit 0.29.x
    → toml_edit 0.25
      → indexmap 2.14
        requires Cargo feature `edition2024` (Rust 1.85+)
```

Our workspace `rust-toolchain.toml` pins **Rust 1.75** (the version
installed on the development host as `/usr/bin/rustc`). edition2024
stabilised in Rust 1.85. Even after pinning `winit = 0.29.4` and
`toml_edit = 0.22.22` the chain still reaches `indexmap-2.14` through
other paths.

## Decision

We **defer the interactive egui GUI** until a separate change bumps the
workspace MSRV. The Phase 5 acceptance criterion ("`cargo run -p
dottir-gui` reproduces the interactive behaviour") is not met by this
release.

The `dottir-gui` binary is repurposed as a thin headless wrapper around
`dottir-core` and `dottir-io::png_export` (functionally a subset of
`dottir-cli batch`), so the crate still does something useful.

## Consequences

* Phase 5 deliverables (pan/zoom canvas, greyramp panel, crosshair,
  alignment view, settings dialog) are NOT done.
* The headless `dottir-cli` path stays the supported way to compute
  dotplots until the MSRV bump.
* When MSRV is bumped (proposed: target Rust 1.85 stable), this ADR is
  superseded and a follow-up ADR documents the new floor. At that point
  the workspace's transitive-dep pins (clap, toml, proptest, indexmap,
  tempfile, rayon-core) can be relaxed too.

## Alternatives considered

1. **Pin every transitive dep through `cargo update --precise`.** This
   was attempted and abandoned: the chain forks at multiple points (each
   `egui-*` sub-crate, plus winit/wgpu/glow), and each new release of
   one of them re-introduces an edition2024 dep that needs pinning. The
   maintenance cost outweighs the benefit while a simple toolchain bump
   would solve it cleanly.

2. **Use iced or a different GUI framework.** Reverses ADR 0001 just to
   avoid an MSRV bump, which is the wrong trade.

3. **Bump the MSRV now as part of this change.** Tempting, but it
   touches every CI pipeline and is a decision orthogonal to dottir's
   algorithmic work. Better as its own change with its own review.
