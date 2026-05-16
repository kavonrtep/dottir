# Installation

dottir is a Rust project. There are no system-library dependencies at
runtime; everything needed to compute a dotplot ships in the binaries.

## From source

```sh
git clone https://github.com/petr/dottir
cd dottir
cargo build --release
# binaries land at target/release/dottir and target/release/dottir-gui
```

The workspace pins **Rust 1.75** as MSRV at the time of writing. A
future release will bump this past 1.85 to unblock the egui GUI; see
ADR 0003.

## From a release binary

Pending Phase 9 of `IMPLEMENTATION_PLAN.md`. The GitHub Actions
release workflow lives at `.github/workflows/release.yml` and will
produce static Linux + Windows artifacts on tag push.
