# Repository Guidelines

## Project Structure & Module Organization
This is a Rust workspace. Core code lives under `crates/`:
- `crates/dottir-core/` holds the I/O-free plotting algorithms and unit/property tests.
- `crates/dottir-io/` contains FASTA, alignment, parameter, and export helpers.
- `crates/dottir-cli/` builds the `dottir` batch binary.
- `crates/dottir-gui/` builds the `dottir-gui` egui/eframe app.

Supporting material lives in `docs/book/` and `docs/adr/`. Golden regression fixtures are under `tests/golden/`, with regeneration helpers in `tests/golden_gen/`.

## Build, Test, and Development Commands
- `cargo build --workspace` compiles all crates.
- `cargo run -p dottir-cli -- --help` runs the CLI binary.
- `cargo run -p dottir-gui` launches the GUI.
- `cargo test --workspace` runs the full test suite, including crate tests.
- `cargo clippy --workspace --all-targets -- -D warnings` enforces lint cleanliness.
- `cargo fmt --check` verifies formatting against `rustfmt.toml`.

## Coding Style & Naming Conventions
Use standard Rust formatting: 4-space indentation, `snake_case` for functions/modules, `CamelCase` for types, and `SCREAMING_SNAKE_CASE` for constants. Keep `dottir-core` free of I/O and GUI dependencies. Prefer explicit, short comments only when the algorithmic intent is not obvious. Avoid unordered iteration in code that affects pixelmap output.

## Testing Guidelines
Add or update tests alongside behavior changes. Use crate-local tests for fast checks and golden tests for output-sensitive algorithm work. Changes to the sliding-window recurrence, anti-diagonal suppression, or pixel merge rules should include golden coverage. If a change affects the scientific contract, regenerate fixtures intentionally and document why.

## Commit & Pull Request Guidelines
Recent history uses short, imperative subjects with scoped prefixes such as `Phase 5: ...` or `chore: ...`. Follow that style: keep the subject specific and include context in the body when needed. Pull requests should explain the change, list validation performed, and call out any algorithmic impact. GUI changes should include screenshots or a brief screen recording. Core algorithm changes should mention any ADR updates and golden-file regeneration.

## Agent-Specific Instructions
Before changing algorithmic behavior, read `docs/dottir_specification.md` and `CLAUDE.md`. For `dottir-core`, prioritize determinism and scientific correctness over convenience. Do not introduce new dependencies lightly; prefer the workspace’s existing crates unless there is a clear justification.
