#!/bin/bash
set -euo pipefail

# The repo pins `rust-toolchain.toml` to a specific channel for dev,
# which makes cargo try to invoke rustup. The conda build env has
# the `rust` package on PATH instead — so move the override aside
# for the duration of the build.
[ -f rust-toolchain.toml ] && mv rust-toolchain.toml rust-toolchain.toml.bak

cargo build --release --workspace --locked

install -m 0755 target/release/dottir     "${PREFIX}/bin/dottir"
install -m 0755 target/release/dottir-gui "${PREFIX}/bin/dottir-gui"
