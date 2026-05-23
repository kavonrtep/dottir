#!/bin/bash
set -euo pipefail

# The repo pins `rust-toolchain.toml` to a specific channel for dev,
# which makes cargo try to invoke rustup. The conda build env has
# the `rust` package on PATH instead — so move the override aside
# for the duration of the build.
[ -f rust-toolchain.toml ] && mv rust-toolchain.toml rust-toolchain.toml.bak

# Conda-forge's rustc + Cargo.toml's `lto = "thin"` miscompiles
# x11-dl's static lookup table: at runtime dlopen() gets called with
# garbage strings (substrings of the symbol table) instead of
# "libX11.so.6", and dottir-gui dies with "Failed to load one of
# xlib's shared libraries". Force-disabling LTO for the conda build
# avoids the miscompile without changing the workspace profile.
export CARGO_PROFILE_RELEASE_LTO=off
export CARGO_PROFILE_RELEASE_CODEGEN_UNITS=16

cargo build --release --workspace --locked

install -m 0755 target/release/dottir     "${PREFIX}/bin/dottir"
install -m 0755 target/release/dottir-gui "${PREFIX}/bin/dottir-gui"
