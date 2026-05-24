#!/bin/bash
set -euo pipefail

# The repo pins `rust-toolchain.toml` to a specific channel for dev,
# which makes cargo try to invoke rustup. The conda build env has
# the `rust` package on PATH instead — so move the override aside
# for the duration of the build.
[ -f rust-toolchain.toml ] && mv rust-toolchain.toml rust-toolchain.toml.bak

# Bake RPATH=$ORIGIN/../lib at link time so dottir-gui finds the
# conda env's X11/Wayland libs without needing conda-build's
# install-time path-relocation step (which is disabled in
# meta.yaml — see the comment on detect_binary_files_with_prefix
# for why). The literal `$ORIGIN` token must survive into the
# linker invocation unchanged — `\$ORIGIN` at the bash level
# gives us `$ORIGIN` literally.
export RUSTFLAGS="-C link-args=-Wl,-rpath,\$ORIGIN/../lib"

# Smaller package, no debug info to ship.
export CARGO_PROFILE_RELEASE_DEBUG=0
export CARGO_PROFILE_RELEASE_STRIP=symbols

cargo build --release --workspace --locked

install -m 0755 target/release/dottir     "${PREFIX}/bin/dottir"
install -m 0755 target/release/dottir-gui "${PREFIX}/bin/dottir-gui"
