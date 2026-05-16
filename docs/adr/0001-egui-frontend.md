# ADR 0001 — Use egui / eframe as the GUI framework

* Status: Accepted
* Date: 2026-05-16
* Deciders: petr
* Spec reference: §4.3, §10 question #6

## Context

The original Dotter is a GTK2 application. GTK2 is end-of-life, its bindings
in Rust are deprecated, and packaging it on Windows or as a single binary is
painful. We need a cross-platform GUI framework that ships as part of the
binary (no system GTK), supports Linux/Windows/macOS, and ideally allows a
later WebAssembly target without rewriting the UI layer.

Candidates considered:

1. **egui / eframe** — immediate-mode UI, pure Rust, single static binary,
   WebAssembly support out of the box. Used by 1Password and Rerun. Pixmap
   rendering via `TextureHandle` is fast and ergonomic.
2. **iced** — retained-mode (Elm-style) UI. Polished, but the WASM story is
   weaker and immediate-mode is a better fit for a real-time crosshair +
   live greyramp.
3. **gtk4-rs** — modern GTK; would mirror the upstream architecture more
   closely. Trades portability and binary-size for native look-and-feel on
   Linux. Windows builds require shipping the GTK runtime.
4. **Tauri** — HTML/JS frontend in a webview. Excellent for app UI; poor
   fit for a high-FPS interactive plot with a custom canvas.

## Decision

We use **egui + eframe**.

## Consequences

* No dependency on system GTK; `cargo build --release` produces a self-
  contained binary on all major desktop platforms.
* Live greyramp updates are cheap: rebuild the texture from the existing
  pixelmap via the LUT, no recomputation.
* WebAssembly target is essentially free (Phase 8). Multithreading on WASM
  is limited; the parallel kernel falls back to single-thread on WASM.
* Look-and-feel will differ from native widgets. For a domain tool, this
  is acceptable.
* `egui_kittest` enables limited GUI smoke tests (spec §4.5.8).
