# Code Review

This review is based on the current source tree, `IMPLEMENTATION_PLAN.md`, `dottir_specification.md`, and the docs under `docs/book/`. I validated the workspace with Rust 1.95.0 (`conda run -n rust_1.95.0 cargo check --workspace` and `cargo test --workspace`), so the notes below focus on product/architecture gaps rather than basic build breakage.

## Summary

The core algorithmic work is in good shape: Karlin/Altschul statistics, the sliding-window kernel, anti-diagonal suppression, deterministic parallelism, and the main test suite are all implemented and passing. The main risks are now around scalability, completeness relative to the plan, and consistency between docs and code.

## Findings

1. Parallel execution can blow past the memory limit. In [`plot.rs`](/home/petr/PycharmProjects/dottir/crates/dottir-core/src/plot.rs#L342), each rayon chunk allocates a full-size `PixelMap`, and all chunk maps are collected before the final merge. The limit check only guards each individual map, not aggregate peak memory. On large inputs this can multiply memory use by the number of chunks and defeat the advertised cap.
2. Custom matrix parsing is too permissive and only handles protein matrices. [`matrix.rs`](/home/petr/PycharmProjects/dottir/crates/dottir-core/src/matrix.rs#L279) silently drops unsupported header letters and fills missing cells with `-4`. That makes malformed matrix files easy to accept accidentally, and it prevents the same parser from supporting DNA matrices or stricter validation for user-provided files.
3. Multi-record FASTA boundaries are discarded too early. [`fasta.rs`](/home/petr/PycharmProjects/dottir/crates/dottir-io/src/fasta.rs#L128) concatenates records into a single byte vector, so the rest of the pipeline loses the offsets needed for breaklines, boundary-aware exports, and track rendering on concatenated assemblies.
4. The GUI is still an MVP, not the feature set described in the plan/spec. [`app.rs`](/home/petr/PycharmProjects/dottir/crates/dottir-gui/src/app.rs#L331) currently exposes file open, save PNG, greyramp, a basic settings dialog, and pan/zoom. The planned alignment pane, GFF3/PAF tracks, sub-dotter spawning, session save/load, export actions, and breakline rendering are still absent.
5. The GUI hardcodes a 1 GiB memory cap and does not expose the spec’s configurable limit. [`app.rs`](/home/petr/PycharmProjects/dottir/crates/dottir-gui/src/app.rs#L197) sets `memory_limit_bytes` directly instead of surfacing the setting in the UI or reusing the CLI’s default.
6. FASTA loading does redundant work. [`main.rs`](/home/petr/PycharmProjects/dottir/crates/dottir-cli/src/main.rs#L172) and [`fasta.rs`](/home/petr/PycharmProjects/dottir/crates/dottir-io/src/fasta.rs#L49) both read whole files eagerly, and the CLI path reads each file twice in the common case. That is acceptable for small inputs but wasteful for the larger assemblies the spec explicitly targets.
7. The docs are stale relative to the actual codebase. [`intro.md`](/home/petr/PycharmProjects/dottir/docs/book/src/intro.md#L38) and [`install.md`](/home/petr/PycharmProjects/dottir/docs/book/src/install.md#L15) still describe GUI support as deferred and claim an MSRV of Rust 1.75, but the workspace now pins Rust 1.95.0 and the GUI crate builds. The documentation should be updated to match the shipped state.

## Missing Features

Compared with `IMPLEMENTATION_PLAN.md` and `dottir_specification.md`, the remaining gaps are substantial:

- BLASTX support is still `NotImplemented`.
- The CLI only writes PNG plus a params sidecar; it does not yet emit SVG, PDF, or `.dot`.
- The CLI cannot load/save `.dot` or accept GFF3/PAF overlays.
- The GUI lacks the annotation track panel, alignment dock, selection/export flow, and self-comparison inverted-repeat visualization.
- The sequence model does not preserve record boundaries for breaklines or coordinate-aware multi-record views.
- There is no shared session file format for round-tripping full GUI state.

## How To Improve It

- Fix the parallel memory model first. Prefer bounded chunking that reuses buffers or merges incrementally instead of materializing one full pixelmap per chunk.
- Tighten matrix loading. Separate protein and DNA parsing, validate all expected cells, and fail loudly on malformed input instead of silently substituting `-4`.
- Preserve sequence metadata through the I/O layer. Keep record boundaries, offsets, and source labels so breaklines and exports can be implemented without re-parsing.
- Consolidate CLI and GUI configuration around a shared state struct so the supported feature surface stays in sync.
- Add an explicit export layer for PNG/SVG/PDF/.dot and a separate persistence layer for session state and provenance.
- Refresh the docs and phase notes after the code changes land. The current book is useful, but it no longer matches the implementation status.

## Zoom Quality

The current GUI zoom is a pure viewport transform over one computed pixelmap. That is useful for smooth navigation, but it is not equivalent to the original Dotter behavior or to a recomputed higher-resolution view. The result is that some zoom levels look artificially blocky, while others look overly smeared, because the underlying data never changes.

The original workflow and this repo’s spec point in a different direction:
- zoom should still support continuous pan/zoom for inspection
- finer detail should appear when the computation zoom changes
- sub-dotter behavior should spawn a new, higher-resolution view for the selected region

Practical improvements:
- Add a recompute threshold after wheel zoom settles, so a user gets a high-quality refreshed plot instead of only an enlarged texture.
- Keep a small cache of pixelmaps at nearby zoom levels, or build a multi-resolution pyramid, so zooming does not always start from the coarsest raster.
- Snap the computation zoom to useful tiers for the data density, while leaving viewport zoom continuous for movement.
- Recompute overlays and axis labels from sequence coordinates, not from the current texture scale, so annotations stay crisp when the plot resolution changes.
- If full recomputation is too expensive during interaction, debounce the expensive work and show a fast preview first, then replace it with the refined view on scroll end.

## Overall Assessment

The project is beyond scaffolding and the algorithmic core is solid. The next risks are product completeness and scaling correctness, not basic correctness. If the goal is to reach the plan/spec faithfully, the highest-value work is to harden memory behavior, add the remaining export/annotation features, and make the docs reflect the current build and GUI status.
