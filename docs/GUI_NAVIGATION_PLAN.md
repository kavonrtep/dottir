# dottir — GUI navigation plan: cache, zoom out, exact plot on demand

Status: **implemented** (2026-11-20), interactive desktop testing pending. GUI only; no change to algorithms,
pixelmaps or file formats, so no ADR. Estimates are person-days.

Target release: **v0.8.0**.

## Progress

| ID | Item | Status |
|----|------|--------|
| N1 | Cancellable exact engine, Stop button | ☑ |
| N2 | Exact plot of a view on demand, with cost estimate | ☑ |
| N3 | Zoom out ×2 | ☑ |
| N4 | Cache and history rework | ☑ (one large session, see N4) |
| N5 | Docs, changelog, release | ☑ docs, changelog; release pending |

Order: N1 is a prerequisite of N2. Implemented N1, N4, N2, N3: N2's cache
entries need the engine in the N4 key.

## Problem

* **Recomputation.** View history and the pixelmap cache hold 5 entries
  each (`HISTORY_MAX`, `app.rs`). `invalidate_caches()` clears both on any
  compute-relevant settings change, so toggling a setting back recomputes
  and Back stops working. A cached large-mode entry keeps the whole
  `LargePlot` (8 score + 8 coverage channels, ~16 B/pixel: 55 MB at
  2000×1700, 256 MB at 4000²), so an entry count is the wrong bound.
* **No way out.** Back and Fit exist; there is no step out from the
  current view.
* **No detailed view on demand.** Exact scoring in large mode needs
  Settings → Engine → Exact (global, no estimate). The exact engine cannot
  be cancelled; a running job only stops being wanted when another
  request is dispatched, and it still runs to completion.

## Non-goals

* Deriving a view from another cached zoom (downsampling). Exact-engine
  pixels depend on the zoom block through the anti-diagonal rule
  (invariant 6), so a derived view is not the computed one.
* Disk persistence of pixelmaps.
* Mouse-wheel zoom in dot-plot mode: a zoom-out can cost minutes, so it is
  an explicit action only.
* A "sticky" exact mode. The existing global Engine setting covers that.
* Partial results after Stop.

## N1 — Cancellable exact engine, Stop (≈ 2 d)

Core: move `Control` from `large` to the crate root (re-exported from
`large` for compatibility) and add a cancellable entry point:

```rust
// dottir-core/src/plot.rs
pub fn compute_dotplot(query: &[u8], subject: &[u8], cfg: &PlotConfig)
    -> Result<DotPlot, DottirError>;           // = compute_dotplot_with(.., &Control::default())

pub fn compute_dotplot_with(query: &[u8], subject: &[u8], cfg: &PlotConfig,
    ctl: &Control<'_>) -> Result<DotPlot, DottirError>;   // Err(Cancelled) on cancel
```

* The parallel loop splits the subject into ~4 × threads chunks; at
  minutes per run that is several seconds per chunk. The cancel flag is
  therefore checked inside the sliding-window loop every 2²⁴ query cells
  (~10 ms; a row count would be too coarse for a 100 Mb query), and the
  cells done are added to an atomic counter that drives `ctl.progress`.
* BLASTX: the same, per frame pass.
* Pixels of a completed run are unchanged. Cancel only aborts.

GUI:

* `ComputeWorker::stop()` sets the cancel flag without dispatching. A
  **Stop** button appears in the compute overlay and the status bar while
  a job runs, for both engines.
* Large-engine index preparation becomes cancellable too. A cancelled
  preparation leaves no session; the next large request prepares again.
* After Stop the previous plot stays; the veil and badge go away. Status
  line: `Stopped.` A stopped request is not cached and adds no history.
* The worker reports `Stage::Exact` progress as a percentage, as it does
  for the large engine.

## N2 — Exact plot of a view on demand (≈ 3 d)

One-off, per view. Entry points:

* **Compute exact** button in the large-mode side panel and the status
  bar, for the current slice.
* **Exact plot of selection** in the right-drag selection menu.

The request carries an engine override; Settings → Engine is not changed.
The result is cached under engine `exact` (N4 key), so it sits next to the
large-mode view of the same slice, and Back returns to that. Zooming into
an exact view stays exact (rect-zoom inherits the override); Fit, zoom out
and Back drop it.

Cost estimate, shown before dispatch:

```
cells   = |q| × |s| × strands          (strands = 2 for Both)
          ÷ 2 for a triangle self-comparison
time    = cells / rate
memory  = w × h × (1 + separate reverse pixelmap)          pixelmaps
        + 4 B × alphabet × |q|                             score vector
        + per-thread window buffers (2 × 4 B × |q| × threads)
```

* `rate` (cells/s) is measured: each completed exact compute above 1 s
  updates a running estimate (exponential mean). Until the first
  measurement, a conservative default of 1 × 10⁹ cells/s is used and the
  dialog says "estimate, not yet measured on this machine".
* Under 5 s estimated: dispatch without asking.
* Otherwise a modal dialog:

  ```
  Exact plot of 1.20 Mb × 1.20 Mb (both strands)
  ≈ 2.9 × 10¹² cells · ≈ 4 min · ≈ 120 MB
  [Compute]  [Cancel]
  ```

* Memory above the exact engine's memory cap (`Settings::memory_limit_bytes`,
  the limit the core enforces): the dialog says so and offers only Cancel.
* The running job shows progress and Stop (N1).

## N3 — Zoom out ×2 (≈ 1 d)

* Centre: the crosshair if it lies in the visible area, otherwise the
  centre of the visible area. Click a point, then zoom out around it.
* New slice: the visible extent × 2 on each axis around the centre,
  shifted (not shrunk) to stay inside the sequences; compute zoom × 2, so
  the pixelmap keeps its size. If the slice covers both sequences, the
  action is Fit.
* Pushes history; goes through the cache (N4) like any view. Under engine
  Auto the larger area may switch to the large engine; that is intended.
* Entry points: toolbar button next to Back / Fit and the `-` key
  (unbound in dot-plot mode). Disabled at the fit view.
* Under engine Auto the exact engine only runs below the switch threshold
  (about a minute at most), so no dialog. With the global Engine = Exact,
  any exact dispatch estimated above 30 s (zoom out, Back, Fit) shows the
  N2 dialog first.

## N4 — Cache and history rework (≈ 3 d)

```rust
/// Hash of every setting that changes computed pixels: engine actually
/// used, mode, matrix, window, pixel factor (resolved), strand,
/// self-comparison, triangle, large seed / depth / class / sampling
/// parameters, repeat-annotation source key. Display settings (classes
/// shown, shading, greyramp, coverage ramp) are excluded.
struct ComputeKey(u64);

struct CacheEntry {
    key: ComputeKey,
    slice: ViewSlice,
    zoom: u32,
    plot: DotPlot,
    large: Option<LargePlot>,
    bytes: usize,          // plot + all LargePlot channels, coverage, chains
}

struct PixelmapCache {
    entries: VecDeque<CacheEntry>,   // LRU: touched entries move to the back
    budget: usize,                   // default 1 GiB, Settings → Performance
}

struct ViewSnapshot {                // history entry, unchanged fields
    slice: ViewSlice,
    compute_zoom: u32,
    view_offset: Vec2,
    crosshair: Option<(u32, u32)>,
    engine_override: Option<Engine>, // N2
}
```

* Cache bound: total `bytes` ≤ budget; least recently used entries are
  evicted first. A single entry above the budget is shown but not cached.
* A settings change no longer clears the cache; entries with another
  `ComputeKey` stop matching and age out. Setting a value back hits the
  cache.
* History: up to 50 snapshots, kept across settings changes. Back restores
  slice, zoom and offset under the **current** settings: from the cache if
  present, otherwise recomputed.
* History and cache are cleared only when sequences are (re)loaded.
* `locked_pixel_fac` / `locked_window` stay per settings state; they are
  part of `ComputeKey` through the resolved values.
* Large-engine sessions: the worker keeps **one** prepared session, as
  before. *Deviation from the first draft (two sessions):*
  `LargeSession::prepare_reusing` consumes the previous session to reuse
  its index after a window change or a record sort; keeping it as a second
  cached session would lose that reuse or double the index memory.
* Settings (next to the memory cap) shows the budget and current use
  (`Cache: 412 MB of 1 GB, 7 views`) and a **Clear** button.

## Failure modes

| Case | Behaviour |
|------|-----------|
| Stop during exact or large compute | previous plot stays; nothing cached; `Stopped.` |
| Stop during index preparation | no session kept; next large request re-prepares |
| Exact estimate above memory limit | dialog offers Cancel only |
| Rate not yet measured | default rate, dialog marks it as unmeasured |
| Entry larger than cache budget | displayed, not cached; logged at `info` |
| Back after settings change, not cached | recomputed under current settings |
| Zoom out at the fit view | button disabled, key ignored |
| New request while exact job runs | job cancelled (now effective for exact too) |

Logged at `info`: dispatch with engine and estimate, cancel, cache
hit/miss/evict with sizes, measured exact rate.

## Tests

* Core: `compute_dotplot_with` with a never-set cancel flag gives pixels
  identical to `compute_dotplot` (forward, reverse, both, self, BLASTX).
  A flag set before the call returns `Err(Cancelled)`; a flag set from a
  progress callback mid-run returns `Err(Cancelled)` before completion.
* Golden pixelmaps unchanged (`tests/golden`).
* GUI unit tests (no window):
  * cost estimate: cells for both strands and triangle; memory formula;
    `< 5 s` skips the dialog.
  * zoom out: centre on crosshair vs view centre; clamping shifts the
    slice; full coverage → Fit; zoom doubles.
  * cache: byte-budget eviction order; settings A → B → A hits the cache;
    an oversize entry is not inserted; `ComputeKey` ignores display
    settings.
  * history survives a settings change; Back restores slice and zoom.
* Worker: Stop on an exact job returns no result and the worker accepts
  the next request (60 s deadline, like the large worker test).

## N5 — Docs and release (≈ 1 d)

Book GUI chapter (navigation, Stop, exact on demand, cache setting),
changelog, release v0.8.0.

## Decisions

* Cache budget 1 GiB by default, adjustable.
* History survives settings changes.
* Exact on demand is one-off per view.
* No wheel zoom in dot-plot mode; zoom out is an explicit action.
