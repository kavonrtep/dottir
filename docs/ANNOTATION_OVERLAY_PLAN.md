# Implementation Plan — GFF3/BED annotation overlay

Companion to **ADR 0005**. Feature: overlay annotation intervals on the dotplot
as full-span bands, for self-comparison (one file → both axes) and pairwise
comparison (one file per axis). GFF3 colored by attribute value; BED single
color. Global alpha slider. Loading via File menu + CLI, per axis.

Naming follows the codebase: **query** = X axis (`q_*`), **subject** = Y axis
(`s_*`); existing overlay precedent is the ridge renderer at `app.rs:2577`.

## Data model (`dottir-io`)

New module `dottir-io/src/annotation.rs`:

```rust
pub struct Feature {
    pub record: String,              // GFF3 col1 seqid / BED chrom
    pub range: Range<usize>,         // LOCAL, 0-based half-open within the record
    pub strand: Option<Strand>,      // +/-/. (parsed; unused by bands in v1)
    pub attrs: BTreeMap<String, String>, // GFF3 col9; BED: {"name": col4} if present
    pub feature_type: String,        // GFF3 col3; BED: "" 
}

pub enum AnnotSource { Gff3, Bed }

pub struct AnnotSet {
    pub source_path: PathBuf,
    pub source: AnnotSource,
    pub features: Vec<Feature>,      // local coords; mapped to buffer offsets at bind time
}
```

- `AnnotSet::load(path) -> Result<AnnotSet, AnnotError>` dispatches on extension
  (`.gff/.gff3` → noodles-gff; `.bed` → hand-rolled reader; transparent `.gz`
  via `flate2`, mirroring the FASTA loader).
- `attribute_keys(&self) -> BTreeSet<String>` — union of col9 keys (drives the
  "Color by" dropdown). `(type)` is a synthetic key handled by the GUI, not stored.
- Typed errors via `thiserror` (`AnnotError`), `anyhow` only inside binaries.
- **Coordinate care**: GFF3 is 1-based inclusive `[start, end]`; convert to
  0-based half-open `[start-1, end)` on parse. BED is already 0-based half-open.
  Both stored as half-open to match `Sequence`/`RecordSpan` byte indexing.

## Axis binding (`dottir-gui`, app state)

A loaded `AnnotSet` has *local* per-record coordinates; the dotplot works in the
*concatenated buffer* coordinate space. Bind at load/assign time:

```rust
struct BoundFeature { range: Range<usize>, color_key: String } // buffer coords
struct AxisAnnot {
    set: AnnotSet,
    bound: Vec<BoundFeature>,        // features whose seqid matched a record on this axis
    skipped: usize,                  // seqid-not-found count (surfaced as a warning)
    color_by: String,                // attr key or "(type)" (GFF3); ignored for BED
    palette: BTreeMap<String, Color32>, // value -> color, with user overrides
    hidden: BTreeSet<String>,        // legend visibility toggles
    bed_color: Color32,              // used when source == Bed
}
```

New `DottirApp` fields: `query_annot: Option<AxisAnnot>`,
`subject_annot: Option<AxisAnnot>`, `annot_alpha: f32` (in `Settings`),
`annot_show: bool`. Self-mode: load once, bind against the single sequence for
both axes (clone the `AxisAnnot` or share via `Rc`).

Binding maps each `Feature.range` (local) → buffer coords using
`Sequence::record_at` / `RecordSpan { id, range }` (sequence.rs:45,125). Features
whose `record` matches no `RecordSpan.id` increment `skipped`.

## Coordinate transform (the testable core)

Reuse the established chain (mirrors crosshair at `app.rs:2615` and ridges at
`app.rs:2582`):

```
buffer_coord  --(- slice.start)/zoom-->  pixelmap_pixel
pixelmap_pixel --/ppp + plot_area.left - draw_offset.x--> screen_x   (query/X)
                                          plot_area.top  - draw_offset.y--> screen_y (subject/Y)
```

A query feature `[lo,hi)` → screen rect `x ∈ [f(lo), f(hi)]`,
`y ∈ [plot_area.top, plot_area.bottom]`; subject feature transposed. Clip to the
plot rect before drawing. **Unit test** `feature_to_pixel_band` pins this with a
known slice+zoom (the rule CLAUDE.md flags as "most likely to silently drift").

## Rendering (`dottir-gui`, canvas)

In the canvas draw block after the pixelmap image, before/after ridges:

```rust
if self.settings.annot_show {
    for ax in [&self.query_annot, &self.subject_annot].into_iter().flatten() {
        for bf in &ax.bound {
            if ax.hidden.contains(&bf.color_key) { continue; }
            let mut c = ax.color_for(&bf.color_key);     // palette or bed_color
            c = c.gamma_multiply(self.settings.annot_alpha); // single global alpha
            let rect = /* band rect, clipped to plot */;
            clip_painter.rect_filled(rect, 0.0, c);
        }
    }
}
```

- Only emit features intersecting the visible slice (`current_slice`, app.rs:327).
- Draw order: over the pixelmap so dots show through alpha; under the crosshair
  and axis labels so those stay legible.
- Overlap self-darkens via alpha compositing — no special intersection pass.

## Side-panel "Annotations" section (`dottir-gui`, settings panel)

New collapsing section modeled on the ridge controls (`app.rs:2025`):

```
▾ Annotations
  Query:   <name> [Load…][Clear]      (or "—")
  Subject: <name> [Load…][Clear]
  [✓] Show bands
  Color by: [ Name ▾ ]                 (GFF3 sources only)
  Alpha:    [──●────] 0.35
  ── Legend (per axis, sorted by count desc) ──
  [✓] ALR/Alpha  ■  (812)
  [✓] HSAT       ■  (240)
  [ ] LINE1      ■  (95)
  [✓] (none)     ■  (12)
  ⚠ 142 subject features skipped (seqid not found)
```

- `Color by` dropdown lists `attribute_keys() ∪ {"(type)"}`; changing it rebuilds
  `palette`/`color_key`. Hidden for BED-only sources.
- Legend swatches are `egui::color_picker`; checkbox toggles `hidden`.
- Self-mode shows a single "GFF3/BED" row instead of query/subject pair.

## Loading wiring

- **File menu** (`app.rs` menu bar): `Load GFF3/BED for query…` /
  `…for subject…` via `rfd` (already a dep, filter `gff,gff3,bed,gz`). Self-mode:
  one `Load annotations…`.
- **CLI** (`dottir-gui` clap args + `dottir-cli` if batch export ever needs it):
  `--gff-query`, `--bed-query`, `--gff-subject`, `--bed-subject`, `--gff`/`--bed`
  (self-mode). Parse → `AnnotSet::load` → bind → store. Mismatched/empty binding
  emits a `tracing::warn!` and the panel ⚠ line.

## Dependencies

- `noodles-gff` (spec §6.3 pin) — add to `dottir-io`.
- Categorical colors: small built-in palette in `dottir-gui` (deterministic, no
  dep). `colorgrad` is for continuous colormaps; not needed here.
- `flate2` — already present (gz). `rfd` — already present (dialogs).
- No `noodles-bed` (hand-rolled).

## Test plan

1. `dottir-io`: parse a small GFF3 (1-based→0-based conversion, attrs, gz) and a
   BED; assert `Feature` ranges and `attribute_keys()`.
2. `dottir-io`: seqid filtering — feature on an absent record is excluded and
   counted.
3. `dottir-gui`: `feature_to_pixel_band` transform unit test (known slice+zoom),
   including a feature partly off-screen (clipping).
4. `egui_kittest` smoke: load fixture, toggle Show bands, assert no panic and a
   non-empty legend (spec §4.5.8 style).

## Out of scope (v1)

- Layered multiple files per axis; PAF HSP overlay (separate ADR/feature).
- Strand-aware rendering (strand is parsed and stored, not yet drawn).
- Feature tooltips / click-to-inspect (nice follow-up; not required).

## Suggested commit sequence

1. `dottir-io`: `annotation` module + parsers + tests (no GUI).
2. `dottir-gui`: app state + binding + coordinate transform + test.
3. `dottir-gui`: band renderer.
4. `dottir-gui`: Annotations panel section + legend.
5. `dottir-gui`: File-menu + CLI wiring + self-mode.
6. Docs: CHANGELOG entry; mark spec §5 item done.
