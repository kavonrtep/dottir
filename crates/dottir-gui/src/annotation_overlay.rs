//! GUI-side annotation overlay state (ADR 0005).
//!
//! Loaded GFF3/BED files ([`dottir_io::AnnotSet`]) carry *record-local*
//! coordinates; the dotplot works in concatenated-buffer coordinates. An
//! [`AxisAnnot`] binds a parsed set to one axis (query = X, subject = Y) by
//! matching each feature's record id against the sequence's
//! [`dottir_io::RecordSpan`]s and translating the local interval into buffer
//! space. It also owns the per-axis coloring state (which attribute to color
//! by, the value→color palette, and per-value visibility).
//!
//! The pure [`band_screen_span`] helper maps a buffer-coordinate interval to a
//! clipped screen-space pixel span; it is the one piece of coordinate math the
//! band renderer relies on, and is unit-tested here.

use std::collections::{BTreeMap, BTreeSet};

use dottir_io::{AnnotSet, AnnotSource, Sequence, TYPE_KEY};
use egui::Color32;

/// Categorical palette for distinct color-by values. A compact, high-contrast
/// set (Tableau-10 inspired); cycled when there are more values than colors.
pub const PALETTE: [Color32; 10] = [
    Color32::from_rgb(0x4e, 0x79, 0xa7),
    Color32::from_rgb(0xf2, 0x8e, 0x2b),
    Color32::from_rgb(0xe1, 0x57, 0x59),
    Color32::from_rgb(0x76, 0xb7, 0xb2),
    Color32::from_rgb(0x59, 0xa1, 0x4f),
    Color32::from_rgb(0xed, 0xc9, 0x48),
    Color32::from_rgb(0xb0, 0x7a, 0xa1),
    Color32::from_rgb(0xff, 0x9d, 0xa7),
    Color32::from_rgb(0x9c, 0x75, 0x5f),
    Color32::from_rgb(0xba, 0xb0, 0xac),
];

/// Color for features missing the active color-by attribute.
pub const NONE_COLOR: Color32 = Color32::from_rgb(0x96, 0x96, 0x96);

/// Value bucket for features lacking the active color-by attribute.
pub const NONE_VALUE: &str = "(none)";

/// A feature bound to one axis, in concatenated-buffer coordinates.
#[derive(Debug, Clone)]
pub struct BoundFeature {
    /// `[start, end)` in the axis sequence's concatenated buffer.
    pub range: std::ops::Range<usize>,
    /// Index into [`AxisAnnot::set`]'s `features`, used to re-derive the
    /// color-by value when the active key changes.
    pub src: usize,
}

/// One axis's loaded annotation set plus its binding and coloring state.
#[derive(Debug, Clone)]
pub struct AxisAnnot {
    pub set: AnnotSet,
    /// Features whose record id matched a record on this axis, in buffer
    /// coords. In file order.
    pub bound: Vec<BoundFeature>,
    /// Count of features dropped because their record id matched nothing on
    /// this axis (surfaced as a warning in the panel).
    pub skipped: usize,
    /// Attribute key to color by, or [`TYPE_KEY`]. Empty string for BED
    /// (single color, see [`Self::bed_color`]).
    pub color_by: String,
    /// value → color. User overrides are preserved across `color_by` changes
    /// when the value still exists.
    pub palette: BTreeMap<String, Color32>,
    /// Values toggled off in the legend (not drawn).
    pub hidden: BTreeSet<String>,
    /// Single color for BED sources.
    pub bed_color: Color32,
}

impl AxisAnnot {
    /// Bind a parsed set to `seq`, choosing a sensible default color-by key
    /// and building the initial palette.
    pub fn new(set: AnnotSet, seq: &Sequence) -> Self {
        // record id -> buffer start. Last span wins on duplicate ids (rare);
        // duplicate FASTA ids are already ambiguous upstream.
        let mut starts: BTreeMap<&str, &std::ops::Range<usize>> = BTreeMap::new();
        for rec in &seq.records {
            starts.insert(rec.id.as_str(), &rec.range);
        }

        let mut bound = Vec::new();
        let mut skipped = 0usize;
        for (i, f) in set.features.iter().enumerate() {
            let Some(span) = starts.get(f.record.as_str()) else {
                skipped += 1;
                continue;
            };
            let rec_len = span.end - span.start;
            // Clamp to the record; drop features starting past its end.
            let lo = f.range.start.min(rec_len);
            let hi = f.range.end.min(rec_len);
            if hi <= lo {
                skipped += 1;
                continue;
            }
            bound.push(BoundFeature {
                range: (span.start + lo)..(span.start + hi),
                src: i,
            });
        }

        let color_by = default_color_by(&set);
        let mut ax = AxisAnnot {
            set,
            bound,
            skipped,
            color_by,
            palette: BTreeMap::new(),
            hidden: BTreeSet::new(),
            bed_color: PALETTE[0],
        };
        ax.rebuild_palette();
        ax
    }

    pub fn is_bed(&self) -> bool {
        self.set.source == AnnotSource::Bed
    }

    /// The color-by value of a bound feature. When the active color-by
    /// attribute (e.g. `Name`) is absent for a feature, fall back to its GFF3
    /// type (column 3) so it still gets a meaningful bucket/color rather than
    /// lumping every nameless feature into `(none)`. `(none)` remains only for
    /// features with neither the attribute nor a type. BED features all share
    /// one bucket (BED uses [`Self::bed_color`] directly).
    pub fn value_of(&self, bf: &BoundFeature) -> String {
        let f = &self.set.features[bf.src];
        f.color_value(&self.color_by)
            .or_else(|| f.color_value(TYPE_KEY))
            .unwrap_or(NONE_VALUE)
            .to_string()
    }

    /// Resolve the draw color for a color-by value.
    pub fn color_for(&self, value: &str) -> Color32 {
        if self.is_bed() {
            return self.bed_color;
        }
        *self.palette.get(value).unwrap_or(&NONE_COLOR)
    }

    /// Visible bound features whose buffer range covers `pos`, in file order.
    /// Features whose color-by value is toggled off in the legend are skipped,
    /// matching what the band renderer draws. Nested/overlapping features all
    /// appear, so callers can surface more than one.
    pub fn features_at(&self, pos: usize) -> Vec<&BoundFeature> {
        self.bound
            .iter()
            .filter(|bf| bf.range.contains(&pos) && !self.hidden.contains(&self.value_of(bf)))
            .collect()
    }

    /// Distinct color-by values with feature counts, sorted by descending
    /// count then value (for a stable legend order).
    pub fn value_counts(&self) -> Vec<(String, usize)> {
        let mut counts: BTreeMap<String, usize> = BTreeMap::new();
        for bf in &self.bound {
            *counts.entry(self.value_of(bf)).or_default() += 1;
        }
        let mut v: Vec<(String, usize)> = counts.into_iter().collect();
        v.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        v
    }

    /// Switch the active color-by key and rebuild the palette, preserving any
    /// user color overrides for values that still exist.
    pub fn set_color_by(&mut self, key: String) {
        self.color_by = key;
        self.rebuild_palette();
    }

    /// Assign palette colors to the current distinct values. Existing entries
    /// are kept (user overrides survive); new values get the next palette slot
    /// in sorted order for determinism.
    pub fn rebuild_palette(&mut self) {
        let values: BTreeSet<String> = self.bound.iter().map(|bf| self.value_of(bf)).collect();
        // Drop palette/hidden entries for values that no longer exist.
        self.palette.retain(|k, _| values.contains(k));
        self.hidden.retain(|k| values.contains(k));
        let mut next = self.palette.len();
        for v in &values {
            if v == NONE_VALUE {
                self.palette.entry(v.clone()).or_insert(NONE_COLOR);
                continue;
            }
            if !self.palette.contains_key(v) {
                self.palette
                    .insert(v.clone(), PALETTE[next % PALETTE.len()]);
                next += 1;
            }
        }
    }
}

/// Default color-by key for a freshly loaded set: `Name` when present, else
/// the synthetic `(type)`, else the first attribute key, else empty (BED).
fn default_color_by(set: &AnnotSet) -> String {
    if set.source == AnnotSource::Bed {
        return String::new();
    }
    let keys = set.attribute_keys();
    if keys.contains("Name") {
        "Name".to_string()
    } else if keys.contains(TYPE_KEY) {
        TYPE_KEY.to_string()
    } else {
        keys.into_iter().next().unwrap_or_default()
    }
}

/// Map a buffer-coordinate interval `[lo, hi)` on one axis to a screen-space
/// pixel span `[s0, s1]`, clipped to `[clip_lo, clip_hi]`.
///
/// The transform mirrors the pixelmap/breakline/ridge math in `draw_canvas`:
/// `pixel = (coord - off) / zoom`, then
/// `screen = axis_origin + pixel / ppp - draw_offset`, where `off` is the
/// slice origin on this axis. Returns `None` when the band lies entirely
/// outside the visible axis range.
// Nine scalar args, but they're the irreducible inputs to one coordinate
// transform; bundling them into a struct would add ceremony without clarity.
#[allow(clippy::too_many_arguments)]
pub fn band_screen_span(
    lo: usize,
    hi: usize,
    off: usize,
    zoom: usize,
    axis_origin: f32,
    draw_offset: f32,
    ppp: f32,
    clip_lo: f32,
    clip_hi: f32,
) -> Option<(f32, f32)> {
    let zoom = zoom.max(1) as f32;
    let to_screen = |c: usize| {
        let pixel = (c as f32 - off as f32) / zoom;
        axis_origin + pixel / ppp - draw_offset
    };
    let mut s0 = to_screen(lo);
    let mut s1 = to_screen(hi);
    if s1 < s0 {
        std::mem::swap(&mut s0, &mut s1);
    }
    let c0 = s0.max(clip_lo);
    let c1 = s1.min(clip_hi);
    (c1 > c0).then_some((c0, c1))
}

#[cfg(test)]
mod tests {
    use super::*;
    use dottir_io::{AnnotSource, Feature, RecordSpan};
    use std::collections::BTreeMap as Map;

    fn seq_with_records(recs: &[(&str, usize)]) -> Sequence {
        // Build a Sequence directly with the desired record spans.
        let mut seq = Vec::new();
        let mut records = Vec::new();
        let mut off = 0;
        for (id, len) in recs {
            seq.extend(std::iter::repeat_n(b'A', *len));
            records.push(RecordSpan {
                id: id.to_string(),
                description: None,
                range: off..off + *len,
            });
            off += *len;
        }
        Sequence {
            seq,
            records,
            source_path: None,
        }
    }

    fn feat(record: &str, start: usize, end: usize, name: &str) -> Feature {
        let mut attrs = Map::new();
        if !name.is_empty() {
            attrs.insert("Name".to_string(), name.to_string());
        }
        Feature {
            record: record.to_string(),
            range: start..end,
            strand: dottir_io::Strand::Unknown,
            feature_type: "region".to_string(),
            attrs,
        }
    }

    #[test]
    fn binding_maps_local_to_buffer_and_counts_skips() {
        // Two records: chr1 [0,100), chr2 [100,150).
        let seq = seq_with_records(&[("chr1", 100), ("chr2", 50)]);
        let set = AnnotSet {
            source_path: None,
            source: AnnotSource::Gff3,
            features: vec![
                feat("chr1", 10, 20, "ALR"),  // -> buffer [10,20)
                feat("chr2", 5, 15, "HSAT"),  // -> buffer [105,115)
                feat("chrX", 0, 10, "ghost"), // unknown record -> skipped
                feat("chr2", 60, 70, "oob"),  // past record end -> skipped
            ],
        };
        let ax = AxisAnnot::new(set, &seq);
        assert_eq!(ax.bound.len(), 2);
        assert_eq!(ax.bound[0].range, 10..20);
        assert_eq!(ax.bound[1].range, 105..115);
        assert_eq!(ax.skipped, 2);
        // Default color-by is Name (present).
        assert_eq!(ax.color_by, "Name");
        // Each distinct Name got a palette color.
        assert_eq!(ax.color_for("ALR"), PALETTE[0]);
        assert_ne!(ax.color_for("ALR"), ax.color_for("HSAT"));
    }

    #[test]
    fn binding_clamps_feature_overrunning_record_end() {
        let seq = seq_with_records(&[("chr1", 100)]);
        let set = AnnotSet {
            source_path: None,
            source: AnnotSource::Gff3,
            features: vec![feat("chr1", 90, 120, "tail")], // clamp end to 100
        };
        let ax = AxisAnnot::new(set, &seq);
        assert_eq!(ax.bound.len(), 1);
        assert_eq!(ax.bound[0].range, 90..100);
    }

    #[test]
    fn bed_uses_single_color_and_empty_color_by() {
        let seq = seq_with_records(&[("chr1", 100)]);
        let mut f = feat("chr1", 0, 10, "");
        f.feature_type = String::new(); // BED features carry no type
        let set = AnnotSet {
            source_path: None,
            source: AnnotSource::Bed,
            features: vec![f],
        };
        let ax = AxisAnnot::new(set, &seq);
        assert!(ax.is_bed());
        assert_eq!(ax.color_by, "");
        assert_eq!(ax.color_for("anything"), ax.bed_color);
    }

    #[test]
    fn value_of_falls_back_to_type_when_color_by_attr_absent() {
        // chr1 with a named and an unnamed feature; color-by defaults to Name.
        let seq = seq_with_records(&[("chr1", 100)]);
        let named = feat("chr1", 0, 10, "ALR"); // feature_type "region", Name=ALR
        let mut unnamed = feat("chr1", 20, 30, ""); // no Name
        unnamed.feature_type = "transposable_element".to_string();
        let set = AnnotSet {
            source_path: None,
            source: AnnotSource::Gff3,
            features: vec![named, unnamed],
        };
        let ax = AxisAnnot::new(set, &seq);
        assert_eq!(ax.color_by, "Name");
        // Named feature uses its Name; unnamed falls back to its GFF3 type,
        // not "(none)", and gets its own palette color.
        assert_eq!(ax.value_of(&ax.bound[0]), "ALR");
        assert_eq!(ax.value_of(&ax.bound[1]), "transposable_element");
        assert_ne!(
            ax.color_for("transposable_element"),
            ax.color_for(NONE_VALUE)
        );
    }

    #[test]
    fn features_at_returns_all_overlapping_and_respects_hidden() {
        // One record chr1 [0,100). Three features, two of them nested around
        // position 50.
        let seq = seq_with_records(&[("chr1", 100)]);
        let set = AnnotSet {
            source_path: None,
            source: AnnotSource::Gff3,
            features: vec![
                feat("chr1", 0, 80, "outer"),  // covers 50
                feat("chr1", 40, 60, "inner"), // covers 50
                feat("chr1", 70, 90, "other"), // does not cover 50
            ],
        };
        let mut ax = AxisAnnot::new(set, &seq);

        // Position inside both nested features → both returned, in file order.
        let hits: Vec<String> = ax
            .features_at(50)
            .iter()
            .map(|bf| ax.value_of(bf))
            .collect();
        assert_eq!(hits, vec!["outer".to_string(), "inner".to_string()]);

        // A gap with no coverage.
        assert!(ax.features_at(95).is_empty());

        // End is exclusive: feature [70,90) does not cover 90.
        assert!(ax
            .features_at(90)
            .iter()
            .all(|bf| ax.value_of(bf) != "other"));

        // Hiding a value drops it from the lookup.
        ax.hidden.insert("inner".to_string());
        let hits: Vec<String> = ax
            .features_at(50)
            .iter()
            .map(|bf| ax.value_of(bf))
            .collect();
        assert_eq!(hits, vec!["outer".to_string()]);
    }

    #[test]
    fn band_span_basic_transform() {
        // off=100, zoom=2, origin=50, draw_offset=10, ppp=1, clip [50,250].
        // feature [120,160): pixels [10,30] -> screen [50,70].
        let s = band_screen_span(120, 160, 100, 2, 50.0, 10.0, 1.0, 50.0, 250.0);
        assert_eq!(s, Some((50.0, 70.0)));
    }

    #[test]
    fn band_span_clips_left_edge() {
        // feature starts before the slice/view: [100,140) -> screen [40,60],
        // clipped to clip_lo=50 -> [50,60].
        let s = band_screen_span(100, 140, 100, 2, 50.0, 10.0, 1.0, 50.0, 250.0);
        assert_eq!(s, Some((50.0, 60.0)));
    }

    #[test]
    fn band_span_fully_offscreen_is_none() {
        // [100,110) -> screen [40,45], entirely left of clip_lo=50.
        let s = band_screen_span(100, 110, 100, 2, 50.0, 10.0, 1.0, 50.0, 250.0);
        assert_eq!(s, None);
    }

    #[test]
    fn band_span_respects_ppp() {
        // ppp=2 halves the on-screen extent. [0,100) at zoom 1 -> 100 px ->
        // 50 logical points.
        let s = band_screen_span(0, 100, 0, 1, 0.0, 0.0, 2.0, 0.0, 1000.0);
        assert_eq!(s, Some((0.0, 50.0)));
    }
}
