//! GFF3 carrying its own sequences after `##FASTA` as a positional
//! sequence input.
//!
//! Corpus: `tests/corpora/gff3_with_sequence/` (regenerate with
//! `scripts/make_gff3_with_sequence_fixtures.py`). The `*.with_seq.*`
//! fixtures are byte-preserving merges of the `annotation_overlay/`
//! pairs, so the central invariant here is equivalence: merged file
//! ≡ separate FASTA + GFF3.

use std::path::{Path, PathBuf};

use dottir_io::{load_sequence_input, AnnotSet, InputError, InputKind, Sequence};

fn corpus(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/corpora")
        .join(rel)
}

/// Merged GFF3 and the separate FASTA + GFF3 pair must agree on every
/// residue, record span and feature.
fn assert_equivalent(merged: &str, fasta: &str, gff: &str) {
    let input = load_sequence_input(corpus(merged)).unwrap();
    assert_eq!(input.kind, InputKind::Gff3WithFasta, "{merged}");

    let separate_seq = Sequence::load(corpus(fasta)).unwrap();
    assert_eq!(input.sequence.seq, separate_seq.seq, "{merged}: residues");
    assert_eq!(
        input.sequence.records.len(),
        separate_seq.records.len(),
        "{merged}: record count"
    );
    for (a, b) in input.sequence.records.iter().zip(&separate_seq.records) {
        assert_eq!(a.id, b.id, "{merged}: record id");
        assert_eq!(a.range, b.range, "{merged}: record span");
    }

    let separate_annot = AnnotSet::load(corpus(gff)).unwrap();
    let embedded = input.annotations.expect("embedded features");
    assert_eq!(
        embedded.features, separate_annot.features,
        "{merged}: features"
    );
    assert!(!embedded.features.is_empty(), "{merged}: no features");
}

#[test]
fn tir_simple_merged_matches_separate_pair() {
    assert_equivalent(
        "gff3_with_sequence/tir_simple.with_seq.gff3",
        "annotation_overlay/tir_simple.fasta",
        "annotation_overlay/tir_simple.gff3",
    );
}

#[test]
fn tir_elements_merged_matches_separate_pair() {
    assert_equivalent(
        "gff3_with_sequence/tir_elements.with_seq.gff3",
        "annotation_overlay/tir_elements.fasta",
        "annotation_overlay/tir_elements.gff3",
    );
}

#[test]
fn gzipped_merged_matches_separate_pair() {
    assert_equivalent(
        "gff3_with_sequence/ltr_angela.with_seq.gff3.gz",
        "annotation_overlay/ltr_angela.fasta",
        "annotation_overlay/ltr_angela.gff3",
    );
}

/// Every feature must name a record present in the embedded FASTA —
/// otherwise the GUI would silently drop it at bind time.
#[test]
fn embedded_features_reference_embedded_records() {
    for fixture in [
        "gff3_with_sequence/tir_simple.with_seq.gff3",
        "gff3_with_sequence/tir_elements.with_seq.gff3",
        "gff3_with_sequence/ltr_angela.with_seq.gff3.gz",
        "gff3_with_sequence/pair_a.gff3",
        "gff3_with_sequence/pair_b.gff3",
    ] {
        let input = load_sequence_input(corpus(fixture)).unwrap();
        let annot = input.annotations.expect("features");
        for f in &annot.features {
            let span = input
                .sequence
                .records
                .iter()
                .find(|r| r.id == f.record)
                .unwrap_or_else(|| panic!("{fixture}: feature on unknown record {}", f.record));
            assert!(
                f.range.end <= span.len(),
                "{fixture}: feature {:?} runs past the end of {} ({} > {})",
                f.range,
                f.record,
                f.range.end,
                span.len()
            );
        }
    }
}

#[test]
fn synthetic_pair_has_the_documented_layout() {
    let a = load_sequence_input(corpus("gff3_with_sequence/pair_a.gff3")).unwrap();
    assert_eq!(a.sequence.records.len(), 2);
    for rec in &a.sequence.records {
        assert_eq!(rec.len(), 600);
    }
    assert_eq!(a.sequence.len(), 1200);

    let annot = a.annotations.expect("features");
    assert_eq!(annot.features.len(), 4);
    let elements: Vec<_> = annot
        .features
        .iter()
        .filter(|f| f.attrs.get("Name").map(String::as_str) == Some("SharedElement"))
        .collect();
    assert_eq!(elements.len(), 2);
    // 1-based inclusive 151..450 → 0-based half-open 150..450.
    for e in elements {
        assert_eq!(e.range, 150..450);
        assert_eq!(e.feature_type, "repeat_region");
    }

    // The shared element really is shared: pair_b's copy is the same
    // length and ~8 % diverged, which is what makes the fixture useful
    // for an end-to-end dotplot.
    let b = load_sequence_input(corpus("gff3_with_sequence/pair_b.gff3")).unwrap();
    let ea = &a.sequence.seq[150..450];
    let eb = &b.sequence.seq[150..450];
    let diffs = ea.iter().zip(eb).filter(|(x, y)| x != y).count();
    assert!(
        (10..=70).contains(&diffs),
        "expected a diverged-but-alignable element, got {diffs}/300 differences"
    );
}

#[test]
fn gff3_without_embedded_fasta_is_rejected_as_a_sequence_input() {
    let err = load_sequence_input(corpus("gff3_with_sequence/no_sequence.gff3")).unwrap_err();
    assert!(
        matches!(err, InputError::NoEmbeddedFasta(_)),
        "expected NoEmbeddedFasta, got {err:?}"
    );
    // The message has to point at the fix, not just the failure.
    let msg = err.to_string();
    assert!(msg.contains("##FASTA"), "{msg}");
    assert!(msg.contains("--gff-query"), "{msg}");

    // The same file is still perfectly loadable as annotations — the
    // separate-file workflow is untouched.
    let annot = AnnotSet::load(corpus("gff3_with_sequence/no_sequence.gff3")).unwrap();
    assert_eq!(annot.features.len(), 1);
}

/// A features-only GFF3 passed to `AnnotSet::load` and the same file
/// with sequences appended must yield identical features — i.e. the
/// `##FASTA` split never leaks sequence lines into the feature parser.
#[test]
fn fasta_section_never_reaches_the_feature_parser() {
    let merged = AnnotSet::load(corpus("gff3_with_sequence/tir_simple.with_seq.gff3")).unwrap();
    let plain = AnnotSet::load(corpus("annotation_overlay/tir_simple.gff3")).unwrap();
    assert_eq!(merged.features, plain.features);
}
