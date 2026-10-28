//! Session save/load (C5).
//!
//! A *session* captures everything needed to reproduce the current GUI
//! state from a single TOML file: the loaded sequence file paths, the
//! plot settings, the greyramp, the view transform, and the crosshair.
//! It deliberately does NOT include the sequence bytes — re-loading
//! the sessions opens the same FASTA files from disk.
//!
//! TOML schema is stable across minor versions; new fields may appear
//! but existing field names won't change meaning.

use std::path::PathBuf;

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Session {
    /// Schema version. Bumped only on breaking changes; readers should
    /// refuse versions they don't recognise.
    pub version: u32,
    pub query: Option<PathBuf>,
    pub subject: Option<PathBuf>,
    pub plot: SessionPlot,
    pub greyramp: SessionGreyramp,
    pub view: SessionView,
    /// Large-sequence mode settings (ADR 0007). Absent in older sessions.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub large: Option<SessionLarge>,
    /// Genome comparison (ADR 0008). Absent in older sessions.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome: Option<SessionGenome>,
    /// Record sorts applied after loading, in order: "vertical" (subject
    /// records along the query) or "horizontal".
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub sorted_records: Vec<String>,
    /// HSP overlay: PAF path and mode (`line`, `score`, `greyscale`, `off`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub hsps: Option<SessionHsps>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SessionHsps {
    pub paf: PathBuf,
    pub mode: String,
}

/// A genome comparison and its view settings. The alignments are not
/// stored: re-opening runs the aligner again (from its cache) or reads
/// the PAF.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SessionGenome {
    pub inputs: crate::genome_view::GenomeInputs,
    /// Whether the genome view (rather than the dotplot) was shown.
    pub active: bool,
    pub min_aln_len: u64,
    pub min_identity: f64,
    pub primary_only: bool,
    pub remove_noise: bool,
    pub sort_query: bool,
    /// Reference records sorted along the query (exclusive with
    /// `sort_query`).
    #[serde(default)]
    pub sort_reference: bool,
    /// Vertical records flipped by hand, by name.
    pub flipped: Vec<String>,
    pub line_width: u32,
    pub drill_limit_mb: f64,
}

impl SessionGenome {
    pub fn from_view(v: &crate::genome_view::GenomeView, active: bool) -> Self {
        Self {
            inputs: v.inputs.clone(),
            active,
            min_aln_len: v.filter.min_len,
            min_identity: v.filter.min_identity,
            primary_only: v.filter.primary_only,
            remove_noise: v.filter.remove_noise,
            sort_query: v.sort == Some(dottir_core::genome::SortAxis::Vertical),
            sort_reference: v.sort == Some(dottir_core::genome::SortAxis::Horizontal),
            flipped: v
                .flips
                .iter()
                .zip(&v.data.y.records)
                .filter(|(&f, _)| f)
                .map(|(_, r)| r.name.clone())
                .collect(),
            line_width: v.line_width,
            drill_limit_mb: v.drill_limit_mb,
        }
    }

    /// Restore the view settings; unknown record names are ignored.
    pub fn apply(&self, v: &mut crate::genome_view::GenomeView) {
        v.filter.min_len = self.min_aln_len;
        v.filter.min_identity = self.min_identity.clamp(0.0, 1.0);
        v.filter.primary_only = self.primary_only;
        v.filter.remove_noise = self.remove_noise;
        use dottir_core::genome::SortAxis;
        v.sort = match (self.sort_query, self.sort_reference) {
            (true, _) => Some(SortAxis::Vertical),
            (false, true) => Some(SortAxis::Horizontal),
            _ => None,
        };
        v.line_width = self.line_width.clamp(1, 5);
        v.drill_limit_mb = self.drill_limit_mb.max(1.0);
        for name in &self.flipped {
            if let Some(i) = v.data.y.index.find(name) {
                v.flips[i] = true;
            }
        }
        v.rearrange();
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SessionLarge {
    /// "auto" / "exact" / "large".
    pub engine: String,
    pub threshold_cells: f64,
    pub seed_pattern: String,
    pub two_hit_dist: u32,
    pub max_occ: u32,
    /// "seed" / "kmer" / "none".
    pub depth: String,
    pub depth_k: u8,
    /// `[low, mid, high]` copy-number thresholds.
    pub class_thresholds: [u16; 3],
    pub softmask_as_repeat: bool,
    pub memory_limit_mib: u64,
    pub shown_classes: Vec<String>,
    pub colour_by_class: bool,
    pub show_chains: bool,
    pub chain_min_score: f32,
    // Added after the first large-mode sessions; absent = default.
    /// "coverage" / "score".
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub shade: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coverage_ramp_bp: Option<[f32; 2]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coverage_auto: Option<bool>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub axis_tracks: Option<bool>,
    /// Seed sampling rates; 0 = auto.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seed_sample: Option<u32>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub depth_sample: Option<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SessionPlot {
    /// Lowercase: "blastn" / "blastp" / "blastx".
    pub mode: String,
    pub matrix_name: String,
    /// `None` ↔ Karlin/Altschul auto.
    pub window_size: Option<u32>,
    pub zoom: u32,
    pub pixel_fac: u32,
    /// "forward" / "reverse" / "both".
    pub strand: String,
    pub self_comparison: bool,
    /// "both" / "upper" / "lower".
    pub triangle: String,
    pub memory_limit_mib: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SessionGreyramp {
    pub white: u8,
    pub black: u8,
    pub swap: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SessionView {
    /// Top-left of canvas in pixelmap coords.
    pub offset_x: f32,
    pub offset_y: f32,
    /// Display (viewport) zoom — distinct from `plot.zoom`, which is
    /// the *computation* zoom.
    pub display_zoom: f32,
    /// Crosshair residue coord `[q_seq, s_seq]` (0-indexed) since
    /// schema v2. Schema v1 stored *full-sequence pixelmap* coords
    /// here; the GUI migrates those on load (multiplies by the
    /// recorded compute zoom + adds half a block) so old sessions
    /// round-trip without visible drift.
    pub crosshair: Option<[u32; 2]>,
    pub light_theme: bool,
    /// H2: whether the alignment-view dock is shown. Default-true
    /// via `#[serde(default = …)]` so older session files
    /// (schema v1 pre-H2) round-trip cleanly.
    #[serde(default = "default_true")]
    pub align_dock_visible: bool,
    /// H2: window size (residues) of the alignment dock. Default 100
    /// via `#[serde(default = …)]` for pre-H2 sessions.
    #[serde(default = "default_align_window")]
    pub align_window_size: u32,
}

fn default_true() -> bool {
    true
}
fn default_align_window() -> u32 {
    100
}

pub const SESSION_VERSION: u32 = 2;

#[derive(Debug, thiserror::Error)]
pub enum SessionError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("TOML parse error: {0}")]
    Parse(#[from] toml::de::Error),
    #[error("TOML serialise error: {0}")]
    Serialize(#[from] toml::ser::Error),
    #[error("unsupported session version {0}; this dottir reads up to {1}")]
    UnsupportedVersion(u32, u32),
}

impl Session {
    pub fn save<P: AsRef<std::path::Path>>(&self, path: P) -> Result<(), SessionError> {
        let toml_text = toml::to_string_pretty(self)?;
        std::fs::write(path, toml_text)?;
        Ok(())
    }

    pub fn load<P: AsRef<std::path::Path>>(path: P) -> Result<Self, SessionError> {
        let text = std::fs::read_to_string(path)?;
        let s: Session = toml::from_str(&text)?;
        if s.version > SESSION_VERSION {
            return Err(SessionError::UnsupportedVersion(s.version, SESSION_VERSION));
        }
        Ok(s)
    }
}

/// String<->enum helpers; the GUI imports these via wildcard.
pub mod codec {
    use dottir_core::{BlastMode, Strand, Triangle};

    pub fn mode_to_str(m: BlastMode) -> &'static str {
        match m {
            BlastMode::Blastn => "blastn",
            BlastMode::Blastp => "blastp",
            BlastMode::Blastx => "blastx",
        }
    }
    pub fn mode_from_str(s: &str) -> Option<BlastMode> {
        match s.to_ascii_lowercase().as_str() {
            "blastn" => Some(BlastMode::Blastn),
            "blastp" => Some(BlastMode::Blastp),
            "blastx" => Some(BlastMode::Blastx),
            _ => None,
        }
    }

    pub fn strand_to_str(s: Strand) -> &'static str {
        match s {
            Strand::Forward => "forward",
            Strand::Reverse => "reverse",
            Strand::Both => "both",
        }
    }
    pub fn strand_from_str(s: &str) -> Option<Strand> {
        match s.to_ascii_lowercase().as_str() {
            "forward" => Some(Strand::Forward),
            "reverse" => Some(Strand::Reverse),
            "both" => Some(Strand::Both),
            _ => None,
        }
    }

    pub fn triangle_to_str(t: Triangle) -> &'static str {
        match t {
            Triangle::Both => "both",
            Triangle::Upper => "upper",
            Triangle::Lower => "lower",
        }
    }
    pub fn triangle_from_str(s: &str) -> Option<Triangle> {
        match s.to_ascii_lowercase().as_str() {
            "both" => Some(Triangle::Both),
            "upper" => Some(Triangle::Upper),
            "lower" => Some(Triangle::Lower),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> Session {
        Session {
            version: SESSION_VERSION,
            query: Some(PathBuf::from("/tmp/q.fa")),
            subject: Some(PathBuf::from("/tmp/s.fa")),
            plot: SessionPlot {
                mode: "blastn".into(),
                matrix_name: "DNA+5/-4".into(),
                window_size: Some(25),
                zoom: 1,
                pixel_fac: 50,
                strand: "both".into(),
                self_comparison: false,
                triangle: "both".into(),
                memory_limit_mib: 512,
            },
            greyramp: SessionGreyramp {
                white: 40,
                black: 100,
                swap: false,
            },
            view: SessionView {
                offset_x: 0.0,
                offset_y: 0.0,
                display_zoom: 1.0,
                crosshair: Some([10, 20]),
                light_theme: true,
                align_dock_visible: true,
                align_window_size: 100,
            },
            large: Some(SessionLarge {
                engine: "auto".into(),
                threshold_cells: 2e10,
                seed_pattern: "1110100110010101111".into(),
                two_hit_dist: 100,
                max_occ: 1000,
                depth: "seed".into(),
                depth_k: 21,
                class_thresholds: [2, 5, 51],
                softmask_as_repeat: true,
                memory_limit_mib: 6144,
                shown_classes: vec!["unique".into(), "low".into()],
                colour_by_class: true,
                show_chains: true,
                chain_min_score: 1000.0,
                shade: Some("coverage".into()),
                coverage_ramp_bp: Some([50.0, 3000.0]),
                coverage_auto: Some(false),
                axis_tracks: Some(true),
                seed_sample: Some(0),
                depth_sample: Some(4),
            }),
            sorted_records: vec!["vertical".into()],
            hsps: Some(SessionHsps {
                paf: "/tmp/a.paf".into(),
                mode: "score".into(),
            }),
            genome: Some(SessionGenome {
                inputs: crate::genome_view::GenomeInputs {
                    reference: PathBuf::from("/tmp/ref.fa.gz"),
                    query: None,
                    paf: Some(PathBuf::from("/tmp/a.paf")),
                    ..Default::default()
                },
                active: true,
                min_aln_len: 2000,
                min_identity: 0.5,
                primary_only: true,
                remove_noise: false,
                sort_query: false,
                sort_reference: true,
                flipped: vec!["ctg2".into()],
                line_width: 2,
                drill_limit_mb: 150.0,
            }),
        }
    }

    #[test]
    fn round_trip_through_toml() {
        let s = fixture();
        let toml_text = toml::to_string_pretty(&s).unwrap();
        let back: Session = toml::from_str(&toml_text).unwrap();
        assert_eq!(back.version, s.version);
        assert_eq!(back.plot.matrix_name, s.plot.matrix_name);
        assert_eq!(back.view.crosshair, s.view.crosshair);
        assert_eq!(back.genome, s.genome);
        assert_eq!(back.sorted_records, s.sorted_records);
        assert_eq!(back.hsps, s.hsps);
        assert_eq!(
            back.large.unwrap().shown_classes,
            s.large.unwrap().shown_classes
        );
    }

    #[test]
    fn sessions_without_large_section_still_load() {
        let mut s = fixture();
        s.large = None;
        let toml_text = toml::to_string_pretty(&s).unwrap();
        assert!(!toml_text.contains("[large]"));
        let back: Session = toml::from_str(&toml_text).unwrap();
        assert!(back.large.is_none());
    }

    #[test]
    fn version_in_future_is_rejected() {
        let mut s = fixture();
        s.version = SESSION_VERSION + 1;
        let toml_text = toml::to_string_pretty(&s).unwrap();
        let dir = std::env::temp_dir();
        let path = dir.join("dottir_session_future.toml");
        std::fs::write(&path, toml_text).unwrap();
        let err = Session::load(&path).unwrap_err();
        assert!(matches!(err, SessionError::UnsupportedVersion(_, _)));
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn codec_round_trips() {
        use codec::*;
        use dottir_core::{BlastMode, Strand, Triangle};
        for s in [Strand::Forward, Strand::Reverse, Strand::Both] {
            assert_eq!(strand_from_str(strand_to_str(s)), Some(s));
        }
        for t in [Triangle::Both, Triangle::Upper, Triangle::Lower] {
            assert_eq!(triangle_from_str(triangle_to_str(t)), Some(t));
        }
        for m in [BlastMode::Blastn, BlastMode::Blastp, BlastMode::Blastx] {
            assert_eq!(mode_from_str(mode_to_str(m)), Some(m));
        }
    }
}
