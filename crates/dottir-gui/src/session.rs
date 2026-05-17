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
    /// Crosshair pixelmap coord, if set.
    pub crosshair: Option<[u32; 2]>,
    pub light_theme: bool,
}

pub const SESSION_VERSION: u32 = 1;

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
            },
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
