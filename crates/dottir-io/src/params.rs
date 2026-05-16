//! Parameter-reproducibility sidecar (spec §4.4.7).
//!
//! Every CLI batch export writes a `<output>.params.toml` next to the
//! main artifact containing:
//!
//! * Input file paths and their SHA-256 hashes.
//! * All resolved plot parameters (window size, zoom, pixel factor,
//!   matrix name, mode, strand, …).
//! * dottir version + (when available) git SHA.
//! * Host info: hostname, OS, ISO-8601 timestamp.
//!
//! The TOML structure is stable across minor versions; new fields may be
//! added but existing ones won't change meaning.

use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParamsSidecar {
    pub dottir: DottirInfo,
    pub query: InputInfo,
    pub subject: InputInfo,
    pub plot: PlotParamsInfo,
    pub host: HostInfo,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DottirInfo {
    pub version: String,
    pub git_sha: Option<String>,
    pub pixelmap_format_version: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InputInfo {
    pub path: String,
    pub sha256: String,
    pub size_bytes: u64,
    pub n_records: usize,
    pub total_residues: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PlotParamsInfo {
    pub mode: String,
    pub matrix: String,
    pub strand: String,
    pub window_size: u32,
    pub zoom: u32,
    pub pixel_fac: u32,
    pub self_comparison: bool,
    pub karlin: Option<KarlinInfo>,
    pub width: u32,
    pub height: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KarlinInfo {
    pub lambda: f64,
    pub k: f64,
    pub h: f64,
    pub expected_residue_score: f64,
    pub expected_msp_score: f64,
    pub predicted_msp_length: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HostInfo {
    pub hostname: String,
    pub os: String,
    pub timestamp_utc: String,
}

impl ParamsSidecar {
    /// Serialize to a pretty TOML string ready for writing.
    pub fn to_toml(&self) -> Result<String> {
        toml::to_string_pretty(self).context("serializing params sidecar to TOML")
    }
}

/// Compute SHA-256 of a file as a lowercase hex string.
pub fn sha256_file<P: AsRef<Path>>(path: P) -> Result<String> {
    let bytes = std::fs::read(path.as_ref())
        .with_context(|| format!("reading {}", path.as_ref().display()))?;
    Ok(sha256_bytes(&bytes))
}

pub fn sha256_bytes(data: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(data);
    let hash = hasher.finalize();
    let mut hex = String::with_capacity(64);
    for byte in hash {
        hex.push_str(&format!("{byte:02x}"));
    }
    hex
}

/// Best-effort hostname lookup. Falls back to "unknown" on platforms
/// where the environment variable isn't set.
pub fn hostname() -> String {
    std::env::var("HOSTNAME")
        .or_else(|_| std::env::var("COMPUTERNAME"))
        .unwrap_or_else(|_| {
            std::process::Command::new("hostname")
                .output()
                .ok()
                .and_then(|o| String::from_utf8(o.stdout).ok())
                .map(|s| s.trim().to_string())
                .unwrap_or_else(|| "unknown".to_string())
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sha256_known_vectors() {
        // Empty input.
        assert_eq!(
            sha256_bytes(b""),
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
        // "abc" — the canonical SHA-256 test vector.
        assert_eq!(
            sha256_bytes(b"abc"),
            "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
        );
    }

    #[test]
    fn params_round_trip_through_toml() {
        let s = ParamsSidecar {
            dottir: DottirInfo {
                version: "0.1.0".into(),
                git_sha: Some("abc1234".into()),
                pixelmap_format_version: 0,
            },
            query: InputInfo {
                path: "/tmp/q.fa".into(),
                sha256: "0".repeat(64),
                size_bytes: 100,
                n_records: 1,
                total_residues: 100,
            },
            subject: InputInfo {
                path: "/tmp/s.fa".into(),
                sha256: "1".repeat(64),
                size_bytes: 200,
                n_records: 1,
                total_residues: 200,
            },
            plot: PlotParamsInfo {
                mode: "Blastn".into(),
                matrix: "DNA+5/-4".into(),
                strand: "Both".into(),
                window_size: 25,
                zoom: 1,
                pixel_fac: 50,
                self_comparison: false,
                karlin: None,
                width: 100,
                height: 200,
            },
            host: HostInfo {
                hostname: "test-host".into(),
                os: "linux".into(),
                timestamp_utc: "2026-05-16T00:00:00Z".into(),
            },
        };
        let toml_text = s.to_toml().unwrap();
        let back: ParamsSidecar = toml::from_str(&toml_text).unwrap();
        assert_eq!(back.plot.window_size, s.plot.window_size);
        assert_eq!(back.query.sha256, s.query.sha256);
    }
}
