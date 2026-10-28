//! The per-user cache directory (genome-mode alignments, FASTA indexes
//! for read-only inputs, decompressed copies of gzipped FASTA).

use std::path::PathBuf;

/// `$DOTTIR_CACHE_DIR`, else `$XDG_CACHE_HOME/dottir`, else
/// `~/.cache/dottir` (Unix) or `%LOCALAPPDATA%\dottir\cache` (Windows),
/// else a `dottir-cache` folder in the system temp directory.
pub fn cache_dir() -> PathBuf {
    let env = |k: &str| {
        std::env::var_os(k)
            .filter(|v| !v.is_empty())
            .map(PathBuf::from)
    };
    if let Some(d) = env("DOTTIR_CACHE_DIR") {
        return d;
    }
    if let Some(d) = env("XDG_CACHE_HOME") {
        return d.join("dottir");
    }
    if cfg!(windows) {
        if let Some(d) = env("LOCALAPPDATA") {
            return d.join("dottir").join("cache");
        }
    } else if let Some(h) = env("HOME") {
        return h.join(".cache").join("dottir");
    }
    std::env::temp_dir().join("dottir-cache")
}
