//! Public error type for `dottir-core`.
//!
//! Per spec §4.5.6 we use `thiserror` for typed errors crossing the public
//! API boundary; binaries wrap these with `anyhow`.

use thiserror::Error;

#[derive(Debug, Error)]
pub enum DottirError {
    #[error("not implemented yet: {0}")]
    NotImplemented(&'static str),

    #[error("Karlin/Altschul statistics failed: {0}")]
    KarlinFailure(String),

    #[error("invalid score matrix: {0}")]
    InvalidMatrix(String),

    #[error("sequence is empty")]
    EmptySequence,

    #[error("no scorable residues in sequence (alphabet mismatch?)")]
    NoScorableResidues,

    #[error(
        "pixelmap allocation of {requested} bytes ({channels} channel(s) × {per_channel} bytes) \
         exceeds memory_limit_bytes = {limit}; try a larger zoom factor or raise the limit"
    )]
    OutOfMemory {
        /// Total bytes the compute would allocate (`channels × per_channel`).
        requested: u64,
        /// Bytes per channel — `width × height`.
        per_channel: u64,
        /// Number of distinct pixelmaps the compute will hold at peak.
        /// `1` for the typical case, `2` for forward + separate reverse.
        channels: u32,
        /// The active cap (typically [`PlotConfig::memory_limit_bytes`]).
        limit: u64,
    },

    #[error("invalid configuration: {0}")]
    InvalidConfig(String),
}
