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
        "pixelmap allocation of {requested} bytes exceeds memory_limit_bytes = {limit}; \
         try a larger zoom factor"
    )]
    OutOfMemory { requested: u64, limit: u64 },

    #[error("invalid configuration: {0}")]
    InvalidConfig(String),
}
