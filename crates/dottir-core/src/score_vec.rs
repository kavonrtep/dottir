//! Precomputed score vector: for each (subject_residue, query_position)
//! pair, the substitution score from the matrix.
//!
//! ## Layout
//!
//! `scores[row * qlen + col]`, where `row` is the encoded subject residue
//! and `col` is the query position. This is the cache-friendly flat layout
//! called out in spec §4.1.3 as a permitted deviation from C dotter's
//! row-of-pointers (`int **scoreVec`).
//!
//! ## Unknown residues
//!
//! Subject bytes that don't encode to a valid alphabet index (e.g. 'N' in
//! DNA, ambiguity codes that aren't in the score matrix) map to a synthetic
//! "unknown" row whose entries are all zero. This means the sliding-window
//! recurrence gracefully tolerates unscorable subject residues — they
//! neither contribute to the sum nor subtract from it. The C dotter
//! achieves the same outcome implicitly via row 5 (DNA) or row 24 (protein)
//! being left uninitialised; we make it explicit.

use crate::alphabet::{AlphabetKind, SENTINEL};
use crate::matrix::ScoreMatrix;

/// Flat precomputed score vector. The last row (index `unknown_row()`)
/// is all zeros and is used for any subject residue that didn't encode
/// to a valid alphabet index.
#[derive(Debug, Clone)]
pub struct ScoreVec {
    pub kind: AlphabetKind,
    /// `(n + 1) * qlen` entries, row-major. Row n is the "unknown" zero row.
    pub scores: Vec<i32>,
    pub qlen: usize,
    pub n_alphabet: usize,
}

impl ScoreVec {
    /// Number of scorable alphabet rows. Subject residues with index `>= this`
    /// hit the zero row.
    #[inline]
    pub fn n_alphabet(&self) -> usize {
        self.n_alphabet
    }

    /// Index of the synthetic "unknown" row.
    #[inline]
    pub fn unknown_row(&self) -> u8 {
        self.n_alphabet as u8
    }

    /// Return a slice covering `scores[row, ..]`.
    #[inline]
    pub fn row(&self, row: usize) -> &[i32] {
        let q = self.qlen;
        let start = row * q;
        &self.scores[start..start + q]
    }

    /// Build the score vector for a given query sequence (already encoded as
    /// alphabet indices, with [`SENTINEL`] for unknown bytes).
    pub fn build(matrix: &ScoreMatrix, query_encoded: &[u8]) -> Self {
        let n = matrix.size();
        let qlen = query_encoded.len();
        let mut scores = vec![0_i32; (n + 1) * qlen];

        for row in 0..n {
            let row_offset = row * qlen;
            for (q, &qb) in query_encoded.iter().enumerate() {
                // Query residues that don't have a column in the matrix
                // (SENTINEL bytes, or alphabet indices outside the
                // matrix's Karlin range — e.g. 'N' which encodes to 4
                // for DNA but the matrix is only 4×4 over ACGT)
                // contribute 0. Symmetric with the unknown-subject row.
                let v = if qb == SENTINEL || (qb as usize) >= n {
                    0
                } else {
                    matrix.get(row, qb as usize)
                };
                scores[row_offset + q] = v;
            }
        }
        // Row `n` is intentionally all zeros — the "unknown subject" row.
        ScoreVec {
            kind: matrix.kind,
            scores,
            qlen,
            n_alphabet: n,
        }
    }

    /// Map an encoded subject residue to a row index that is always valid
    /// (in `0..=n_alphabet`). SENTINEL maps to the zero row.
    #[inline]
    pub fn subject_row(&self, subject_residue: u8) -> u8 {
        if subject_residue == SENTINEL || (subject_residue as usize) >= self.n_alphabet {
            self.n_alphabet as u8
        } else {
            subject_residue
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::{encode, AlphabetKind};

    #[test]
    fn dna_score_vec_basic() {
        // Identity DNA matrix gives +5 on the diagonal, -4 off.
        let m = ScoreMatrix::dna_identity();
        let q = encode(b"ACGT", AlphabetKind::Dna);
        let sv = ScoreVec::build(&m, &q);
        assert_eq!(sv.n_alphabet(), 4);
        assert_eq!(sv.qlen, 4);
        // Row 0 = scores against 'A' subject.
        // q='A' → +5, q='C' → -4, q='G' → -4, q='T' → -4.
        assert_eq!(sv.row(0), &[5, -4, -4, -4]);
        assert_eq!(sv.row(1), &[-4, 5, -4, -4]); // subject C
        assert_eq!(sv.row(2), &[-4, -4, 5, -4]); // subject G
        assert_eq!(sv.row(3), &[-4, -4, -4, 5]); // subject T
        // Row 4 (unknown) is all zero.
        assert_eq!(sv.row(4), &[0, 0, 0, 0]);
    }

    #[test]
    fn unknown_query_residues_score_zero() {
        let m = ScoreMatrix::dna_identity();
        let q = encode(b"ANGN", AlphabetKind::Dna); // positions 1 and 3 are N
        let sv = ScoreVec::build(&m, &q);
        // Subject 'A' (row 0): q[0]=A→5, q[1]=N→0, q[2]=G→-4, q[3]=N→0.
        assert_eq!(sv.row(0), &[5, 0, -4, 0]);
    }

    #[test]
    fn subject_row_maps_unknown_to_zero_row() {
        let m = ScoreMatrix::dna_identity();
        let q = encode(b"ACGT", AlphabetKind::Dna);
        let sv = ScoreVec::build(&m, &q);
        assert_eq!(sv.subject_row(0), 0);
        assert_eq!(sv.subject_row(3), 3);
        // N (encoded as 4) is outside the 4-letter Karlin alphabet.
        assert_eq!(sv.subject_row(4), 4);
        assert_eq!(sv.subject_row(SENTINEL), 4);
    }

    #[test]
    fn protein_score_vec_shape() {
        let m = ScoreMatrix::blosum62();
        let q = encode(b"MKTAYIAKQ", AlphabetKind::Protein);
        let sv = ScoreVec::build(&m, &q);
        assert_eq!(sv.n_alphabet(), 24);
        assert_eq!(sv.qlen, 9);
        assert_eq!(sv.scores.len(), 25 * 9);
        // Last row (unknown) is all zero.
        assert!(sv.row(24).iter().all(|&v| v == 0));
    }
}
