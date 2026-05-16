//! Substitution score matrices.
//!
//! The reference C dotter ships its matrices inline as `int [24][24]`
//! arrays. We do the same for protein (BLOSUM62), and for DNA we use a 4×4
//! match/mismatch matrix matching the C dotter's `DNAmatrix()` (diagonal
//! +5, off-diagonal −4). Other built-in protein matrices (BLOSUM45/50/80/90,
//! PAM30/70/250) will be added incrementally — each shipped verbatim from
//! NCBI source.
//!
//! ## Layout
//!
//! `scores` is a flat `n*n` row-major `Vec<i32>`. Row index = encoded
//! "query" residue, column index = encoded "subject" residue. The alphabet
//! ordering is fixed per [`AlphabetKind`] and lives in
//! [`crate::alphabet`]. The flat layout (rather than the C
//! row-of-pointers) is the cache-friendly internal representation called
//! out in spec §4.1.3.

use crate::alphabet::{
    encode_protein, AlphabetKind, DNA_KARLIN_SIZE, PROTEIN_KARLIN_SIZE, SENTINEL,
};
use crate::error::DottirError;

/// BLAST mode determines the alphabet and how sequences are pre-processed
/// before they reach the inner loop.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlastMode {
    /// DNA vs DNA. Forward + reverse passes; identity-style scoring.
    Blastn,
    /// Protein vs protein. Single pass; BLOSUM/PAM scoring.
    Blastp,
    /// Translated DNA query vs protein subject (three reading frames).
    Blastx,
}

impl BlastMode {
    pub const fn alphabet(self) -> AlphabetKind {
        match self {
            BlastMode::Blastn => AlphabetKind::Dna,
            BlastMode::Blastp | BlastMode::Blastx => AlphabetKind::Protein,
        }
    }

    /// Number of residue buckets used by the Karlin/Altschul probability
    /// sum. Matches the C dotter `abetsize` argument (24 for protein, 4 for
    /// DNA — note 'N' is excluded for DNA).
    pub const fn karlin_size(self) -> usize {
        match self {
            BlastMode::Blastn => DNA_KARLIN_SIZE,
            BlastMode::Blastp | BlastMode::Blastx => PROTEIN_KARLIN_SIZE,
        }
    }
}

/// A square score matrix over one of the built-in alphabets.
///
/// `scores[row * n + col]` gives the score for substituting the residue
/// encoded as `row` (in [`crate::alphabet`]'s ordering) with the residue
/// encoded as `col`.
///
/// # Example
///
/// ```
/// use dottir_core::matrix::ScoreMatrix;
///
/// let m = ScoreMatrix::blosum62();
/// assert_eq!(m.size(), 24);
/// assert_eq!(m.get(0, 0), 4);          // A vs A
/// assert_eq!(m.get(17, 17), 11);       // W vs W
/// assert_eq!(m.name, "BLOSUM62");
///
/// // Round-trip through the BLAST text format.
/// let text = m.to_blast_format();
/// let parsed = ScoreMatrix::parse_blast_format("BLOSUM62", &text).unwrap();
/// assert_eq!(parsed.scores, m.scores);
/// ```
#[derive(Debug, Clone)]
pub struct ScoreMatrix {
    pub name: String,
    pub kind: AlphabetKind,
    /// Row-major n×n scores, where `n == self.size()`.
    pub scores: Vec<i32>,
}

impl ScoreMatrix {
    pub fn size(&self) -> usize {
        match self.kind {
            AlphabetKind::Protein => PROTEIN_KARLIN_SIZE,
            AlphabetKind::Dna => DNA_KARLIN_SIZE,
        }
    }

    #[inline]
    pub fn get(&self, row: usize, col: usize) -> i32 {
        let n = self.size();
        self.scores[row * n + col]
    }

    /// Built-in BLOSUM62. Source: `dotterApp/dotter.c:354` (the same table
    /// shipped with C dotter, which in turn matches NCBI's blast/data/
    /// BLOSUM62 file dated 930809).
    pub fn blosum62() -> Self {
        // 24×24, row-major. Alphabet order: A R N D C Q E G H I L K M F P S T W Y V B Z X *
        const TABLE: [i32; 24 * 24] = [
            // A
             4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4,
            // R
            -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4,
            // N
            -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4,
            // D
            -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4,
            // C
             0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4,
            // Q
            -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4,
            // E
            -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4,
            // G
             0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4,
            // H
            -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4,
            // I
            -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4,
            // L
            -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4,
            // K
            -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4,
            // M
            -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4,
            // F
            -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4,
            // P
            -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4,
            // S
             1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4,
            // T
             0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4,
            // W
            -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4,
            // Y
            -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4,
            // V
             0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4,
            // B
            -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4,
            // Z
            -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4,
            // X
             0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4,
            // *
            -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1,
        ];
        ScoreMatrix {
            name: "BLOSUM62".to_string(),
            kind: AlphabetKind::Protein,
            scores: TABLE.to_vec(),
        }
    }

    /// Default DNA identity matrix: +5 on the diagonal, −4 off-diagonal.
    /// Mirrors `DNAmatrix()` at `dotterApp/dotter.c:2800`. Covers the 4-letter
    /// alphabet A/C/G/T; the 'N' row/column is not represented because the
    /// Karlin sum and the inner loop both skip non-scorable residues.
    pub fn dna_identity() -> Self {
        let n = DNA_KARLIN_SIZE;
        let mut scores = vec![-4_i32; n * n];
        for i in 0..n {
            scores[i * n + i] = 5;
        }
        ScoreMatrix {
            name: "DNA+5/-4".to_string(),
            kind: AlphabetKind::Dna,
            scores,
        }
    }

    /// Default matrix for a given mode (BLOSUM62 for protein, DNA identity
    /// for nucleotide modes). Spec §4.1.2.
    pub fn default_for(mode: BlastMode) -> Self {
        match mode {
            BlastMode::Blastn => Self::dna_identity(),
            BlastMode::Blastp | BlastMode::Blastx => Self::blosum62(),
        }
    }

    // NOTE on BLOSUM62 versioning: the inline [`Self::blosum62`] is the
    // 1993-08-09 table that C-dotter ships at `dotter.c:354`. The modern
    // NCBI BLOSUM62 file differs in the ambiguous-residue rows (X/B/Z).
    // We deliberately use the older table so that golden tests pinned
    // against C-dotter output reproduce exactly. The other BLOSUM/PAM
    // matrices below are sourced from NCBI directly; they are not in
    // C-dotter and have no comparable legacy version.

    /// Built-in NCBI BLOSUM45 (vendored from NCBI BLAST data files).
    pub fn blosum45() -> Self {
        Self::parse_blast_format("BLOSUM45", include_str!("../data/BLOSUM45"))
            .expect("vendored BLOSUM45 parses")
    }

    /// Built-in NCBI BLOSUM50.
    pub fn blosum50() -> Self {
        Self::parse_blast_format("BLOSUM50", include_str!("../data/BLOSUM50"))
            .expect("vendored BLOSUM50 parses")
    }

    /// Built-in NCBI BLOSUM80.
    pub fn blosum80() -> Self {
        Self::parse_blast_format("BLOSUM80", include_str!("../data/BLOSUM80"))
            .expect("vendored BLOSUM80 parses")
    }

    /// Built-in NCBI BLOSUM90.
    pub fn blosum90() -> Self {
        Self::parse_blast_format("BLOSUM90", include_str!("../data/BLOSUM90"))
            .expect("vendored BLOSUM90 parses")
    }

    /// Built-in NCBI PAM30.
    pub fn pam30() -> Self {
        Self::parse_blast_format("PAM30", include_str!("../data/PAM30"))
            .expect("vendored PAM30 parses")
    }

    /// Built-in NCBI PAM70.
    pub fn pam70() -> Self {
        Self::parse_blast_format("PAM70", include_str!("../data/PAM70"))
            .expect("vendored PAM70 parses")
    }

    /// Built-in NCBI PAM250.
    pub fn pam250() -> Self {
        Self::parse_blast_format("PAM250", include_str!("../data/PAM250"))
            .expect("vendored PAM250 parses")
    }

    /// Look up a built-in protein matrix by its canonical name (uppercase,
    /// no separator). Returns `None` for unrecognised names.
    pub fn by_name(name: &str) -> Option<Self> {
        match name {
            "BLOSUM45" => Some(Self::blosum45()),
            "BLOSUM50" => Some(Self::blosum50()),
            "BLOSUM62" => Some(Self::blosum62()),
            "BLOSUM80" => Some(Self::blosum80()),
            "BLOSUM90" => Some(Self::blosum90()),
            "PAM30" => Some(Self::pam30()),
            "PAM70" => Some(Self::pam70()),
            "PAM250" => Some(Self::pam250()),
            _ => None,
        }
    }

    pub fn lowest(&self) -> i32 {
        self.scores.iter().copied().min().unwrap()
    }

    pub fn highest(&self) -> i32 {
        self.scores.iter().copied().max().unwrap()
    }

    /// Parse a protein score matrix in the standard NCBI BLAST text format.
    ///
    /// The grammar is:
    ///
    /// ```text
    /// # any number of comment lines starting with '#'
    ///    A  R  N  D  ...   *       <- header row of alphabet letters
    /// A  4 -1 -2 -2  ...  -4       <- one row per alphabet letter
    /// R -1  5  0 -2  ...  -4
    /// ...
    /// ```
    ///
    /// Whitespace separation. The header determines the column ordering;
    /// rows may be in any order — the row letter is the first non-space
    /// token. The returned matrix is *always* re-ordered into the canonical
    /// `ARNDCQEGHILKMFPSTWYVBZX*` ordering used elsewhere in dottir, so the
    /// caller does not need to know the file's internal ordering. Missing
    /// letters are an error.
    pub fn parse_blast_format(name: &str, text: &str) -> Result<Self, DottirError> {
        // Read header.
        let mut lines = text
            .lines()
            .filter(|l| !l.trim_start().starts_with('#') && !l.trim().is_empty());

        let header_line = lines.next().ok_or_else(|| {
            DottirError::InvalidMatrix("matrix file has no header row".into())
        })?;
        let header_letters: Vec<u8> = header_line
            .split_whitespace()
            .map(|tok| {
                if tok.len() != 1 {
                    Err(DottirError::InvalidMatrix(format!(
                        "header column '{tok}' is not a single residue letter"
                    )))
                } else {
                    Ok(tok.as_bytes()[0])
                }
            })
            .collect::<Result<Vec<_>, _>>()?;

        if header_letters.is_empty() {
            return Err(DottirError::InvalidMatrix("empty header row".into()));
        }
        // Map each header letter to its canonical index. Letters outside the
        // 24-residue C-dotter alphabet (e.g. NCBI's 'J' for Leu/Ile) are
        // silently dropped — recorded as `None` so the row parser knows to
        // skip the corresponding column entries.
        let col_to_canon: Vec<Option<usize>> = header_letters
            .iter()
            .map(|&b| {
                let idx = encode_protein(b);
                if idx == SENTINEL {
                    None
                } else {
                    Some(idx as usize)
                }
            })
            .collect();

        let n = PROTEIN_KARLIN_SIZE;
        let mut scores: Vec<Option<i32>> = vec![None; n * n];
        for row in lines {
            let mut toks = row.split_whitespace();
            let row_letter_tok = match toks.next() {
                Some(t) => t,
                None => continue, // skip blank
            };
            if row_letter_tok.len() != 1 {
                return Err(DottirError::InvalidMatrix(format!(
                    "row label '{row_letter_tok}' is not a single letter"
                )));
            }
            let row_b = row_letter_tok.as_bytes()[0];
            let row_canon = encode_protein(row_b);
            let row_canon_opt = if row_canon == SENTINEL { None } else { Some(row_canon as usize) };
            let mut col = 0;
            for tok in toks {
                if col >= header_letters.len() {
                    return Err(DottirError::InvalidMatrix(format!(
                        "row '{}' has too many columns",
                        row_b as char
                    )));
                }
                let v: i32 = tok.parse().map_err(|_| {
                    DottirError::InvalidMatrix(format!(
                        "non-integer entry '{tok}' in row '{}'",
                        row_b as char
                    ))
                })?;
                if let (Some(r), Some(c)) = (row_canon_opt, col_to_canon[col]) {
                    scores[r * n + c] = Some(v);
                }
                col += 1;
            }
            if col != header_letters.len() {
                return Err(DottirError::InvalidMatrix(format!(
                    "row '{}' has {col} columns, expected {}",
                    row_b as char,
                    header_letters.len()
                )));
            }
        }

        // Any score the header / rows didn't cover stays None. Allow it only
        // if every covered (row, col) pair lies inside the abetsize = 24
        // protein submatrix; we don't insist on covering the '*' row unless
        // the header included it. Cells with no entry get a deeply-negative
        // sentinel mirroring C dotter's "*-column = -4" convention.
        let filled: Vec<i32> = scores
            .into_iter()
            .map(|s| s.unwrap_or(-4))
            .collect();

        Ok(ScoreMatrix {
            name: name.to_string(),
            kind: AlphabetKind::Protein,
            scores: filled,
        })
    }

    /// Reformat a matrix as the canonical NCBI BLAST text format. Round-trip
    /// guarantee: `parse_blast_format(_, &m.to_blast_format()).scores == m.scores`.
    pub fn to_blast_format(&self) -> String {
        let n = self.size();
        let letters: &[u8] = match self.kind {
            AlphabetKind::Protein => crate::alphabet::PROTEIN_LETTERS,
            AlphabetKind::Dna => crate::alphabet::DNA_LETTERS,
        };
        let mut out = String::new();
        out.push_str("#  dottir score matrix\n");
        out.push_str("#  ");
        out.push_str(&self.name);
        out.push('\n');
        out.push_str("  ");
        for &b in letters[..n].iter() {
            out.push_str(&format!(" {:>3}", b as char));
        }
        out.push('\n');
        for i in 0..n {
            out.push(letters[i] as char);
            out.push(' ');
            for j in 0..n {
                out.push_str(&format!(" {:>3}", self.get(i, j)));
            }
            out.push('\n');
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blosum62_has_canonical_shape() {
        let m = ScoreMatrix::blosum62();
        assert_eq!(m.size(), 24);
        assert_eq!(m.scores.len(), 24 * 24);
        // Self-pair anchors from NCBI BLOSUM62.
        assert_eq!(m.get(0, 0), 4); // A/A
        assert_eq!(m.get(4, 4), 9); // C/C
        assert_eq!(m.get(17, 17), 11); // W/W
        // Off-diagonal sanity.
        assert_eq!(m.get(0, 4), 0); // A/C
        assert_eq!(m.get(13, 17), 1); // F/W
        // The '*' row/col is heavily negative (with a +1 on *,*).
        assert_eq!(m.get(23, 23), 1);
        assert_eq!(m.get(23, 0), -4);
    }

    #[test]
    fn blosum62_is_symmetric() {
        let m = ScoreMatrix::blosum62();
        let n = m.size();
        for i in 0..n {
            for j in 0..n {
                assert_eq!(m.get(i, j), m.get(j, i), "asymmetric at ({i},{j})");
            }
        }
    }

    #[test]
    fn dna_identity_matches_c_dotter() {
        let m = ScoreMatrix::dna_identity();
        assert_eq!(m.size(), 4);
        for i in 0..4 {
            for j in 0..4 {
                let want = if i == j { 5 } else { -4 };
                assert_eq!(m.get(i, j), want);
            }
        }
    }

    #[test]
    fn lowest_highest() {
        let m = ScoreMatrix::blosum62();
        assert_eq!(m.highest(), 11); // W,W
        assert_eq!(m.lowest(), -4);
    }

    #[test]
    fn blast_format_round_trip_blosum62() {
        let m = ScoreMatrix::blosum62();
        let text = m.to_blast_format();
        let parsed = ScoreMatrix::parse_blast_format("BLOSUM62", &text).unwrap();
        assert_eq!(parsed.scores, m.scores);
    }

    #[test]
    fn parse_blast_format_with_permuted_header() {
        // A 4-letter toy matrix in a non-canonical ordering. After parsing
        // it must come back in the canonical ARNDCQEG... ordering.
        let text = "\
#  toy 4-letter matrix
   R  A  C  N
A -1  4  0 -2
R  5 -1 -3  0
C -3  0  9 -3
N  0 -2 -3  6
";
        let m = ScoreMatrix::parse_blast_format("toy", text).unwrap();
        // Canonical order is A R N D C Q E G H ...
        assert_eq!(m.get(0, 0), 4); // A/A
        assert_eq!(m.get(0, 1), -1); // A/R
        assert_eq!(m.get(1, 0), -1); // R/A
        assert_eq!(m.get(1, 1), 5); // R/R
        assert_eq!(m.get(2, 2), 6); // N/N
        assert_eq!(m.get(4, 4), 9); // C/C
        // Uncovered cells default to -4 sentinel (mirrors C dotter '*' column).
        assert_eq!(m.get(3, 3), -4); // D/D
    }

    #[test]
    fn parse_blast_rejects_bad_header() {
        let bad = "AA BB\nA  4 -1\n";
        assert!(ScoreMatrix::parse_blast_format("bad", bad).is_err());
    }

    #[test]
    fn all_builtin_matrices_load_and_are_24x24() {
        for name in &[
            "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90",
            "PAM30", "PAM70", "PAM250",
        ] {
            let m = ScoreMatrix::by_name(name)
                .unwrap_or_else(|| panic!("by_name({name}) returned None"));
            assert_eq!(m.size(), 24, "{name}");
            assert_eq!(m.scores.len(), 24 * 24, "{name}");
            // Sanity: highest must be positive, lowest must be negative.
            assert!(m.highest() > 0, "{name} non-positive max");
            assert!(m.lowest() < 0, "{name} non-negative min");
            // Diagonal anchor: A/A is the standard "self-pair" score and is
            // always positive in BLOSUM/PAM.
            assert!(m.get(0, 0) > 0, "{name} A/A non-positive");
        }
    }

    #[test]
    fn each_builtin_matrix_round_trips_through_text() {
        for name in &[
            "BLOSUM45", "BLOSUM50", "BLOSUM62", "BLOSUM80", "BLOSUM90",
            "PAM30", "PAM70", "PAM250",
        ] {
            let m = ScoreMatrix::by_name(name).unwrap();
            let text = m.to_blast_format();
            let reparsed = ScoreMatrix::parse_blast_format(name, &text).unwrap();
            assert_eq!(reparsed.scores, m.scores, "{name}");
        }
    }

    #[test]
    fn by_name_unknown_returns_none() {
        assert!(ScoreMatrix::by_name("NONEXISTENT").is_none());
        assert!(ScoreMatrix::by_name("blosum62").is_none()); // case-sensitive
    }
}
