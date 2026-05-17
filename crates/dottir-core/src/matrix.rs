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

    /// Build a custom DNA matrix with the given match / mismatch scores.
    /// Convenient for users who want `+1 / -1`, `+2 / -3`, etc. without
    /// hand-writing a BLAST-format file.
    ///
    /// Both diagonals are filled with `match_score`, all off-diagonals
    /// with `mismatch_score`. The matrix is 4×4 over the canonical ACGT
    /// ordering (matches [`crate::alphabet::DNA_LETTERS`]).
    pub fn custom_dna(match_score: i32, mismatch_score: i32) -> Self {
        let n = DNA_KARLIN_SIZE;
        let mut scores = vec![mismatch_score; n * n];
        for i in 0..n {
            scores[i * n + i] = match_score;
        }
        ScoreMatrix {
            name: format!("DNA+{match_score}/{mismatch_score}"),
            kind: AlphabetKind::Dna,
            scores,
        }
    }

    /// Validate the matrix against the Karlin/Altschul prerequisites:
    /// at least one negative score, at least one positive score, and
    /// dimensions consistent with the declared alphabet. Returns an
    /// error rather than a corrupted result so callers can fail loudly
    /// on hand-written matrices that wouldn't produce useful goldens.
    pub fn validate(&self) -> Result<(), DottirError> {
        let n = self.size();
        if self.scores.len() != n * n {
            return Err(DottirError::InvalidMatrix(format!(
                "matrix '{}': scores len {} != n×n = {}",
                self.name,
                self.scores.len(),
                n * n
            )));
        }
        let any_neg = self.scores.iter().any(|&s| s < 0);
        let any_pos = self.scores.iter().any(|&s| s > 0);
        if !any_neg {
            return Err(DottirError::InvalidMatrix(format!(
                "matrix '{}' has no negative scores — Karlin/Altschul \
                 statistics require at least one (otherwise λ doesn't \
                 converge). Most BLAST matrices have mismatch ≤ -1.",
                self.name
            )));
        }
        if !any_pos {
            return Err(DottirError::InvalidMatrix(format!(
                "matrix '{}' has no positive scores — there's no way to \
                 score a match. The diagonal should be > 0.",
                self.name
            )));
        }
        Ok(())
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

    /// Parse a protein score matrix in the standard NCBI BLAST text
    /// format. **Strict** — every cell of the 24×24 protein submatrix
    /// must be covered, and any header letter not in dottir's 24-letter
    /// alphabet causes a parse error unless it appears in
    /// `extra_letters`.
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
    /// token. The returned matrix is *always* re-ordered into the
    /// canonical `ARNDCQEGHILKMFPSTWYVBZX*` ordering used elsewhere in
    /// dottir.
    ///
    /// `extra_letters` lets the caller opt into header columns that are
    /// outside the 24-letter dottir alphabet — pass `&["J"]` for current
    /// NCBI BLOSUM/PAM files, which include 'J' (Leu/Ile ambiguity).
    /// Anything in `extra_letters` is silently dropped from the row/col
    /// data; an unexpected letter not in `extra_letters` is an error.
    ///
    /// After parsing the matrix is checked with [`Self::validate`].
    pub fn parse_blast_protein(
        name: &str,
        text: &str,
        extra_letters: &[&str],
    ) -> Result<Self, DottirError> {
        let n = PROTEIN_KARLIN_SIZE;
        let (header_letters, body) = read_header_and_body(text)?;
        let col_to_canon = map_header(&header_letters, extra_letters, encode_protein, n)?;
        let scores =
            parse_rows(body, &header_letters, &col_to_canon, n, extra_letters, encode_protein)?;
        let m = ScoreMatrix {
            name: name.to_string(),
            kind: AlphabetKind::Protein,
            scores,
        };
        m.validate()?;
        Ok(m)
    }

    /// Parse a DNA score matrix in BLAST text format over the 4-letter
    /// `ACGT` alphabet. Strict in the same way as
    /// [`Self::parse_blast_protein`]: every 4×4 cell must be covered;
    /// extra header letters (e.g. `N`) must be passed in
    /// `extra_letters` and are silently dropped.
    pub fn parse_blast_dna(
        name: &str,
        text: &str,
        extra_letters: &[&str],
    ) -> Result<Self, DottirError> {
        let n = DNA_KARLIN_SIZE;
        let (header_letters, body) = read_header_and_body(text)?;
        let col_to_canon = map_header(
            &header_letters,
            extra_letters,
            |b| crate::alphabet::encode_dna(b),
            n,
        )?;
        let scores = parse_rows(
            body,
            &header_letters,
            &col_to_canon,
            n,
            extra_letters,
            |b| crate::alphabet::encode_dna(b),
        )?;
        let m = ScoreMatrix {
            name: name.to_string(),
            kind: AlphabetKind::Dna,
            scores,
        };
        m.validate()?;
        Ok(m)
    }

    /// Back-compat alias: parse a protein matrix accepting NCBI's `J`
    /// column. Pre-existing callers (the built-in BLOSUM/PAM
    /// constructors) used this signature.
    pub fn parse_blast_format(name: &str, text: &str) -> Result<Self, DottirError> {
        Self::parse_blast_protein(name, text, &["J"])
    }

}

/// Helpers shared by the protein and DNA parsers.
type EncodeFn = fn(u8) -> u8;

/// Split a BLAST text matrix into its header letters and the remainder
/// lines (comments stripped). The body lines are owned `String`s
/// because the iterator passes through the input twice in
/// [`parse_rows`] otherwise — owned is simpler than juggling lifetimes.
fn read_header_and_body(text: &str) -> Result<(Vec<u8>, Vec<String>), DottirError> {
    let mut lines = text
        .lines()
        .filter(|l| !l.trim_start().starts_with('#') && !l.trim().is_empty());
    let header_line = lines
        .next()
        .ok_or_else(|| DottirError::InvalidMatrix("matrix file has no header row".into()))?;
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
    Ok((header_letters, lines.map(|l| l.to_string()).collect()))
}

/// Map each header column to its canonical alphabet index, or `None`
/// for letters that the caller explicitly allowlisted as "extra"
/// (e.g. NCBI's J for protein). Letters that don't encode AND aren't
/// in `extra_letters` are a parse error.
fn map_header(
    header_letters: &[u8],
    extra_letters: &[&str],
    encode: EncodeFn,
    n: usize,
) -> Result<Vec<Option<usize>>, DottirError> {
    let extras = extras_set(extra_letters);
    header_letters
        .iter()
        .map(|&b| {
            let idx = encode(b);
            // "In the Karlin alphabet" = encoded AND inside the n×n
            // submatrix we score against. DNA's 'N' encodes to 4 but
            // the matrix is only 4×4 (ACGT), so it falls in the
            // `extra_letters` path.
            if idx != SENTINEL && (idx as usize) < n {
                Ok(Some(idx as usize))
            } else if extras.contains(&b) {
                Ok(None)
            } else {
                Err(DottirError::InvalidMatrix(format!(
                    "header has unknown residue letter '{}'. \
                     Pass it in `extra_letters` to ignore.",
                    b as char
                )))
            }
        })
        .collect()
}

fn extras_set(extra_letters: &[&str]) -> std::collections::HashSet<u8> {
    extra_letters
        .iter()
        .filter_map(|s| {
            if s.len() == 1 {
                Some(s.as_bytes()[0])
            } else {
                None
            }
        })
        .collect()
}

/// Parse the body of a BLAST matrix file into a `n × n` row-major
/// score vector. Every cell of the canonical alphabet must be covered,
/// otherwise an [`DottirError::InvalidMatrix`] names the missing
/// (row, col).
fn parse_rows(
    body: Vec<String>,
    header_letters: &[u8],
    col_to_canon: &[Option<usize>],
    n: usize,
    extra_letters: &[&str],
    encode: EncodeFn,
) -> Result<Vec<i32>, DottirError> {
    let extras = extras_set(extra_letters);
    let mut scores: Vec<Option<i32>> = vec![None; n * n];
    for row in body {
        let mut toks = row.split_whitespace();
        let row_letter_tok = match toks.next() {
            Some(t) => t,
            None => continue,
        };
        if row_letter_tok.len() != 1 {
            return Err(DottirError::InvalidMatrix(format!(
                "row label '{row_letter_tok}' is not a single letter"
            )));
        }
        let row_b = row_letter_tok.as_bytes()[0];
        let row_canon = encode(row_b);
        let row_canon_opt = if row_canon != SENTINEL && (row_canon as usize) < n {
            Some(row_canon as usize)
        } else if extras.contains(&row_b) {
            None
        } else {
            return Err(DottirError::InvalidMatrix(format!(
                "row label '{}' is not in the alphabet. Pass it in \
                 `extra_letters` to ignore.",
                row_b as char
            )));
        };
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
    // Every cell of the canonical n×n submatrix must be present —
    // silent fill-with-default was the bug REVIEW.md called out.
    for r in 0..n {
        for c in 0..n {
            if scores[r * n + c].is_none() {
                return Err(DottirError::InvalidMatrix(format!(
                    "matrix has no entry for ({}, {}) — every cell of the \
                     {}×{} submatrix must be covered",
                    canonical_letter(r, n) as char,
                    canonical_letter(c, n) as char,
                    n,
                    n
                )));
            }
        }
    }
    Ok(scores.into_iter().map(|s| s.unwrap()).collect())
}

/// Map a canonical alphabet index back to its ASCII letter for error
/// messages. The protein alphabet (n=24) uses [`crate::alphabet::PROTEIN_LETTERS`];
/// DNA (n=4) uses [`crate::alphabet::DNA_LETTERS`].
fn canonical_letter(idx: usize, n: usize) -> u8 {
    if n == PROTEIN_KARLIN_SIZE {
        crate::alphabet::PROTEIN_LETTERS[idx]
    } else {
        crate::alphabet::DNA_LETTERS[idx]
    }
}

impl ScoreMatrix {
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
    fn parse_blast_dna_with_permuted_header() {
        // 4-letter DNA matrix in non-canonical column ordering.
        // After parsing it must come back in the canonical ACGT order.
        let text = "\
#  toy DNA matrix in non-canonical column ordering
   T  A  C  G
A -4  5 -4 -4
C -4 -4  5 -4
G -4 -4 -4  5
T  5 -4 -4 -4
";
        let m = ScoreMatrix::parse_blast_dna("toy_dna", text, &[]).unwrap();
        // Canonical order is A C G T.
        assert_eq!(m.get(0, 0), 5); // A/A
        assert_eq!(m.get(1, 1), 5); // C/C
        assert_eq!(m.get(2, 2), 5); // G/G
        assert_eq!(m.get(3, 3), 5); // T/T
        // Off-diagonal sample.
        assert_eq!(m.get(0, 3), -4); // A/T
        assert_eq!(m.get(3, 0), -4); // T/A
    }

    #[test]
    fn parse_blast_dna_with_extra_n_column() {
        // The 'N' ambiguity column must be passed via extra_letters
        // or it's an error. With extra_letters = &["N"] it's dropped.
        let text = "\
   A  C  G  T  N
A  5 -4 -4 -4  0
C -4  5 -4 -4  0
G -4 -4  5 -4  0
T -4 -4 -4  5  0
N  0  0  0  0  0
";
        // Without the N allow-list, the parser refuses.
        let err = ScoreMatrix::parse_blast_dna("dna_with_n", text, &[]).unwrap_err();
        assert!(format!("{err}").contains("unknown residue letter"));
        // With the allow-list, the matrix loads cleanly.
        let m = ScoreMatrix::parse_blast_dna("dna_with_n", text, &["N"]).unwrap();
        assert_eq!(m.size(), 4);
        assert_eq!(m.get(0, 0), 5);
    }

    #[test]
    fn parse_blast_protein_rejects_unknown_header_letter() {
        // 'J' is in NCBI matrices but not in dottir's 24-letter alphabet;
        // strict mode rejects unless allowlisted.
        let bad = "\
   A  R  J
A  4 -1  0
R -1  5  0
J  0  0  0
";
        let err = ScoreMatrix::parse_blast_protein("bad", bad, &[]).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("unknown residue letter"), "got: {msg}");
    }

    #[test]
    fn parse_blast_protein_rejects_missing_cell() {
        // A protein matrix that omits 'D' from the header → D row/column
        // entries aren't covered → strict parse refuses.
        let bad = "\
   A  R
A  4 -1
R -1  5
";
        let err = ScoreMatrix::parse_blast_protein("bad", bad, &[]).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("no entry for"), "got: {msg}");
    }

    #[test]
    fn parse_blast_rejects_bad_header() {
        let bad = "AA BB\nA  4 -1\n";
        assert!(ScoreMatrix::parse_blast_format("bad", bad).is_err());
    }

    #[test]
    fn custom_dna_diagonal_and_off_diagonal() {
        let m = ScoreMatrix::custom_dna(2, -3);
        assert_eq!(m.size(), 4);
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(m.get(i, j), if i == j { 2 } else { -3 });
            }
        }
        assert!(m.validate().is_ok());
    }

    #[test]
    fn validate_rejects_all_positive_or_all_negative() {
        let mut m = ScoreMatrix::dna_identity();
        // Make everything positive (no mismatch penalty) → validate fails.
        for s in &mut m.scores {
            if *s < 0 {
                *s = 1;
            }
        }
        let err = m.validate().unwrap_err();
        assert!(format!("{err}").contains("no negative scores"));

        // All-negative: no way to score a match.
        let mut m = ScoreMatrix::dna_identity();
        for s in &mut m.scores {
            if *s > 0 {
                *s = -1;
            }
        }
        let err = m.validate().unwrap_err();
        assert!(format!("{err}").contains("no positive scores"));
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
