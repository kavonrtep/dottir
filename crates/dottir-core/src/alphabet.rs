//! Residue alphabets and ASCII → numeric encoding tables.
//!
//! Mirrors the C dotter encoding (see `dotterApp/dotplot.c` near
//! `atob_0[]` and `ntob[]`). Two alphabets are supported:
//!
//! * Protein, 24 letters: `ARNDCQEGHILKMFPSTWYVBZX*` (index 0..23).
//! * Nucleotide, 5 letters: `ACGTN` (index 0..4). For Karlin/Altschul
//!   purposes only `ACGT` (indices 0..3) are counted; `N` and any other
//!   character are ignored.
//!
//! The encoded sequence type is `Vec<u8>` of small indices. Out-of-alphabet
//! bytes are encoded as [`Alphabet::na`] (a single sentinel index past the
//! end of the scorable range).

/// Marker for an unknown residue: any input byte not in the alphabet maps
/// here. By construction this is always `alphabet.size()` so the test
/// `tob[byte] < alphabet.scorable_size()` excludes it.
pub const SENTINEL: u8 = u8::MAX;

/// Which alphabet to use when encoding a sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlphabetKind {
    Protein,
    Dna,
}

/// The protein alphabet in canonical C-dotter order. Index = encoded residue,
/// value = ASCII letter.
pub const PROTEIN_LETTERS: &[u8; 24] = b"ARNDCQEGHILKMFPSTWYVBZX*";

/// The nucleotide alphabet. Index = encoded residue, value = ASCII letter.
pub const DNA_LETTERS: &[u8; 5] = b"ACGTN";

/// Alphabet size used when summing Karlin/Altschul residue probabilities.
/// For DNA we count only A/C/G/T (4); the 'N' bucket is *not* included.
/// For protein all 24 buckets are included (matching C dotter abetsize=24).
pub const PROTEIN_KARLIN_SIZE: usize = 24;
pub const DNA_KARLIN_SIZE: usize = 4;

/// Look up the encoded residue for an ASCII byte. Returns [`SENTINEL`] for
/// any unknown character (including whitespace, digits, '-').
#[inline]
pub fn encode_protein(b: u8) -> u8 {
    PROTEIN_TABLE[b as usize]
}

/// Look up the encoded residue for a DNA ASCII byte (case-insensitive).
/// Returns [`SENTINEL`] for unknown bytes; 4 for 'N'/'n'.
#[inline]
pub fn encode_dna(b: u8) -> u8 {
    DNA_TABLE[b as usize]
}

/// Convenience: encode a whole sequence using the chosen alphabet.
pub fn encode(seq: &[u8], kind: AlphabetKind) -> Vec<u8> {
    let lut: &[u8; 256] = match kind {
        AlphabetKind::Protein => &PROTEIN_TABLE,
        AlphabetKind::Dna => &DNA_TABLE,
    };
    seq.iter().map(|&b| lut[b as usize]).collect()
}

/// Build the protein table at compile time. Mirrors `atob_0[]` from
/// `dotplot.c:55`. The C table uses `NR = 23` (= the `*` index) for unknown
/// characters; we use [`SENTINEL`] instead so that callers can filter
/// unscorable residues explicitly. The Karlin port keeps the bucket-23 (`*`)
/// behaviour by encoding `*` to 23 directly (it is part of the alphabet).
const PROTEIN_TABLE: [u8; 256] = build_protein_table();

const fn build_protein_table() -> [u8; 256] {
    let mut t = [SENTINEL; 256];
    let mut i = 0;
    // PROTEIN_LETTERS lists the residues *in encoded order*.
    while i < PROTEIN_LETTERS.len() {
        let letter = PROTEIN_LETTERS[i];
        t[letter as usize] = i as u8;
        // Case insensitivity for A..Z (does nothing for '*').
        if letter >= b'A' && letter <= b'Z' {
            t[(letter | 0x20) as usize] = i as u8;
        }
        i += 1;
    }
    t
}

const DNA_TABLE: [u8; 256] = build_dna_table();

const fn build_dna_table() -> [u8; 256] {
    let mut t = [SENTINEL; 256];
    let mut i = 0;
    while i < DNA_LETTERS.len() {
        let letter = DNA_LETTERS[i];
        t[letter as usize] = i as u8;
        if letter >= b'A' && letter <= b'Z' {
            t[(letter | 0x20) as usize] = i as u8;
        }
        i += 1;
    }
    t
}

/// Reverse complement of a *raw ASCII* DNA byte (case preserved is not
/// required; output is uppercase). Unknown bytes map to `N`.
#[inline]
pub fn complement_dna_byte(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

/// Reverse-complement a DNA sequence.
pub fn reverse_complement_dna(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement_dna_byte(b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protein_table_round_trip() {
        for (idx, &letter) in PROTEIN_LETTERS.iter().enumerate() {
            assert_eq!(encode_protein(letter), idx as u8);
            if letter.is_ascii_uppercase() {
                assert_eq!(encode_protein(letter | 0x20), idx as u8);
            }
        }
        assert_eq!(encode_protein(b' '), SENTINEL);
        assert_eq!(encode_protein(b'-'), SENTINEL);
        assert_eq!(encode_protein(b'\n'), SENTINEL);
    }

    #[test]
    fn dna_table_matches_c_dotter_ntob() {
        // ntob in dotplot.c:103 — A=0, C=1, G=2, T=3, N=4; everything else NN=5.
        assert_eq!(encode_dna(b'A'), 0);
        assert_eq!(encode_dna(b'a'), 0);
        assert_eq!(encode_dna(b'C'), 1);
        assert_eq!(encode_dna(b'G'), 2);
        assert_eq!(encode_dna(b'T'), 3);
        assert_eq!(encode_dna(b'N'), 4);
        assert_eq!(encode_dna(b'X'), SENTINEL);
    }

    #[test]
    fn reverse_complement_basic() {
        assert_eq!(reverse_complement_dna(b"ACGT"), b"ACGT".to_vec());
        assert_eq!(reverse_complement_dna(b"AAAACCCGG"), b"CCGGGTTTT".to_vec());
        assert_eq!(reverse_complement_dna(b"NXa"), b"TNN".to_vec());
    }

    #[test]
    fn encode_filters_whitespace() {
        let enc = encode(b"AC GT\nN", AlphabetKind::Dna);
        assert_eq!(enc, vec![0, 1, SENTINEL, 2, 3, SENTINEL, 4]);
    }
}
