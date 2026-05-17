//! DNA → protein translation for BLASTX (Phase F1).
//!
//! Standard NCBI genetic code (translation table 1). Stop codons are
//! translated to `*`, which dottir's protein alphabet handles as
//! index 23. Codons containing any non-ACGT base translate to `X` so
//! they score against the matrix's ambiguity column rather than
//! breaking the kernel.

/// Translate a DNA codon (three uppercased ASCII bases) to its
/// single-letter amino-acid code via the standard NCBI genetic code.
/// Any codon containing a base outside ACGT returns `b'X'`; stop
/// codons return `b'*'`.
#[inline]
pub fn translate_codon(c: [u8; 3]) -> u8 {
    let i = match codon_index(c) {
        Some(i) => i,
        None => return b'X',
    };
    NCBI_GENETIC_CODE_1[i]
}

/// Map a codon to an index `0..64` for table lookup. Returns `None`
/// if any base isn't ACGT.
#[inline]
fn codon_index(c: [u8; 3]) -> Option<usize> {
    let a = base_index(c[0])?;
    let b = base_index(c[1])?;
    let z = base_index(c[2])?;
    // Standard NCBI ordering: T C A G — base 0 is T.
    Some(a * 16 + b * 4 + z)
}

#[inline]
fn base_index(b: u8) -> Option<usize> {
    match b {
        b'T' | b't' | b'U' | b'u' => Some(0),
        b'C' | b'c' => Some(1),
        b'A' | b'a' => Some(2),
        b'G' | b'g' => Some(3),
        _ => None,
    }
}

/// Standard NCBI translation table 1, indexed by `T·16 + C·4 + A`.
/// Source: NCBI Taxonomy translation tables, ID 1 ("Standard Code").
const NCBI_GENETIC_CODE_1: [u8; 64] = *b"\
FFLLSSSSYY**CC*W\
LLLLPPPPHHQQRRRR\
IIIMTTTTNNKKSSRR\
VVVVAAAADDEEGGGG";

/// Translate a full DNA sequence in one reading frame. `frame_offset`
/// is `0`, `1`, or `2`; output length is `(seq.len() - frame_offset) / 3`.
/// Trailing partial codons are dropped.
pub fn translate_frame(seq: &[u8], frame_offset: usize) -> Vec<u8> {
    if frame_offset >= 3 || seq.len() < frame_offset + 3 {
        return Vec::new();
    }
    let tail = &seq[frame_offset..];
    let n_codons = tail.len() / 3;
    let mut out = Vec::with_capacity(n_codons);
    for i in 0..n_codons {
        let off = i * 3;
        let codon = [tail[off], tail[off + 1], tail[off + 2]];
        out.push(translate_codon(codon));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Canonical codon spot-checks against NCBI table 1.
    #[test]
    fn known_codons() {
        // Start codon.
        assert_eq!(translate_codon(*b"ATG"), b'M');
        // Stop codons.
        assert_eq!(translate_codon(*b"TAA"), b'*');
        assert_eq!(translate_codon(*b"TAG"), b'*');
        assert_eq!(translate_codon(*b"TGA"), b'*');
        // Sample amino acids.
        assert_eq!(translate_codon(*b"GAA"), b'E');
        assert_eq!(translate_codon(*b"GGG"), b'G');
        assert_eq!(translate_codon(*b"TTT"), b'F');
        assert_eq!(translate_codon(*b"CGT"), b'R');
        assert_eq!(translate_codon(*b"AGA"), b'R');
        assert_eq!(translate_codon(*b"TGG"), b'W');
        // Case insensitivity.
        assert_eq!(translate_codon(*b"atg"), b'M');
    }

    #[test]
    fn unknown_bases_translate_to_x() {
        assert_eq!(translate_codon(*b"NTG"), b'X');
        assert_eq!(translate_codon(*b"ANG"), b'X');
        assert_eq!(translate_codon(*b"AT-"), b'X');
    }

    #[test]
    fn translate_frame_drops_partial_codon() {
        // 7 bases → 2 codons in frame 0.
        let p = translate_frame(b"ATGGCATG", 0);
        assert_eq!(p, b"MA"); // ATG | GCA | "TG" dropped
                              // Frame 1: shift by 1, 6 usable bases → 2 codons.
        let p = translate_frame(b"ATGGCATG", 1);
        assert_eq!(p, b"WH"); // TGG | CAT, last base dropped
    }

    #[test]
    fn translate_frame_empty_for_too_short() {
        assert!(translate_frame(b"", 0).is_empty());
        assert!(translate_frame(b"AT", 0).is_empty());
        assert!(translate_frame(b"ATG", 1).is_empty());
        assert!(translate_frame(b"ATGG", 2).is_empty());
        assert!(translate_frame(b"ATG", 3).is_empty());
    }
}
