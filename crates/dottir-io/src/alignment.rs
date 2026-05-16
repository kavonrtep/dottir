//! Alignment-export helper (spec §4.4.5).
//!
//! Slices a window of `±N` residues around a (q, s) coordinate and emits
//! the pair in one of three formats: a FASTA pair, a Stockholm 1-row
//! block, or plain text. No alignment is done here — that's deferred to
//! Phase 6-extra (which will optionally shell out to mafft/muscle, or
//! call into the `bio` crate on windows ≤ 1 kb).
//!
//! For the dottir use case (interactive crosshair → "show me what's
//! aligning here") the ungapped slice is usually all the user needs:
//! the dotplot already certifies the diagonal alignment, and the user
//! visually inspects the match in the surrounding context.

use std::fmt::Write as _;

/// Output format for [`slice_pair`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SliceFormat {
    /// Two FASTA records, one per sequence.
    FastaPair,
    /// Single Stockholm block with two `# STOCKHOLM 1.0` lines.
    Stockholm,
    /// Plain `q: ...\ns: ...\n` text — good for clipboards.
    Plain,
}

/// Slice ±`window` residues around `(q_idx, s_idx)` from `q_full` and
/// `s_full`, clamped to sequence boundaries, and serialise to `format`.
///
/// The slice always pairs the *same* number of residues from each
/// sequence: the slice length is
/// `min(window*2 + 1, qlen - q_lo, slen - s_lo)`, where
/// `lo = max(0, idx - window)`. The caller can recover the alignment of
/// the slice in the original sequence coordinates from `q_lo` / `s_lo`.
pub fn slice_pair(
    q_full: &[u8],
    s_full: &[u8],
    q_idx: usize,
    s_idx: usize,
    window: usize,
    format: SliceFormat,
    q_name: &str,
    s_name: &str,
) -> SlicedAlignment {
    let q_lo = q_idx.saturating_sub(window);
    let s_lo = s_idx.saturating_sub(window);
    let q_hi = (q_idx + window + 1).min(q_full.len());
    let s_hi = (s_idx + window + 1).min(s_full.len());
    let pair_len = (q_hi - q_lo).min(s_hi - s_lo);
    let q = &q_full[q_lo..q_lo + pair_len];
    let s = &s_full[s_lo..s_lo + pair_len];

    let text = match format {
        SliceFormat::FastaPair => format_fasta_pair(q_name, q, s_name, s, q_lo, s_lo),
        SliceFormat::Stockholm => format_stockholm(q_name, q, s_name, s),
        SliceFormat::Plain => format_plain(q_name, q, s_name, s),
    };

    SlicedAlignment {
        q_range: q_lo..q_lo + pair_len,
        s_range: s_lo..s_lo + pair_len,
        text,
    }
}

/// Result of [`slice_pair`].
#[derive(Debug, Clone)]
pub struct SlicedAlignment {
    pub q_range: std::ops::Range<usize>,
    pub s_range: std::ops::Range<usize>,
    pub text: String,
}

fn format_fasta_pair(q_name: &str, q: &[u8], s_name: &str, s: &[u8], q_lo: usize, s_lo: usize) -> String {
    let mut out = String::new();
    let _ = writeln!(out, ">{q_name}:{}..{} (length {})", q_lo + 1, q_lo + q.len(), q.len());
    out.push_str(&wrap_fasta_seq(q, 60));
    let _ = writeln!(out, ">{s_name}:{}..{} (length {})", s_lo + 1, s_lo + s.len(), s.len());
    out.push_str(&wrap_fasta_seq(s, 60));
    out
}

fn wrap_fasta_seq(seq: &[u8], width: usize) -> String {
    let mut out = String::with_capacity(seq.len() + seq.len() / width + 1);
    for chunk in seq.chunks(width) {
        out.push_str(&String::from_utf8_lossy(chunk));
        out.push('\n');
    }
    out
}

fn format_stockholm(q_name: &str, q: &[u8], s_name: &str, s: &[u8]) -> String {
    let mut out = String::new();
    out.push_str("# STOCKHOLM 1.0\n");
    let name_w = q_name.len().max(s_name.len()) + 2;
    let _ = writeln!(
        out,
        "{:width$}{}",
        q_name,
        String::from_utf8_lossy(q),
        width = name_w
    );
    let _ = writeln!(
        out,
        "{:width$}{}",
        s_name,
        String::from_utf8_lossy(s),
        width = name_w
    );
    out.push_str("//\n");
    out
}

fn format_plain(q_name: &str, q: &[u8], s_name: &str, s: &[u8]) -> String {
    let name_w = q_name.len().max(s_name.len()) + 2;
    let mut out = String::new();
    let _ = writeln!(
        out,
        "{:width$}{}",
        format!("{q_name}:"),
        String::from_utf8_lossy(q),
        width = name_w
    );
    let _ = writeln!(
        out,
        "{:width$}{}",
        format!("{s_name}:"),
        String::from_utf8_lossy(s),
        width = name_w
    );
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slice_centered_inside_sequences() {
        let q = b"AAAAACCCCCGGGGGTTTTT".to_vec();
        let s = b"NNNNNAAAAACCCCCGGGGG".to_vec();
        let r = slice_pair(&q, &s, 10, 10, 3, SliceFormat::Plain, "q", "s");
        // window of 3 → 7 residues centred at idx 10 on both.
        // q[7..14] = "CCCGGGG", s[7..14] = "AAACCCC"
        assert_eq!(r.q_range, 7..14);
        assert_eq!(r.s_range, 7..14);
        assert!(r.text.contains("CCCGGGG"));
        assert!(r.text.contains("AAACCCC"));
    }

    #[test]
    fn slice_near_left_boundary_clamps() {
        let q = b"AAAAACCCCC".to_vec();
        let s = b"GGGGGAAAAA".to_vec();
        let r = slice_pair(&q, &s, 0, 0, 3, SliceFormat::Plain, "q", "s");
        // q[0..4], s[0..4]
        assert_eq!(r.q_range, 0..4);
        assert_eq!(r.s_range, 0..4);
        assert!(r.text.contains("AAAA"));
        assert!(r.text.contains("GGGG"));
    }

    #[test]
    fn slice_near_right_boundary_clamps_to_shorter() {
        let q = b"AAAA".to_vec();
        let s = b"GGGGGGGGGGG".to_vec();
        let r = slice_pair(&q, &s, 3, 8, 5, SliceFormat::Plain, "q", "s");
        // Pair length is min(qlen - q_lo, slen - s_lo) starting at the
        // clamped lo. q_lo = 3 - 5 → 0 (saturating). s_lo = 8 - 5 = 3.
        // q_hi = min(qlen, 3+5+1) = 4. s_hi = min(11, 8+5+1) = 11.
        // pair_len = min(4 - 0, 11 - 3) = 4.
        assert_eq!(r.q_range, 0..4);
        assert_eq!(r.s_range, 3..7);
    }

    #[test]
    fn fasta_pair_format_has_two_records() {
        let q = b"ACGTACGTAC".to_vec();
        let s = b"GTACGTACGT".to_vec();
        let r = slice_pair(&q, &s, 5, 5, 3, SliceFormat::FastaPair, "Q", "S");
        let header_count = r.text.matches('>').count();
        assert_eq!(header_count, 2);
        assert!(r.text.contains(">Q:"));
        assert!(r.text.contains(">S:"));
    }

    #[test]
    fn stockholm_format_has_terminator() {
        let q = b"ACGT".to_vec();
        let s = b"ACGT".to_vec();
        let r = slice_pair(&q, &s, 1, 1, 5, SliceFormat::Stockholm, "alpha", "beta");
        assert!(r.text.starts_with("# STOCKHOLM 1.0"));
        assert!(r.text.contains("//"));
        assert!(r.text.contains("alpha"));
        assert!(r.text.contains("beta"));
    }
}
