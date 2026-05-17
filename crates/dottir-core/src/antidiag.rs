//! Anti-diagonal suppression rule (spec §4.1.6, C source
//! `dotterApp/dotplot.c:1405`).
//!
//! Each output pixel covers a `zoom × zoom` block of the underlying score
//! matrix. Multiple diagonals project into the same pixel; without
//! suppression the dotplot acquires diagonal-direction noise that obscures
//! true alignments. The rule keeps only sub-pixel positions whose local
//! coordinates are consistent with the diagonal we are scanning:
//!
//! * Forward strand: `s_local >= q_local`.
//! * Reverse strand: `(zoom - 1 - s_local) >= q_local`.
//!
//! This rule is load-bearing: per spec §4.1.6 ("Add an explicit unit test")
//! we keep [`keep_dot`] testable in isolation. The Phase 1 inner loop will
//! call this for every candidate pixel before the max-merge step.

/// True iff a candidate dot at sub-pixel local position `(q_local, s_local)`
/// should be merged into the output pixel for the given strand.
///
/// `zoom` is the linear pixel size (1 means no zoom; a single block = one
/// matrix cell, so the suppression rule degenerates to `s >= q`).
#[inline]
pub fn keep_dot(zoom: u32, q_local: u32, s_local: u32, reverse: bool) -> bool {
    debug_assert!(q_local < zoom, "q_local {q_local} >= zoom {zoom}");
    debug_assert!(s_local < zoom, "s_local {s_local} >= zoom {zoom}");
    let s = if reverse { zoom - 1 - s_local } else { s_local };
    s >= q_local
}

#[cfg(test)]
mod tests {
    use super::keep_dot;

    /// At zoom = 1 the only sub-pixel is (0, 0); both strands keep it.
    #[test]
    fn zoom_one_always_keeps() {
        assert!(keep_dot(1, 0, 0, false));
        assert!(keep_dot(1, 0, 0, true));
    }

    /// Forward rule keeps the upper-right triangle of the sub-pixel
    /// (s_local >= q_local), including the diagonal.
    #[test]
    fn forward_keeps_upper_right_triangle() {
        let zoom = 4;
        // Expected mask: true iff s_local >= q_local.
        let expected: [[bool; 4]; 4] = [
            // s = 0    1     2     3
            [true, true, true, true],    // q = 0
            [false, true, true, true],   // q = 1
            [false, false, true, true],  // q = 2
            [false, false, false, true], // q = 3
        ];
        for q in 0..zoom {
            for s in 0..zoom {
                let got = keep_dot(zoom, q, s, false);
                let want = expected[q as usize][s as usize];
                assert_eq!(got, want, "q={q} s={s}");
            }
        }
    }

    /// Reverse rule keeps the lower-right triangle of the sub-pixel:
    /// `(zoom - 1 - s_local) >= q_local`  ⇔  `s_local <= zoom - 1 - q_local`.
    /// Note this is the C-dotter behaviour: the "axis flip" is on s, not on q.
    #[test]
    fn reverse_keeps_anti_diagonal_triangle() {
        let zoom = 4;
        let expected: [[bool; 4]; 4] = [
            // s = 0    1     2     3
            [true, true, true, true],    // q = 0
            [true, true, true, false],   // q = 1
            [true, true, false, false],  // q = 2
            [true, false, false, false], // q = 3
        ];
        for q in 0..zoom {
            for s in 0..zoom {
                let got = keep_dot(zoom, q, s, true);
                let want = expected[q as usize][s as usize];
                assert_eq!(got, want, "q={q} s={s}");
            }
        }
    }

    /// Forward and reverse rules pick disjoint *strict* triangles (off the
    /// boundary): if `s > q` *and* `s < zoom - 1 - q`, both keep; otherwise
    /// at most one keeps. Above all, the rules are pure functions of the
    /// (q_local, s_local, zoom) triple, so behaviour is fully described by
    /// [`forward_keeps_upper_right_triangle`] and
    /// [`reverse_keeps_anti_diagonal_triangle`].
    #[test]
    fn rules_are_pure_functions_of_inputs() {
        // 100 calls with the same input return the same answer.
        for _ in 0..100 {
            for zoom in [2_u32, 4, 8, 16] {
                for q in 0..zoom {
                    for s in 0..zoom {
                        assert_eq!(keep_dot(zoom, q, s, false), s >= q);
                        assert_eq!(keep_dot(zoom, q, s, true), zoom - 1 - s >= q);
                    }
                }
            }
        }
    }
}
