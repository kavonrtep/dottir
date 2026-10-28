# ADR 0009 — BLASTX Karlin estimate over three frames

* Status: Accepted
* Date: 2026-10-28
* Deciders: petr
* Spec reference: §4.1.1 (window size), §4.1.7 (strand and frame handling)

## Context

BLASTX translates the horizontal DNA sequence in its three forward reading
frames and max-merges the three protein-kernel passes (spec §4.1.7). The
automatic window size and pixel factor come from a Karlin/Altschul
estimate. dottir computed it from frame 0 only, with a comment claiming
that matched C dotter. It does not: `dotplot.c:1107` runs
`winsizeFromlambdak` on each frame's translation and uses

* window = `(win1 + win2 + win3) / 3.0` stored in an `int` (truncated),
  each `winN` the unclamped MSP length; the usual clamp applies after;
* `expResScore` = the mean of the three, which sets the auto pixel factor.

With frames of different composition the frame-0 window and pixel factor
differ from C dotter's, so BLASTX pixelmaps with an automatic window or
pixel factor did not match.

## Decision

Compute the estimate per frame and combine it as C dotter does
(`dottir_core::plot::blastx_karlin`). λ, K and H reported for BLASTX are
frame 0's (C prints all three; only window and residue score are used).
Frames too short to translate (queries under 5 bases) are left out of
the mean instead of failing.

`PIXELMAP_FORMAT_VERSION` goes from 0 to 1: BLASTX pixelmaps with an
automatic window or pixel factor change. BLASTN and BLASTP are
unaffected. No golden pixelmaps existed for BLASTX, so none are
regenerated.

## Known remaining differences from C dotter (not changed here)

* Unknown protein characters (and the missing last codon of a short
  frame) score 0 in dottir; C maps them to the `*` row/column of the
  matrix (`NR = 23`). This is a documented choice in `alphabet.rs`.
* A Karlin window below 3 is clamped to 3 (spec §4.1.1); C dotter
  substitutes 10.

Both affect the zero-tolerance golden contract and need a decision of
their own once C dotter goldens can be generated.
