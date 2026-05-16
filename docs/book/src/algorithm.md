# Algorithm overview

The dottir kernel is faithful to the C dotter
(`dotterApp/dotplot.c:1308` onward). The high-level shape:

1. **Encode** the query and subject into alphabet indices using the
   tables in [`dottir_core::alphabet`] (the same ASCII→index mapping
   the C `atob_0[]` and `ntob[]` tables use).
2. **Build a flat score vector** `scoreVec[row * qlen + col]` where
   `row` is an encoded subject residue and `col` is a query position.
   This is the cache-friendly internal layout called out in spec
   §4.1.3 as a permitted deviation from the C row-of-pointers.
3. **Walk the subject axis** in either Forward or Reverse direction.
   For each subject position `s`, ping-pong two row buffers
   (`sum1`, `sum2`) and compute the recurrence
   ```text
   newsum[q] = oldsum[q-1] + scoreVec[s_idx][q] - scoreVec[s_idx - W][q - W]
   ```
   with a `delrow = zeros` warm-up for the first W subject positions
   (gated on the chunk's `iter_idx` so parallel chunking preserves
   determinism — see [crate layout](./crates.md)).
4. **Emit a pixel** at the integer-divided coordinates
   `(q/zoom, s/zoom)` if and only if `newsum > 0`, the absolute
   coordinate is in the valid range, **and** the anti-diagonal
   suppression rule passes:
   ```text
   forward: s_local >= q_local
   reverse: (zoom - 1 - s_local) >= q_local
   ```
   The pixel value is `min(255, newsum * pixel_fac / W)`, max-merged
   into the pixelmap.

## Karlin/Altschul window-size estimation

Given the score matrix and the residue composition of both inputs,
the kernel computes λ, K, and H per Karlin & Altschul 1990, then
derives an expected MSP length for a nominal 100×100 matrix, rounded
to integer and clamped to `[3, 50]`. The `-W` flag bypasses this
estimate.

The implementation in `dottir-core::karlin` is a structural port of
`dotterApp/dotterKarlin.c:144`. Its λ/K/H output is bit-identical to
the C reference on the four pinned fixtures in
`tests/golden/karlin/values.tsv`.

## Reverse-strand pass (BLASTN)

For dual-strand BLASTN, dottir runs a second pass against a score
vector built from the **complement** of the query (mirroring C's
`ntob_compl[]` table), with the subject iterated backwards. A hit
at `(q, s)` then represents `subject[s..s+W]` (read backwards) versus
`complement(query)[q..q+W]`, which is the BLASTN reverse-strand match
semantics.

## Self-comparison

When `--self-comparison` is set, the inner loop caps `qmax = s + 1`
so only the lower triangle is filled. The mirror step then
populates the other half per the `Triangle` mode (`Both` mirrors
symmetrically, `Upper` mirrors then zeros the lower, `Lower` leaves
the lower as computed). `--disable-mirror` short-circuits the
post-process.

## Determinism

Per spec §4.1.11 dottir produces byte-identical pixelmaps across runs
and across thread counts. This is verified by
`crates/dottir-core/tests/parallel_determinism.rs` at `n_threads = 1,
2, 4, 8`.
