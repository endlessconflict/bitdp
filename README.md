# bitdp

Bit-parallel sequence alignment kernels that Zig derives at compile time from the scoring scheme you give it.

Fast bit-vector aligners exist, but each one was worked out by hand for one cost model. Myers did it for unit edit costs in 1999, and BitPAl later built a family of constructions for linear-gap integer weights. However, if your scheme falls outside what someone already derived, you are back to the scalar dynamic program. bitdp takes the recurrence itself and produces the bit-parallel program with `comptime`, so a new scheme costs a recompile instead of a paper.

## How it works

Read the DP matrix one text column at a time, with the pattern rows packed into the bits of a 64-bit word. Down a column, each cell hands the next one its horizontal score difference. When differences are bounded, that makes the column a finite-state machine whose state is the difference, and for min-plus recurrences every step of that machine is monotone in its state.

Monotonicity is the key fact. Write the state in thermometer code (one bit per threshold, `[d >= t]`). Each threshold bit of the next state is then either a constant or a copy of one threshold bit of the previous state. A threshold that copies itself behaves like a carry chain with generate, propagate and kill positions, and a carry chain across a machine word is exactly what integer addition computes. A threshold that copies a different threshold is a shift and some bitwise logic. When the copies between different thresholds form no cycle, the whole column becomes a short cascade of additions and boolean operations. Every boolean function in that cascade is synthesized from its truth table at compile time.

Myers' algorithm comes out of this as the special case with one carry chain. Nothing about it is hard-coded here.

Wide score ranges need one more idea. Most planes ask whether a difference of two thermometer-coded numbers clears a threshold, and on thermometer codes, addition is merging, since the merged sequence of two sorted bit strings is the thermometer code of their sum. A Batcher odd-even merging network (one OR and one AND per comparator) therefore delivers every threshold at once. The algebra also fixes a single cut-off, mismatch cost minus gap cost, above which no level needs a carry chain. Levels below it keep their chains and the rest of the column comes out of two merges. For small score ranges the compiler also tries exhaustive truth-table synthesis, and it keeps whichever program is shorter.

## Status

Early research code. What it does today, measured on one machine:

| Scheme (match, mismatch, gap) | Distinct differences | Additions | Word ops per column | Check |
|---|---|---|---|---|
| 0, 1, 1 (edit distance) | 3 | 1 | 16 | 10^6 random pairs, bit-exact |
| 0, 2, 1 (indel distance) | 2 | 1 | 8 | 2 x 10^5 pairs, bit-exact |
| 0, 1, 2 | 5 | 1 | 30 | 2 x 10^5 pairs, bit-exact |
| 0, 3, 2 | 5 | 3 | 63 | 2 x 10^5 pairs, bit-exact |
| -2, 3, 5 | 13 | 5 | 190 | 10^5 pairs, bit-exact |
| -3, 4, 6 | 16 | 7 | 267 | 10^5 pairs, bit-exact |
| -4, 5, 9 | 23 | 9 | 414 | 10^5 pairs, bit-exact |
| -4, 7, 11 | 27 | 11 | 528 | 10^5 pairs, bit-exact |

The last four rows are the weight sets that BitPAl (Loving, Hernandez and Benson, 2014) benchmarks, written as costs. A score scheme (M, I, G) becomes costs (-M, -I, -G), which gives the same optimal alignments. For (2, -3, -5), the one set whose operation counts the BitPAl paper reports, BitPAl needs 265 operations per 64-bit word and its packed variant 166. The derived kernel needs 190, fewer than the first and still more than the second. Operation counts here include every AND, OR, XOR, NOT, shift and addition, and a + b + 1 counts as two.

For comparison, our own implementation of Myers' hand-derived kernel uses 15 operations per column. On a 64-base pattern against a 50 Mbp random text, the derived edit-distance kernel ran at 15.6 GCUPS and our Myers implementation at 12.7 GCUPS. That is a single run on one CPU with both kernels written by us, so read it as a sanity check rather than a benchmark.

Current limits:

- Global alignment with linear gap costs only.
- Pattern length up to 64, one machine word.
- Two substitution classes (match and mismatch), at most 32 distinct score differences.
- The carry-chain part still grows quadratically with the gap between match and mismatch scores, which is where BitPAl's packed variant stays ahead.

## Usage

Requires Zig 0.16.0.

```zig
const bitdp = @import("bitdp");

const Edit = bitdp.Kernel(.{ .match = 0, .mismatch = 1, .gap = 1 });
const d = Edit.distance("ACGTTGCA", "ACGTGCA"); // 1
// Edit.ops is the derived number of word operations per column.
```

```sh
zig build test
zig build mve -Doptimize=ReleaseFast   # prints derived programs, verifies, times
```

## Next steps

Getting the carry-chain block below quadratic, affine gaps, local alignment, full substitution matrices, and patterns longer than one word.

## References

See [REFERENCES.md](REFERENCES.md).

## License

MIT
