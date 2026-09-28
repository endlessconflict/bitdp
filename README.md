# bitdp

Bit-parallel sequence alignment kernels that Zig derives at compile time from the scoring scheme you give it.

Fast bit-vector aligners exist, but each one was worked out by hand for one cost model. Myers did it for unit edit costs in 1999, and BitPAl later built a family of constructions for linear-gap integer weights. However, if your scheme falls outside what someone already derived, you are back to the scalar dynamic program. bitdp takes the recurrence itself and produces the bit-parallel program with `comptime`, so a new scheme costs a recompile instead of a paper.

## How it works

Read the DP matrix one text column at a time, with the pattern rows packed into the bits of a 64-bit word. Down a column, each cell hands the next one its horizontal score difference. When differences are bounded, that makes the column a finite-state machine whose state is the difference, and for min-plus recurrences every step of that machine is monotone in its state.

Monotonicity is the key fact. Write the state in thermometer code (one bit per threshold, `[d >= t]`). Each threshold bit of the next state is then either a constant or a copy of one threshold bit of the previous state. A threshold that copies itself behaves like a carry chain with generate, propagate and kill positions, and a carry chain across a machine word is exactly what integer addition computes. A threshold that copies a different threshold is a shift and some bitwise logic. When the copies between different thresholds form no cycle, which holds for every linear-gap scheme, the whole column becomes a short cascade of additions and boolean operations.

Myers' algorithm comes out of this as the special case with one carry chain. Nothing about it is hard-coded here, and the derived kernel uses 14 word operations per column where Myers' uses 15.

Wide score ranges need one more idea. Most planes ask whether a difference of two thermometer-coded numbers clears a threshold, and on thermometer codes, addition is merging, since the merged sequence of two sorted bit strings is the thermometer code of their sum. A Batcher odd-even merging network (one OR and one AND per comparator) therefore delivers every threshold at once. The algebra also fixes a single cut-off, the highest substitution cost minus the gap cost, above which no level needs a carry chain. Levels below it keep their chains and the rest of the column comes out of two merges.

Two passes finish the program. For small score ranges the compiler also tries exhaustive truth-table synthesis and keeps the shorter result. Then it removes NOT gates by choosing, for every AND and OR and for every input plane, whether to build the value or its complement, since De Morgan makes both equally cheap and the kernel can store its inputs either way.

In addition, any integer substitution matrix works, not only match and mismatch, because the cost itself becomes one more thermometer-coded input and the same rules apply.

## Status

Research code, not yet stable. Word operations per column, as derived:

| Scheme (match, mismatch, gap) | Distinct differences | Additions | Word ops | Check |
|---|---|---|---|---|
| 0, 1, 1 (edit distance) | 3 | 1 | 14 | 10^6 random pairs, bit-exact |
| 0, 2, 1 (indel distance) | 2 | 1 | 6 | 2 x 10^5 pairs, bit-exact |
| 0, 1, 2 | 5 | 1 | 28 | 2 x 10^5 pairs, bit-exact |
| 0, 3, 2 | 5 | 3 | 57 | 2 x 10^5 pairs, bit-exact |
| -2, 3, 5 | 13 | 5 | 181 | 10^5 pairs, bit-exact |
| -3, 4, 6 | 16 | 7 | 253 | 10^5 pairs, bit-exact |
| -4, 5, 9 | 23 | 9 | 399 | 10^5 pairs, bit-exact |
| -4, 7, 11 | 27 | 11 | 509 | 10^5 pairs, bit-exact |
| DNA: transition 1, transversion 2, gap 2 | 5 | 2 | 51 | 10^5 pairs, bit-exact |
| BLOSUM62 costs, gap 4 | 20 | 15 | 799 | 10^5 protein pairs, bit-exact |

Rows five to eight are the weight sets that BitPAl (Loving, Hernandez and Benson, 2014) benchmarks, written as costs. A score scheme (M, I, G) becomes costs (-M, -I, -G), which gives the same optimal alignments. For (2, -3, -5), the one set whose operation counts the BitPAl paper reports, BitPAl needs 265 operations per 64-bit word and its packed variant 166. The derived kernel needs 181, fewer than the first and still more than the second. Operation counts include every AND, OR, XOR, NOT, shift and addition, and a + b + 1 counts as two.

## Speed

Measured against edlib and parasail on one core of a Ryzen 7 8840U, all tools on the same inputs, with every tool's total cost checked against the others. Billions of DP cells per second (median of three runs):

| Workload | Scheme | bitdp, 1 pair | bitdp, 4 lanes AVX2 | bitdp, 8 lanes AVX-512 | edlib | parasail (fastest correct kernel) |
|---|---|---|---|---|---|---|
| 150 bp, 10 % divergence | edit | 7.43 | 10.47 | 9.67 | 6.19 | 2.16 |
| | bitpal | 0.91 | 4.35 | 5.51 | | 2.18 |
| | tstv | 2.34 | 3.59 | 3.69 | | 2.21 |
| 1 kbp, 15 % divergence | edit | 18.23 | 30.70 | 50.44 | 21.73 | 5.19 |
| | bitpal | 1.16 | 7.29 | 11.59 | | 5.18 |
| | tstv | 5.52 | 14.04 | 16.82 | | 5.18 |
| 300 aa, 30 % divergence | blosum | 0.33 | 0.59 | 0.63 | | 2.87 |

The lane columns run one alignment per SIMD lane through `Kernel.distances`. The Ubuntu parasail build has no AVX-512 kernels, so the AVX2 column is the like-for-like comparison. On DNA, batched bitdp is faster than every correct baseline in the table. On BLOSUM62 it is several times slower than parasail: fifteen cost classes make the column program too long to pay off.

parasail's striped global kernels, which would otherwise be its fastest, returned a different total cost than every other tool on every workload, so they are left out of the comparison. [bench/RESULTS.md](bench/RESULTS.md) has all numbers, including them, and [bench/README.md](bench/README.md) explains how to reproduce the run.

## Limits

- Global alignment with linear gap costs only. With affine gaps the value carried down a column can shift both up and down between thresholds, the copy graph acquires cycles, and the cascade above no longer applies.
- At most 32 distinct score differences and 24 distinct substitution costs.
- The carry-chain part still grows quadratically with the spread between substitution costs. This is where BitPAl's packed variant stays ahead, and why protein matrices are slow.
- Deriving a scheme happens inside the Zig compiler. Small schemes take seconds, while BLOSUM62 takes about 20 seconds and 1.5 GB of compiler memory.

## Usage

Requires Zig 0.16.0.

```zig
const bitdp = @import("bitdp");

const Edit = bitdp.Kernel(.{ .match = 0, .mismatch = 1, .gap = 1 });

// Patterns up to 63 characters, no allocation.
const d = Edit.distance("ACGTTGCA", "ACGTGCA"); // 1

// Any pattern length: prepare it once, align it against many texts.
var a = try Edit.Aligner.init(gpa, long_pattern);
defer a.deinit(gpa);
const d2 = a.distance(text);

// Many pairs, one per SIMD lane.
try Edit.distances(gpa, patterns, texts, out);

// Any substitution cost function over a declared alphabet.
const TsTv = bitdp.Kernel(bitdp.schemes.ts_tv);
const Blosum = bitdp.Kernel(bitdp.schemes.blosum62Linear(4));
```

```sh
zig build test
zig build mve -Doptimize=ReleaseFast   # prints derived programs, verifies, times
```

## References

See [REFERENCES.md](REFERENCES.md).

## License

MIT
