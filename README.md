# bitdp

Bit-parallel sequence alignment kernels that Zig derives at compile time from the scoring scheme you give it.

Fast bit-vector aligners exist, but each one was worked out by hand for one cost model. Myers did it for unit edit costs in 1999, and BitPAl later built a family of constructions for linear-gap integer weights. However, if your scheme falls outside what someone already derived, you are back to the scalar dynamic program. bitdp takes the recurrence itself and produces the bit-parallel program with `comptime`, so a new scheme costs a recompile instead of a paper.

## How it works

Read the DP matrix one text column at a time, with the pattern rows packed into the bits of a 64-bit word. Down a column, each cell hands the next one its horizontal score difference. When differences are bounded, that makes the column a finite-state machine whose state is the difference, and for min-plus recurrences every step of that machine is monotone in its state.

Monotonicity is the key fact. Write the state in thermometer code (one bit per threshold, `[d >= t]`). Each threshold bit of the next state is then either a constant or a copy of one threshold bit of the previous state. A threshold that copies itself behaves like a carry chain with generate, propagate and kill positions, and a carry chain across a machine word is exactly what integer addition computes. A threshold that copies a different threshold is a shift and some bitwise logic. When the copies between different thresholds form no cycle, which holds for every linear-gap scheme, the whole column becomes a short cascade of additions and boolean operations.

Myers' algorithm comes out of this as the special case with one carry chain. Nothing about it is hard-coded here, and the derived kernel uses 14 word operations per column where Myers' uses 15.

Wide score ranges need one more idea. Most planes ask whether a difference of two thermometer-coded numbers clears a threshold, and on thermometer codes, addition is merging, since the merged sequence of two sorted bit strings is the thermometer code of their sum. A Batcher odd-even merging network (one OR and one AND per comparator) therefore delivers every threshold at once. The algebra also fixes a single cut-off, the highest substitution cost minus the gap cost, above which no level needs a carry chain. Levels below it keep their chains and the rest of the column comes out of two merges.

For small score ranges the compiler also tries exhaustive truth-table synthesis and keeps whichever program is shorter.

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

All numbers are from one core of a Ryzen 7 8840U, with every tool's total cost checked against the others; [bench/RESULTS.md](bench/RESULTS.md) has the full tables and [bench/README.md](bench/README.md) the way to reproduce them.

The closest tool is BGSA, which runs Myers' and BitPAl's hand-derived kernels with one alignment per SIMD lane. In its own setting (all queries against all subjects of equal length) and at the same vector width, the derived kernels land close to it. For edit distance, bitdp reaches 0.73 to 0.93 of BGSA's speed with AVX2 and 0.90 to 1.20 with AVX-512. On BitPAl's (2, -3, -5) weights with AVX2, bitdp ties or beats BGSA's standard BitPAl (8.6 against 8.6 and 11.7 against 9.2 GCUPS) and trails its packed variant (8.6 against 13.0 and 11.7 against 14.0). Schemes BGSA cannot express still run at the same order of speed, for example 45.8 GCUPS with AVX2 for transition/transversion costs on 1 kbp sequences.

On independent pairs, where each lane has its own text, batched bitdp is faster than edlib and than parasail's correct kernels on every DNA workload measured (for instance 35.2 against 21.7 and 5.2 GCUPS for edit distance on 1 kbp pairs with AVX2). For BLOSUM62 it is several times slower than parasail, since fifteen cost classes make the column program too long to pay off.

parasail's striped global kernels, which would otherwise be its fastest, scored some pairs below their optimum (426 of 100 000 on one workload), so the comparison leaves them out.

## Limits

- Global alignment with linear gap costs only. With affine gaps the value carried down a column can shift both up and down between thresholds, the copy graph acquires cycles, and the cascade above no longer applies.
- At most 32 distinct score differences and 24 distinct substitution costs.
- The carry-chain part still grows quadratically with the spread between substitution costs. This is where BitPAl's packed variant stays ahead, and why protein matrices are slow.
- Bytes outside the scheme's alphabet (N, for instance) count as the costliest substitution.
- Deriving a scheme happens inside the Zig compiler. Small schemes take seconds, while BLOSUM62 takes about 20 seconds and 1.5 GB of compiler memory.

## Usage

Requires Zig 0.16.0.

```zig
const bitdp = @import("bitdp");

const Edit = bitdp.Kernel(.{ .match = 0, .mismatch = 1, .gap = 1 });

// Patterns up to 62 characters, no allocation.
const d = Edit.distance("ACGTTGCA", "ACGTGCA"); // 1

// Any pattern length: prepare it once, align it against many texts.
var a = try Edit.Aligner.init(gpa, long_pattern);
defer a.deinit(gpa);
const d2 = a.distance(text);

// Many pairs, one per SIMD lane.
try Edit.distances(gpa, patterns, texts, out);

// Up to Edit.lanes patterns against the same text, one per lane.
var g = try Edit.Group.init(gpa, subjects);
defer g.deinit(gpa);
const costs = g.distances(query);

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
