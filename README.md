# bitdp

Bit-parallel sequence alignment kernels that Zig derives at compile time from the scoring scheme you give it.

Fast, implemented bit-vector kernels exist for a few cost models, each worked out by hand: Myers' for unit edit costs (1999), kernels for indel distance via the longest common subsequence, and BitPAl's construction for linear-gap match, mismatch and gap weights, which its generator instantiates per weight set. Bergeron and Hamel showed in 2002 and 2004 that the DP column of weighted edit distance with a single indel cost is an automaton that decomposes into threshold levels, each computed with one addition, but they gave no implementation or word-level operation count. bitdp builds that cascade at the bit level, without arithmetic on vectors of multi-bit integers, and derives the program from the scoring scheme with `comptime`, so a new linear-gap scheme costs a recompile.

## How it works

Read the DP matrix one text column at a time, with the pattern rows packed into the bits of a 64-bit word. Down a column, each cell hands the next one its horizontal score difference. When differences are bounded, that makes the column a finite-state machine whose state is the difference (Bergeron and Hamel's automaton), and for the linear-gap recurrence every step of that machine is monotone in its state.

Monotonicity is the key fact. Write the state in thermometer code (one bit per threshold, `[d >= t]`). Each threshold bit of the next state is then either a constant or a copy of one threshold bit of the previous state. A threshold that copies itself behaves like a carry chain with generate, propagate and kill positions, and a carry chain across a machine word is exactly what integer addition computes. A threshold that copies a different threshold is a shift and some bitwise logic. When the copies between different thresholds form no cycle, which holds for every linear-gap scheme, the whole column becomes a short cascade of additions and boolean operations: Bergeron and Hamel's cascade decomposition, here written out bit by bit.

Myers' algorithm comes out of this as the special case with one carry chain. Nothing about it is hard-coded here. In global mode the derived kernel uses 14 word operations per column where Myers' kernel, as written in `reference.zig`, uses 16, and in search mode both use 15. The 14 relies on the truth-table synthesis described below; the cascade alone gives 19.

Wide score ranges need one more idea. Most planes ask whether a difference of two thermometer-coded numbers clears a threshold, and on thermometer codes, addition is merging, since the merged sequence of two sorted bit strings is the thermometer code of their sum. A Batcher odd-even merging network (one OR and one AND per comparator) therefore delivers every threshold at once. The algebra also fixes a single cut-off, the highest substitution cost minus the gap cost, above which no level needs a carry chain. Levels below it keep their chains and the rest of the column comes out of two merges.

For two cost classes and at most five distinct differences, the compiler also synthesizes the program from its truth table (prime implicants with a greedy cover) and keeps whichever program is shorter.

In addition, integer substitution matrices work (up to the limits below), not only match and mismatch, because the cost itself becomes one more thermometer-coded input and the same rules apply.

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

Rows five to eight are the weight sets that BitPAl (Loving, Hernandez and Benson, 2014) benchmarks, written as costs. A score scheme (M, I, G) becomes costs (-M, -I, -G), which gives the same optimal alignments. For (2, -3, -5), the set for which the BitPAl paper states both counts in its text, BitPAl needs 265 operations per 64-bit word and its packed variant 166, by BitPAl's own count. The derived kernel needs 181, fewer than the first and still more than the second. Bergeron and Hamel's listing, counted with one operation per integer-vector operation (a lower bound for any implementation that executes it operation for operation), needs 175; it grows more slowly than bitdp when substitution costs span a wide range, and is lower on six of the ten schemes above: (0, 3, 2), the four BitPAl sets and BLOSUM62. Operation counts include every AND, OR, XOR, NOT, shift and addition, and a + b + 1 and (x << 1) | 1 count as two.

## Speed

All throughput numbers come from one session of [bench/reproduce.sh](bench/README.md) on one core of a Ryzen 7 8840U under Ubuntu 26.04 in WSL2, in the Windows power mode "Best power efficiency": every tool is built from a pinned commit, every cell runs ten times after a warm-up with the tools alternating, and every run's cost of every pair is checked against an independent scalar DP. [bench/raw/r2/summary.csv](bench/raw/r2/summary.csv) has every median, quartile and paired ratio below.

The closest tool is BGSA, which runs Myers' and BitPAl's hand-derived kernels with one alignment per SIMD lane. In its own setting (all queries against all subjects of equal length) and at the same vector width, the derived edit-distance kernel reaches 0.78 (150 bp) and 1.18 (1 kbp) of BGSA's speed with AVX2, and 0.82 and 1.27 with AVX-512. On BitPAl's (2, -3, -5) weights it runs at 1.11 and 1.49 of BGSA's standard BitPAl with AVX2, and at 0.79 and 1.04 of its packed variant. BGSA's generator emits AVX-512 BitPAl code that does not compile, but the repository ships a hand-written AVX-512 packed BitPAl kernel for these weights; against it, bitdp with AVX-512 reaches 0.77 at 150 bp and 1.01 at 1 kbp. Schemes BGSA cannot express still run at the same order of speed, for example 49.8 GCUPS with AVX2 for transition/transversion costs on 1 kbp sequences.

On independent pairs, where each lane has its own text, batched bitdp is faster on every DNA workload than the aligners that fill the whole DP matrix: SeqAn's vectorized alignment, which also runs one pair per lane, parasail's correct kernels and ksw2 (for instance 62.1 against 3.35, 5.54 and 2.46 GCUPS for edit distance on 1 kbp pairs with AVX2). Aligners whose time grows with the alignment cost win where that cost is small: WFA2-lib on near-identical pairs (101 against 11.9 GCUPS on real Illumina reads checked against the reference span a mapper placed them on), and edlib, WFA2-lib and A*PA2 on long pairs of low divergence (at 100 kbp and 5 % divergence A*PA2 reaches 1240 GCUPS against 133 for bitdp). For BLOSUM62 parasail's scan kernel is faster (2.87 against 1.30 GCUPS with AVX2), since the spread of its costs puts 15 carry chains and 799 operations in every column.

On a real genome, 32 CRISPR-style 20 nt guides scanned along *E. coli* K-12 (4.6 Mbp) with transition/transversion costs return the same 78 sites as the plain DP, at 12.6 GCUPS with AVX-512. A 20 nt guide leaves most of a 64-bit word empty, so `Packed` puts three guides in each word and a spacer row after each guide to stop carries and shifts. On 96 guides that raises the scan from 12.5 to 31.1 GCUPS with AVX-512 and from 8.45 to 22.7 with AVX2, with hits identical to the unpacked kernel and the plain DP. Hyyrö, Fredriksson and Navarro packed patterns this way for unit costs in 2005; here the same spacer rule works for every derived kernel in search mode.

parasail's striped global kernels, which would otherwise be its fastest, return wrong costs on every workload (on the 150 bp edit set, 426 of 100 000 pairs score below their optimum with 16-bit lanes and 231 with 32-bit lanes, and the saturation flag is never set), so the comparison leaves them out. The wrong scores appear only when the gap opening penalty equals the extension penalty, which is how linear gaps are expressed in parasail; `bench/parasail_check.c` reproduces them.

## Limits

- Global and search (semi-global) modes, linear gap costs only. With affine gaps the value carried down a column can shift both up and down between thresholds, the copy graph acquires cycles, and the cascade above no longer applies.
- At most 32 distinct score differences and 24 distinct substitution costs.
- The derived program is not always the shortest. For indel distance it needs 6 operations, while the bit-vector LCS step (V + (V & M)) | (V - (V & M)), in the form Hyyrö gives, does the same column in 4.
- The carry-chain part still grows quadratically with the spread between substitution costs, which is why protein matrices are slow; Bergeron and Hamel's construction grows as O(c log c) in the gap cost c and the count of its listing is lower on six of the ten schemes in the table above.
- `Packed` handles search mode only, with patterns of at most 63 characters.
- Bytes outside the scheme's alphabet (N, for instance) count as the costliest substitution.
- Deriving a scheme happens inside the Zig compiler. Under Windows 11 with Zig 0.16.0, a debug build of a program that derives one kernel compiles in 2 seconds for unit costs and in 6 seconds and 487 MB of compiler memory for BLOSUM62; an optimized build takes 8 and 14 seconds.

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

// Search mode: the pattern against its best-matching stretch of the text,
// or every end position with cost at most k.
const Find = bitdp.Kernel(.{ .sub = &bitdp.schemes.tsTvCost, .gap = 2, .mode = .search });
var guides = try Find.Group.init(gpa, guide_list);
defer guides.deinit(gpa);
try guides.scan(gpa, genome, 4, &hits);

// Short patterns, several per lane word. `count` says how many it took.
const pk = Find.Packed.init(guide_list);
try pk.scan(gpa, genome, 4, &hits); // Hit.lane is the pattern index

// Any substitution cost function over a declared alphabet.
const TsTv = bitdp.Kernel(bitdp.schemes.ts_tv);
const Blosum = bitdp.Kernel(bitdp.schemes.blosum62Linear(4));
```

```sh
zig build test
zig build mve -Doptimize=ReleaseFast   # prints derived programs, verifies, times
```

## Using from C

Precompiled kernels are exported through a C interface ([include/bitdp.h](include/bitdp.h)) for edit distance and transition/transversion costs (global and search mode) and for BitPAl's (2, -3, -5) weights (global mode). Custom schemes stay Zig-only, because the derivation runs at compile time.

```c
#include <stdio.h>
#include "bitdp.h"

int main(void) {
    int64_t cost;
    if (bitdp_edit_global_distance("ACGTTGCA", 8, "ACGTGCA", 7, &cost) != BITDP_OK) return 1;
    printf("%lld\n", (long long)cost); /* 1 */
    return 0;
}
```

```sh
zig build lib -Doptimize=ReleaseFast   # zig-out/lib: libbitdp.a and libbitdp.so (bitdp.lib, bitdp.dll on Windows)
cc -O2 example.c -Izig-out/include zig-out/lib/libbitdp.a -o example
zig build test-c                       # the C tests, built as C99 and as C++
```

Each scheme has `_distance` for one pair with a pattern of up to 62 characters, an `_aligner_new` / `_aligner_distance` / `_aligner_free` trio for a pattern of any length aligned against many texts, and `_distances` for many independent pairs, one per SIMD lane. Search mode adds `_scan`, which returns every end position where a pattern costs at most k; the hit array is released with `bitdp_hits_free`. Every function returns 0 on success and writes its results through pointers.

On x86-64, `zig build lib` targets `x86_64_v3` (AVX2) unless you pass `-Dcpu`, so the library runs on most x86-64 machines from the last decade. Pass `-Dcpu=native` to build for the machine at hand, `-Dcpu=x86_64_v4` for AVX-512, or `-Dcpu=x86_64_v2` for older CPUs. The CPU level sets the number of SIMD lanes, so it changes speed but never results.

## References

See [REFERENCES.md](REFERENCES.md).

## License

MIT
