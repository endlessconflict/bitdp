# Results

The machine is an AMD Ryzen 7 8840U (Zen 4, AVX-512), running Ubuntu under WSL2 (kernel 6.6.114.1) on one core, no threading. bitdp was built with Zig 0.16.0 (`-Doptimize=ReleaseFast`) twice, once with AVX-512 (`-Dcpu=znver4`, 8 lanes of 64 bits) and once with AVX2 only (`-Dcpu=x86_64_v3`, 4 lanes). Baselines were built with gcc 15.2.0 (`-O3 -march=native`): parasail 2.6.2 and edlib 1.2.7 from the Ubuntu packages, and BGSA from its repository (commit dffcf9e, 2022-02-24) with kernels from its own generator.

Numbers are medians of three runs, in billions of DP cells per second, where a pair contributes pattern length times text length cells. Every run's total cost was checked against the other tools, and all totals agree except where noted.

## All against all (BGSA's setting)

Every query against every subject, all of the same length: 64 × 12 500 random 150 bp sequences, and 16 × 1 250 random 1 kbp sequences. bitdp prepares subjects in groups of one per lane and aligns each query against each group (`Kernel.Group`), and its time includes that preparation. BGSA's figure is its own "cal GCUPS", which leaves out its preprocessing.

| Scheme | Tool | 150 bp | 1 kbp |
|---|---|---|---|
| edit | BGSA Myers, AVX2 | 94.3 | 116.9 |
| | bitdp, AVX2 | 68.9 | 108.9 |
| | BGSA Myers, AVX-512 | 102.4 | 120.7 |
| | bitdp, AVX-512 | 92.4 | 145.1 |
| bitpal (2, -3, -5) | BGSA BitPAl packed, AVX2 | 13.0 | 14.0 |
| | BGSA BitPAl non-packed, AVX2 | 8.6 | 9.2 |
| | bitdp, AVX2 | 8.6 | 11.7 |
| | bitdp, AVX-512 | 12.2 | 16.5 |
| tstv | bitdp, AVX2 | 32.5 | 45.8 |
| | bitdp, AVX-512 | 41.5 | 59.8 |

BGSA's generator produced AVX-512 BitPAl sources that do not compile (they use SSE type names), so there is no AVX-512 BitPAl row. BGSA supports only match/mismatch weights, so it has no `tstv` row.

## Pairwise (distinct pattern and text per pair)

Pairs from `gen.py`: random patterns, and texts derived from them by substitutions, insertions and deletions at the given divergence. With `batch`, bitdp runs one pair per SIMD lane (`Kernel.distances`), so each lane reads its own text and planes are gathered lane by lane.

| Workload | Scheme | bitdp, 1 pair | bitdp, 4 lanes AVX2 | bitdp, 8 lanes AVX-512 | edlib | parasail scan16 | parasail diag16 | parasail striped16 * |
|---|---|---|---|---|---|---|---|---|
| 100 000 pairs, 150 bp, 10 % | edit | 5.69 | 7.75 | 7.11 | 6.19 | 2.16 | 1.66 | 5.07 |
| | bitpal | 1.36 | 3.80 | 4.28 | | 2.18 | 1.74 | 5.11 |
| | tstv | 3.15 | 4.70 | 4.45 | | 2.21 | 1.72 | 5.04 |
| 5 000 pairs, 1 kbp, 15 % | edit | 18.07 | 35.20 | 34.89 | 21.73 | 5.19 | 1.98 | 8.62 |
| | bitpal | 2.14 | 8.93 | 11.30 | | 5.18 | 2.10 | 8.76 |
| | tstv | 8.19 | 19.23 | 20.83 | | 5.18 | 1.98 | 8.83 |
| 10 000 pairs, 300 aa, 30 % | blosum | 0.40 | 0.77 | 0.82 | | 2.87 | 2.19 | 5.00 |

\* parasail's striped global kernels returned a different total cost than every other tool on every workload, so they are listed but not counted as a correct baseline. On the 150 bp edit set, `parasail_nw_striped_16` scored 426 of the 100 000 pairs below their optimum. On the first such pair it returned -17 while the scan kernel, edlib, bitdp and an independent textbook DP all give -16.

The Ubuntu parasail build has no AVX-512 kernels, so its fair comparison is bitdp's AVX2 column. In pairwise mode the AVX-512 build gains little over AVX2 because each lane gathers its planes separately; the previous version of the kernel, which detected carries with a 64-bit overflow compare, ran pairwise AVX-512 edit about 17 % faster (45.6 against 38 GCUPS on the 1 kbp set) and AVX2 slower.

## Scaling of the derived programs

Word operations per column as the score range grows (`zig build ablation -Dscheme=100..123`). In the first family the chained block stays at one carry chain while the number k of distinct differences grows. In the second, the spread between match and mismatch grows too, and so does the number of carry chains. The chain count equals the number of difference values at or below mismatch cost minus gap cost, as the derivation predicts.

| Costs (match, mismatch, gap) | k | chains | ops | | Costs | k | chains | ops |
|---|---|---|---|---|---|---|---|---|
| 0, 1, 1 | 3 | 1 | 14 | | -1, 1, 2 | 6 | 2 | 56 |
| 0, 1, 2 | 5 | 1 | 28 | | -1, 2, 3 | 8 | 3 | 89 |
| 0, 1, 3 | 7 | 1 | 55 | | -1, 3, 4 | 10 | 4 | 126 |
| 0, 1, 4 | 9 | 1 | 71 | | -1, 4, 5 | 12 | 5 | 165 |
| 0, 1, 5 | 11 | 1 | 95 | | -1, 5, 6 | 14 | 6 | 210 |
| 0, 1, 6 | 13 | 1 | 119 | | -1, 6, 7 | 16 | 7 | 253 |
| 0, 1, 7 | 15 | 1 | 139 | | -1, 7, 8 | 18 | 8 | 302 |
| 0, 1, 8 | 17 | 1 | 155 | | -1, 8, 9 | 20 | 9 | 351 |
| 0, 1, 9 | 19 | 1 | 183 | | -1, 9, 10 | 22 | 10 | 408 |
| 0, 1, 10 | 21 | 1 | 215 | | -1, 10, 11 | 24 | 11 | 461 |
| 0, 1, 11 | 23 | 1 | 243 | | -1, 11, 12 | 26 | 12 | 522 |
| 0, 1, 12 | 25 | 1 | 267 | | -1, 12, 13 | 28 | 13 | 579 |

## Ablation

Word operations per column with parts of the compiler switched off (`zig build ablation -Dscheme=0..9`).

| Scheme | all | no merge networks | no truth tables | symbolic direct only |
|---|---|---|---|---|
| edit (0,1,1) | 14 | 14 | 17 | 23 |
| indel (0,2,1) | 6 | 6 | 11 | 11 |
| (0,1,2) | 28 | 28 | 33 | 47 |
| (0,3,2) | 57 | 57 | 57 | 68 |
| bitpal (2,-3,-5) | 181 | 301 | 181 | 301 |
| bitpal (3,-4,-6) | 253 | 470 | 253 | 470 |
| bitpal (4,-5,-9) | 399 | 851 | 399 | 851 |
| bitpal (4,-7,-11) | 509 | 1176 | 509 | 1176 |
| tstv | 51 | 73 | 51 | 73 |
| blosum62, gap 4 | 799 | 1742 | 799 | 1742 |

## Genome scan (search mode)

`bitdp-scan` takes the *E. coli* K-12 MG1655 genome (NCBI NC_000913.3, 4 641 652 bp, forward strand), 32 guides of 20 nt (genome 20-mers with two random substitutions each), transition/transversion costs with gap 2, and reports every text position where a guide ends with cost at most 4 (`Kernel.Group.scan`, one guide per lane).

| Build | Hits | Time | GCUPS |
|---|---|---|---|
| bitdp, 8 lanes AVX-512 | 78 | 0.30 s | 9.85 |
| bitdp, 4 lanes AVX2 | 78 | 0.36 s | 8.29 |
| plain scalar DP (the correctness check) | 78 | 13.3 s | 0.22 |

The hit lists agree guide by guide. A 20 nt guide fills only 20 of the 63 rows of a word, which is why the rate is well below the long-pattern numbers; packing several guides into one word is the obvious improvement. The scalar DP here is the plain reference loop, not an optimized aligner.

## Schemes

The `edit` scheme uses unit costs, `bitpal` scores match 2, mismatch -3 and gap -5, `tstv` charges 1 for a transition, 2 for a transversion and 2 per gap position, and `blosum` is BLOSUM62 with gap 4. All gaps are linear.
