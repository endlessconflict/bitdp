# Results

Every throughput number comes from the measurement session of `bench/reproduce.sh` of 2026-11-10 (`raw/r2/`): `raw/r2/summary.csv` gives, per cell and tool, the median, quartiles and range of ten runs, and the paired ratio of every two tools of a cell (`summary-bestperf.csv`: the cells of steps `t4`, `t5` and `t6` again, five runs in the second power mode); the paper's Tables 2 and 3 print them through `paper/generated/numbers.tex`, and `bench/README.md` describes the harness. The session of 2026-10-23, measured at commit ac80e44, is superseded and kept in `raw/session-2026-10-23/`: it timed some tools around each pair and others around a whole set, and its shortest timed region was 0.02 s. The sessions of 2026-09-28 to 2026-10-12 that this file reported before (three runs per cell, measured by hand, the Ubuntu packages of edlib and parasail, BGSA's AVX2 builds with `-march=native`) are superseded; their logs stay in `raw/`, and their tables are in this file at commit 39cfd5b. `raw/r2/env.txt` records commit 1575ba7, the commit of the unpublished paper branch on which the session ran, and the hashes of its `src/` and `bench/` trees. The `src/` tree is that of the published commit eba1660: `git rev-parse eba1660:src` prints 09801db90058607ea25d351fea8a9244948d13ad. The recorded `bench/` hash, ad4f55c12b99f7eaa2e8cc800d53d5cba61ed1d1, is that of the tree before the machine's host name was replaced by `<host>` in the `uname` line of `raw/r2/env.txt` and `raw/session-2026-10-23/env.txt` on 2026-11-11, before publication (each file ends with a note that says so); the `bench/` tree of eba1660 differs from the recorded one only in that line and the note of `raw/session-2026-10-23/env.txt`.

## Genome scan (search mode)

`bitdp-scan` takes the *E. coli* K-12 MG1655 genome (NCBI NC_000913.3, 4 641 652 bp, forward strand), 32, 48 or 96 guides of 20 nt (genome 20-mers with up to two random substitutions each: the two random draws can keep a base or hit the same position; for the 32 guides, 20 have two substitutions and 12 have one), transition/transversion costs with gap 2, and reports every text position where a guide ends with cost at most 4. In every configuration of the session (AVX2 and AVX-512, one or three guides per lane), the hit sets (guide, end position, cost) of the packed kernel, the unpacked kernel and the scalar DP are all identical (78, 128 and 306 hits). Its rate counts 20 cells per guide and genome position, and its time includes preparing the guides.

## Scaling of the derived programs

Word operations per column as the score range grows (`zig build ablation -Dscheme=100..123`). In the first family the chained block stays at one carry chain while the number k of distinct differences grows. In the second, the spread between match and mismatch grows too, and so does the number of carry chains. The chain count equals the number of difference values above the smallest one and at or below mismatch cost minus gap cost, as the derivation predicts.

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

Word operations per column with parts of the compiler switched off (`zig build ablation -Dscheme=0..9`), re-run at commit 0af2499. An earlier version of this table had 17 and 33 in the "no truth tables" column, from an older commit.

| Scheme | all | no merge networks | no truth tables | symbolic direct only |
|---|---|---|---|---|
| edit (0,1,1) | 14 | 14 | 19 | 23 |
| indel (0,2,1) | 6 | 6 | 11 | 11 |
| (0,1,2) | 28 | 28 | 35 | 47 |
| (0,3,2) | 57 | 57 | 57 | 68 |
| bitpal (2,-3,-5) | 181 | 301 | 181 | 301 |
| bitpal (3,-4,-6) | 253 | 470 | 253 | 470 |
| bitpal (4,-5,-9) | 399 | 851 | 399 | 851 |
| bitpal (4,-7,-11) | 509 | 1176 | 509 | 1176 |
| tstv | 51 | 73 | 51 | 73 |
| blosum62, gap 4 | 799 | 1742 | 799 | 1742 |

## Bergeron and Hamel's algorithm on the same back end

`Options.bh_only` makes the compiler emit Bergeron and Hamel's vector algorithm with every integer vector held as bit planes and added by ripple carry; the same kernels run it (`bitdp-bench bh-edit`, `bh-bitpal`, `bh-tstv`, `bh-blosum`). Word operations per column (`raw/r2/ablation.txt`): the derived program, their algorithm as compiled (shared subexpressions, constant planes folded, without the factor not K in N), and their listing with every printed operation emitted on bit planes.

| Scheme | derived | B&H compiled | B&H listing on bit planes |
|---|---|---|---|
| edit (0,1,1) | 14 | 45 | 52 |
| indel (0,2,1) | 6 | 48 | 52 |
| (0,1,2) | 28 | 115 | 129 |
| (0,3,2) | 57 | 118 | 129 |
| bitpal (2,-3,-5) | 181 | 393 | 414 |
| bitpal (3,-4,-6) | 253 | 588 | 615 |
| bitpal (4,-5,-9) | 399 | 798 | 831 |
| bitpal (4,-7,-11) | 509 | 934 | 971 |
| tstv | 51 | 118 | 129 |
| blosum62, gap 4 | 799 | 732 | 760 |

In the all-against-all setting the derived kernels run 1.86 to 3.12 times as fast as the compiled Bergeron-Hamel programs (unit costs, BitPAl's (2, -3, -5) weights and transition/transversion costs, 100 bp to 3 kbp, both register widths), and 1.16 to 2.76 times as fast on the independent DNA pairs of 150 bp and 1 kbp. On BLOSUM62, the one scheme where their program is shorter, it is the faster one in batches (1.17 times the derived kernel's rate with AVX2, 1.13 with AVX-512) and the slower one for one pair at a time. Rates and ratios: `raw/r2/summary.csv`, tools `bitdp-bh-*`.

## Schemes

The `edit` scheme uses unit costs, `bitpal` scores match 2, mismatch -3 and gap -5, `tstv` charges 1 for a transition, 2 for a transversion and 2 per gap position, and `blosum` is BLOSUM62 with gap 4. All gaps are linear.
