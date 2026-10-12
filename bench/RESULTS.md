# Results

The machine is an AMD Ryzen 7 8840U (Zen 4, AVX-512), running Ubuntu under WSL2 (kernel 6.6.114.1) on one core, no threading. bitdp was built with Zig 0.16.0 (`-Doptimize=ReleaseFast`) twice, once with AVX-512 (`-Dcpu=znver4`, 8 lanes of 64 bits) and once with AVX2 only (`-Dcpu=x86_64_v3`, 4 lanes). parasail 2.6.2 and edlib 1.2.7 are the Ubuntu library packages (not rebuilt), called from `baselines.c` compiled with gcc 15.2.0 (`-O3 -march=native`). BGSA is built from its repository (commit dffcf9e, 2022-02-24) with its own makefile and gcc 15.2.0, with `-O3 -march=native` for every build (for AVX-512, `make cc=gcc CFLAGS="-march=native -w"` replaces the makefile's `-mfma -mavx512f -mavx512cd`). On this Zen 4 machine `-march=native` also allows AVX-512VL instructions in BGSA's 256-bit code.

Numbers are medians of three runs, in billions of DP cells per second, where a pair contributes pattern length times text length cells. Every run's total cost was checked against the other tools, and all totals agree except where noted.

## All against all (BGSA's setting)

Every query against every subject, all of the same length: 64 × 12 500 random 150 bp sequences, and 16 × 1 250 random 1 kbp sequences. bitdp prepares subjects in groups of one per lane and aligns each query against each group (`Kernel.Group`), and its time includes that preparation. BGSA's figure is its own "cal GCUPS", which leaves out its preprocessing.

| Scheme | Tool | 150 bp | 1 kbp |
|---|---|---|---|
| edit | BGSA Myers, AVX2 | 110.8 | 114.3 |
| | bitdp, AVX2 | 73.5 | 114.9 |
| | BGSA Myers, AVX-512 | 116.8 | 125.9 |
| | bitdp, AVX-512 | 102.5 | 156.2 |
| bitpal (2, -3, -5) | BGSA BitPAl packed, AVX2 | 12.7 | 13.7 |
| | BGSA BitPAl non-packed, AVX2 | 8.3 | 8.8 |
| | bitdp, AVX2 | 8.1 | 11.4 |
| | BGSA BitPAl packed, AVX-512 (hand-written) | 15.2 | 16.4 |
| | bitdp, AVX-512 | 12.0 | 16.3 |
| tstv | bitdp, AVX2 | 28.6 | 48.3 |
| | bitdp, AVX-512 | 42.6 | 62.1 |

BGSA's generator produces AVX-512 BitPAl sources that do not compile: they declare `sse_read_t` variables, which the AVX-512 sources don't define. BGSA's repository, however, ships a hand-written AVX-512 packed BitPAl kernel for exactly (2, -3, -5) in `original/BGSA_AVX512/align_core.c`. Built from source with the same command as the AVX-512 Myers build (prebuilt binaries removed), its scores matched an independent DP on all 800 000 + 20 000 pairs, and it ran at 15.05, 15.09 and 15.23 GCUPS at 150 bp and 15.59, 15.35 and 15.87 at 1 kbp (`median3.sh`, 2026-10-09). BGSA supports only match/mismatch weights, so it has no `tstv` row.

The all-against-all, pairwise and real-read tables on this page were measured plugged in, with the Windows power plan Balanced, the Windows power mode "Best power efficiency" (`bench/raw/power_mode_events.txt`) and the tools of each row alternating run by run. The AVX-512 rows for edit and (2, -3, -5) come from a session on 2026-10-12 (`bench/raw/t4_avx512_interleaved.txt`): BGSA's two AVX-512 builds, compiled with `-march=native`, and bitdp (built from 0af2499 with `-Dcpu=znver4`, the binary of the other session) alternate run by run, and the three runs of a cell differ by at most 5 % of its median. The rebuilt BGSA binaries write result files byte-identical to those of the earlier builds, and their totals equal bitdp's on all four inputs (`bench/raw/bgsa_avx512_native.txt`, which also holds an earlier measurement of these rows with each tool's three runs back to back; it is no longer used). The power mode matters: with "Best performance" the same binaries ran 1.38 to 1.47 times as fast in the same hour (`bench/raw/t4_avx512_bestperf.txt`), so rates from sessions in different modes must not be combined. Between sessions the machine's rates also drift: the earlier AVX-512 Myers build ran at 116.0 GCUPS at 150 bp on 2026-10-10, against 102.4 twelve days before (the value in this file at commit cb765b7, 2026-09-28; its individual runs were not kept). The other rows of this table, and the pairwise and real-read tables, come from a session on 2026-10-10 (`bench/raw/e10_session.txt`), with bitdp, `baselines.c` and `ksw2_baseline.c` built from commit 0af2499 and every tool's totals equal. Earlier versions of these tables (2026-09-28) combined numbers from different sessions, and their pairwise and real-read bitdp columns were measured on code from before commit 0af2499, which runs those workloads up to about twice as fast (same-session comparison of the builds in `bench/raw/e10_versions.txt` and `e10_versions2.txt`). Within a session the three runs of a cell differed by a median of 4 % of their median and by at most 28 %.

## Pairwise (distinct pattern and text per pair)

Pairs from `gen.py`: random patterns, and texts derived from them by substitutions, insertions and deletions at the given divergence. With `batch`, bitdp runs one pair per SIMD lane (`Kernel.distances`), so each lane reads its own text and planes are gathered lane by lane.

| Workload | Scheme | bitdp, 1 pair | bitdp, 4 lanes AVX2 | bitdp, 8 lanes AVX-512 | edlib | ksw2 extz2_sse | parasail scan16 | parasail diag16 | parasail striped16 * |
|---|---|---|---|---|---|---|---|---|---|
| 100 000 pairs, 150 bp, 10 % | edit | 9.83 | 15.53 | 14.25 | 7.15 | 1.50 | 2.22 | 1.91 | 5.14 |
| | bitpal | 1.82 | 5.44 | 6.48 |  | 1.30 | 2.19 | 1.82 | 5.20 |
| | tstv | 5.12 | 10.71 | 10.74 |  | 0.72 | 2.21 | 1.88 | 5.23 |
| 5 000 pairs, 1 kbp, 15 % | edit | 23.55 | 60.88 | 60.44 | 25.49 | 2.20 | 5.28 | 2.18 | 8.95 |
| | bitpal | 2.44 | 9.68 | 13.08 |  | 2.28 | 5.41 | 2.17 | 8.90 |
| | tstv | 8.98 | 31.02 | 32.87 |  | 0.96 | 5.42 | 2.13 | 8.05 |
| 10 000 pairs, 300 aa, 30 % | blosum | 0.49 | 1.20 | 1.39 |  |  | 2.82 | 1.92 | 4.70 |

\* parasail's striped global kernels returned a different total cost than every other tool on every workload, so they are listed but not counted as a correct baseline. `parasail_check.c` compares every pair with a plain DP:

| Workload | striped16 wrong | striped32 wrong | scan16 / scan32 wrong |
|---|---|---|---|
| 150 bp, edit / bitpal / tstv | 426 / 426 / 430 of 100 000 | 231 / 232 / 223 | 0 / 0 |
| 1 kbp, edit / bitpal / tstv | 46 / 46 / 51 of 5 000 | 18 / 18 / 19 | 0 / 0 |
| 300 aa, blosum | 842 of 10 000 | 395 | 0 / 0 |

Every wrong score is below the optimum, and the saturation flag is never set. On short random pairs the striped kernels fail only when the gap opening penalty equals the extension penalty, i.e. for linear gaps, which parasail's README allows ("greater than or equal"); with open > extend they were correct on 4 × 300 000 pairs. A minimal case is `TAA` against `T` with match 0, mismatch -1, open = extend = 1: every striped entry point (SSE2, SSE4.1, AVX2; 16 and 32 bit) returns -3, while scan and the scalar `parasail_nw` return the optimum -2.

ksw2 is the SIMD DP of minimap2 (github.com/lh3/ksw2, commit 289609b), run through `ksw_extz2_sse` in score-only mode with gap open 0 and built by `ksw2_baseline.c` with gcc 15.2.0 (`-O3 -march=native`). It uses 128-bit SSE registers. The `tstv` run needs the `KSW_EZ_GENERIC_SC` flag, because without it ksw2 reads only the match and mismatch entries of the matrix.

The Ubuntu parasail build has no AVX-512 kernels, so its fair comparison is bitdp's AVX2 column. In pairwise mode the AVX-512 build gains little over AVX2 because each lane gathers its planes separately. A previous version of the kernel, which detected carries with a 64-bit overflow compare, ran pairwise AVX-512 edit faster and AVX2 slower; its measurements were not kept.

## Real reads (candidate verification)

Illumina reads of *E. coli* K-12 MG1655 from ENA run ERR022075 (study PRJEB2323, Genome Analyzer IIx): the first 200 000 reads of mate 1, mapped to NC_000913.3 with bwa mem 0.7.19. Each primary alignment without clipping and without N gives one pair, the read against the reference span bwa placed it on (`reads.py`), which is the verification step of a read mapper. That leaves 193 660 pairs, all 100 bp. The reads are close to the reference: the sum of unit edit distances is 38 308, about 0.2 per read, and the sum of bwa's NM tags is 38 353. NM counts the edits of the alignment bwa chose, so it can only be at or above the edit distance.

| Scheme | bitdp, 1 pair | bitdp, 4 lanes AVX2 | bitdp, 8 lanes AVX-512 | edlib | ksw2 extz2_sse | parasail scan16 |
|---|---|---|---|---|---|---|
| edit | 7.70 | 10.98 | 10.47 | 4.79 | 1.55 | 1.59 |
| bitpal | 1.64 | 4.55 | 5.36 |  | 1.52 | 1.57 |
| tstv | 4.42 | 8.04 | 7.90 |  | 0.82 | 1.58 |

All tools give the same totals (edit 38 308, bitpal -38 540 294, tstv 65 477). Pairs this short spend a large share of their time outside the column loop (pattern preparation, reading the pair), which is why every rate is lower than on the 1 kbp set. edlib bands its DP by the edit distance it finds, so near-identical pairs favour it; batched bitdp is still ahead, by a factor of 2.3.

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

## Genome scan (search mode)

`bitdp-scan` takes the *E. coli* K-12 MG1655 genome (NCBI NC_000913.3, 4 641 652 bp, forward strand), 32 guides of 20 nt (genome 20-mers with up to two random substitutions each: the two random draws can keep a base or hit the same position; for these 32 guides, 20 have two substitutions and 12 have one), transition/transversion costs with gap 2, and reports every text position where a guide ends with cost at most 4 (`Kernel.Group.scan`, one guide per lane).

| Build | Hits | Time | GCUPS |
|---|---|---|---|
| bitdp, 8 lanes AVX-512 | 78 | 0.30 s | 9.85 |
| bitdp, 4 lanes AVX2 | 78 | 0.36 s | 8.29 |
| plain scalar DP (the correctness check) | 78 | 13.3 s | 0.22 |

The hit lists agree guide by guide. The scalar DP here is the plain reference loop, not an optimized aligner.

A 20 nt guide fills only 20 of the 63 rows of a word. `Kernel.Packed` puts three guides in each lane word, each followed by a spacer row where additions and shifts see forced bits, so every guide starts from the row-0 boundary; the scores live in bit fields of one accumulator word per lane. Same genome, scheme and K, medians of three runs (`bitdp-scan ref.fa N 4 [packed]`):

| Build | Guides | Hits | One guide per lane | Three guides per lane |
|---|---|---|---|---|
| AVX-512, 8 lanes | 48 | 128 | 0.374 s, 11.9 GCUPS | 0.165 s, 27.0 GCUPS |
| | 96 | 306 | 0.748 s, 11.9 GCUPS | 0.329 s, 27.1 GCUPS |
| AVX2, 4 lanes | 48 | 128 | 0.547 s, 8.1 GCUPS | 0.227 s, 19.7 GCUPS |
| | 96 | 306 | 1.119 s, 8.0 GCUPS | 0.464 s, 19.2 GCUPS |

The packed hits match the unpacked ones, and on 48 guides the scalar DP (128 hits, 20.6 s). A later check from the same code (2026-10-09) compared the full hit sets, guide, end position and cost, of the packed kernel, the unpacked kernel and the scalar DP on 32, 48 and 96 guides: all identical (78, 128 and 306 hits). Packing is 2.3 to 2.4 times faster rather than 3 times because each column also pays for the spacer masks (one or two operations per addition and shift) and for reading three score fields.

## Schemes

The `edit` scheme uses unit costs, `bitpal` scores match 2, mismatch -3 and gap -5, `tstv` charges 1 for a transition, 2 for a transversion and 2 per gap position, and `blosum` is BLOSUM62 with gap 4. All gaps are linear.
