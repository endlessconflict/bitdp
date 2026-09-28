# Results

The machine is an AMD Ryzen 7 8840U (Zen 4, AVX-512), running Ubuntu under WSL2 (kernel 6.6.114.1) on one core, no threading. bitdp built with Zig 0.16.0 (`-Doptimize=ReleaseFast`); baselines with gcc 15.2.0 (`-O3 -march=native`) against the Ubuntu packages parasail 2.6.2 and edlib 1.2.7. That parasail build ships SSE and AVX2 kernels but no AVX-512 ones, so bitdp was measured both with AVX-512 (`-Dcpu=znver4`, 8 lanes) and with AVX2 only (`-Dcpu=x86_64_v3`, 4 lanes).

Numbers are the median of three runs, in billions of DP cells per second (pattern length times text length, summed over pairs). Time includes pattern preprocessing and excludes file reading.

Inputs come from `gen.py`: random patterns and texts derived from them by substitutions, insertions and deletions at the given divergence.

| Workload | Scheme | bitdp, 1 pair | bitdp, 4 lanes AVX2 | bitdp, 8 lanes AVX-512 | edlib | parasail scan16 | parasail diag16 | parasail striped16 * |
|---|---|---|---|---|---|---|---|---|
| 100 000 pairs, 150 bp, 10 % | edit | 7.43 | 10.47 | 9.67 | 6.19 | 2.16 | 1.66 | 5.07 |
| | bitpal | 0.91 | 4.35 | 5.51 | | 2.18 | 1.74 | 5.11 |
| | tstv | 2.34 | 3.59 | 3.69 | | 2.21 | 1.72 | 5.04 |
| 5 000 pairs, 1 kbp, 15 % | edit | 18.23 | 30.70 | 50.44 | 21.73 | 5.19 | 1.98 | 8.62 |
| | bitpal | 1.16 | 7.29 | 11.59 | | 5.18 | 2.10 | 8.76 |
| | tstv | 5.52 | 14.04 | 16.82 | | 5.18 | 1.98 | 8.83 |
| 10 000 pairs, 300 aa, 30 % | blosum | 0.33 | 0.59 | 0.63 | | 2.87 | 2.19 | 5.00 |

\* parasail's striped global kernels returned a different total cost than every other tool on every workload (for example 1 316 488 against 1 316 060 for edit on the 150 bp set), so they are listed but not counted as a correct baseline. bitdp, edlib and parasail's scan and diagonal kernels agree exactly on every total.

The `edit` scheme uses unit costs, `bitpal` scores match 2, mismatch -3 and gap -5, `tstv` charges 1 for a transition, 2 for a transversion and 2 per gap position, and `blosum` is BLOSUM62 with gap 4. All gaps are linear.
