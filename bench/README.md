# Benchmark

Global alignment cost of every pair in a file, computed by bitdp and by the other aligners on the same machine, the same operating system and the same input. Every driver prints its rate, the sum of all alignment costs and a hash of the per-pair costs, and a run only counts if that hash equals the one of an independent scalar DP.

## One command (WSL2 / Ubuntu)

```sh
bench/reproduce.sh --dry-run     # list every step and command
bench/reproduce.sh --quick       # tiny inputs, 2 repetitions, output under ~/bitdp-r2-quick: a test of the harness
bench/reproduce.sh               # the measurement session (about 1-2 h)
bench/reproduce.sh t5 summary    # selected steps: deps build env inputs ref t4 t5 t6 long scan robust summary
```

- `deps`: `deps.sh` fetches every competitor at a pinned commit into `~/bitdp-deps` and builds the libraries (edlib v1.2.7, parasail v2.6.2, WFA2-lib v2.3.6, A\*PA2 bf2e14e with its pinned Rust nightly; ksw2 289609b, BGSA dffcf9e and SeqAn 2.5.3 are compiled with their drivers); `versions.txt` records the commits.
- `build`: bitdp (Zig 0.16.0, `ReleaseFast`) for AVX-512 (`-Dcpu=znver4`, 8 lanes) and AVX2 (`-Dcpu=x86_64_v3`, 4 lanes); the C and C++ drivers with gcc 15.2.0 and `-O3 -march=native`, except the AVX2 rows (BGSA's AVX2 kernels, SeqAn's AVX2 build), which use `-march=x86-64-v3`; a scripted `objdump` check confirms that the AVX2 binaries hold no EVEX-encoded (AVX-512) instruction (`raw/r2/evex.txt`).
- `inputs`: `gen.py` with recorded seeds; NCBI RefSeq NC_000913.3 (*E. coli* K-12 MG1655); the first 200 000 reads of mate 1 of ENA run ERR022075, mapped with `bwa mem -K 10000000` and turned into read/reference pairs by `reads.py`. `raw/r2/inputs.sha256` lists every file.
- `ref`: `dp_check.c`, a plain scalar DP that shares no code with the aligners, gives the reference cost of every pair.
- `t4`, `t5`, `t6`, `long`: the cells of the tables. Each cell runs every tool once as a discarded warm-up, then `REPS` (10) times in round-robin order, the starting tool rotating by one per repetition, each run pinned with `taskset -c 2`. After every run the hash is compared with the reference; on a mismatch the harness reruns the tool with `BENCH_DUMP` set, prints the first differing pair and stops. parasail's striped kernels are known to return wrong costs for linear gaps (below): their runs are kept and marked `expected-mismatch`.
- `scan`: the genome scan, hit sets checked once against the scalar DP per configuration, then timed runs.
- `robust`: Table 4 again in the Windows power mode "Best performance", as a robustness check of its ratios.
- `summary`: `summary.py` writes `raw/r2/summary.csv`: per cell and tool the median, quartiles and range of the rate, and for every pair of tools the median and quartiles of their ratio within the same repetition.

The measurement steps set the Windows power mode "Best power efficiency" (`powermode.ps1`, through WSL interop), log processor performance once per second (`raw/r2/power-samples.txt`) and restore the mode found at the start (`raw/r2/power.txt`). Rates from sessions in different power modes must not be combined (RESULTS.md).

## Drivers

| File | Tool | Notes |
|---|---|---|
| `bench.zig` | bitdp | one pair at a time (`Kernel.Aligner`), batched one pair per lane (`batch`), or all against all (`ava`, `Kernel.Group`) |
| `baselines.c` | parasail (`scan16`, `diag16`, `striped16`, ...), edlib | parasail 2.6.2 has no AVX-512 alignment kernels; its dispatcher runs the AVX2 ones |
| `ksw2_baseline.c` | ksw2 `ksw_extz2_sse` | score only, gap open 0; `KSW_EZ_GENERIC_SC` for the transition/transversion matrix |
| `wfa2_baseline.c` | WFA2-lib | score only, end to end, no heuristic (its default wf-adaptive pruning is switched off); edit and gap-linear metrics only |
| `astarpa2_baseline.c` | A\*PA2 `simple` and `full` | unit costs only; no score-only mode, so every call also builds a CIGAR |
| `seqan_baseline.cpp` | SeqAn 2 `globalAlignmentScore` with `ExecutionPolicy<Serial, Vectorial>` | 16-bit scores when no pair can overflow them, else 32-bit; the line says which |
| `bgsa_check.py` | BGSA (`aligner`) | reads BGSA's result file; BGSA's rate is its own "cal GCUPS", which leaves out its preprocessing |
| `dp_check.c` | reference | plain scalar DP |

All drivers share the per-pair hash of `pairhash.h`: FNV-1a (64 bit) over each cost in decimal followed by a newline, in input order (all against all: query-major), which equals the FNV-1a of the file `BENCH_DUMP` names. A cost is the negated score for score-based tools. Time covers pattern preprocessing plus alignment and excludes file reading (BGSA: its own figure).

BGSA kernels: `java -jar generator.jar -m 0 -a avx2` (and `-a avx512`) for Myers; `-M 2 -I -3 -G -5 -t packed -a avx2` and `-t non-packed -a avx2` for BitPAl. The generator's AVX-512 BitPAl output does not compile, so the AVX-512 packed BitPAl row uses the hand-written `original/BGSA_AVX512/align_core.c` that BGSA ships for exactly (2, -3, -5). BGSA supports only match/mismatch weights, so it has no transition/transversion row.

## Schemes

`edit` (unit costs), `bitpal` (score 2, -3, gap -5), `tstv` (transition 1, transversion 2, gap 2), `blosum` (BLOSUM62, gap 4). All gaps are linear.

## parasail's striped kernels

`parasail_check.c` compares every pair of parasail's striped and scan kernels with a plain DP (counts, saturation flags, first differing pair; add `list` to print every differing pair):

```sh
gcc -O3 -march=native parasail_check.c -o parasail_check -I$HOME/bitdp-deps/parasail/install/include $HOME/bitdp-deps/parasail/install/lib/libparasail.a -lm
./parasail_check edit dna150.tsv
```

The striped global kernels fail only when the gap opening penalty equals the extension penalty, that is, for linear gaps, which parasail's README allows; RESULTS.md has the counts. A minimal case is `TAA` against `T` with match 0, mismatch -1, open = extend = 1: every striped entry point returns -3, while scan and the scalar `parasail_nw` return the optimum -2.

## Other programs

- `scan.zig` (`zig build scan`): `bitdp-scan GENOME.fa GUIDES K [verify] [packed]`, guides drawn from the genome with a fixed seed.
- `ablation.zig` (`zig build ablation -Dscheme=N`): operation counts per builder configuration.
- `median3.sh`: the median of three runs, as RESULTS.md used before `reproduce.sh`.

RESULTS.md and `raw/` (outside `raw/r2/`) hold the sessions of 2026-09-28 to 2026-10-12, measured by hand with three runs per cell; `raw/r2/` holds the sessions of `reproduce.sh`.
