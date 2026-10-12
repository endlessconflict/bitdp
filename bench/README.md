# Benchmark

Global alignment cost of every pair in a file, computed by bitdp, parasail, edlib and ksw2 on the same machine, the same operating system and the same input. Every tool prints the sum of all alignment costs, and a run only counts if that sum agrees with the other tools.

## Reproduce (Linux)

```sh
# Inputs: pattern<TAB>text per line.
python3 gen.py dna150.tsv 100000 150 0.10 dna 1
python3 gen.py dna1000.tsv 5000 1000 0.15 dna 2
python3 gen.py prot300.tsv 10000 300 0.30 protein 3

# Baselines (Debian/Ubuntu: apt install libparasail-dev libedlib-dev).
gcc -O3 -march=native baselines.c -o baselines -lparasail -ledlib -lstdc++
./baselines edlib edit dna150.tsv
./baselines scan16 bitpal dna150.tsv      # also striped16, diag16, striped_sat, scan_sat

# Per-pair check of parasail's striped and scan kernels against a plain DP (counts, saturation flags,
# first differing pair; add "list" to print every differing pair).
gcc -O3 -march=native parasail_check.c -o parasail_check -lparasail
./parasail_check edit dna150.tsv

# ksw2 (inside a checkout of https://github.com/lh3/ksw2).
gcc -O3 -march=native ksw2_baseline.c ksw2_extz2_sse.c kalloc.c -o ksw2_baseline
./ksw2_baseline tstv dna150.tsv           # also edit, bitpal

# bitdp (from the repository root). The paper's builds: -Dcpu=znver4 (AVX-512, 8 lanes) and
# -Dcpu=x86_64_v3 (AVX2, 4 lanes); -Dcpu=native equals znver4 on the benchmark machine.
zig build bench -Doptimize=ReleaseFast -Dcpu=znver4
./zig-out/bin/bitdp-bench bitpal bench/dna150.tsv
```

All against all, as BGSA runs it (one sequence per line, all of the same length):

```sh
./zig-out/bin/bitdp-bench edit ava queries.txt subjects.txt
# BGSA (commit dffcf9e): generate a kernel, copy it into the matching original/BGSA_* folder, remove the
# prebuilt aligner/convert binaries, build with make cc=gcc (AVX-512 folders: make cc=gcc CFLAGS="-march=native -w"), then
./aligner -q queries.txt -d subjects.txt -f result.txt -N 1

# Inputs with a fixed seed (one sequence per line), and the total checks used for Table 4:
python3 gen.py ava queries.txt subjects.txt 16 1250 1000 2
python3 bgsa_check.py result.txt queries.txt subjects.txt COST_SUM   # BGSA's scores summed vs bitdp's cost_sum
gcc -O3 -march=native dp_check.c -o dp_check
./dp_check 2 -3 -5 queries.txt subjects.txt                          # independent DP total (scores)
```

BGSA kernels used in RESULTS.md: `java -jar generator.jar -m 0 -a avx2` (and `-a avx512`) for Myers; `-M 2 -I -3 -G -5 -t packed -a avx2` and `-t non-packed -a avx2` for BitPAl. The generator's AVX-512 BitPAl output does not compile, so the AVX-512 packed BitPAl row uses the hand-written `original/BGSA_AVX512/align_core.c` that BGSA ships. The inputs were uniform random DNA, 64 queries and 12 500 subjects of 150 bp, and 16 queries and 1 250 subjects of 1 kbp; the script that drew them was not kept, so new random files (`gen.py ava`) give the same setting but not the same bytes.

Genome scan (search mode, every hit with cost at most K):

```sh
zig build scan -Doptimize=ReleaseFast -Dcpu=znver4
./zig-out/bin/bitdp-scan genome.fa 32 4 verify   # 32 guides, K = 4, hit sets (guide, end, cost) checked against the scalar DP
./zig-out/bin/bitdp-scan genome.fa 96 4 verify packed
```

Real reads (the first 200 000 reads of ENA run ERR022075 against *E. coli* K-12 MG1655, NC_000913.3):

```sh
curl -s ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022075/ERR022075_1.fastq.gz | gunzip -c | head -n 800000 > r1.fq
bwa index ref.fa && bwa mem -t 8 ref.fa r1.fq | samtools view -F 0x904 - | python3 reads.py ref.fa > reads.tsv
./zig-out/bin/bitdp-bench edit bench/reads.tsv batch
```

RESULTS.md reports the median of three runs; `bench/median3.sh COMMAND ARGS...` runs any of the commands above three times and prints the median.

Schemes: `edit` (unit costs), `bitpal` (score 2, -3, gap -5), `tstv` (transition 1, transversion 2, gap 2), `blosum` (BLOSUM62, gap 4). Time covers pattern preprocessing plus alignment, and excludes file reading.

Environment of the measurements (`bench/raw/env.txt`, `env2.txt`): Ubuntu 26.04 LTS under WSL2, gcc 15.2.0, the Ubuntu packages of parasail 2.6.2 and edlib 1.2.7, bwa 0.7.19 (`-t 8`, no `-K`) and samtools 1.22.1. The parasail package has SSE2, SSE4.1 and AVX2 kernels but no AVX-512 ones, so on this CPU its dispatcher runs the AVX2 kernels.
