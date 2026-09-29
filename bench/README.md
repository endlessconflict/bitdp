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

# ksw2 (inside a checkout of https://github.com/lh3/ksw2).
gcc -O3 -march=native ksw2_baseline.c ksw2_extz2_sse.c kalloc.c -o ksw2_baseline
./ksw2_baseline tstv dna150.tsv           # also edit, bitpal

# bitdp (from the repository root).
zig build bench -Doptimize=ReleaseFast -Dcpu=native
./zig-out/bin/bitdp-bench bitpal bench/dna150.tsv
```

All against all, as BGSA runs it (one sequence per line, all of the same length):

```sh
./zig-out/bin/bitdp-bench edit ava queries.txt subjects.txt
# BGSA: generate a kernel (java -jar generator.jar -m 0 -a avx2), build with make cc=gcc, then
./aligner -q queries.txt -d subjects.txt -f result.txt -N 1
```

Genome scan (search mode, every hit with cost at most K):

```sh
zig build scan -Doptimize=ReleaseFast -Dcpu=native
./zig-out/bin/bitdp-scan genome.fa 32 4 verify   # 32 guides, K = 4, check against the scalar DP
```

Real reads (the first 200 000 reads of ENA run ERR022075 against *E. coli* K-12 MG1655, NC_000913.3):

```sh
curl -s ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022075/ERR022075_1.fastq.gz | gunzip -c | head -n 800000 > r1.fq
bwa index ref.fa && bwa mem ref.fa r1.fq | samtools view -F 0x904 - | python3 reads.py ref.fa > reads.tsv
./zig-out/bin/bitdp-bench edit bench/reads.tsv batch
```

Schemes: `edit` (unit costs), `bitpal` (score 2, -3, gap -5), `tstv` (transition 1, transversion 2, gap 2), `blosum` (BLOSUM62, gap 4). Time covers pattern preprocessing plus alignment, and excludes file reading.
