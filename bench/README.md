# Benchmark

Global alignment cost of every pair in a file, computed by bitdp, parasail and edlib on the same machine, the same operating system and the same input. Every tool prints the sum of all alignment costs, and a run only counts if that sum agrees with the other tools.

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

# bitdp (from the repository root).
zig build bench -Doptimize=ReleaseFast -Dcpu=native
./zig-out/bin/bitdp-bench bitpal bench/dna150.tsv
```

Schemes: `edit` (unit costs), `bitpal` (score 2, -3, gap -5), `tstv` (transition 1, transversion 2, gap 2), `blosum` (BLOSUM62, gap 4). Time covers pattern preprocessing plus alignment, and excludes file reading.
