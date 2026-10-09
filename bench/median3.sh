#!/bin/sh
# Run a benchmark command three times and report the median throughput, as in
# RESULTS.md. Reads "gcups=X" (bitdp-bench, baselines, ksw2_baseline) or BGSA's
# "cal GCUPS is X" from each run's output and keeps every run's line.
#
# usage: bench/median3.sh COMMAND [ARGS...]
#   e.g. bench/median3.sh ./zig-out/bin/bitdp-bench edit bench/dna150.tsv batch
#        bench/median3.sh ./aligner -q queries.txt -d subjects.txt -f result.txt -N 1
set -eu
[ $# -ge 1 ] || { echo "usage: $0 COMMAND [ARGS...]" >&2; exit 2; }
runs=""
for i in 1 2 3; do
    out=$("$@" 2>&1)
    rate=$(printf '%s\n' "$out" | sed -n -e 's/.*gcups=\([0-9.]*\).*/\1/p' -e 's/.*cal GCUPS is \([0-9.]*\).*/\1/p' | tail -n 1)
    [ -n "$rate" ] || { printf 'run %s: no GCUPS figure in output:\n%s\n' "$i" "$out" >&2; exit 1; }
    printf 'run %s: %s\n' "$i" "$(printf '%s\n' "$out" | grep -iE 'gcups' | tail -n 1)"
    runs="$runs $rate"
done
printf '%s\n' $runs | sort -g | sed -n 2p | sed 's/^/median gcups=/'
