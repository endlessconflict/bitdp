#!/bin/bash
# One command for every timed measurement of the paper. It builds the pinned competitors (bench/deps.sh)
# and the drivers, draws every input from recorded seeds or fixed public data, computes the reference
# cost of every pair with an independent scalar DP (dp_check.c), runs every cell REPS times in
# round-robin order across its tools (the start rotates by one tool per repetition) after one discarded
# warm-up run per tool, pinned to one core, and checks every run pair by pair against the reference
# (per-pair cost hash, pairhash.h); on a mismatch it reports the first differing pair and stops. Raw
# logs go to bench/raw/r2/, the summary to bench/raw/r2/summary.csv (bench/summary.py).
# Runs in WSL2 (Ubuntu). Windows side: bench/powermode.ps1 sets the power mode "Best power efficiency"
# for the measurement steps, samples the processor performance once per second, and restores the mode
# found at the start; step "robust" re-runs Table 4 in "Best performance"; step "compile" times one
# derivation with Zig in WSL and, through WSL interop, with the Windows Zig (ZIG_WIN, else zig.exe on PATH).
#
# usage: bench/reproduce.sh [--dry-run] [--quick] [STEP ...]
#   steps (default: all, in this order): deps build ops env inputs ref t4 t5 t6 long scan compile robust insn summary
#   --dry-run  print every command instead of running it
#   --quick    tiny inputs, 2 repetitions, output under $WORK-quick: a test of the harness, not a measurement
# environment: DEPS ($HOME/bitdp-deps), WORK ($HOME/bitdp-r2), REPS (10), CORE (2), RREPS (5, step robust),
#   CREPS (5, step compile), ZIG_WIN
set -euo pipefail
R=$(cd "$(dirname "$0")/.." && pwd)
DEPS=${DEPS:-$HOME/bitdp-deps}
WORK=${WORK:-$HOME/bitdp-r2}
REPS=${REPS:-10}
RREPS=${RREPS:-5}
CREPS=${CREPS:-5}
CORE=${CORE:-2}
DRY=0 QUICK=0 STEPS=()
for a in "$@"; do
    case $a in --dry-run) DRY=1 ;; --quick) QUICK=1 ;; *) STEPS+=("$a") ;; esac
done
[ ${#STEPS[@]} -gt 0 ] || STEPS=(deps build ops env inputs ref t4 t5 t6 long scan compile robust insn summary)
RAW="$R/bench/raw/r2" GEN="$R/paper/generated"
if [ $QUICK = 1 ]; then WORK=$WORK-quick REPS=2 RREPS=2 CREPS=1 RAW=$WORK/raw GEN=$WORK/generated; fi
B=$WORK/bin IN=$WORK/in S=$WORK/src LOG=runs.txt
ZIG=$DEPS/zig/zig
PMW=$(wslpath -w "$R/bench/powermode.ps1")
pm() { powershell.exe -NoProfile -ExecutionPolicy Bypass -File "$PMW" "$@" | tr -d '\r'; }
x() { # print a command, run it unless --dry-run
    echo "+ $*"
    [ $DRY = 1 ] || eval "$@"
}
mkdir -p "$RAW" "$B" "$IN" "$S"

# ---- deps: pinned competitor sources and builds (skipped when $DEPS/versions.txt exists)
deps() {
    if [ -f "$DEPS/versions.txt" ]; then echo "deps: $DEPS/versions.txt exists, keeping the builds"; else x "sh '$R/bench/deps.sh' all"; fi
    x "cp '$DEPS/versions.txt' '$RAW/versions.txt'"
}

# ---- build: drivers; AVX2 rows with -march=x86-64-v3 (no AVX-512 instruction), the rest -march=native
bgsa_build() { # NAME FOLDER MARCH [GENERATOR ARGS]: no generator args = the folder's own align_core.c
    local name=$1 folder=$2 march=$3 d=$WORK/bgsa-$1
    shift 3
    x "rm -rf $d && cp -r $DEPS/BGSA/original/$folder $d && rm -f $d/aligner $d/convert"
    [ $# -eq 0 ] || x "(cd $d && java -jar $DEPS/BGSA/generator/generator.jar $* -d gen > gen.log && mv gen/align_core.c .)"
    x "make -s -C $d cc=gcc CFLAGS='-march=$march -w' > $d/build.log 2>&1 && cp $d/aligner $B/bgsa-$name"
}
evex() { # FILE: count EVEX-encoded instructions (first opcode byte 0x62 in 64-bit code); fail unless 0
    local n
    n=$(objdump -d "$1" | grep -c -P '^\s+[0-9a-f]+:\t62 ' || true)
    echo "$(basename "$1") evex=$n" | tee -a "$RAW/evex.txt"
    [ "$n" = 0 ]
}
build() {
    local PI="-I$DEPS/parasail/install/include" PL="$DEPS/parasail/install/lib/libparasail.a"
    x "cp '$R'/bench/*.c '$R'/bench/*.cpp '$R'/bench/*.h '$R'/bench/*.py $S/"
    for cpu in znver4 x86_64_v3; do
        x "(cd '$R' && $ZIG build bench scan -Doptimize=ReleaseFast -Dcpu=$cpu --cache-dir $HOME/.cache/bitdp-zig --prefix $WORK/zig-$cpu)"
        x "cp $WORK/zig-$cpu/bin/bitdp-bench $B/bitdp-bench-$cpu && cp $WORK/zig-$cpu/bin/bitdp-scan $B/bitdp-scan-$cpu"
    done
    x "gcc -O3 -march=native $S/dp_check.c -o $B/dp_check $PI $PL -lm"
    x "gcc -O3 -march=native $S/baselines.c -o $B/baselines $PI -I$DEPS/edlib/edlib/include $DEPS/edlib/edlib.o $PL -lstdc++ -lm"
    x "gcc -O3 -march=native $S/ksw2_baseline.c $DEPS/ksw2/ksw2_extz2_sse.c $DEPS/ksw2/kalloc.c -I$DEPS/ksw2 -o $B/ksw2_baseline"
    x "gcc -O3 -march=native $S/wfa2_baseline.c -I$DEPS/WFA2-lib -o $B/wfa2_baseline $DEPS/WFA2-lib/lib/libwfa.a -lm -fopenmp"
    x "gcc -O3 -march=native $S/astarpa2_baseline.c -I$DEPS/astar-pairwise-aligner/astarpa-c -o $B/astarpa2_baseline $DEPS/astar-pairwise-aligner/target/release/libastarpa_c.a -lpthread -ldl -lm"
    x "g++ -O3 -march=native -std=c++17 $S/seqan_baseline.cpp -I$DEPS/seqan/include -o $B/seqan_avx512"
    x "g++ -O3 -march=x86-64-v3 -std=c++17 $S/seqan_baseline.cpp -I$DEPS/seqan/include -o $B/seqan_avx2"
    bgsa_build myers-avx2 BGSA_AVX2 x86-64-v3 -m 0 -a avx2
    bgsa_build bitpal-packed-avx2 BGSA_AVX2 x86-64-v3 -M 2 -I -3 -G -5 -t packed -a avx2
    bgsa_build bitpal-nonpacked-avx2 BGSA_AVX2 x86-64-v3 -M 2 -I -3 -G -5 -t non-packed -a avx2
    bgsa_build myers-avx512 BGSA_AVX512 native -m 0 -a avx512
    bgsa_build bitpal-packed-avx512 BGSA_AVX512 native
    x ": > '$RAW/evex.txt'"
    for f in bitdp-bench-x86_64_v3 bitdp-scan-x86_64_v3 seqan_avx2 bgsa-myers-avx2 bgsa-bitpal-packed-avx2 bgsa-bitpal-nonpacked-avx2; do
        x "evex $B/$f"
    done
}

# ---- ops: operation counts (deterministic, not timed): Table 1 and Fig. 2 from `zig build mve`, Tables 2-3 and
# Fig. 3 from `zig build ablation` (schemes 0-9: ablation, 100-123: scaling sweep)
ops() {
    local schemes="$(seq 0 9) $(seq 100 123)"
    [ $QUICK = 0 ] || schemes="0 100"
    x "(cd '$R' && $ZIG build mve -Doptimize=ReleaseFast --cache-dir $HOME/.cache/bitdp-zig) > '$RAW/mve.txt' 2>&1"
    x ": > '$RAW/ablation.txt'"
    for s in $schemes; do
        x "(cd '$R' && $ZIG build ablation -Dscheme=$s --cache-dir $HOME/.cache/bitdp-zig) >> '$RAW/ablation.txt' 2>&1"
    done
}

# ---- env: machine, toolchain and repository state
env_log() {
    x "{ date -Iseconds; git -C '$R' rev-parse HEAD; git -C '$R' status --short | head -20; uname -a; grep PRETTY /etc/os-release;
        lscpu | grep -E 'Model name|^CPU\\(s\\)|Thread|MHz|Flags' | cut -c1-300; free -h | head -2; gcc --version | head -1;
        g++ --version | head -1; java -version 2>&1 | head -1; python3 --version; $ZIG version; bwa 2>&1 | grep Version;
        samtools --version | head -1; cat '$DEPS/versions.txt'; pm get; } > '$RAW/env.txt' 2>&1 || true"
}

# ---- inputs: generated from seeds, plus ERR022075 reads mapped to NC_000913.3
inputs() {
    local G="python3 -I $S/gen.py" nreads=200000
    if [ $QUICK = 1 ]; then
        nreads=2000
        x "cd $IN && $G ava q150.txt s150.txt 4 50 150 4 && $G ava q1000.txt s1000.txt 2 20 1000 5"
        x "cd $IN && $G dna150.tsv 500 150 0.10 dna 1 && $G dna1000.tsv 50 1000 0.15 dna 2 && $G prot300.tsv 50 300 0.30 protein 3"
        x "cd $IN && $G dna10k5.tsv 4 10000 0.05 dna 6 && $G dna10k15.tsv 4 10000 0.15 dna 7 && $G dna100k5.tsv 1 30000 0.05 dna 8"
    else
        x "cd $IN && $G ava q150.txt s150.txt 64 12500 150 4 && $G ava q1000.txt s1000.txt 16 1250 1000 5"
        x "cd $IN && $G dna150.tsv 100000 150 0.10 dna 1 && $G dna1000.tsv 5000 1000 0.15 dna 2 && $G prot300.tsv 10000 300 0.30 protein 3"
        x "cd $IN && $G dna10k5.tsv 128 10000 0.05 dna 6 && $G dna10k15.tsv 128 10000 0.15 dna 7 && $G dna100k5.tsv 8 100000 0.05 dna 8"
    fi
    x "cd $IN && curl -sSf 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=fasta&retmode=text' > ref.fa"
    x "cd $IN && { curl -sS https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR022/ERR022075/ERR022075_1.fastq.gz | gunzip -c | head -n $((4 * nreads)) > r1.fq || true; }"
    x "cd $IN && bwa index ref.fa 2> bwa-index.log && bwa mem -t 8 -K 10000000 ref.fa r1.fq 2> bwa-mem.log | samtools view -F 0x904 - | python3 -I $S/reads.py ref.fa > reads.tsv 2> reads.log"
    x "cd $IN && sha256sum *.tsv q*.txt s*.txt ref.fa r1.fq > '$RAW/inputs.sha256' && cat reads.log >> '$RAW/inputs.sha256'"
}

# ---- ref: reference costs (sum, hash and per-pair dump) for every input and scheme a cell uses
ref() {
    x ": > '$RAW/ref.txt'"
    for key in edit:dna150.tsv bitpal:dna150.tsv tstv:dna150.tsv edit:dna1000.tsv bitpal:dna1000.tsv tstv:dna1000.tsv \
        blosum:prot300.tsv edit:reads.tsv bitpal:reads.tsv tstv:reads.tsv edit:dna10k5.tsv edit:dna10k15.tsv edit:dna100k5.tsv; do
        x "cd $IN && BENCH_DUMP=ref-${key%%:*}-${key#*:}.txt $B/dp_check ${key%%:*} ${key#*:} | sed 's/^/key=$key /' >> '$RAW/ref.txt'"
    done
    for s in edit bitpal tstv; do for n in 150 1000; do
        x "cd $IN && BENCH_DUMP=ref-$s-ava$n.txt $B/dp_check $s ava q$n.txt s$n.txt | sed 's/^/key=$s:ava$n /' >> '$RAW/ref.txt'"
    done; done
}

# ---- one cell: "cell ID REFKEY NAME=COMMAND ..." (COMMAND runs in $IN). A "!" before NAME marks a tool
# known to return wrong costs: its runs are kept and marked instead of stopping the session.
one() { # REP NAME=COMMAND, with id, key and want of the calling cell
    local rep=$1 name=${2%%=*} cmd=${2#*=} out got
    echo "+ [$id rep $rep] taskset -c $CORE $cmd"
    [ $DRY = 1 ] && return
    out=$(cd $IN && eval "taskset -c $CORE $cmd" 2>&1 | tail -n 1) || true
    echo "cell=$id rep=$rep tool=${name#!} $out" >> "$RAW/$LOG"
    got=$(sed -n 's/.* hash=\([0-9a-f]*\).*/\1/p' <<< "$out")
    [ "$got" != "$want" ] || return 0
    if [ "${name:0:1}" = "!" ]; then echo "cell=$id rep=$rep tool=${name#!} expected-mismatch" >> "$RAW/$LOG"; return 0; fi
    (cd $IN && eval "BENCH_DUMP=$WORK/dump.txt $cmd" > /dev/null 2>&1) || true
    echo "MISMATCH in $id, tool $name: hash '$got', reference $want ($key); last line: $out" >&2
    awk 'NR == FNR { r[FNR] = $0; next } $0 != r[FNR] { print "first differing pair: " FNR ", reference " r[FNR] ", tool " $0; exit }' \
        "$IN/ref-${key%%:*}-${key#*:}.txt" "$WORK/dump.txt" >&2 || true
    exit 1
}
cell() {
    local id=$1 key=$2 want
    shift 2
    local tools=("$@") n=$#
    want=$(grep "^key=$key " "$RAW/ref.txt" 2> /dev/null | sed -n 's/.* hash=\([0-9a-f]*\).*/\1/p') || true
    [ $DRY = 1 ] || [ -n "$want" ] || { echo "no reference for $key" >&2; exit 1; }
    for t in "${tools[@]}"; do one warm "$t"; done
    for ((rep = 1; rep <= REPS; rep++)); do
        for ((i = 0; i < n; i++)); do one $rep "${tools[(i + rep - 1) % n]}"; done
    done
}

bd512=$B/bitdp-bench-znver4 bd2=$B/bitdp-bench-x86_64_v3
bgsa_cmd() { # VARIANT {distance|score} LENGTH
    echo "$B/bgsa-$1 -q q$3.txt -d s$3.txt -f bgsa-$1.bin -N 1 > bgsa-$1.log && python3 -I $S/bgsa_check.py bgsa-$1 $2 bgsa-$1.bin q$3.txt s$3.txt bgsa-$1.log"
}

t4() { # all against all, BGSA's setting
    for n in 150 1000; do
        cell t4-edit-$n edit:ava$n "bgsa-myers-avx2=$(bgsa_cmd myers-avx2 distance $n)" "bitdp-avx2=$bd2 edit ava q$n.txt s$n.txt" \
            "bgsa-myers-avx512=$(bgsa_cmd myers-avx512 distance $n)" "bitdp-avx512=$bd512 edit ava q$n.txt s$n.txt"
        cell t4-bitpal-$n bitpal:ava$n "bgsa-bitpal-packed-avx2=$(bgsa_cmd bitpal-packed-avx2 score $n)" \
            "bgsa-bitpal-nonpacked-avx2=$(bgsa_cmd bitpal-nonpacked-avx2 score $n)" "bitdp-avx2=$bd2 bitpal ava q$n.txt s$n.txt" \
            "bgsa-bitpal-packed-avx512=$(bgsa_cmd bitpal-packed-avx512 score $n)" "bitdp-avx512=$bd512 bitpal ava q$n.txt s$n.txt"
        cell t4-tstv-$n tstv:ava$n "bitdp-avx2=$bd2 tstv ava q$n.txt s$n.txt" "bitdp-avx512=$bd512 tstv ava q$n.txt s$n.txt"
    done
}

pairs_cell() { # ID SCHEME FILE [EXTRA TOOLS]: the tools that run on independent pairs of this scheme
    local id=$1 s=$2 f=$3
    shift 3
    local t=("bitdp-1pair=$bd512 $s $f" "bitdp-avx2=$bd2 $s $f batch" "bitdp-avx512=$bd512 $s $f batch"
        "seqan-avx2=$B/seqan_avx2 $s $f" "seqan-avx512=$B/seqan_avx512 $s $f" "parasail-scan16=$B/baselines scan16 $s $f")
    [ $s = blosum ] || t+=("ksw2=$B/ksw2_baseline $s $f")
    [ $s != edit ] || t+=("edlib=$B/baselines edlib edit $f" "astarpa2-simple=$B/astarpa2_baseline simple $f"
        "astarpa2-full=$B/astarpa2_baseline full $f")
    [ $s != edit ] && [ $s != bitpal ] || t+=("wfa2=$B/wfa2_baseline $s $f")
    cell $id $s:$f "${t[@]}" "$@"
}

t5() { # independent pairs; parasail's diag16 and striped16 as before (striped16 is wrong for linear gaps)
    for s in edit bitpal tstv; do for f in dna150.tsv dna1000.tsv; do
        pairs_cell t5-$s-${f%.tsv} $s $f "parasail-diag16=$B/baselines diag16 $s $f" "!parasail-striped16=$B/baselines striped16 $s $f"
    done; done
    pairs_cell t5-blosum-prot300 blosum prot300.tsv "parasail-diag16=$B/baselines diag16 blosum prot300.tsv" \
        "!parasail-striped16=$B/baselines striped16 blosum prot300.tsv"
}

t6() { for s in edit bitpal tstv; do pairs_cell t6-$s $s reads.tsv; done; }

long() { # long pairs, unit costs; at 100 kbp bitdp and the tools that are not quadratic per pair
    pairs_cell long-edit-10k5 edit dna10k5.tsv
    pairs_cell long-edit-10k15 edit dna10k15.tsv
    local f=dna100k5.tsv
    cell long-edit-100k5 edit:$f "bitdp-1pair=$bd512 edit $f" "bitdp-avx2=$bd2 edit $f batch" "bitdp-avx512=$bd512 edit $f batch" \
        "edlib=$B/baselines edlib edit $f" "wfa2=$B/wfa2_baseline edit $f" "astarpa2-simple=$B/astarpa2_baseline simple $f" \
        "astarpa2-full=$B/astarpa2_baseline full $f"
}

scan() { # genome scan: hit sets checked once against the scalar DP per configuration, then timed runs
    local gs="32 48 96"
    [ $QUICK = 0 ] || gs=32
    for cpu in znver4 x86_64_v3; do for g in $gs; do for p in "" packed; do
        x "cd $IN && taskset -c $CORE $B/bitdp-scan-$cpu ref.fa $g 4 verify $p 2>&1 | sed 's/^/cell=scan-$cpu-$g$p rep=verify /' >> '$RAW/scan.txt'"
        for ((rep = 1; rep <= REPS; rep++)); do
            x "cd $IN && taskset -c $CORE $B/bitdp-scan-$cpu ref.fa $g 4 $p 2>&1 | tail -n 1 | sed 's/^/cell=scan-$cpu-$g$p rep=$rep /' >> '$RAW/scan.txt'"
        done
    done; done; done
}

cc_project() { # DIR SCHEME ROOT: a Zig project with one derivation, Kernel(SCHEME) of the bitdp module at ROOT
    local d=$1
    echo "+ write $d/build.zig and $d/cc.zig for $2"
    [ $DRY = 1 ] && return
    mkdir -p "$d"
    cat > "$d/build.zig" << 'EOF'
const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});
    const bitdp = b.createModule(.{ .root_source_file = .{ .cwd_relative = "@ROOT@" }, .target = target, .optimize = optimize });
    const exe = b.addExecutable(.{ .name = "cc", .root_module = b.createModule(.{
        .root_source_file = b.path("cc.zig"),
        .target = target,
        .optimize = optimize,
        .imports = &.{.{ .name = "bitdp", .module = bitdp }},
    }) });
    b.installArtifact(exe);
}
EOF
    cat > "$d/cc.zig" << 'EOF'
const std = @import("std");
const bitdp = @import("bitdp");

pub fn main() void {
    std.debug.print("ops={d}\n", .{bitdp.Kernel(bitdp.schemes.@SCHEME@).ops});
}
EOF
    sed -i "s|@ROOT@|$3|" "$d/build.zig"
    sed -i "s|@SCHEME@|$2|" "$d/cc.zig"
}

compile() { # compile cost of one derivation, Debug (mostly the derivation) and ReleaseFast (what a user builds):
    # CREPS builds per scheme, OS and mode, each with a fresh local cache, after one warm-up per OS and mode
    local zw=${ZIG_WIN:-$(command -v zig.exe || true)} wd one
    wd=$(wslpath "$(cmd.exe /c echo %TEMP% 2> /dev/null | tr -d '\r')")/bitdp-cc
    x "rm -rf $WORK/srczig && cp -r '$R/src' $WORK/srczig" # WSL builds read the sources from the Linux file system
    for mode in Debug ReleaseFast; do
        for s in edit ts_tv 'blosum62Linear(4)'; do
            cc_project $WORK/cc "$s" "$WORK/srczig/root.zig"
            [ -z "$zw" ] || cc_project "$wd" "$s" "$(wslpath -m "$R/src/root.zig")"
            one=$(seq 1 $CREPS)
            [ "$s" != edit ] || one="warm $one"
            for rep in $one; do
                x "(cd $WORK/cc && rm -rf .zig-cache zig-out && $ZIG build -Doptimize=$mode --summary all 2>&1 | grep 'compile exe' | sed 's/^/scheme=$s os=wsl rep=$rep /' >> '$RAW/compile.txt')"
                [ -z "$zw" ] || x "(cd '$wd' && rm -rf .zig-cache zig-out && '$zw' build -Doptimize=$mode --summary all 2>&1 | tr -d '\\r' | grep 'compile exe' | sed 's/^/scheme=$s os=windows rep=$rep /' >> '$RAW/compile.txt')"
            done
        done
    done
}

session() { # power mode for the measurement steps; the mode found at the start is restored on exit
    if [ $DRY = 1 ]; then echo "+ powermode.ps1 set efficiency; sampler on; restore on exit"; return; fi
    local now
    now=$(pm get)
    # every table was measured on mains power; on battery Windows applies its DC settings
    case $now in *ac=True*) ;; *) [ "${ALLOW_DC:-0}" = 1 ] || { echo "on battery ($now): plug in, or set ALLOW_DC=1" >&2; exit 1; } ;; esac
    ORIG=$(sed -n 's/.*overlay=\([^ ]*\).*/\1/p' <<< "$now")
    echo "$(date -Iseconds) start $(pm get)" >> "$RAW/power.txt"
    echo "$(date -Iseconds) set $(pm set efficiency)" >> "$RAW/power.txt"
    rm -f "$WORK/stop"
    pm sample "$(wslpath -w "$RAW/power-samples.txt")" "$(wslpath -w "$WORK/stop")" &
    SAMPLER=$!
    trap 'touch "$WORK/stop"; wait $SAMPLER || true; echo "$(date -Iseconds) restore $(pm set $ORIG)" >> "$RAW/power.txt"' EXIT
    trap 'exit 130' INT TERM
}

# ---- insn: retired instructions and cycles of one run per bitdp configuration (perfcount.c, hardware
# counters), with the column words the run executes (wordcols.py); counts, not timings, so outside the session
icount() { # CELL SCHEME CPU ARGS...: one counted run of bitdp-bench-CPU SCHEME ARGS
    local cell=$1 s=$2 cpu=$3 lanes=8 tool=bitdp-avx512
    shift 3
    [ $cpu = znver4 ] || { lanes=4 tool=bitdp-avx2; }
    x "cd $IN && echo \"cell=$cell-$s tool=$tool \$(python3 -I $S/wordcols.py $lanes $*) \$(taskset -c $CORE $B/perfcount $B/bitdp-bench-$cpu $s $* 2>&1 | tr '\n' ' ')\" >> '$RAW/insn.txt'"
}
insn() {
    x "cp '$R/bench/perfcount.c' '$R/bench/wordcols.py' $S/ && gcc -O2 -Wall $S/perfcount.c -o $B/perfcount"
    x ": > '$RAW/insn.txt'"
    for cpu in znver4 x86_64_v3; do
        for s in edit bitpal tstv; do
            icount insn-ava1000 $s $cpu ava q1000.txt s1000.txt
            icount insn-dna1000 $s $cpu dna1000.tsv batch
            icount insn-dna150 $s $cpu dna150.tsv batch
        done
        icount insn-prot300 blosum $cpu prot300.tsv batch
    done
    x "objdump -d $B/bitdp-bench-znver4 $B/bitdp-bench-x86_64_v3 | grep -c -P '\tvp?gather' | sed 's/^/gather instructions in both drivers: /' >> '$RAW/insn.txt' || true"
    x "objdump -d $B/bitdp-bench-znver4 | grep -c -P '\tvpternlog' | sed 's/^/vpternlog instructions in the AVX-512 driver: /' >> '$RAW/insn.txt' || true"
}

started=0
for step in "${STEPS[@]}"; do
    case $step in t4 | t5 | t6 | long | scan | compile | robust) [ $started = 1 ] || { session; started=1; } ;; esac
    echo "== $step $(date -Iseconds)"
    case $step in
    deps) deps ;; build) build ;; ops) ops ;; env) env_log ;; inputs) inputs ;; ref) ref ;;
    t4) t4 ;; t5) t5 ;; t6) t6 ;; long) long ;; scan) scan ;; compile) compile ;; insn) insn ;;
    robust) # Table 4 in "Best performance", RREPS repetitions: a robustness check of its ratios
        [ $DRY = 1 ] || echo "$(date -Iseconds) set $(pm set performance)" >> "$RAW/power.txt"
        REPS=$RREPS LOG=runs-bestperf.txt t4
        [ $DRY = 1 ] || echo "$(date -Iseconds) set $(pm set efficiency)" >> "$RAW/power.txt" ;;
    summary) # optional inputs quoted here: the repository path may hold spaces
        opt=""
        [ ! -f "$RAW/insn.txt" ] || opt="'$RAW/insn.txt'"
        x "python3 -I '$R/bench/summary.py' '$RAW/runs.txt' '$RAW/scan.txt' '$RAW/compile.txt' $opt > '$RAW/summary.csv'"
        opt=""
        [ ! -f "$RAW/runs-bestperf.txt" ] || { x "python3 -I '$R/bench/summary.py' '$RAW/runs-bestperf.txt' > '$RAW/summary-bestperf.csv'"; opt="'$RAW/summary-bestperf.csv'"; }
        x "mkdir -p '$GEN' && python3 -I '$R/bench/texmacros.py' '$RAW/summary.csv' $opt > '$GEN/numbers.tex'" ;;
    *) echo "unknown step $step" >&2; exit 2 ;;
    esac
done
echo "== done $(date -Iseconds)"
