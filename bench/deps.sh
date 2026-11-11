#!/bin/sh
# Fetch the pinned competitors and toolchain into $DEPS (default $HOME/bitdp-deps) and build them.
# Needs git, curl, make, autoconf, automake, libtool, gcc/g++ (the paper used 15.2.0), java (BGSA's
# generator) and rustup (A*PA2 pins its own nightly). Every checkout is a fixed commit;
# $DEPS/versions.txt records them, and each build log lands next to its sources.
#
# usage: bench/deps.sh [fetch|build|all]   (default all)
set -eu
command -v cargo > /dev/null || PATH="$HOME/.cargo/bin:$PATH"
DEPS=${DEPS:-$HOME/bitdp-deps}
STEP=${1:-all}
mkdir -p "$DEPS"
cd "$DEPS"

pin() { # dir url rev
    [ -d "$1/.git" ] || git clone -q "$2" "$1"
    git -C "$1" fetch -q --tags origin
    git -C "$1" -c advice.detachedHead=false checkout -q "$3"
    git -C "$1" submodule -q update --init --recursive
    printf '%s %s %s\n' "$1" "$3" "$(git -C "$1" rev-parse HEAD)" >> versions.txt
}

if [ "$STEP" = fetch ] || [ "$STEP" = all ]; then
    : > versions.txt
    pin edlib https://github.com/Martinsos/edlib v1.2.7
    pin parasail https://github.com/jeffdaily/parasail v2.6.2
    pin ksw2 https://github.com/lh3/ksw2 289609b
    pin BGSA https://github.com/sdu-hpcl/BGSA dffcf9e
    pin WFA2-lib https://github.com/smarco/WFA2-lib v2.3.6
    pin seqan https://github.com/seqan/seqan seqan-v2.5.3
    pin astar-pairwise-aligner https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner bf2e14e0cbc3a9a03600dcda0641d7f89e401e63
    # Zig 0.16.0 (2026-04-13), checksum from ziglang.org/download/index.json
    if [ ! -x zig/zig ]; then
        curl -sSfL -o zig.tar.xz https://ziglang.org/download/0.16.0/zig-x86_64-linux-0.16.0.tar.xz
        echo "70e49664a74374b48b51e6f3fdfbf437f6395d42509050588bd49abe52ba3d00  zig.tar.xz" | sha256sum -c -
        mkdir -p zig && tar -xJf zig.tar.xz -C zig --strip-components=1 && rm zig.tar.xz
    fi
    printf 'zig %s\n' "$(zig/zig version)" >> versions.txt
fi

if [ "$STEP" = build ] || [ "$STEP" = all ]; then
    J=$(nproc)
    # edlib: one translation unit, linked into the driver.
    g++ -O3 -march=native -c edlib/edlib/src/edlib.cpp -I edlib/edlib/include -o edlib/edlib.o
    # parasail: autotools (its configure script selects the instruction sets; this version has no AVX-512 alignment kernels).
    (cd parasail && { [ -x configure ] || autoreconf -fi; } &&
        ./configure -q --enable-static --disable-shared --prefix="$DEPS/parasail/install" CFLAGS="-O3 -march=native" &&
        make -j"$J" && make install) > parasail/build.log 2>&1
    # WFA2-lib: its `all` target adds -O3 -march=native.
    make -C WFA2-lib -j"$J" all BUILD_TOOLS=0 BUILD_EXAMPLES=0 > WFA2-lib/build.log 2>&1
    # A*PA2: the C library of its repository, with the nightly its rust-toolchain.toml pins.
    (cd astar-pairwise-aligner && cargo build --release -p astarpa-c) > astar-pairwise-aligner/build.log 2>&1
    # ksw2, BGSA and SeqAn (header-only) are compiled with their drivers by bench/reproduce.sh.
    printf '%s\n' "gcc $(gcc -dumpfullversion)" "cargo $(cargo --version)" >> versions.txt
fi
