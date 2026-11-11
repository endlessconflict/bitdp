// A*PA2 baseline: exact global edit distance of every pair in a "pattern<TAB>text" file through the C API
// of A*PA2 (astar-pairwise-aligner, crate astarpa-c), in its "simple" variant (band doubling with
// bit-packed blocks) or its "full" variant (with the A* heuristics). Unit costs only. A*PA2 has no
// score-only mode: every call also builds a CIGAR, which is freed inside the timed region.
//
// build: gcc -O3 -march=native astarpa2_baseline.c -I$DEPS/astar-pairwise-aligner/astarpa-c -o astarpa2_baseline \
//            $DEPS/astar-pairwise-aligner/target/release/libastarpa_c.a -lpthread -ldl -lm
// usage: astarpa2_baseline {simple|full} FILE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "astarpa.h"
#include "pairhash.h"

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
    if (argc != 3 || (strcmp(argv[1], "simple") && strcmp(argv[1], "full"))) {
        fprintf(stderr, "usage: %s {simple|full} FILE\n", argv[0]);
        return 2;
    }
    int full = !strcmp(argv[1], "full");

    pairset in = pairs_load(argv[2]);
    if (!in.p) return 1;
    long long *cost = malloc((in.count + 1) * sizeof *cost);
    long long total = 0, rounds = 0;
    double sec, min_sec = bench_min_seconds(), t0 = now();
    do {
        for (long long i = 0; i < in.count; i++) {
            uint8_t *cigar;
            uintptr_t cigar_len;
            cost[i] = full ? astarpa2_full((const uint8_t *)in.p[i], in.m[i], (const uint8_t *)in.t[i], in.n[i], &cigar, &cigar_len)
                           : astarpa2_simple((const uint8_t *)in.p[i], in.m[i], (const uint8_t *)in.t[i], in.n[i], &cigar, &cigar_len);
            astarpa_free_cigar(cigar);
        }
        rounds++;
    } while ((sec = now() - t0) < min_sec);
    pairhash ph;
    pairhash_init(&ph);
    for (long long i = 0; i < in.count; i++) total += cost[i], pairhash_add(&ph, cost[i]);
    printf("astarpa2-%s edit pairs=%lld cost_sum=%lld hash=%016llx seconds=%.4f gcups=%.3f rounds=%lld\n", argv[1], in.count, total,
           pairhash_end(&ph), sec, in.cells * rounds / sec * 1e-9, rounds);
    return 0;
}
