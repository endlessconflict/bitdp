// ksw2 baseline for the bitdp benchmark: global alignment score of every pair
// in a "pattern<TAB>text" file with ksw_extz2_sse (the SIMD DP of minimap2).
// Linear gaps are affine gaps with open cost 0.
//
// build (inside a ksw2 checkout, https://github.com/lh3/ksw2):
//   gcc -O3 -march=native ksw2_baseline.c ksw2_extz2_sse.c kalloc.c -o ksw2_baseline
// usage: ksw2_baseline {edit|bitpal|tstv} FILE

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "ksw2.h"
#include "pairhash.h"

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

static uint8_t code(char c) {
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default: return 4;
    }
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: %s {edit|bitpal|tstv} FILE\n", argv[0]);
        return 2;
    }
    int8_t mat[25];
    int match = 0, mismatch = -1, gap = 1;
    if (!strcmp(argv[1], "bitpal")) match = 2, mismatch = -3, gap = 5;
    if (!strcmp(argv[1], "tstv")) mismatch = -2, gap = 2;
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++) {
            int s = (i == j && i < 4) ? match : mismatch;
            // transitions A<->G (0,2) and C<->T (1,3) cost 1 in tstv
            if (!strcmp(argv[1], "tstv") && i != j && i < 4 && j < 4 && (i ^ j) == 2) s = -1;
            mat[i * 5 + j] = (int8_t)s;
        }

    // without KSW_EZ_GENERIC_SC ksw2 reads only mat[0] (match) and mat[1] (mismatch)
    int flag = KSW_EZ_SCORE_ONLY | (!strcmp(argv[1], "tstv") ? KSW_EZ_GENERIC_SC : 0);
    ksw_extz_t ez;
    memset(&ez, 0, sizeof ez);
    pairset in = pairs_load(argv[2]);
    if (!in.p) return 1;
    for (long long i = 0; i < in.count; i++) { // nucleotide codes, in place
        for (int j = 0; j < in.m[i]; j++) in.p[i][j] = (char)code(in.p[i][j]);
        for (int j = 0; j < in.n[i]; j++) in.t[i][j] = (char)code(in.t[i][j]);
    }
    long long *cost = malloc((in.count + 1) * sizeof *cost);
    long long total = 0, rounds = 0;
    double sec, min_sec = bench_min_seconds(), t0 = now();
    do {
        for (long long i = 0; i < in.count; i++) {
            ksw_extz2_sse(0, in.m[i], (uint8_t *)in.p[i], in.n[i], (uint8_t *)in.t[i], 5, mat, 0, (int8_t)gap, -1, -1, 0,
                          flag, &ez);
            cost[i] = -(long long)ez.score;
        }
        rounds++;
    } while ((sec = now() - t0) < min_sec);
    pairhash ph;
    pairhash_init(&ph);
    for (long long i = 0; i < in.count; i++) total += cost[i], pairhash_add(&ph, cost[i]);
    printf("ksw2_extz2_sse %s pairs=%lld cost_sum=%lld hash=%016llx seconds=%.4f gcups=%.3f rounds=%lld\n", argv[1], in.count, total,
           pairhash_end(&ph), sec, in.cells * rounds / sec * 1e-9, rounds);
    return 0;
}
