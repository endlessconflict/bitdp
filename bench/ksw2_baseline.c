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

    FILE *f = fopen(argv[2], "rb");
    if (!f) return 1;
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *buf = malloc(size + 1);
    if (fread(buf, 1, size, f) != (size_t)size) return 1;
    buf[size] = 0;
    fclose(f);

    uint8_t *q = malloc(1 << 20), *t = malloc(1 << 20);
    long long total = 0, cells = 0, pairs = 0;
    double busy = 0;
    // without KSW_EZ_GENERIC_SC ksw2 reads only mat[0] (match) and mat[1] (mismatch)
    int flag = KSW_EZ_SCORE_ONLY | (!strcmp(argv[1], "tstv") ? KSW_EZ_GENERIC_SC : 0);
    ksw_extz_t ez;
    memset(&ez, 0, sizeof ez);
    for (char *line = strtok(buf, "\n"); line; line = strtok(NULL, "\n")) {
        char *tab = strchr(line, '\t');
        if (!tab) continue;
        *tab = 0;
        int m = strlen(line), n = strlen(tab + 1);
        for (int i = 0; i < m; i++) q[i] = code(line[i]);
        for (int i = 0; i < n; i++) t[i] = code(tab[1 + i]);
        double t0 = now();
        ksw_extz2_sse(0, m, q, n, t, 5, mat, 0, (int8_t)gap, -1, -1, 0, flag, &ez);
        busy += now() - t0;
        total -= ez.score;
        cells += (long long)m * n;
        pairs++;
    }
    printf("ksw2_extz2_sse %s pairs=%lld cost_sum=%lld seconds=%.4f gcups=%.3f\n", argv[1], pairs, total, busy, cells / busy * 1e-9);
    return 0;
}
