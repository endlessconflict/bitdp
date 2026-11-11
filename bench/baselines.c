// Baselines for the bitdp benchmark: global alignment score of every pair in
// a "pattern<TAB>text" file, with parasail (SIMD) or edlib (Myers).
//
// build: gcc -O3 -march=native baselines.c -o baselines -lparasail -ledlib -lstdc++
// usage: baselines {edlib|striped16|scan16|diag16|striped_sat|scan_sat} {edit|bitpal|tstv|blosum} FILE [K]
// K (edlib only): distance bound passed to edlib; pairs above it count as failures.
//
// Prints the sum of alignment costs (negated parasail scores) and their per-pair hash (pairhash.h), so
// the result can be checked against bitdp's pair by pair, and the time of the whole loop over the pairs (pairhash.h).

#include <edlib.h>
#include <parasail.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "pairhash.h"

typedef parasail_result_t *(*nw_fn)(const char *, int, const char *, int, int, int, const parasail_matrix_t *);

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

static int is_purine(char c) { return c == 'A' || c == 'G'; }

int main(int argc, char **argv) {
    if (argc != 4 && argc != 5) {
        fprintf(stderr, "usage: %s TOOL SCHEME FILE [K]\n", argv[0]);
        return 2;
    }
    const char *tool = argv[1], *scheme = argv[2];
    int k = argc == 5 ? atoi(argv[4]) : -1;

    // Scores and linear gap penalty (parasail: gap of length k costs open + (k-1) * extend).
    parasail_matrix_t *matrix = NULL;
    int gap = 1;
    if (!strcmp(scheme, "edit")) {
        matrix = parasail_matrix_create("ACGT", 0, -1);
        gap = 1;
    } else if (!strcmp(scheme, "bitpal")) {
        matrix = parasail_matrix_create("ACGT", 2, -3);
        gap = 5;
    } else if (!strcmp(scheme, "tstv")) {
        matrix = parasail_matrix_create("ACGT", 0, -2);
        const char *a = "ACGT";
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                if (i != j && is_purine(a[i]) == is_purine(a[j])) parasail_matrix_set_value(matrix, i, j, -1);
        gap = 2;
    } else if (!strcmp(scheme, "blosum")) {
        matrix = parasail_matrix_copy(parasail_matrix_lookup("blosum62"));
        gap = 4;
    } else {
        fprintf(stderr, "unknown scheme %s\n", scheme);
        return 2;
    }

    nw_fn fn = NULL;
    if (!strcmp(tool, "striped16")) fn = parasail_nw_striped_16;
    else if (!strcmp(tool, "scan16")) fn = parasail_nw_scan_16;
    else if (!strcmp(tool, "diag16")) fn = parasail_nw_diag_16;
    else if (!strcmp(tool, "striped_sat")) fn = parasail_nw_striped_sat;
    else if (!strcmp(tool, "scan_sat")) fn = parasail_nw_scan_sat;
    else if (strcmp(tool, "edlib")) {
        fprintf(stderr, "unknown tool %s\n", tool);
        return 2;
    }

    pairset in = pairs_load(argv[3]);
    if (!in.p) return 1;
    long long *cost = malloc((in.count + 1) * sizeof *cost);
    long long total = 0, rounds = 0;
    double sec, min_sec = bench_min_seconds(), t0 = now();
    do {
        for (long long i = 0; i < in.count; i++) {
            if (fn) {
                parasail_result_t *r = fn(in.p[i], in.m[i], in.t[i], in.n[i], gap, gap, matrix);
                cost[i] = -parasail_result_get_score(r);
                parasail_result_free(r);
            } else {
                EdlibAlignResult r = edlibAlign(in.p[i], in.m[i], in.t[i], in.n[i],
                                                edlibNewAlignConfig(k, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
                cost[i] = r.editDistance; // -1 above K
                edlibFreeAlignResult(r);
            }
        }
        rounds++;
    } while ((sec = now() - t0) < min_sec);
    pairhash ph;
    pairhash_init(&ph);
    long long over_k = 0;
    for (long long i = 0; i < in.count; i++) {
        if (cost[i] < 0 && !fn) over_k++;
        else total += cost[i];
        pairhash_add(&ph, cost[i]);
    }
    printf("%s %s pairs=%lld cost_sum=%lld hash=%016llx seconds=%.4f gcups=%.3f rounds=%lld", tool, scheme, in.count, total,
           pairhash_end(&ph), sec, in.cells * rounds / sec * 1e-9, rounds);
    if (k >= 0) printf(" k=%d over_k=%lld", k, over_k);
    printf("\n");
    parasail_matrix_free(matrix);
    return 0;
}
