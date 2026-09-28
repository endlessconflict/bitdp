// Baselines for the bitdp benchmark: global alignment score of every pair in
// a "pattern<TAB>text" file, with parasail (SIMD) or edlib (Myers).
//
// build: gcc -O3 -march=native baselines.c -o baselines -lparasail -ledlib -lstdc++
// usage: baselines {edlib|striped16|scan16|diag16|striped_sat|scan_sat} {edit|bitpal|tstv|blosum} FILE
//
// Prints the sum of alignment costs (negated parasail scores), so the result
// can be checked against bitdp's, and the time spent aligning.

#include <edlib.h>
#include <parasail.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef parasail_result_t *(*nw_fn)(const char *, int, const char *, int, int, int, const parasail_matrix_t *);

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

static int is_purine(char c) { return c == 'A' || c == 'G'; }

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr, "usage: %s TOOL SCHEME FILE\n", argv[0]);
        return 2;
    }
    const char *tool = argv[1], *scheme = argv[2];

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

    FILE *f = fopen(argv[3], "rb");
    if (!f) return 1;
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *buf = malloc(size + 1);
    if (fread(buf, 1, size, f) != (size_t)size) return 1;
    buf[size] = 0;
    fclose(f);

    long long total = 0, cells = 0, pairs = 0;
    double busy = 0;
    for (char *line = strtok(buf, "\n"); line; line = strtok(NULL, "\n")) {
        char *tab = strchr(line, '\t');
        if (!tab) continue;
        *tab = 0;
        const char *p = line, *t = tab + 1;
        int m = strlen(p), n = strlen(t);
        double t0 = now();
        if (fn) {
            parasail_result_t *r = fn(p, m, t, n, gap, gap, matrix);
            total -= parasail_result_get_score(r);
            parasail_result_free(r);
        } else {
            EdlibAlignResult r = edlibAlign(p, m, t, n, edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
            total += r.editDistance;
            edlibFreeAlignResult(r);
        }
        busy += now() - t0;
        cells += (long long)m * n;
        pairs++;
    }
    printf("%s %s pairs=%lld cost_sum=%lld seconds=%.4f gcups=%.3f\n", tool, scheme, pairs, total, busy, cells / busy * 1e-9);
    parasail_matrix_free(matrix);
    free(buf);
    return 0;
}
