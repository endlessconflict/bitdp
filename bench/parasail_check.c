// Per-pair check of parasail's global (nw) kernels against a plain DP, for the
// bitdp benchmark files: reports, per kernel, how many pairs differ from the
// optimum and how many results carry parasail's saturation flag, and prints the
// first differing pair. Same schemes and gap arguments as baselines.c.
//
// build: gcc -O3 -march=native parasail_check.c -o parasail_check -lparasail
// usage: parasail_check {edit|bitpal|tstv|blosum} FILE [list]
// With "list", every differing pair is printed as: diff KERNEL PAIR M N PARASAIL DP

#include <parasail.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef parasail_result_t *(*nw_fn)(const char *, int, const char *, int, int, int, const parasail_matrix_t *);

static const struct { const char *name; nw_fn fn; } kernels[] = {
    {"striped16", parasail_nw_striped_16},
    {"striped32", parasail_nw_striped_32},
    {"scan16", parasail_nw_scan_16},
    {"scan32", parasail_nw_scan_32},
};
#define NK (int)(sizeof kernels / sizeof kernels[0])

static int is_purine(char c) { return c == 'A' || c == 'G'; }

// Score of the best global alignment with a linear gap: every gap position costs `gap`.
static long dp(const char *p, int m, const char *t, int n, const parasail_matrix_t *mat, int gap, long *row) {
    for (int j = 0; j <= n; j++) row[j] = -(long)j * gap;
    for (int i = 1; i <= m; i++) {
        long diag = row[0];
        row[0] = -(long)i * gap;
        for (int j = 1; j <= n; j++) {
            long s = diag + mat->matrix[mat->mapper[(unsigned char)p[i - 1]] * mat->size + mat->mapper[(unsigned char)t[j - 1]]];
            long up = row[j] - gap, left = row[j - 1] - gap;
            if (up > s) s = up;
            if (left > s) s = left;
            diag = row[j];
            row[j] = s;
        }
    }
    return row[n];
}

int main(int argc, char **argv) {
    if (argc != 3 && argc != 4) {
        fprintf(stderr, "usage: %s {edit|bitpal|tstv|blosum} FILE [list]\n", argv[0]);
        return 2;
    }
    int list = argc == 4 && !strcmp(argv[3], "list");
    const char *scheme = argv[1];
    parasail_matrix_t *matrix;
    int gap;
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

    // Gap semantics probe: parasail is called with open = extend = gap, as in baselines.c.
    // A pattern four longer than the text needs one gap of length 4; a linear gap scores -4*gap.
    long row[8192];
    {
        const char *p = "ACGTACGT", *t = "ACGT";
        parasail_result_t *r = parasail_nw_scan_32(p, 8, t, 4, gap, gap, matrix);
        printf("gap probe: parasail scan32 %d, linear-gap DP %ld\n", parasail_result_get_score(r), dp(p, 8, t, 4, matrix, gap, row));
        parasail_result_free(r);
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

    long pairs = 0, wrong[NK] = {0}, saturated[NK] = {0}, below[NK] = {0};
    int shown[NK] = {0};
    for (char *line = strtok(buf, "\n"); line; line = strtok(NULL, "\n")) {
        char *tab = strchr(line, '\t');
        if (!tab) continue;
        *tab = 0;
        const char *p = line, *t = tab + 1;
        int m = strlen(p), n = strlen(t);
        if (n + 1 > (int)(sizeof row / sizeof row[0])) return 3;
        long want = dp(p, m, t, n, matrix, gap, row);
        pairs++;
        for (int k = 0; k < NK; k++) {
            parasail_result_t *r = kernels[k].fn(p, m, t, n, gap, gap, matrix);
            long got = parasail_result_get_score(r);
            if (parasail_result_is_saturated(r)) saturated[k]++;
            if (got != want) {
                wrong[k]++;
                if (got < want) below[k]++;
                if (list) printf("diff %s %ld %d %d %ld %ld\n", kernels[k].name, pairs, m, n, got, want);
                if (!shown[k]++)
                    printf("first %s difference: pair %ld (m=%d, n=%d): %s %ld, DP %ld, saturated %d\n  %s\n  %s\n",
                           kernels[k].name, pairs, m, n, kernels[k].name, got, want, parasail_result_is_saturated(r), p, t);
            }
            parasail_result_free(r);
        }
    }
    for (int k = 0; k < NK; k++)
        printf("%-9s %s %s: %ld of %ld pairs differ from the DP (%ld below it), %ld saturated\n",
               kernels[k].name, scheme, argv[2], wrong[k], pairs, below[k], saturated[k]);
    parasail_matrix_free(matrix);
    free(buf);
    return 0;
}
