// Independent scalar DP, the reference of every benchmark run: global alignment cost (the negated score)
// of every pair with linear gaps, summed and hashed as in pairhash.h. Shares no code with bitdp or the
// aligners under test; the BLOSUM62 values are read from parasail's table (data only).
//
// usage: dp_check {edit|bitpal|tstv|blosum} FILE                 pairs, one "pattern<TAB>text" per line
//        dp_check {edit|bitpal|tstv|blosum} ava QUERIES SUBJECTS  every query against every subject
// build: gcc -O3 -march=native dp_check.c -o dp_check -I$DEPS/parasail/install/include \
//            $DEPS/parasail/install/lib/libparasail.a -lm
#include <parasail.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pairhash.h"

static int S[256][256], gap; // substitution scores and the (negative) gap score

static char **read_lines(const char *path, int *n) {
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(1); }
    char **v = NULL, *buf = malloc(1 << 22);
    *n = 0;
    while (fgets(buf, 1 << 22, f)) {
        size_t l = strcspn(buf, "\r\n");
        if (!l) continue;
        buf[l] = 0;
        v = realloc(v, (*n + 1) * sizeof *v);
        v[(*n)++] = strdup(buf);
    }
    fclose(f);
    free(buf);
    return v;
}

static long long cost(const char *p, const char *t) {
    static int *row;
    static size_t cap;
    size_t m = strlen(p), n = strlen(t);
    if (m + 1 > cap) row = realloc(row, (cap = m + 1) * sizeof *row);
    for (size_t i = 0; i <= m; i++) row[i] = (int)i * gap; // column 0
    for (size_t j = 1; j <= n; j++) {
        int diag = row[0];
        row[0] = (int)j * gap;
        for (size_t i = 1; i <= m; i++) {
            int best = diag + S[(unsigned char)p[i - 1]][(unsigned char)t[j - 1]];
            if (row[i] + gap > best) best = row[i] + gap;         // gap in the pattern
            if (row[i - 1] + gap > best) best = row[i - 1] + gap; // gap in the text
            diag = row[i];
            row[i] = best;
        }
    }
    return -(long long)row[m];
}

static int purine(int c) { return c == 'A' || c == 'G'; }

int main(int argc, char **argv) {
    if (argc != 3 && !(argc == 5 && !strcmp(argv[2], "ava"))) {
        fprintf(stderr, "usage: %s SCHEME FILE | %s SCHEME ava QUERIES SUBJECTS\n", argv[0], argv[0]);
        return 2;
    }
    const char *scheme = argv[1];
    for (int a = 0; a < 256; a++)
        for (int b = 0; b < 256; b++) {
            if (!strcmp(scheme, "edit")) S[a][b] = a == b ? 0 : -1, gap = -1;
            else if (!strcmp(scheme, "bitpal")) S[a][b] = a == b ? 2 : -3, gap = -5;
            else if (!strcmp(scheme, "tstv")) S[a][b] = a == b ? 0 : purine(a) == purine(b) ? -1 : -2, gap = -2;
            else if (!strcmp(scheme, "blosum")) gap = -4;
            else { fprintf(stderr, "unknown scheme %s\n", scheme); return 2; }
        }
    if (!strcmp(scheme, "blosum")) {
        const parasail_matrix_t *mx = parasail_matrix_lookup("blosum62");
        for (int i = 0; i < mx->size; i++)
            for (int j = 0; j < mx->size; j++)
                S[(unsigned char)mx->alphabet[i]][(unsigned char)mx->alphabet[j]] = mx->matrix[i * mx->size + j];
    }

    pairhash ph;
    pairhash_init(&ph);
    long long total = 0, pairs = 0;
    if (argc == 5) {
        int nq, ns;
        char **q = read_lines(argv[3], &nq), **s = read_lines(argv[4], &ns);
        for (int a = 0; a < nq; a++)
            for (int b = 0; b < ns; b++, pairs++) {
                long long c = cost(q[a], s[b]);
                total += c;
                pairhash_add(&ph, c);
            }
    } else {
        int n;
        char **lines = read_lines(argv[2], &n);
        for (int i = 0; i < n; i++) {
            char *tab = strchr(lines[i], '\t');
            if (!tab) continue;
            *tab = 0;
            long long c = cost(lines[i], tab + 1);
            total += c;
            pairhash_add(&ph, c);
            pairs++;
        }
    }
    printf("dp_check %s pairs=%lld cost_sum=%lld hash=%016llx\n", scheme, pairs, total, pairhash_end(&ph));
    return 0;
}
