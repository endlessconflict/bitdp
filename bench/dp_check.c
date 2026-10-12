// Independent scalar DP for the all-against-all check: global alignment score with linear gaps
// for every query x subject pair, summed. Shares no code with bitdp or BGSA.
//
// usage: dp_check MATCH MISMATCH GAP QUERIES SUBJECTS
//   unit costs as scores: dp_check 0 -1 -1 q.txt s.txt;  BitPAl's weights: dp_check 2 -3 -5 q.txt s.txt
// build: gcc -O3 -march=native dp_check.c -o dp_check
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char **read_lines(const char *path, int *n) {
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(1); }
    char **v = NULL, buf[1 << 16];
    *n = 0;
    while (fgets(buf, sizeof buf, f)) {
        size_t l = strcspn(buf, "\r\n");
        if (!l) continue;
        buf[l] = 0;
        v = realloc(v, (*n + 1) * sizeof *v);
        v[(*n)++] = strdup(buf);
    }
    fclose(f);
    return v;
}

int main(int argc, char **argv) {
    if (argc != 6) { fprintf(stderr, "usage: %s MATCH MISMATCH GAP QUERIES SUBJECTS\n", argv[0]); return 2; }
    int ma = atoi(argv[1]), mi = atoi(argv[2]), g = atoi(argv[3]), nq, ns;
    char **q = read_lines(argv[4], &nq), **s = read_lines(argv[5], &ns);
    long long total = 0;
    int *row = NULL;
    for (int a = 0; a < nq; a++) {
        int m = strlen(q[a]);
        row = realloc(row, (m + 1) * sizeof *row);
        for (int b = 0; b < ns; b++) {
            int n = strlen(s[b]);
            for (int i = 0; i <= m; i++) row[i] = i * g;              // column 0
            for (int j = 1; j <= n; j++) {
                int diag = row[0];
                row[0] = j * g;
                for (int i = 1; i <= m; i++) {
                    int best = diag + (q[a][i - 1] == s[b][j - 1] ? ma : mi);
                    if (row[i] + g > best) best = row[i] + g;               // gap in the query
                    if (row[i - 1] + g > best) best = row[i - 1] + g;       // gap in the subject
                    diag = row[i];
                    row[i] = best;
                }
            }
            total += row[m];
        }
    }
    printf("dp_check: %d x %d pairs, score sum %lld\n", nq, ns, total);
    return 0;
}
