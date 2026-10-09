/* Tests for the C interface (include/bitdp.h) against a plain scalar DP written here.
 * Compiles as C99 and as C++ (zig build test-c builds and runs both). */
#include "bitdp.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;

#define CHECK(cond, ...)                          \
    do {                                          \
        if (!(cond)) {                            \
            if (failures < 20) {                  \
                printf("FAIL %s:%d: ", __FILE__, __LINE__); \
                printf(__VA_ARGS__);              \
                printf("\n");                     \
            }                                     \
            failures++;                           \
        }                                         \
    } while (0)

static uint64_t rng_state = 0x9e3779b97f4a7c15ull;
static uint64_t rnd(void) {
    rng_state ^= rng_state << 13;
    rng_state ^= rng_state >> 7;
    rng_state ^= rng_state << 17;
    return rng_state;
}
static size_t rnd_range(size_t lo, size_t hi) { return lo + (size_t)(rnd() % (hi - lo + 1)); }
static void rnd_dna(char *s, size_t n) {
    size_t i;
    for (i = 0; i < n; i++) s[i] = "ACGT"[rnd() & 3];
}

/* ---- the schemes, as costs ---- */
static int cost_edit(char a, char b) { return a == b ? 0 : 1; }
static int cost_bitpal(char a, char b) { return a == b ? -2 : 3; }
static int purine(char c) { return c == 'A' || c == 'G'; }
static int cost_tstv(char a, char b) {
    if (a == b) return 0;
    return purine(a) == purine(b) ? 1 : 2;
}

/* Plain DP. search = 0: D(m, n) with D(i, 0) = i*gap, D(0, j) = j*gap.
 * search = 1: D(0, j) = 0, result min_j D(m, j); col_out (if not null) receives D(m, j), j = 1..n. */
static int64_t dp(const char *p, size_t m, const char *t, size_t n, int (*cost)(char, char), int64_t gap, int search, int64_t *col_out) {
    int64_t *col = (int64_t *)malloc((m + 1) * sizeof(int64_t));
    int64_t best, diag, cell, up;
    size_t i, j;
    for (i = 0; i <= m; i++) col[i] = (int64_t)i * gap;
    best = col[m];
    for (j = 1; j <= n; j++) {
        diag = col[0];
        col[0] = search ? 0 : col[0] + gap;
        for (i = 1; i <= m; i++) {
            cell = diag + cost(p[i - 1], t[j - 1]);
            up = col[i] + gap;
            if (up < cell) cell = up;
            if (col[i - 1] + gap < cell) cell = col[i - 1] + gap;
            diag = col[i];
            col[i] = cell;
        }
        if (col[m] < best) best = col[m];
        if (col_out) col_out[j - 1] = col[m];
    }
    best = search ? best : col[m];
    free(col);
    return best;
}

typedef struct { size_t pattern, end; int64_t cost; } hit_t;
static int hit_cmp(const void *a, const void *b) {
    const hit_t *x = (const hit_t *)a, *y = (const hit_t *)b;
    if (x->pattern != y->pattern) return x->pattern < y->pattern ? -1 : 1;
    if (x->end != y->end) return x->end < y->end ? -1 : 1;
    return 0;
}

#define NPAIRS 10000
#define MAXN 400

/* One test function per exported scheme and mode; PRE is the symbol prefix. */
#define TEST_SCHEME(PRE, COST, GAP, SEARCH)                                                               \
    static void test_##PRE(void) {                                                                        \
        static char pbuf[NPAIRS][62], tbuf[NPAIRS][200];                                                  \
        static const char *ps[NPAIRS], *ts[NPAIRS];                                                       \
        static size_t ms[NPAIRS], ns[NPAIRS];                                                             \
        static int64_t want[NPAIRS], got[NPAIRS];                                                         \
        char lp[300], lt[MAXN];                                                                           \
        int64_t c = 0;                                                                                    \
        size_t k, i, lens[2] = {100, 300};                                                                \
        PRE##_aligner *a = NULL;                                                                          \
        int before = failures;                                                                            \
        for (k = 0; k < NPAIRS; k++) {                                                                    \
            ms[k] = rnd_range(1, 62);                                                                     \
            ns[k] = rnd_range(0, 200);                                                                    \
            rnd_dna(pbuf[k], ms[k]);                                                                      \
            rnd_dna(tbuf[k], ns[k]);                                                                      \
            ps[k] = pbuf[k];                                                                              \
            ts[k] = tbuf[k];                                                                              \
            want[k] = dp(ps[k], ms[k], ts[k], ns[k], COST, GAP, SEARCH, NULL);                           \
            CHECK(PRE##_distance(ps[k], ms[k], ts[k], ns[k], &c) == BITDP_OK && c == want[k],             \
                  #PRE "_distance pair %zu: got %lld want %lld", k, (long long)c, (long long)want[k]);    \
            CHECK(PRE##_aligner_new(ps[k], ms[k], &a) == BITDP_OK, #PRE "_aligner_new");                  \
            CHECK(PRE##_aligner_distance(a, ts[k], ns[k], &c) == BITDP_OK && c == want[k],                \
                  #PRE "_aligner_distance pair %zu", k);                                                  \
            PRE##_aligner_free(a);                                                                        \
        }                                                                                                 \
        CHECK(PRE##_distances(ps, ms, ts, ns, NPAIRS, got) == BITDP_OK, #PRE "_distances status");        \
        for (k = 0; k < NPAIRS; k++)                                                                      \
            CHECK(got[k] == want[k], #PRE "_distances pair %zu: got %lld want %lld", k, (long long)got[k], \
                  (long long)want[k]);                                                                    \
        /* edge cases: empty text (null pointer allowed), length 1, exactly 62, 63 rejected */            \
        CHECK(PRE##_distance("A", 1, NULL, 0, &c) == BITDP_OK && c == dp("A", 1, "", 0, COST, GAP, SEARCH, NULL), \
              #PRE " empty text");                                                                        \
        CHECK(PRE##_distance("G", 1, "ACGT", 4, &c) == BITDP_OK && c == dp("G", 1, "ACGT", 4, COST, GAP, SEARCH, NULL), \
              #PRE " length-1 pattern");                                                                  \
        rnd_dna(lp, 63);                                                                                  \
        rnd_dna(lt, 150);                                                                                 \
        CHECK(PRE##_distance(lp, 62, lt, 150, &c) == BITDP_OK && c == dp(lp, 62, lt, 150, COST, GAP, SEARCH, NULL), \
              #PRE " pattern of exactly 62");                                                             \
        CHECK(PRE##_distance(lp, 63, lt, 150, &c) == BITDP_EINVAL, #PRE " 63 rejected by _distance");     \
        CHECK(PRE##_distance(lp, 0, lt, 150, &c) == BITDP_EINVAL, #PRE " empty pattern rejected");        \
        CHECK(PRE##_distance(lp, 10, lt, 150, NULL) == BITDP_EINVAL, #PRE " null output rejected");       \
        CHECK(PRE##_distance(NULL, 10, lt, 150, &c) == BITDP_EINVAL, #PRE " null pattern rejected");      \
        CHECK(PRE##_aligner_new(lp, 0, &a) == BITDP_EINVAL, #PRE " empty pattern rejected by aligner");   \
        /* multi-word patterns: 100 and 300 rows, aligner and batched */                                 \
        for (i = 0; i < 2; i++) {                                                                         \
            for (k = 0; k < 100; k++) {                                                                   \
                size_t n = rnd_range(0, MAXN);                                                            \
                const char *pp[1], *tt[1];                                                                \
                size_t mm[1], nn[1];                                                                      \
                int64_t w, g = 0;                                                                         \
                rnd_dna(lp, lens[i]);                                                                     \
                rnd_dna(lt, n);                                                                           \
                w = dp(lp, lens[i], lt, n, COST, GAP, SEARCH, NULL);                                      \
                CHECK(PRE##_aligner_new(lp, lens[i], &a) == BITDP_OK, #PRE "_aligner_new %zu", lens[i]);  \
                CHECK(PRE##_aligner_distance(a, lt, n, &c) == BITDP_OK && c == w,                         \
                      #PRE " aligner, pattern %zu: got %lld want %lld", lens[i], (long long)c, (long long)w); \
                PRE##_aligner_free(a);                                                                    \
                pp[0] = lp; tt[0] = lt; mm[0] = lens[i]; nn[0] = n;                                       \
                CHECK(PRE##_distances(pp, mm, tt, nn, 1, &g) == BITDP_OK && g == w,                       \
                      #PRE " batched, pattern %zu", lens[i]);                                             \
            }                                                                                             \
        }                                                                                                 \
        PRE##_aligner_free(NULL);                                                                         \
        printf("%-22s %s (%d pairs + edge cases + multi-word 100/300)\n", #PRE, failures == before ? "ok" : "FAILED", NPAIRS); \
    }

TEST_SCHEME(bitdp_edit_global, cost_edit, 1, 0)
TEST_SCHEME(bitdp_edit_search, cost_edit, 1, 1)
TEST_SCHEME(bitdp_tstv_global, cost_tstv, 2, 0)
TEST_SCHEME(bitdp_tstv_search, cost_tstv, 2, 1)
TEST_SCHEME(bitdp_bitpal_global, cost_bitpal, 5, 0)

/* Search scan: every (pattern, end, cost) with cost <= K, against the scalar DP's last row. */
#define TEST_SCAN(PRE, COST, GAP)                                                                         \
    static void test_##PRE##_scan(void) {                                                                 \
        char pat[20][40], text[3000];                                                                     \
        const char *ps[20];                                                                               \
        size_t ms[20], trial, i, j, nwant, nhits;                                                         \
        int64_t col[3000];                                                                                \
        hit_t *want = (hit_t *)malloc(20 * 3000 * sizeof(hit_t));                                         \
        bitdp_hit *hits;                                                                                  \
        int before = failures;                                                                            \
        size_t total = 0;                                                                                 \
        for (trial = 0; trial < 60; trial++) {                                                            \
            size_t count = rnd_range(1, 20), n = rnd_range(0, 3000);                                      \
            int64_t k = (int64_t)rnd_range(0, 6);                                                         \
            rnd_dna(text, n);                                                                             \
            nwant = 0;                                                                                    \
            for (i = 0; i < count; i++) {                                                                 \
                ms[i] = rnd_range(1, 40);                                                                 \
                rnd_dna(pat[i], ms[i]);                                                                   \
                if (i > 0 && rnd() % 3 == 0) memcpy(pat[i], text + (n > 40 ? rnd_range(0, n - 41) : 0), n > 40 ? ms[i] : 0); \
                ps[i] = pat[i];                                                                           \
                dp(pat[i], ms[i], text, n, COST, GAP, 1, col);                                            \
                for (j = 0; j < n; j++)                                                                   \
                    if (col[j] <= k) { want[nwant].pattern = i; want[nwant].end = j + 1; want[nwant].cost = col[j]; nwant++; } \
            }                                                                                             \
            CHECK(PRE##_scan(ps, ms, count, text, n, k, &hits, &nhits) == BITDP_OK, #PRE "_scan status"); \
            CHECK(nhits == nwant, #PRE "_scan trial %zu: %zu hits, want %zu", trial, nhits, nwant);       \
            total += nwant;                                                                               \
            if (nhits == nwant && nhits > 0) {                                                            \
                hit_t *got = (hit_t *)malloc(nhits * sizeof(hit_t));                                      \
                for (i = 0; i < nhits; i++) { got[i].pattern = hits[i].pattern; got[i].end = hits[i].end; got[i].cost = hits[i].cost; } \
                qsort(got, nhits, sizeof(hit_t), hit_cmp);                                                \
                qsort(want, nwant, sizeof(hit_t), hit_cmp);                                               \
                for (i = 0; i < nhits; i++)                                                               \
                    CHECK(got[i].pattern == want[i].pattern && got[i].end == want[i].end && got[i].cost == want[i].cost, \
                          #PRE "_scan trial %zu hit %zu differs", trial, i);                              \
                free(got);                                                                                \
            }                                                                                             \
            bitdp_hits_free(hits, nhits);                                                                 \
        }                                                                                                 \
        CHECK(PRE##_scan(ps, ms, 0, text, 10, 3, &hits, &nhits) == BITDP_OK && hits == NULL && nhits == 0, #PRE "_scan no patterns"); \
        CHECK(PRE##_scan(ps, ms, 1, text, 10, 3, NULL, &nhits) == BITDP_EINVAL, #PRE "_scan null output"); \
        free(want);                                                                                       \
        printf("%-22s %s (60 random scans, up to 20 patterns, %zu hits checked against the scalar DP)\n", \
               #PRE "_scan", failures == before ? "ok" : "FAILED", total);                                 \
    }

TEST_SCAN(bitdp_edit_search, cost_edit, 1)
TEST_SCAN(bitdp_tstv_search, cost_tstv, 2)

int main(void) {
    test_bitdp_edit_global();
    test_bitdp_edit_search();
    test_bitdp_tstv_global();
    test_bitdp_tstv_search();
    test_bitdp_bitpal_global();
    test_bitdp_edit_search_scan();
    test_bitdp_tstv_search_scan();
    printf("%s: %d failure(s)\n", failures ? "FAILED" : "PASSED", failures);
    return failures ? 1 : 0;
}
