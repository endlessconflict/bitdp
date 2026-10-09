/* bitdp: bit-parallel alignment kernels, C interface (C99, usable from C++).
 *
 * Precompiled schemes, all as costs to minimize with linear gaps:
 *   edit    unit costs (match 0, mismatch 1, gap 1)
 *   tstv    transition 1 (A<->G, C<->T), transversion 2, gap 2
 *   bitpal  BitPAl's score weights (2, -3, -5) as costs (match -2, mismatch 3, gap 5)
 * Modes: global (whole pattern against whole text) and search (whole pattern
 * against the best-matching substring of the text; cost = minimum over end positions).
 * Sequences are bytes; DNA letters are upper-case ACGT, and any other byte counts as
 * the costliest substitution. Custom schemes are Zig-only (they are derived at compile time).
 *
 * Every function returns BITDP_OK (0) or an error code; results go through pointers.
 * Memory the library returns is released with the matching bitdp_*_free function only.
 * Functions are reentrant; a prepared aligner must not be used by two threads at once.
 */
#ifndef BITDP_H
#define BITDP_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define BITDP_OK 0     /* success */
#define BITDP_EINVAL 1 /* bad argument: null pointer, empty pattern, or pattern too long */
#define BITDP_ENOMEM 2 /* allocation failed */

/* One search hit: pattern index, end position in the text (exclusive, 1-based), cost. */
typedef struct bitdp_hit {
    size_t pattern;
    size_t end;
    int64_t cost;
} bitdp_hit;

typedef struct bitdp_edit_global_aligner bitdp_edit_global_aligner;
typedef struct bitdp_edit_search_aligner bitdp_edit_search_aligner;
typedef struct bitdp_tstv_global_aligner bitdp_tstv_global_aligner;
typedef struct bitdp_tstv_search_aligner bitdp_tstv_search_aligner;
typedef struct bitdp_bitpal_global_aligner bitdp_bitpal_global_aligner;

/* ---- edit, global ---- */
/* Cost of one pair; pattern length 1..62, no allocation. */
int bitdp_edit_global_distance(const char *pattern, size_t m, const char *text, size_t n, int64_t *cost);
/* Prepare a pattern of any length >= 1 for many texts. */
int bitdp_edit_global_aligner_new(const char *pattern, size_t m, bitdp_edit_global_aligner **out);
/* Cost of the prepared pattern against one text. */
int bitdp_edit_global_aligner_distance(bitdp_edit_global_aligner *a, const char *text, size_t n, int64_t *cost);
/* Release a prepared pattern (null is ignored). */
void bitdp_edit_global_aligner_free(bitdp_edit_global_aligner *a);
/* Costs of count independent pairs, one pair per SIMD lane; cost[i] for pair i. */
int bitdp_edit_global_distances(const char *const *patterns, const size_t *m, const char *const *texts, const size_t *n, size_t count, int64_t *cost);

/* ---- edit, search ---- */
/* Best cost of the pattern (length 1..62) against any substring of the text. */
int bitdp_edit_search_distance(const char *pattern, size_t m, const char *text, size_t n, int64_t *cost);
/* Prepare a pattern of any length >= 1 for many texts. */
int bitdp_edit_search_aligner_new(const char *pattern, size_t m, bitdp_edit_search_aligner **out);
/* Best substring cost of the prepared pattern against one text. */
int bitdp_edit_search_aligner_distance(bitdp_edit_search_aligner *a, const char *text, size_t n, int64_t *cost);
/* Release a prepared pattern (null is ignored). */
void bitdp_edit_search_aligner_free(bitdp_edit_search_aligner *a);
/* Best substring costs of count independent pairs, one pair per SIMD lane. */
int bitdp_edit_search_distances(const char *const *patterns, const size_t *m, const char *const *texts, const size_t *n, size_t count, int64_t *cost);
/* Every (pattern, end) with cost <= max_cost; *hits is freed with bitdp_hits_free (null when *nhits is 0). */
int bitdp_edit_search_scan(const char *const *patterns, const size_t *m, size_t count, const char *text, size_t n, int64_t max_cost, bitdp_hit **hits, size_t *nhits);

/* ---- tstv, global ---- */
/* Cost of one pair; pattern length 1..62, no allocation. */
int bitdp_tstv_global_distance(const char *pattern, size_t m, const char *text, size_t n, int64_t *cost);
/* Prepare a pattern of any length >= 1 for many texts. */
int bitdp_tstv_global_aligner_new(const char *pattern, size_t m, bitdp_tstv_global_aligner **out);
/* Cost of the prepared pattern against one text. */
int bitdp_tstv_global_aligner_distance(bitdp_tstv_global_aligner *a, const char *text, size_t n, int64_t *cost);
/* Release a prepared pattern (null is ignored). */
void bitdp_tstv_global_aligner_free(bitdp_tstv_global_aligner *a);
/* Costs of count independent pairs, one pair per SIMD lane. */
int bitdp_tstv_global_distances(const char *const *patterns, const size_t *m, const char *const *texts, const size_t *n, size_t count, int64_t *cost);

/* ---- tstv, search ---- */
/* Best cost of the pattern (length 1..62) against any substring of the text. */
int bitdp_tstv_search_distance(const char *pattern, size_t m, const char *text, size_t n, int64_t *cost);
/* Prepare a pattern of any length >= 1 for many texts. */
int bitdp_tstv_search_aligner_new(const char *pattern, size_t m, bitdp_tstv_search_aligner **out);
/* Best substring cost of the prepared pattern against one text. */
int bitdp_tstv_search_aligner_distance(bitdp_tstv_search_aligner *a, const char *text, size_t n, int64_t *cost);
/* Release a prepared pattern (null is ignored). */
void bitdp_tstv_search_aligner_free(bitdp_tstv_search_aligner *a);
/* Best substring costs of count independent pairs, one pair per SIMD lane. */
int bitdp_tstv_search_distances(const char *const *patterns, const size_t *m, const char *const *texts, const size_t *n, size_t count, int64_t *cost);
/* Every (pattern, end) with cost <= max_cost; *hits is freed with bitdp_hits_free (null when *nhits is 0). */
int bitdp_tstv_search_scan(const char *const *patterns, const size_t *m, size_t count, const char *text, size_t n, int64_t max_cost, bitdp_hit **hits, size_t *nhits);

/* ---- bitpal, global ---- */
/* Cost of one pair; pattern length 1..62, no allocation. */
int bitdp_bitpal_global_distance(const char *pattern, size_t m, const char *text, size_t n, int64_t *cost);
/* Prepare a pattern of any length >= 1 for many texts. */
int bitdp_bitpal_global_aligner_new(const char *pattern, size_t m, bitdp_bitpal_global_aligner **out);
/* Cost of the prepared pattern against one text. */
int bitdp_bitpal_global_aligner_distance(bitdp_bitpal_global_aligner *a, const char *text, size_t n, int64_t *cost);
/* Release a prepared pattern (null is ignored). */
void bitdp_bitpal_global_aligner_free(bitdp_bitpal_global_aligner *a);
/* Costs of count independent pairs, one pair per SIMD lane. */
int bitdp_bitpal_global_distances(const char *const *patterns, const size_t *m, const char *const *texts, const size_t *n, size_t count, int64_t *cost);

/* ---- shared ---- */
/* Release a hit array returned by a *_scan function (pass the nhits it returned). */
void bitdp_hits_free(bitdp_hit *hits, size_t nhits);

#ifdef __cplusplus
}
#endif

#endif /* BITDP_H */
