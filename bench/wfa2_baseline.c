// WFA2-lib baseline: exact global (end-to-end) alignment cost of every pair in a "pattern<TAB>text" file,
// score only. WFA2's default wf-adaptive pruning is a heuristic, so it is switched off. Unit costs use
// WFA2's edit metric; BitPAl's weights use its gap-linear metric with the weights as costs (match -2,
// mismatch 3, indel 5), since WFA2 accepts a negative match score (Eizenga's formulation). WFA2 takes no
// substitution matrix, so it has no transition/transversion or BLOSUM62 run.
//
// build: gcc -O3 -march=native wfa2_baseline.c -I$DEPS/WFA2-lib -o wfa2_baseline $DEPS/WFA2-lib/lib/libwfa.a -lm
// usage: wfa2_baseline {edit|bitpal} FILE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "wavefront/wavefront_align.h"
#include "pairhash.h"

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
    if (argc != 3 || (strcmp(argv[1], "edit") && strcmp(argv[1], "bitpal"))) {
        fprintf(stderr, "usage: %s {edit|bitpal} FILE\n", argv[0]);
        return 2;
    }
    wavefront_aligner_attr_t attr = wavefront_aligner_attr_default;
    attr.alignment_scope = compute_score;
    attr.alignment_form.span = alignment_end2end;
    attr.heuristic.strategy = wf_heuristic_none;
    if (!strcmp(argv[1], "edit")) {
        attr.distance_metric = edit;
    } else {
        attr.distance_metric = gap_linear;
        attr.linear_penalties.match = -2;
        attr.linear_penalties.mismatch = 3;
        attr.linear_penalties.indel = 5;
    }
    wavefront_aligner_t *wf = wavefront_aligner_new(&attr);

    pairset in = pairs_load(argv[2]);
    if (!in.p) return 1;
    long long *cost = malloc((in.count + 1) * sizeof *cost);
    long long total = 0, rounds = 0;
    double sec, min_sec = bench_min_seconds(), t0 = now();
    do {
        for (long long i = 0; i < in.count; i++) {
            wavefront_align(wf, in.p[i], in.m[i], in.t[i], in.n[i]);
            // WFA2 reports the edit metric as the distance and gap-linear as the negated penalty
            cost[i] = attr.distance_metric == edit ? wf->cigar->score : -(long long)wf->cigar->score;
        }
        rounds++;
    } while ((sec = now() - t0) < min_sec);
    pairhash ph;
    pairhash_init(&ph);
    for (long long i = 0; i < in.count; i++) total += cost[i], pairhash_add(&ph, cost[i]);
    printf("wfa2 %s pairs=%lld cost_sum=%lld hash=%016llx seconds=%.4f gcups=%.3f rounds=%lld\n", argv[1], in.count, total,
           pairhash_end(&ph), sec, in.cells * rounds / sec * 1e-9, rounds);
    wavefront_aligner_delete(wf);
    return 0;
}
