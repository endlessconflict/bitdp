// Per-pair cost hash of the benchmark drivers: FNV-1a (64 bit) over each pair's cost, written in decimal
// and followed by a newline, in input order (all-against-all: query-major). It equals the FNV-1a of the
// file that the environment variable BENCH_DUMP names, where the driver writes those lines, so the
// harness can find the first pair on which two tools differ. bench.zig computes the same hash.
#include <stdio.h>
#include <stdlib.h>

typedef struct {
    unsigned long long h;
    FILE *dump;
} pairhash;

static void pairhash_init(pairhash *ph) {
    const char *path = getenv("BENCH_DUMP");
    ph->h = 0xcbf29ce484222325ULL;
    ph->dump = path ? fopen(path, "w") : NULL;
}

static void pairhash_add(pairhash *ph, long long cost) {
    char buf[32];
    int n = snprintf(buf, sizeof buf, "%lld\n", cost);
    for (int i = 0; i < n; i++) ph->h = (ph->h ^ (unsigned char)buf[i]) * 0x100000001b3ULL;
    if (ph->dump) fputs(buf, ph->dump);
}

static unsigned long long pairhash_end(pairhash *ph) {
    if (ph->dump) fclose(ph->dump);
    return ph->h;
}
