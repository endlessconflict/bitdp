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

// The timed region of every driver is the whole loop over the input, repeated until BENCH_MIN_SECONDS
// (default 0.25) have passed; the rate is over all rounds. Costs are stored in the loop and hashed after it.
static inline double bench_min_seconds(void) {
    const char *s = getenv("BENCH_MIN_SECONDS");
    return s ? atof(s) : 0.25;
}

// A "pattern<TAB>text" file, split in place before the timed region.
typedef struct {
    char **p, **t;
    int *m, *n;
    long long count, cells;
} pairset;

static inline pairset pairs_load(const char *path) {
    pairset s = {0};
    FILE *f = fopen(path, "rb");
    if (!f) return s;
    fseek(f, 0, SEEK_END);
    long size = ftell(f);
    fseek(f, 0, SEEK_SET);
    char *buf = (char *)malloc(size + 1);
    if (fread(buf, 1, size, f) != (size_t)size) size = 0;
    buf[size] = 0;
    fclose(f);
    long long lines = 1;
    for (long i = 0; i < size; i++) lines += buf[i] == '\n';
    s.p = (char **)malloc(lines * sizeof *s.p), s.t = (char **)malloc(lines * sizeof *s.t);
    s.m = (int *)malloc(lines * sizeof *s.m), s.n = (int *)malloc(lines * sizeof *s.n);
    for (char *line = buf, *end; *line; line = end) {
        char *stop = line;
        while (*stop && *stop != '\n') stop++;
        end = *stop ? stop + 1 : stop;
        *stop = 0;
        char *tab = line;
        while (*tab && *tab != '\t') tab++;
        if (!*tab) continue;
        *tab = 0;
        s.p[s.count] = line, s.t[s.count] = tab + 1;
        s.m[s.count] = (int)(tab - line), s.n[s.count] = (int)(stop - tab - 1);
        s.cells += (long long)s.m[s.count] * s.n[s.count];
        s.count++;
    }
    return s;
}
