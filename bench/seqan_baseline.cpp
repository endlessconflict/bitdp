// SeqAn 2 baseline: global alignment cost of every pair in a "pattern<TAB>text" file with SeqAn's
// inter-sequence vectorized alignment, globalAlignmentScore(ExecutionPolicy<Serial, Vectorial>, ...),
// score only, linear gaps (gap open = gap extension). Scores are 16-bit when no pair can leave the 16-bit
// range (|score| <= (m + n) * largest weight), else 32-bit; the line names the width used.
//
// build: g++ -O3 -march=native -std=c++17 seqan_baseline.cpp -I$DEPS/seqan/include -o seqan_baseline
// usage: seqan_baseline {edit|bitpal|tstv|blosum} FILE
#include <seqan/align_parallel.h>
#include <seqan/score.h>
#include <seqan/sequence.h>

#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
extern "C" {
#include "pairhash.h"
}

using namespace seqan2;

template <typename TAlphabet, typename TScore>
static int run(char const *scheme, std::vector<std::string> const &ps, std::vector<std::string> const &ts, TScore const &sc,
               int width) {
    StringSet<String<TAlphabet>> setH, setV;
    long long cells = 0;
    for (size_t i = 0; i < ps.size(); i++) {
        appendValue(setH, ps[i].c_str());
        appendValue(setV, ts[i].c_str());
        cells += (long long)ps[i].size() * ts[i].size();
    }
    ExecutionPolicy<Serial, Vectorial> exec;
    auto t0 = std::chrono::steady_clock::now();
    auto scores = globalAlignmentScore(exec, setH, setV, sc);
    auto since = [&] { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
    long long rounds = 1; // the whole call, repeated until the minimum region of pairhash.h
    double sec, min_sec = bench_min_seconds();
    for (; (sec = since()) < min_sec; rounds++) scores = globalAlignmentScore(exec, setH, setV, sc);
    pairhash ph;
    pairhash_init(&ph);
    long long total = 0;
    for (size_t i = 0; i < length(scores); i++) {
        total -= scores[i];
        pairhash_add(&ph, -(long long)scores[i]);
    }
#ifdef SEQAN_SIMD_ENABLED
    char const *mode = "simd"; // the Vectorial policy runs SIMD kernels only when SeqAn detects SIMD support
#else
    char const *mode = "novec";
#endif
    std::printf("seqan-%s%d %s pairs=%zu cost_sum=%lld hash=%016llx seconds=%.4f gcups=%.3f rounds=%lld\n", mode, width, scheme,
                (size_t)length(scores), total, pairhash_end(&ph), sec, cells * rounds / sec * 1e-9, rounds);
    return 0;
}

template <typename TValue>
static int dispatch(char const *scheme, std::vector<std::string> const &ps, std::vector<std::string> const &ts, int width) {
    if (!std::strcmp(scheme, "edit")) return run<Dna>(scheme, ps, ts, Score<TValue, Simple>(0, -1, -1, -1), width);
    if (!std::strcmp(scheme, "bitpal")) return run<Dna>(scheme, ps, ts, Score<TValue, Simple>(2, -3, -5, -5), width);
    if (!std::strcmp(scheme, "tstv")) {
        Score<TValue, ScoreMatrix<Dna, Default>> sc(-2, -2); // gap extension, gap open
        char const *a = "ACGT";
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++) {
                bool pi = a[i] == 'A' || a[i] == 'G', pj = a[j] == 'A' || a[j] == 'G';
                setScore(sc, Dna(a[i]), Dna(a[j]), i == j ? 0 : pi == pj ? -1 : -2);
            }
        return run<Dna>(scheme, ps, ts, sc, width);
    }
    if (!std::strcmp(scheme, "blosum"))
        return run<AminoAcid>(scheme, ps, ts, Score<TValue, ScoreMatrix<AminoAcid, ScoreSpecBlosum62>>(-4, -4), width);
    std::fprintf(stderr, "unknown scheme %s\n", scheme);
    return 2;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        std::fprintf(stderr, "usage: %s {edit|bitpal|tstv|blosum} FILE\n", argv[0]);
        return 2;
    }
    std::ifstream in(argv[2]);
    std::vector<std::string> ps, ts;
    size_t longest = 0;
    for (std::string line; std::getline(in, line);) {
        size_t tab = line.find('\t');
        if (tab == std::string::npos) continue;
        ps.push_back(line.substr(0, tab));
        ts.push_back(line.substr(tab + 1));
        longest = std::max(longest, line.size() - 1);
    }
    int weight = !std::strcmp(argv[1], "bitpal") ? 5 : !std::strcmp(argv[1], "blosum") ? 11 : 2; // largest |weight|
    if ((long long)longest * weight <= 32767) return dispatch<int16_t>(argv[1], ps, ts, 16);
    return dispatch<int32_t>(argv[1], ps, ts, 32);
}
