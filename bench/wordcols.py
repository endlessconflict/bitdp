"""Column words of one bitdp-bench run: over every group of LANES alignments the kernel runs together,
the number of text columns times the number of 63-row words of the longest pattern, as in
src/root.zig (columns): a group runs every word of every column until its longest text ends.
Also the DP cells of the run, as bench.zig counts them (pattern length times text length per pair).

usage: wordcols.py LANES ava QUERIES SUBJECTS   (subjects in groups of LANES against each query)
       wordcols.py LANES PAIRS batch             (consecutive pairs in groups of LANES)
Prints "lanes=LANES wordcols=N cells=M".
"""
import sys

lanes = int(sys.argv[1])
words = lambda m: m // 63 + 1  # one spare row for the bottom-row state
if sys.argv[2] == "ava":
    queries = open(sys.argv[3]).read().split()
    subjects = open(sys.argv[4]).read().split()
    per_query = sum(words(max(map(len, subjects[i:i + lanes]))) for i in range(0, len(subjects), lanes))
    n = sum(len(q) for q in queries) * per_query
    cells = sum(len(q) for q in queries) * sum(len(s) for s in subjects)
else:
    pairs = [line.split("\t", 1) for line in open(sys.argv[2]).read().split("\n") if "\t" in line]
    n = sum(words(max(len(p) for p, _ in g)) * max(len(t) for _, t in g)
            for g in (pairs[i:i + lanes] for i in range(0, len(pairs), lanes)))
    cells = sum(len(p) * len(t) for p, t in pairs)
print(f"lanes={lanes} wordcols={n} cells={cells}")
