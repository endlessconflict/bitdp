"""Summary of a reproduce.sh session. From the run log: per cell and tool the median, first and third
quartile and range of the rate (GCUPS) over the timed repetitions (the warm-up run is left out), and for
every pair of tools of a cell the paired ratio of their rates in the same repetition (median, quartiles,
range of those ratios). From the scan log: the same for the scan rates (tool "bitdp") and the hit count of
each configuration. From the compile log: the compile step's seconds and MaxRSS (MB) per scheme, operating
system and optimization mode (cell "compile-OS-MODE"), warm-up builds left out. From the instruction log
(step insn): retired instructions per column word and per DP cell of each counted run.
Values are read as decimals and never pass through binary floats, so later rounding (texmacros.py) is exact.
Tools known to return wrong costs keep their rates and are flagged in the last column.

usage: python3 summary.py runs.txt [scan.txt] [compile.txt] [insn.txt] > summary.csv
"""
import re
import statistics
import sys
from collections import defaultdict
from decimal import Decimal


def quartiles(v):
    if len(v) == 1:
        return v[0], v[0], v[0]
    q1, med, q3 = statistics.quantiles(v, n=4, method="inclusive")
    return q1, med, q3


def row(kind, cell, a, b, v, wrong=""):
    q1, med, q3 = quartiles(v)
    print(f"{kind},{cell},{a},{b},{len(v)},{med},{q1},{q3},{min(v)},{max(v)},{wrong}")


def main():
    sys.stdout.reconfigure(newline="\n")  # LF on Windows too
    rate = defaultdict(dict)  # (cell, tool) -> {rep: gcups}
    order, wrong = [], set()
    hits, sec, mb, insn = {}, defaultdict(list), defaultdict(list), []
    for path in sys.argv[1:]:
        for line in open(path):
            m = re.match(r"cell=(\S+) tool=(\S+) lanes=\d+ wordcols=(\d+) cells=(\d+) .* instructions=(\d+) ", line)
            if m:
                insn.append((m[1], m[2], Decimal(m[5]) / Decimal(m[3]), Decimal(m[5]) / Decimal(m[4])))
                continue
            m = re.match(r"scheme=(\S+) os=(\S+) rep=(\d+) .*compile exe cc (\w+) \S+ success ([0-9.]+)(ms|s) MaxRSS:(\d+)M", line)
            if m:
                key = (f"compile-{m[2]}-{m[4]}", m[1])
                sec[key].append(Decimal(m[5]) / (1000 if m[6] == "ms" else 1))
                mb[key].append(Decimal(m[7]))
                continue
            m = re.match(r"cell=(\S+) rep=(\S+) (?:tool=(\S+) (.*)|genome .*: (\d+) hits in \S+ s, ([0-9.]+) GCUPS)", line)
            if not m or m[2] in ("warm", "verify"):
                continue
            cell, rep = m[1], int(m[2])
            if m[3] is None:  # scan line
                tool, value = "bitdp", m[6]
                hits[cell] = m[5]
            elif m[4].strip() == "expected-mismatch":
                wrong.add((cell, m[3]))
                continue
            else:
                g = re.search(r"gcups=([0-9.]+)", m[4])
                if not g:
                    continue
                tool, value = m[3], g[1]
            if (cell, tool) not in rate:
                order.append((cell, tool))
            rate[cell, tool][rep] = Decimal(value)
    print("kind,cell,tool,versus,n,median,q1,q3,min,max,wrong_costs")
    for cell in dict.fromkeys(c for c, _ in order):
        tools = [t for c, t in order if c == cell]
        for t in tools:
            row("rate", cell, t, "", list(rate[cell, t].values()), int((cell, t) in wrong))
        for a in tools:
            for b in tools:
                reps = sorted(set(rate[cell, a]) & set(rate[cell, b]))
                if a != b and reps:
                    row("ratio", cell, a, b, [rate[cell, a][k] / rate[cell, b][k] for k in reps])
    for cell, n in hits.items():
        row("hits", cell, "bitdp", "", [Decimal(n)])
    for key in sec:
        row("seconds", *key, "", sec[key])
        row("megabytes", *key, "", mb[key])
    for cell, tool, per_word, per_cell in insn:
        row("insnword", cell, tool, "", [per_word])
        row("insncell", cell, tool, "", [per_cell])


if __name__ == "__main__":
    main()
