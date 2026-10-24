"""Summary of a reproduce.sh session: per cell and tool the median, first and third quartile and range of
the rate (GCUPS) over the timed repetitions (the warm-up run is left out), and for every pair of tools of
a cell the paired ratio of their rates in the same repetition (median, quartiles, range of those ratios).
Rates are read as decimals and never pass through binary floats, so later rounding (texmacros.py) is exact.
Tools known to return wrong costs keep their rates and are flagged in the last column.

usage: python3 summary.py runs.txt > summary.csv
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
    rate = defaultdict(dict)  # (cell, tool) -> {rep: gcups}
    order = []
    wrong = set()
    for line in open(sys.argv[1]):
        m = re.match(r"cell=(\S+) rep=(\S+) tool=(\S+) (.*)", line)
        if not m or m.group(2) == "warm":
            continue
        cell, rep, tool, rest = m.groups()
        if rest.strip() == "expected-mismatch":
            wrong.add((cell, tool))
            continue
        g = re.search(r"gcups=([0-9.]+)", rest)
        if g:
            if (cell, tool) not in rate:
                order.append((cell, tool))
            rate[cell, tool][int(rep)] = Decimal(g.group(1))
    print("kind,cell,tool,versus,n,median,q1,q3,min,max,wrong_costs")
    cells = list(dict.fromkeys(c for c, _ in order))
    for cell in cells:
        tools = [t for c, t in order if c == cell]
        for t in tools:
            row("rate", cell, t, "", list(rate[cell, t].values()), int((cell, t) in wrong))
        for a in tools:
            for b in tools:
                reps = sorted(set(rate[cell, a]) & set(rate[cell, b]))
                if a != b and reps:
                    row("ratio", cell, a, b, [rate[cell, a][k] / rate[cell, b][k] for k in reps])


if __name__ == "__main__":
    main()
