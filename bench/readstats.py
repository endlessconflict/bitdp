"""Length statistics of read pair files ("read<TAB>reference span" per line), one line per file:
pairs, the shortest, median and longest read, and how many reads have the most common length.

usage: python3 readstats.py FILE.tsv ...
"""
import statistics
import sys
from collections import Counter

for path in sys.argv[1:]:
    lengths = [len(line.split("\t", 1)[0]) for line in open(path) if "\t" in line]
    mode, count = Counter(lengths).most_common(1)[0]
    print(f"{path}: pairs={len(lengths)} read_min={min(lengths)} read_median={statistics.median(lengths):g} "
          f"read_max={max(lengths)} mode={mode} at_mode={count}")
