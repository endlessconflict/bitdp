"""Turn mapped reads into pattern<TAB>text pairs: each read against the reference span it maps to.

usage: samtools view -F 0x904 aln.bam | python3 reads.py ref.fa > reads.tsv

Keeps primary alignments without clipping and without N. Prints the sum of the NM tags
(bwa's edit count) to stderr, an upper bound on the sum of unit edit distances.
"""
import re
import sys

ref = "".join(l.strip() for l in open(sys.argv[1]) if not l.startswith(">")).upper()
kept = nm = 0
for line in sys.stdin:
    f = line.split("\t")
    cigar, seq = f[5], f[9]
    if "S" in cigar or "H" in cigar or "N" in seq:
        continue
    span = sum(int(n) for n, op in re.findall(r"(\d+)([MIDN=X])", cigar) if op in "MDN=X")
    pos = int(f[3]) - 1
    print(f"{seq}\t{ref[pos:pos + span]}")
    kept += 1
    nm += next(int(t[5:]) for t in f[11:] if t.startswith("NM:i:"))
print(f"pairs={kept} nm_sum={nm}", file=sys.stderr)
