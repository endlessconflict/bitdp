"""Turn one BGSA run into the line the other drivers print: per-pair costs from its result file, their
sum and hash (pairhash.h), and BGSA's own "cal GCUPS" from its output.

BGSA (`aligner -f RESULT`) writes one int16 value per pair, query-major, with the subject count padded
up to its lane multiple; the padding columns are ignored. Its Myers kernels write the edit distance
("distance": the cost itself), its BitPAl kernels the score ("score": the cost is its negation). If
BENCH_DUMP is set, the costs are written there, one per line, as the other drivers do.

usage: python3 bgsa_check.py LABEL {distance|score} RESULT QUERIES SUBJECTS ALIGNER_OUTPUT
"""
import os
import re
import struct
import sys


def main():
    label, kind, res, qf, sf, log = sys.argv[1:7]
    sign = {"distance": 1, "score": -1}[kind]
    nq = sum(1 for l in open(qf) if l.strip())
    ns = sum(1 for l in open(sf) if l.strip())
    data = open(res, "rb").read()
    v = struct.unpack("<%dh" % (len(data) // 2), data)
    assert len(v) % nq == 0 and len(v) // nq >= ns, (len(v), nq, ns)
    padded = len(v) // nq
    costs = [sign * v[q * padded + s] for q in range(nq) for s in range(ns)]
    text = "".join(f"{c}\n" for c in costs)
    h = 0xCBF29CE484222325
    for b in text.encode():
        h = ((h ^ b) * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    if os.environ.get("BENCH_DUMP"):
        open(os.environ["BENCH_DUMP"], "w", newline="\n").write(text)
    total = sum(costs)
    rates = re.findall(r"[Cc]al GCUPS is ([0-9.]+)", open(log).read())
    assert rates, f"no GCUPS figure in {log}"
    print(f"{label} pairs={nq * ns} cost_sum={total} hash={h:016x} gcups={rates[-1]}")


if __name__ == "__main__":
    main()
