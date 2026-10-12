"""Sum BGSA's result file over the real query x subject pairs and compare it with a reference total.

BGSA (`aligner -f RESULT`) writes one int16 score per pair, query-major, with the subject count
padded up to its lane multiple; the padding columns are ignored.

usage: python bgsa_check.py RESULT QUERIES SUBJECTS [TOTAL]
  TOTAL: bitdp-bench's cost_sum, or dp_check's total (scores and costs differ only in sign).
"""
import struct
import sys


def main():
    res, qf, sf = sys.argv[1:4]
    nq = sum(1 for l in open(qf) if l.strip())
    ns = sum(1 for l in open(sf) if l.strip())
    data = open(res, "rb").read()
    v = struct.unpack("<%dh" % (len(data) // 2), data)
    assert len(v) % nq == 0 and len(v) // nq >= ns, (len(v), nq, ns)
    padded = len(v) // nq
    total = sum(v[q * padded + s] for q in range(nq) for s in range(ns))
    print(f"BGSA: {nq} x {ns} pairs (subjects padded to {padded}), score sum {total}")
    if len(sys.argv) > 4:
        ref = int(sys.argv[4])
        ok = abs(ref) == abs(total)
        print("MATCH" if ok else f"MISMATCH: reference {ref}")
        sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
