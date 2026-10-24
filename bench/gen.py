"""Generate sequence pairs for the benchmark: one "pattern<TAB>text" per line.

Each text is its pattern after random substitutions, insertions and
deletions, each at rate div/3 per position.

usage: python gen.py OUT PAIRS LENGTH DIVERGENCE {dna|protein} [SEED]
       python gen.py ava QUERIES SUBJECTS NQ NS LENGTH [SEED]
The second form writes uniform random DNA, one sequence per line, for the all-against-all
setting (bench/reproduce.sh draws Table 4's 64 x 12 500 of 150 bp and 16 x 1 250 of 1 kbp
with seeds 4 and 5).
"""
import random
import sys

DNA = "ACGT"
PROTEIN = "ARNDCQEGHILKMFPSTWYV"


def mutate(seq, div, alphabet, rnd):
    out = []
    for c in seq:
        r = rnd.random()
        if r < div / 3:
            out.append(rnd.choice(alphabet))  # substitution
        elif r < 2 * div / 3:
            out.append(c)
            out.append(rnd.choice(alphabet))  # insertion
        elif r < div:
            pass  # deletion
        else:
            out.append(c)
    return "".join(out) or rnd.choice(alphabet)


def main():
    if sys.argv[1] == "ava":
        qf, sf, nq, ns, length = sys.argv[2:7]
        rnd = random.Random(int(sys.argv[7]) if len(sys.argv) > 7 else 1)
        for path, count in ((qf, nq), (sf, ns)):
            with open(path, "w", newline="\n") as f:
                for _ in range(int(count)):
                    f.write("".join(rnd.choice(DNA) for _ in range(int(length))) + "\n")
        return
    out, pairs, length, div, kind = sys.argv[1:6]
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    alphabet = DNA if kind == "dna" else PROTEIN
    rnd = random.Random(seed)
    with open(out, "w", newline="\n") as f:
        for _ in range(int(pairs)):
            p = "".join(rnd.choice(alphabet) for _ in range(int(length)))
            f.write(p + "\t" + mutate(p, float(div), alphabet, rnd) + "\n")


if __name__ == "__main__":
    main()
