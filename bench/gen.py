"""Generate sequence pairs for the benchmark: one "pattern<TAB>text" per line.

Each text is its pattern after random substitutions, insertions and
deletions, each at rate div/3 per position.

usage: python gen.py OUT PAIRS LENGTH DIVERGENCE {dna|protein} [SEED]
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
