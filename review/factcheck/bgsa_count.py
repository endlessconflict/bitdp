"""E2 (round 3): count the word operations of BGSA's generated AVX2 BitPAl kernels under Table 1's rule.

Usage: python bgsa_count.py ALIGN_CORE.C [...]   (the generated file of each build; bench/raw/e2_bgsa_count.txt is its output)

Input: align_core.c as BGSA's generator (commit dffcf9e) emitted it for match 2, mismatch -3,
gap -5, standard ("nonpacked") and packed. The files are not shipped; the script only counts.

Counting rule (Table 1): one operation per AND, OR, XOR, NOT, shift and addition (and
subtraction) on a vector word, per column and per word of the pattern. NOT is XOR with all
ones, one operation. Not counted: loading or storing the state and the match masks,
constants (set1), plain copies, and the score update after the column loop, as for bitdp
and Myers' kernel. The counted region is the body of the per-word loop inside the
per-column loop (the first `for(j = 0; j < word_num; j++)` after `for(i = 0, itr = ref`).

BGSA runs every pattern on its multi-word path (31 rows per 32-bit lane, bit 31 as
carry catcher). Statements that pass a carry or a top bit between words are counted as
"transfer": they read or write overflowN, *_prevbit, tmpbitN, or the inter-word bitN
words. total - transfer is the count of the same kernel for a pattern of one word,
where those statements have nothing to pass; it is the counterpart of Table 1.
"""
import hashlib
import re
import sys

COUNTED = {"and_si256", "or_si256", "xor_si256", "andnot_si256", "add_epi32", "sub_epi32", "slli_epi32", "srli_epi32"}
IGNORED = {"load_si256", "store_si256", "set1_epi32"}
TRANSFER = re.compile(r"\b(overflow\d+|\w+_prevbit|tmpbit\d+|bit(?:1|2|4|8|16))\b")


def loop_body(src):
    col = src.index("for(i = 0, itr = ref")
    start = src.index("for(j = 0; j < word_num; j++)", col)
    i = src.index("{", start)
    depth = 0
    for k in range(i, len(src)):
        depth += {"{": 1, "}": -1}.get(src[k], 0)
        if depth == 0:
            return src[i + 1:k]
    raise ValueError("unbalanced loop")


def count(path):
    with open(path, encoding="utf-8") as f:
        src = f.read()
    weights = re.search(r"int match_score = (-?\d+);\s*int mismatch_score = (-?\d+);\s*int gap_score = (-?\d+);", src)
    total = transfer = 0
    unknown = set()
    for stmt in loop_body(src).split(";"):
        ops = re.findall(r"_mm256_(\w+)\(", stmt)
        for op in ops:
            if op in COUNTED:
                total += 1
                transfer += bool(TRANSFER.search(stmt))
            elif op not in IGNORED:
                unknown.add(op)
    return dict(md5=hashlib.md5(src.encode()).hexdigest(), weights=tuple(int(x) for x in weights.groups()),
                total=total, transfer=transfer, single=total - transfer, unknown=sorted(unknown))


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for p in sys.argv[1:]:
        c = count(p)
        kind = "packed" if "packed-" in p and "nonpacked" not in p else "standard"
        print(f"BGSA generated BitPAl {kind} AVX2 {c['weights']}: md5 {c['md5']}; per word and column: "
              f"total {c['total']}, inter-word transfer {c['transfer']}, single-word {c['single']}; unknown intrinsics {c['unknown']}")


if __name__ == "__main__":
    main()
