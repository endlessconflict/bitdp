"""Bergeron & Hamel's vector algorithm (IJFCS 13(1):53-65, 2002, Section 3.2, p. 64)
checked against the scalar DP, with word-operation counts under our counting rule.

Independent of the Zig implementation: nothing here reads or imports src/.

The paper gives the algorithm at the level of vectors of small integers. It does not
spell out the integer operations it relies on (assigning, comparing, and adding small
bounded integers, which it costs at O(log c) bit operations each, pp. 64-65). We complete
it in the most direct way: an integer vector is a list of bit planes, addition is a
ripple-carry adder over the planes, equality with a constant is a product of plane
literals with shared prefixes. Two tallies are kept:

  vec  - every integer-vector operation counts 1, every boolean operation as in our rule.
         This is the paper's own level of abstraction and a lower bound on any
         bit-level implementation.
  word - every word operation of our completion, under our counting rule
         (NOT counts one, (x<<1)|1 counts two, a+b counts one).

Their recurrence (7) has one indel cost c and substitution costs 0 <= delta < 2c
(p. 61-62), in search mode. Everything outside that is our extension, labelled in
the output: global mode (initial state 2c instead of c), delta = 2c (no change to
the algorithm), and negative costs, mapped into their range by the shift
sigma' = sigma + 2x, g' = g + x (see bh_cost; it holds in both modes).

The listing on p. 64 is run with two implementation choices of ours (see bh_column):
  (a) the last two assignments of the loop are read as simultaneous, so V is advanced
      before K is updated;
  (b) the loop test given to Ind, in two variants that give the same distances:
      restricted, (Δv = 0) ∧ ¬K, so that the vectors (r = k) are disjoint and Δh_j is
      assembled from them; or as printed, (Δv = 0), with Δh_j read from the cumulative
      vectors K = (r ≤ k), which are the same in both variants.
The final update of V, which nothing reads, is dropped. Both variants are checked
against the DP and counted under the same rule; the counts reported are those of the
printed test with Δh_j read from K, which never needs more operations than the other.

These counts are for the listing with every printed operation executed. They are not a
lower bound for the algorithm as such: with Δh_j read from K, the factor ¬K in N can be
left out ("K-" below), which saves a NOT and an AND per level at the vector level. That
variant is checked against the DP too, and its count is printed next to the others.
"""
import math
import os
import re
import random

from common import (PAPER, Report, dp_cost, match_mismatch, num, paper_text, parse_cost_triple,
                    table1_schemes, table_rows)

Z, O = "ZERO", "ONES"   # constant planes, folded at no cost


class Words:
    def __init__(self, m):
        self.mask = (1 << m) - 1
        self.word = self.vec = 0
        self.depth = 0

    def _n(self, k=1):
        self.word += k
        if self.depth == 0:
            self.vec += k

    def val(self, a):
        return 0 if a is Z else self.mask if a is O else a

    # boolean word operations with constant folding
    def AND(self, a, b):
        if a is Z or b is Z:
            return Z
        if a is O:
            return b
        if b is O:
            return a
        self._n()
        return a & b

    def OR(self, a, b):
        if a is O or b is O:
            return O
        if a is Z:
            return b
        if b is Z:
            return a
        self._n()
        return a | b

    def NOT(self, a):
        if a is Z:
            return O
        if a is O:
            return Z
        self._n()
        return ~a & self.mask

    def XOR(self, a, b):
        if a is Z:
            return b
        if b is Z:
            return a
        if a is O:
            return self.NOT(b)
        if b is O:
            return self.NOT(a)
        self._n()
        return a ^ b

    def SHL(self, a, fill):
        """Their up-arrow: row i moves to row i+1, row 1 gets `fill`."""
        self._n(2 if fill else 1)
        return ((self.val(a) << 1) | fill) & self.mask

    def ADD(self, a, b):
        """Their +_b: binary addition, carries run from row 1 upward."""
        self._n()
        return (self.val(a) + self.val(b)) & self.mask

    # integer-vector operations (planes LSB first); each call is one vector operation
    def kit(self):
        self.vec += 1
        self.depth += 1

    def unkit(self):
        self.depth -= 1

    def add(self, A, B, width, cin=Z):
        self.kit()
        out = []
        for k in range(width):
            a = A[k] if k < len(A) else Z
            b = B[k] if k < len(B) else Z
            t = self.XOR(a, b)
            out.append(self.XOR(t, cin))
            if k + 1 < width:
                cin = self.OR(self.AND(a, b), self.AND(cin, t))
        self.unkit()
        return out

    def add_bit(self, A, x, width):
        self.kit()
        out = []
        for k in range(width):
            a = A[k] if k < len(A) else Z
            out.append(self.XOR(a, x))
            if k + 1 < width:
                x = self.AND(a, x)
        self.unkit()
        return out

    def sub(self, A, B, width):
        self.kit()
        nb = [self.NOT(B[k] if k < len(B) else Z) for k in range(width)]
        self.unkit()
        self.vec -= 1  # the NOTs belong to this subtraction
        out = self.add(A, nb, width, cin=O)
        return out

    def shift(self, A, fill_value):
        self.kit()
        out = [self.SHL(a, (fill_value >> k) & 1) for k, a in enumerate(A)]
        self.unkit()
        return out

    def eq(self, A, ks):
        """[A == k] for each k, one vector operation per k; literals and prefixes shared."""
        lit, pre, out = {}, {}, {}
        for k in ks:
            self.kit()
            acc = O
            for t in reversed(range(len(A))):
                bit = (k >> t) & 1
                key = (t, bit)
                if key not in lit:
                    lit[key] = A[t] if bit else self.NOT(A[t])
                pk = (t, k >> t)
                if pk not in pre:
                    pre[pk] = self.AND(acc, lit[key])
                acc = pre[pk]
            out[k] = acc
            self.unkit()
        return out

    def encode(self, ind, width):
        """Integer planes from the indicator vectors (r = k)."""
        self.kit()
        out = []
        for t in range(width):
            acc = Z
            for k, v in ind.items():
                if (k >> t) & 1:
                    acc = self.OR(acc, v)
            out.append(acc)
        self.unkit()
        return out

    def encode_k(self, Ks, nKs, width):
        """Integer planes of the number of k with K_k = 0, from the nested vectors
        K_0 <= K_1 <= ... (K_k = (r <= k)). nKs holds the complements already computed."""
        self.kit()
        top = len(Ks)                                      # values 0..top; [value >= j] = not K_{j-1}
        ge = lambda j: nKs[j - 1] if j - 1 < len(nKs) else self.NOT(Ks[j - 1])
        out = []
        for t in range(width):
            acc = Z
            for q in range(1, (top >> t) + 1, 2):          # value in [q 2^t, (q+1) 2^t)
                hi = (q + 1) << t
                acc = self.OR(acc, self.AND(ge(q << t), Ks[hi - 1] if hi <= top else O))
            out.append(acc)
        self.unkit()
        return out


def bits(x):
    return max(1, x.bit_length())


def bh_column(w, dv, dl, c, init, v_mode, y_mode):
    """One text column. dv: planes of Δv_{.,j-1}; dl: planes of δ(p, t_j).
    Returns (planes of Δv_{.,j}, planes of Δh_{.,j}).
    v_mode: where V is advanced in the loop.
      "printed" - after K <- K ∨ (r = k), as listed on p. 64 (uses the new K = (r ≤ k));
      "swap"    - the last two instructions swapped: V advanced before K (uses K = (r < k)).
                  N in step k then reads V_{k-1}, which equals V_k wherever ↑(r < k) holds,
                  the only positions where Proposition 2 reads it;
      "early"   - V_k computed before N, literally as p. 63 defines it.
    y_mode: the loop set passed to Ind.
      "none"    - (Δv = 0), as on p. 63-64;
      "loop"    - (Δv = 0) ∧ ¬K;
      "decode"  - (Δv = 0), but each indicator is masked with ¬K before Δh_j is assembled;
      "K"       - (Δv = 0), and Δh_j is read from the cumulative vectors K = (r ≤ k)
                  instead of the vectors (r = k);
      "K-"      - as "K", and N is computed without its factor ¬K."""
    nh, ns = bits(2 * c), bits(4 * c)
    S = w.add(dv, dl, ns)                                  # Δv + δ (bounded sum, see header)
    E = w.eq(S, range(2 * c))                              # (Δv + δ = k)
    Zv = w.eq(dv, [0])[0]                                  # (Δv = 0)

    def ind(k, X, Y):                                      # p. 58
        if k == init:
            return w.OR(X, w.AND(Y, w.ADD(w.NOT(X), w.NOT(w.OR(X, Y)))))
        return w.OR(X, w.AND(Y, w.NOT(w.ADD(X, w.OR(X, Y)))))

    r = {0: ind(0, E[0], Zv)}
    K, V = r[0], list(dv)
    Ks, nKs = [K], []
    for k in range(1, 2 * c):
        last = k == 2 * c - 1                              # V is dead after the last step
        sK = w.SHL(K, int(init < k))
        if v_mode == "early":                              # V_k = V_{k-1} + ¬(↑(r<k)), p. 63
            V = w.add_bit(V, w.NOT(sK), ns)
        Vk = w.eq(V, [k])[k]
        if y_mode == "K-":
            N = w.OR(w.AND(sK, Vk), E[k])
        else:
            nK = w.NOT(K)
            nKs.append(nK)
            N = w.AND(w.OR(w.AND(sK, Vk), E[k]), nK)
        r[k] = ind(k, N, w.AND(Zv, nK) if y_mode == "loop" else Zv)
        if v_mode == "swap" and not last:
            V = w.add_bit(V, w.NOT(sK), ns)
        K = w.OR(K, r[k])
        Ks.append(K)
        if y_mode == "decode":
            r[k] = w.AND(r[k], nK)
        if v_mode == "printed" and not last:               # p. 64 as printed
            V = w.add_bit(V, w.NOT(w.SHL(K, int(init < k))), ns)
    if y_mode in ("K", "K-"):
        Dh = w.encode_k(Ks, nKs, nh)
    else:
        r[2 * c] = w.NOT(K)
        Dh = w.encode(r, nh)
    up = w.shift(Dh, init)                                 # ↑_c Δh_j (fill: Δh_{0,j})
    dv_new = w.sub(w.add(up, dv, nh), Dh, nh)              # Δv_j = ↑Δh_j + Δv_{j-1} − Δh_j
    return dv_new, Dh


def shifted(s):
    """Their constraints: one gap cost c, 0 <= δ < 2c. Returns (x, sigma', c)."""
    x = max(0, math.ceil(-min(s.costs) / 2))
    return x, (lambda a, b: s.sigma(a, b) + 2 * x), s.gap + x


def bh_cost(s, p, t, search=False, v_mode="swap", y_mode="loop", counts=None):
    """Cost of p against t with their algorithm on the shifted instance (x, sigma', c).
    C'[i,j] = C[i,j] + x(i+j) satisfies recurrence (7) with sigma', c, and boundaries
    C'[i,0] = ic, C'[0,j] = j*c (global) or j*x (search). So only the initial state
    Δh_{0,j} = C'[0,j] - C'[0,j-1] + c changes: 2c in global mode, c + x in search mode
    (c + x <= 2c because the original gap cost is >= 0). Their paper: search, x = 0, state c."""
    x, sig, c = shifted(s)
    m = len(p)
    w = Words(m)
    nh = bits(2 * c)
    planes = {}
    for a in s.alphabet:
        vals = [sig(pi, a) for pi in p]
        planes[a] = [sum(((v >> k) & 1) << i for i, v in enumerate(vals)) for k in range(bits(4 * c))]
    init = c + x if search else 2 * c                       # Δh_{0,j}
    dv = [0] * nh                                           # Δv_{i,0} = c - c = 0
    score = m * c                                           # C'[m,0]
    best = score - x * m
    per_col = set()
    for j, ch in enumerate(t, 1):
        w0, v0 = w.word, w.vec
        dv, Dh = bh_column(w, dv, planes[ch], c, init, v_mode, y_mode)
        per_col.add((w.word - w0, w.vec - v0))
        dh = sum(((w.val(a) >> (m - 1)) & 1) << k for k, a in enumerate(Dh))
        score += dh - c                                     # C'[m,j] = C'[m,j-1] + Δh_{m,j} - c
        best = min(best, score - x * (m + j))
    if counts is not None:
        counts |= per_col
    if search:
        return best
    return score - x * (m + len(t))


def rand_str(rng, alphabet, n):
    return "".join(rng.choice(alphabet) for _ in range(n))


def main():
    r = Report("bh")
    rng = random.Random(20021)
    s_edit = match_mismatch("edit", 0, 1, 1)

    # Readings of the listing: where V is advanced, which loop test, how Δh_j is assembled.
    s_02 = match_mismatch("(0,1,2)", 0, 1, 2)
    s_03 = match_mismatch("(0,1,3)", 0, 1, 3)
    for vm, ym, should_differ in [("printed", "none", True), ("printed", "loop", True), ("printed", "decode", True),
                                  ("printed", "K", True), ("early", "none", True), ("swap", "none", True),
                                  ("early", "loop", False), ("swap", "loop", False),
                                  ("early", "decode", False), ("swap", "decode", False),
                                  ("early", "K", False), ("swap", "K", False), ("swap", "K-", False)]:
        bad = {}
        for s in (s_edit, s_02, s_03):
            for search in (False, True):
                key = f"{s.label} {'search' if search else 'global'}"
                bad[key] = 0
                for _ in range(300):
                    p, t = rand_str(rng, "ACGT", rng.randint(1, 20)), rand_str(rng, "ACGT", rng.randint(0, 20))
                    if bh_cost(s, p, t, search, v_mode=vm, y_mode=ym) != dp_cost(s, p, t, search):
                        bad[key] += 1
        ok = any(bad.values()) if should_differ else not any(bad.values())
        r.check("BH0", ok, f"V {vm}, loop test {ym}: {'differs from the DP on some pairs' if should_differ else 'agrees with the DP on every pair'}",
                f"differing of 300 pairs: {bad}")

    restricted, reduced = {}, {}

    def run(label, s, ours):
        x, sig, c = shifted(s)
        dmax = max(sig(a, b) for a in s.alphabet for b in s.alphabet)
        counts = {(ym, search): set() for ym in ("K", "loop", "K-") for search in (False, True)}
        for search in (False, True):
            mode = "search" if search else "global"
            differ = {"K": 0, "loop": 0, "K-": 0}
            n_pairs = 400 if c <= 6 else 150
            for _ in range(n_pairs):
                p = rand_str(rng, s.alphabet, rng.randint(1, 62))
                t = rand_str(rng, s.alphabet, rng.randint(0, 60))
                want = dp_cost(s, p, t, search)
                for ym in differ:
                    differ[ym] += bh_cost(s, p, t, search, y_mode=ym, counts=counts[ym, search]) != want
            r.check("BH1", not any(differ.values()),
                    f"{label}: B&H = scalar DP with either loop test (printed, Δh from K; restricted) and without ¬K in N, {mode}",
                    f"{n_pairs} pairs, differing {differ}; x={x} c={c} max δ={dmax}")
            r.check("BH2", all(len(counts[ym, search]) == 1 for ym in differ),
                    f"{label}: per-column count is data-independent, {mode}",
                    f"printed/K {counts['K', search]}, restricted {counts['loop', search]}, without ¬K in N {counts['K-', search]}")
            (kw, kv), = counts["K", search]
            (lw, lv), = counts["loop", search]
            r.check("BH5", kv < lv and kw <= lw,
                    f"{label}: the printed test with Δh from K is the cheaper variant (fewer vector ops, no more plane ops), {mode}",
                    f"vec {kv} < {lv}, word {kw} <= {lw}")
            restricted[label, search] = (lv, lw)
            (_, rv), = counts["K-", search]
            r.check("BH6", rv == kv - 2 * (2 * c - 1),
                    f"{label}: leaving ¬K out of N saves two vector-level operations per level, {mode}", f"{kv} -> {rv}")
            reduced[label, search] = rv
        (gw, gv), = counts["K", False]
        (sw, sv), = counts["K", True]
        if x:
            scope = "shift"
        elif dmax >= 2 * c:
            scope = "δ=2c"
        else:
            scope = "in"
        return (label, scope, x, c, ours, gv, gw, sv, sw)

    # Without the paper source (a checkout without paper/), count the Table 1 schemes only.
    tex = paper_text() if os.path.exists(PAPER) else None
    ours = {row[0]: num(row[3]) for row in table_rows(tex, "tab:ops")} if tex else {}
    rows = [run(label, s, ours.get(label, "-")) for label, s in table1_schemes()]
    for row in table_rows(tex, "tab:scaling") if tex else []:
        for costs, ops in ((row[0], row[3]), (row[4], row[7])):
            a, b, g = parse_cost_triple(costs)
            rows.append(run(costs.replace("$", ""), match_mismatch(costs, a, b, g), num(ops)))

    r.info("")
    r.info("Word operations per column. ours: paper Tables 1 and 3 (global mode).")
    r.info("B&H: last two loop assignments simultaneous; loop test (Δv=0) as printed, Δh read from K;")
    r.info("     on the shifted instance. In brackets: the restricted test (Δv=0)∧¬K with (r=k) as indicators.")
    r.info("     Without ¬K in N (not the listing as printed), global vec: "
           + ", ".join(f"{label} {reduced[label, False]}" for label, *_ in rows))
    r.info("scope: in = 0 <= δ < 2c as published; δ=2c = one cost at 2c, outside their stated range;")
    r.info("       shift = negative costs, run on sigma+2x, c=g+x (our mapping, both modes).")
    r.info(f"{'scheme':<38}{'scope':>6}{'x':>3}{'c':>4}{'ours':>6}{'glob vec':>12}{'glob word':>13}"
           f"{'srch vec':>12}{'srch word':>13}")
    for label, scope, x, c, o, gv, gw, sv, sw in rows:
        (rgv, rgw), (rsv, rsw) = restricted[label, False], restricted[label, True]
        cell = lambda a, b: f"{a} [{b}]"
        r.info(f"{label:<38}{scope:>6}{x:>3}{c:>4}{o:>6}{cell(gv, rgv):>12}{cell(gw, rgw):>13}"
               f"{cell(sv, rsv):>12}{cell(sw, rsw):>13}")
    t1 = rows[:len(table1_schemes())]
    if tex is None:
        r.done()
    # P7 (round 3): Table 1's B&H columns = the global counts above; ‡ exactly on the rows outside their range
    printed = [(row[0], row[4], row[5]) for row in table_rows(tex, "tab:ops")]
    want = [(row[0], f"{row[5]}" + ("$^\\ddagger$" if row[1] == "δ=2c" else ""), f"{row[6]}") for row in t1]
    r.check("BH4", printed == want, "Table 1 columns B&H (vector ops) and B&H planes = global counts; ‡ = cost 2c",
            f"{len(printed)} rows" if printed == want else f"printed {printed} want {want}")
    r.info(f"Table 1, global: ours < B&H vec on {sum(row[4] < row[5] for row in t1)}, = on "
           f"{sum(row[4] == row[5] for row in t1)}, > on {sum(row[4] > row[5] for row in t1)}; "
           f"ours < B&H word on {sum(row[4] < row[6] for row in t1)}, > on {sum(row[4] > row[6] for row in t1)}")

    # R04 E2: the bit-plane realization is below the derived program on BLOSUM62 only
    below = [str(row[0]) for row in t1 if row[6] < row[4]]
    r.check("BH4b", len(t1) == 10 and len(below) == 1 and "BLOSUM62" in below[0]
            and "and 760 for BLOSUM62, below our 799 there and above our count on the other nine schemes of Table~\\ref{tab:ops}" in tex,
            "B&H planes below ours on BLOSUM62 only, above on the other nine schemes of Table 1 (Appendix B)", f"below on {below}")

    # §4.3 read back: every number in the Bergeron-Hamel passage = the counts above (global mode).
    sec = re.search(r"Bergeron and Hamel give no operation count.*?not fewer than the lower bound\.", tex, re.S)
    r.check("BH3", sec is not None, "§4.3 has the Bergeron-Hamel count passage")
    r.check("BH3", "in \\texttt{review/factcheck/bh.py} and \\texttt{bh\\_listing.py} of release v0.2.0-paper" in tex
            and all(os.path.exists(os.path.join(os.path.dirname(__file__), f)) for f in ("bh.py", "bh_listing.py")),
            "§4.3 cites the B&H implementation by path and release (E3/R3); both files exist")
    if sec:
        txt = re.sub(r"\s+", " ", sec.group(0))
        by = {row[0]: row for row in t1}
        e, b, B = by["edit (0, 1, 1)"], by["BitPAl (2, $-3$, $-5$)"], by["BLOSUM62, gap 4"]
        want = (f"it needs {e[5]} operations per column for unit costs and {b[5]} for BitPAl's (2, $-3$, $-5$) weights, "
                f"against our {e[4]} and {b[4]}")
        r.check("BH3", want in txt, "§4.3 lower-bound counts for unit costs and BitPAl (2,-3,-5)", want)
        lower = [row for row in t1 if row[5] < row[4]]
        names = [row[0] for row in lower]
        listed = (f"({', '.join(str(row[5]) for row in lower[:-1])} and {lower[-1][5]} against "
                  f"{', '.join(str(row[4]) for row in lower[:-1])} and {lower[-1][4]})")
        words = {4: "four", 6: "six", 10: "ten"}
        share = f"{words.get(len(lower))} of the {words.get(len(t1))} schemes of Table~\\ref{{tab:ops}}"
        r.check("BH3", names == ["(0, 3, 2)", "BitPAl (2, $-3$, $-5$)", "BitPAl (3, $-4$, $-6$)", "BitPAl (4, $-5$, $-9$)",
                                 "BitPAl (4, $-7$, $-11$)", "BLOSUM62, gap 4"]
                and f"it is lower than ours on {share}, namely (0, 3, 2), BitPAl's four weight sets and BLOSUM62 " + listed in txt,
                "§4.3: rows where their lower bound is below ours, and their numbers", f"{names} {listed}")
        r.check("BH3", f"and the count of its listing is below ours on {share} (Section" in re.sub(r"\s+", " ", tex),
                "§5: share of Table 1 schemes where their lower bound is below ours", share)
        tv = "transition 1, transversion 2, gap 2"
        want = (f"the factor $\\neg K$ in its vector $N$ can be left out, which saves two operations per level and gives {reduced[b[0], False]} for "
                f"BitPAl's (2, $-3$, $-5$) weights and {reduced[tv, False]} for transition/transversion costs, below our {by[tv][4]}")
        r.check("BH3", want in txt and reduced[tv, False] < by[tv][4] and "not a lower bound for the algorithm itself" in txt
                and "a lower bound for any implementation that executes it operation for operation" in re.sub(r"\s+", " ", tex),
                "§4.3 and Table 1 caption: the count is of the listing, operation for operation; variant without ¬K in N (BH1, BH6)", want)
        r.check("BH3", "restrict its loop test to rows not yet below the current level" in txt
                and "leaving the test as printed and reading $\\Delta h$ from the cumulative vector $K$ instead gives the same distances "
                    "and never needs more operations, so the counts below are for that variant" in txt,
                "§4.3 states both loop-test variants and which one is counted (BH1, BH5)")
        want = f"needs {e[6]} and {b[6]} operations, and {B[6]} for BLOSUM62"
        r.check("BH3", want in txt and B[6] < B[4], "§4.3 bit-plane realization counts", want)
        abst = re.sub(r"\s+", " ", re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, re.S).group(1))
        want = f"needs {b[4]} operations on BitPAl's (2, $-3$, $-5$) weights"
        want2 = f"{b[5]} for Bergeron and Hamel's listing counted in vector operations"
        r.check("BH3", want in abst and want2 in abst, "abstract: our and their listing's counts on BitPAl (2,-3,-5)", f"{b[4]} / {b[5]}")
        one_chain = [row for row in rows if row[0].startswith("(0,1,")]
        r.check("BH3", all(row[4] < row[5] for row in one_chain),
                "'grows more slowly when substitution costs span a wide range': on the one-chain family (0,1,c) ours stays lower",
                f"{[(row[0], row[4], row[5]) for row in one_chain]}")
    r.done()


if __name__ == "__main__":
    main()
