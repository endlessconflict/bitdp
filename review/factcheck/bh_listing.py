"""The vector algorithm of Bergeron and Hamel, "Vector algorithms for approximate
string matching", IJFCS 13(1):53-65, 2002, Section 3.2 (p. 64), as we implemented it,
checked against the dynamic program (7) on random inputs.

Their setting: search mode, one indel cost c, substitution costs 0 <= delta < 2c.
Vectors are Python lists indexed by row (index 0 = row 1); boolean vectors are
lists of 0/1, and +_b is binary addition with carries running from row 1 upward.

Our implementation choices (marked below):
  (a) The last two assignments of the loop are read as simultaneous: V is advanced
      with K = (r < k), before K absorbs (r = k).
  (b) The loop test given to Ind, in two variants that give the same distances:
      "restricted" - (dv = 0) AND NOT K, so that the vectors (r = k) are disjoint and
                     Dh_j is assembled from them as indicator vectors. The test
                     (dv = 0) AND (dv + delta > k) gives the same vectors (r = k).
      "printed"    - (dv = 0), the test as printed. The cumulative vectors
                     K = (r <= k) are the same as with the restricted test, and
                     Dh_j is read from them: row i has Dh = #{k : K_k[i] = 0}.

Each column of the restricted variant also asserts that the vectors (r = k) are the
indicator vectors of Dh_j, exactly one of them set in each row.

Run: python bh_listing.py   (silent on success)
"""
import random


def to_int(x):
    return sum(b << i for i, b in enumerate(x))


def from_int(n, m):
    return [(n >> i) & 1 for i in range(m)]


def add_b(x, y):                       # x +_b y, last carry dropped
    return from_int(to_int(x) + to_int(y), len(x))


def NOT(x):
    return [1 - b for b in x]


def AND(x, y):
    return [a & b for a, b in zip(x, y)]


def OR(x, y):
    return [a | b for a, b in zip(x, y)]


def up(x, a):                          # their up-arrow: a x1 ... x_{m-1}
    return [a] + x[:-1]


def ind(s, X, Y, initial):             # p. 58
    if s == initial:
        return OR(X, AND(Y, add_b(NOT(X), NOT(OR(X, Y)))))
    return OR(X, AND(Y, NOT(add_b(X, OR(X, Y)))))


def column(dv, delta, c, test="restricted"):
    """One text column: (Dv_{j-1}, delta(p, t_j)) -> (Dv_j, Dh_j). Initial state c."""
    m = len(dv)
    s = [min(v + d, 2 * c) for v, d in zip(dv, delta)]           # bounded sum, p. 65
    eq = lambda vec, k: [int(x == k) for x in vec]
    zero = eq(dv, 0)
    r = {0: ind(0, eq(s, 0), zero, c)}
    K, V = r[0], list(dv)
    Ks = [K]
    for k in range(1, 2 * c):
        sK = up(K, int(c < k))
        N = AND(OR(AND(sK, eq(V, k)), eq(s, k)), NOT(K))
        if test == "restricted":                                  # choice (b)
            r[k] = ind(k, N, AND(zero, NOT(K)), c)
            assert r[k] == ind(k, N, AND(zero, [int(x > k) for x in s]), c)
        else:
            r[k] = ind(k, N, zero, c)
        V = [v + 1 - b for v, b in zip(V, sK)]                    # choice (a): before K changes
        K = OR(K, r[k])
        Ks.append(K)
    if test == "restricted":
        r[2 * c] = NOT(K)
        # (r = k) must be the indicator vectors of Dh_j: exactly one of them is 1 in each row.
        assert all(sum(r[k][i] for k in r) == 1 for i in range(m)), "(r = k) vectors overlap"
        dh = [sum(k * r[k][i] for k in r) for i in range(m)]
    else:
        dh = [sum(1 - Kk[i] for Kk in Ks) for i in range(m)]     # Dh read from K
    dv_new = [a + v - h for a, v, h in zip(up(dh, c), dv, dh)]   # Dv_j = up_c Dh_j + Dv_{j-1} - Dh_j
    return dv_new, dh, Ks


def search(p, t, delta, c):
    """min over j of C[m, j], C[0, j] = 0, C[i, 0] = ic. Runs both variants of the loop
    test side by side and asserts that they give the same K and the same Dh_j."""
    dv, score = [0] * len(p), len(p) * c
    best = score
    for ch in t:
        d = [delta(a, ch) for a in p]
        dv_p, dh_p, Ks_p = column(dv, d, c, "printed")
        dv, dh, Ks = column(dv, d, c)
        assert (dv_p, dh_p, Ks_p) == (dv, dh, Ks), "the two loop tests differ"
        score += dh[-1] - c                                       # C[m, j] = C[m, j-1] + Dh_{m,j} - c
        best = min(best, score)
    return best


def dp(p, t, delta, c):
    col = [i * c for i in range(len(p) + 1)]
    best = col[-1]
    for ch in t:
        diag, col[0] = col[0], 0
        for i in range(1, len(p) + 1):
            diag, col[i] = col[i], min(diag + delta(p[i - 1], ch), col[i] + c, col[i - 1] + c)
        best = min(best, col[-1])
    return best


if __name__ == "__main__":
    rng = random.Random(2002)
    for c in range(1, 6):
        for _ in range(200):
            table = {(a, b): 0 if a == b else rng.randrange(2 * c) for a in "ACGT" for b in "ACGT"}
            delta = lambda a, b: table[a, b]
            p = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 30)))
            t = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 40)))
            assert search(p, t, delta, c) == dp(p, t, delta, c), (c, p, t, table)
