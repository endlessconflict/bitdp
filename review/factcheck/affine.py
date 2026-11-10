"""Affine-gap outlook of §5 (stage 07). Every statement in the new paragraph is checked here.

Conventions (paper §5 paragraph): costs to minimize, a gap of length k costs o + k*e.
Gotoh's recurrence, global mode:
  E(i,j) = min(E(i,j-1) + e, H(i,j-1) + o + e)          horizontal gap
  F(i,j) = min(F(i-1,j) + e, H(i-1,j) + o + e)          vertical gap (the column's chain)
  D_i    = min(H(i-1,j-1) + sigma(p_i, t_j), E(i,j))     everything not in the vertical chain
  H(i,j) = min(D_i, F(i,j))
"""
import itertools
import math
import random

from common import Report, has_cycle, level_of

INF = math.inf
SEED = 77


def comp(f, g):
    """(f o g) for f = (a, b), g = (a2, b2), f(x) = min(a, x + b).  [§5: closure formula]"""
    (a, b), (a2, b2) = f, g
    return (min(a, a2 + b), b2 + b)


def apply(f, x):
    return min(f[0], x + f[1])


def gotoh(p, t, sigma, o, e):
    """Full matrices H, E, F (lists of rows) for global affine alignment."""
    m, n = len(p), len(t)
    H = [[INF] * (n + 1) for _ in range(m + 1)]
    E = [[INF] * (n + 1) for _ in range(m + 1)]
    F = [[INF] * (n + 1) for _ in range(m + 1)]
    D = [[INF] * (n + 1) for _ in range(m + 1)]
    H[0][0] = 0
    for j in range(1, n + 1):
        H[0][j] = E[0][j] = o + j * e
    for i in range(1, m + 1):
        H[i][0] = F[i][0] = o + i * e
    for j in range(1, n + 1):
        for i in range(1, m + 1):
            E[i][j] = min(E[i][j - 1] + e, H[i][j - 1] + o + e)
            F[i][j] = min(F[i - 1][j] + e, H[i - 1][j] + o + e)
            D[i][j] = min(H[i - 1][j - 1] + sigma(p[i - 1], t[j - 1]), E[i][j])
            H[i][j] = min(D[i][j], F[i][j])
    return H, E, F, D


def plain_affine(p, t, sigma, o, e):
    """Independent check of gotoh(): minimum over all alignments by memoized recursion on (i, j, state)."""
    from functools import lru_cache

    @lru_cache(None)
    def best(i, j, state):          # cheapest alignment of p[:i], t[:j] whose last column is:
        if i == 0 and j == 0:       # 0 = a substitution (or the empty start), 1 = p_i against a gap,
            return 0 if state == 0 else INF          # 2 = t_j against a gap
        if state == 0:
            if i == 0 or j == 0:
                return INF
            return min(best(i - 1, j - 1, s) for s in (0, 1, 2)) + sigma(p[i - 1], t[j - 1])
        if state == 1:
            return min(best(i - 1, j, 1) + e, min(best(i - 1, j, s) for s in (0, 2)) + o + e) if i > 0 else INF
        return min(best(i, j - 1, 2) + e, min(best(i, j - 1, s) for s in (0, 1)) + o + e) if j > 0 else INF

    return min(best(len(p), len(t), s) for s in (0, 1, 2))


def main():
    r = Report("affine.py")
    rng = random.Random(SEED)
    r.info(f"seed {SEED}")

    # 1. Closure under composition (brute force, a in Z ∪ {+inf}).
    vals, bs, xs = [-5, -1, 0, 2, 7, INF], [-4, -1, 0, 3, 6], range(-12, 13)
    ok = all(apply(comp((a, b), (a2, b2)), x) == apply((a, b), apply((a2, b2), x))
             for a, b, a2, b2, x in itertools.product(vals, bs, vals, bs, xs))
    r.check("7.1", ok, "f_{a,b} o f_{a',b'} = f_{min(a, a'+b), b'+b} (grid incl. a = +inf)")
    # 2. Associativity (follows from function composition; spot-checked).
    fs = [(rng.choice(vals), rng.choice(bs)) for _ in range(300)]
    ok = all(comp(comp(f, g), h) == comp(f, comp(g, h)) for f, g, h in zip(fs, fs[1:], fs[2:]))
    r.check("7.1", ok, "composition is associative (300 random triples); identity is (+inf, 0)")

    # Gotoh itself, against an independent recursion.
    sigma = lambda a, b: 0 if a == b else 1
    o, e = 2, 1
    bad = 0
    for _ in range(300):
        p = "".join(rng.choice("ACG") for _ in range(rng.randint(1, 7)))
        t = "".join(rng.choice("ACG") for _ in range(rng.randint(1, 7)))
        bad += gotoh(p, t, sigma, o, e)[0][-1][-1] != plain_affine(p, t, sigma, o, e)
    r.check("7.4", bad == 0, "Gotoh recurrence = independent recursion over alignments (300 pairs, o=2, e=1)", f"{bad} disagreements")

    # 3. The vertical chain in absolute scores: F(i,j) = min(a_i, F(i-1,j) + e), a_i = D_{i-1} + o + e.  [§5 formula]
    ok, pairs, cols = True, 0, []
    for _ in range(200):
        p = "".join(rng.choice("ACGT") for _ in range(rng.randint(2, 40)))
        t = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 30)))
        H, E, F, D = gotoh(p, t, sigma, o, e)
        for j in range(1, len(t) + 1):
            for i in range(2, len(p) + 1):
                ok &= F[i][j] == min(D[i - 1][j] + o + e, F[i - 1][j] + e)
                pairs += 1
            cols.append(([D[i - 1][j] + o + e for i in range(2, len(p) + 1)], F[1][j], [F[i][j] for i in range(1, len(p) + 1)]))
    r.check("7.3", ok, "F(i,j) = min(D_{i-1} + o + e, F(i-1,j) + e): each row applies x -> min(a_i, x + e)",
            f"{pairs} cells")

    # 4. A prefix scan over the monoid reproduces the chain in ceil(log2 w) rounds (Hillis-Steele).  [§5 O(log w)]
    ok, worst = True, 0
    for a_list, f1, want in cols:
        maps = [(INF, 0)] + [(a, e) for a in a_list]          # row 1 starts from F(1,j); identity first
        w, rounds, step = len(maps), 0, 1
        while step < w:
            maps = [comp(maps[k], maps[k - step]) if k >= step else maps[k] for k in range(w)]
            step, rounds = step * 2, rounds + 1
        ok &= [apply(f, f1) for f in maps] == want and rounds == math.ceil(math.log2(w)) if w > 1 else True
        worst = max(worst, rounds)
    r.check("7.5", ok, "prefix scan of the row maps = sequential chain, in ceil(log2 w) composition rounds",
            f"{len(cols)} columns, up to {worst} rounds")

    # 5. What breaks the cascade: relative to the previous column, level τ of the gap state copies
    #    level τ + v - e, and v can exceed e (copy downward) as well as fall below it (copy upward).
    ident, down, up, levels, shifts = True, None, None, set(), {}
    for _ in range(400):
        p = "".join(rng.choice("ACGT") for _ in range(rng.randint(2, 25)))
        t = "".join(rng.choice("ACGT") for _ in range(rng.randint(2, 25)))
        H, E, F, D = gotoh(p, t, sigma, o, e)
        for j in range(2, len(t) + 1):
            for i in range(2, len(p) + 1):
                f_prev = F[i - 1][j] - H[i - 1][j - 1]               # gap state of the row above
                s = H[i - 1][j] - H[i - 1][j - 1]
                v = H[i][j - 1] - H[i - 1][j - 1]                   # incoming vertical difference
                f_i = F[i][j] - H[i][j - 1]
                ident &= f_i == min(f_prev + e - v, s + o + e - v)
                levels.update((f_prev, f_i))
                shift = v - e                                       # source level = τ + shift
                shifts.setdefault(shift, (p, t, i, j, v))
                if shift > 0 and down is None:
                    down = (p, t, i, j, v, shift)
                if shift < 0 and up is None:
                    up = (p, t, i, j, v, shift)
    r.check("7.4", ident, "relative form: f_i = min(f_{i-1} + e - v, s + o + e - v), so [f_i>=τ] reads [f_{i-1} >= τ + v - e]")
    r.check("7.4", down is not None and up is not None, "copies run both ways: real cells with v > e and with v < e",
            "" if down is None or up is None else
            f"down: p={down[0]} t={down[1]} (i={down[2]}, j={down[3]}): v={down[4]} > e={e}, level τ copies τ+{down[5]}; "
            f"up: p={up[0]} t={up[1]} (i={up[2]}, j={up[3]}): v={up[4]}, level τ copies τ{up[5]}")
    # Copy graph over the gap-state levels that occur, with every shift that occurs (same rule as
    # theorems.py D04 for linear gaps): edge from level [x >= τ + d] to level τ.
    lv = sorted(x for x in levels if x != INF)
    nodes = [t0 for t0 in lv if t0 > lv[0]]
    edges = {(level_of(t0 + d, set(lv)), t0) for t0 in nodes for d in shifts} - {(None, t0) for t0 in nodes}
    r.check("7.4", has_cycle(edges, nodes), "the gap-state copy graph has a cycle (cascade not acyclic)",
            f"shifts v-e seen: {sorted(shifts)}; {len(nodes)} levels, {len(edges)} edges")
    cyc = find_cycle(edges)
    r.info("example cycle: " + " -> ".join(f"[f>={x}]" for x in cyc) +
           "  (each arrow: the target level reads the source level of the row above)")
    for d in sorted({cyc[k + 1] - cyc[k] for k in range(len(cyc) - 1)}):
        sp = shifts[-d]          # edge src -> dst means dst copies dst + shift = src, so shift = src - dst
        r.info(f"  shift {-d:+d} realised at p={sp[0]} t={sp[1]} (i={sp[2]}, j={sp[3]}), v={sp[4]}, e={e}")
    # 6. ISS-0118 (R03): with the horizontal difference of F as state, the affine column has the linear-gap form.
    #    phi_i = F(i,j) - F(i,j-1);  psi = F(i-1,j-1) - H(i-1,j-1) >= 0 (previous column, row i-1):
    #    phi_i = min(alpha_i, phi_{i-1} + beta_i),  beta_i = max(0, psi - o) >= 0,
    #    alpha_i = (D(i-1,j) - H(i-1,j-1)) + max(0, o - psi);  h(i,j) = min(D(i,j) - H(i,j-1), phi_i + F(i,j-1) - H(i,j-1)).
    schemes = [("unit, o=1, e=1", (0, 1), 1, 1), ("(0,2), o=3, e=1", (0, 2), 3, 1), ("(0,10), o=12, e=3", (0, 10), 12, 3),
               ("(0,6), o=4, e=3", (0, 6), 4, 3), ("(0,3), o=0, e=2", (0, 3), 0, 2), ("(0,1,2), o=3, e=1", (0, 1, 2), 3, 1),
               ("(-2,3), o=4, e=1", (-2, 3), 4, 1), ("(0,1), o=2, e=1", (0, 1), 2, 1), ("(0,5), o=6, e=1", (0, 5), 6, 1)]
    ident = hform = nonneg = bounded = True
    cells, most = 0, 0
    for name, costs, o2, e2 in schemes:
        sg = (lambda cs: lambda a, b: cs[0] if a == b else cs[1 + (ord(a) + ord(b)) % (len(cs) - 1)])(costs)
        seen = set()
        for _ in range(250):
            p = "".join(rng.choice("ACGT") for _ in range(rng.randint(2, 22)))
            t = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 22)))
            if rng.random() < 0.5:      # near-copies of the pattern: long gap runs
                t = "".join(ch if rng.random() > 0.15 else rng.choice("ACGT") for ch in (p * 2)[rng.randint(0, 3):][:len(t)])
            H, E, F, D = gotoh(p, t, sg, o2, e2)
            for j in range(1, len(t) + 1):
                for i in range(1, len(p) + 1):
                    phi = F[i][j] - F[i][j - 1]
                    hform &= H[i][j] - H[i][j - 1] == min(D[i][j] - H[i][j - 1], phi + F[i][j - 1] - H[i][j - 1])
                    seen.add(phi)
                    if i >= 2:
                        psi = F[i - 1][j - 1] - H[i - 1][j - 1]
                        beta = max(0, psi - o2)
                        alpha = D[i - 1][j] - H[i - 1][j - 1] + max(0, o2 - psi)
                        ident &= phi == min(alpha, F[i - 1][j] - F[i - 1][j - 1] + beta)
                        nonneg &= psi >= 0 and beta >= 0
                    cells += 1
        bounded &= min(seen) >= min(costs[0], 0) - (o2 + e2) and max(seen) <= o2 + e2
        most = max(most, len(seen))
    r.check("7.6", ident, "state phi_i = F(i,j) - F(i,j-1): phi_i = min(alpha_i, phi_{i-1} + beta_i), beta_i = max(0, psi - o)",
            f"{len(schemes)} schemes, {cells} cells")
    r.check("7.6", nonneg, "the shift beta_i is never negative (psi >= 0): level tau copies a level <= tau, so no cycle")
    r.check("7.6", hform, "h(i,j) = min(D(i,j) - H(i,j-1), phi_i + psi_i): the row's cost follows from phi")
    r.check("7.6", bounded, "phi stays in [min(c_0, 0) - o - e, o + e]", f"at most {most} distinct values in a scheme")
    r.done()


def find_cycle(edges):
    """A shortest cycle among distinct nodes (self-loops ignored), as a node list closing on its start."""
    adj = {}
    for a, b in edges:
        if a != b:
            adj.setdefault(a, set()).add(b)
    best = []
    for start in sorted(adj):                      # BFS from each node back to itself
        prev, frontier = {start: None}, [start]
        while frontier and start not in [n for f in frontier for n in adj.get(f, ())]:
            nxt = []
            for f in frontier:
                for n in sorted(adj.get(f, ())):
                    if n not in prev:
                        prev[n] = f
                        nxt.append(n)
            frontier = nxt
        last = next((f for f in frontier if start in adj.get(f, ())), None)
        if last is None:
            continue
        path = [last]
        while prev[path[-1]] is not None:
            path.append(prev[path[-1]])
        cyc = path[::-1] + [start]
        if not best or len(cyc) < len(best):
            best = cyc
    return best


if __name__ == "__main__":
    main()
