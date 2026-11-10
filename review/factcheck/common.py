"""Shared helpers for the independent fact-check scripts.

Independent of the Zig implementation: nothing here reads or imports src/.
Inputs are the paper source (paper/main.tex) and, where a script says so,
benchmark outputs under bench/ and review/.
"""
import os
import re
import sys

sys.stdout.reconfigure(encoding="utf-8")  # Windows consoles default to cp1252

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PAPER = os.path.join(ROOT, "paper", "main.tex")

# ---------------------------------------------------------------- reporting

class Report:
    """Collects OK/FAIL lines per ledger ID; exit code 1 if any FAIL."""

    def __init__(self, name):
        self.name, self.fails, self.oks = name, [], 0

    def check(self, ledger_id, ok, what, detail=""):
        tag = "OK  " if ok else "FAIL"
        print(f"[{tag}] {ledger_id:<5} {what}" + (f" | {detail}" if detail else ""))
        if ok:
            self.oks += 1
        else:
            self.fails.append(f"{ledger_id} {what}")
        return ok

    def info(self, text):
        print(f"       {text}")

    def done(self):
        status = "PASS" if not self.fails else "FAIL"
        summary = f"{self.name}: {status} ({self.oks} ok, {len(self.fails)} fail)"
        if self.fails:
            summary += " -- " + "; ".join(self.fails)
        print(summary)
        sys.exit(0 if not self.fails else 1)


# ---------------------------------------------------------------- paper text

def paper_text():
    with open(PAPER, encoding="utf-8") as f:
        return f.read()


def table_rows(tex, label):
    """Rows of the tabular in the float carrying \\label{label}, split on '&'."""
    m = re.search(r"\\label\{%s\}(.*?)\\end\{tabular\}" % re.escape(label), tex, re.S)
    if not m:
        raise ValueError(f"table {label} not found")
    body = m.group(1).split("\\midrule", 1)[1].split("\\bottomrule")[0]
    rows = []
    for line in body.split("\\\\"):
        line = line.strip()
        if line:
            rows.append([c.strip() for c in line.split("&")])
    return rows


def num(cell):
    """'$-3$' -> -3, '14' -> 14, '0.73' -> 0.73, '--' -> None."""
    c = cell.replace("$", "").replace("\\,", "").replace("{", "").replace("}", "").strip()
    if c in ("--", "-", ""):
        return None
    return float(c) if "." in c else int(c)


def verbatim_blocks(tex):
    return re.findall(r"\\begin\{verbatim\}\n(.*?)\\end\{verbatim\}", tex, re.S)


# ---------------------------------------------------------------- schemes

DNA = "ACGT"
AA20 = "ARNDCQEGHILKMFPSTWYV"

# BLOSUM62 as published by NCBI: https://ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62
# ("Matrix made by matblas from blosum62.iij", 1/2-bit units). Only the 20
# standard amino acids are used (rows/columns B, Z, X, * dropped).
_BLOSUM62_ROWS = """
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
BLOSUM62 = {}
for _line in _BLOSUM62_ROWS.strip().splitlines():
    _f = _line.split()
    for _b, _s in zip(AA20, _f[1:]):
        BLOSUM62[(_f[0], _b)] = int(_s)
assert len(BLOSUM62) == 400 and all(BLOSUM62[a, b] == BLOSUM62[b, a] for a in AA20 for b in AA20)

PURINES = set("AG")


class Scheme:
    """Costs to minimize: sigma(a, b) and a linear gap g > 0 (paper §3, Setting)."""

    def __init__(self, label, alphabet, sigma, gap):
        self.label, self.alphabet, self.sigma, self.gap = label, alphabet, sigma, gap
        self.costs = sorted({sigma(a, b) for a in alphabet for b in alphabet})

    def __repr__(self):
        return self.label


def match_mismatch(label, match, mismatch, gap):
    return Scheme(label, DNA, lambda a, b: match if a == b else mismatch, gap)


def ts_tv():
    def s(a, b):
        if a == b:
            return 0
        return 1 if (a in PURINES) == (b in PURINES) else 2
    return Scheme("ts/tv (0, ts 1, tv 2; gap 2)", DNA, s, 2)


def blosum62_gap4():
    return Scheme("BLOSUM62, gap 4", AA20, lambda a, b: -BLOSUM62[a, b], 4)


def table1_schemes():
    """(Table 1 row label prefix as printed, scheme). BitPAl rows: score weights negated."""
    return [
        ("edit (0, 1, 1)", match_mismatch("edit (0,1,1)", 0, 1, 1)),
        ("indel (0, 2, 1)", match_mismatch("indel (0,2,1)", 0, 2, 1)),
        ("(0, 1, 2)", match_mismatch("(0,1,2)", 0, 1, 2)),
        ("(0, 3, 2)", match_mismatch("(0,3,2)", 0, 3, 2)),
        ("BitPAl (2, $-3$, $-5$)", match_mismatch("(-2,3,5)", -2, 3, 5)),
        ("BitPAl (3, $-4$, $-6$)", match_mismatch("(-3,4,6)", -3, 4, 6)),
        ("BitPAl (4, $-5$, $-9$)", match_mismatch("(-4,5,9)", -4, 5, 9)),
        ("BitPAl (4, $-7$, $-11$)", match_mismatch("(-4,7,11)", -4, 7, 11)),
        ("transition 1, transversion 2, gap 2", ts_tv()),
        ("BLOSUM62, gap 4", blosum62_gap4()),
    ]


def parse_cost_triple(cell):
    """'(0,1,3)' or '($-1$,3,4)' -> (0, 1, 3)."""
    c = cell.replace("$", "").strip("() ")
    return tuple(int(x) for x in c.split(","))


# ---------------------------------------------------------------- DP and Δ

def dp_differences(s, p, t, out):
    """Global DP (paper §3): add every v(i,j) and h(i,j) to `out`; return D(m,n)."""
    g, sig = s.gap, s.sigma
    prev = [j * g for j in range(len(t) + 1)]
    for j in range(1, len(t) + 1):
        out.add(prev[j] - prev[j - 1])            # h(0, j)
    for i in range(1, len(p) + 1):
        cur = [i * g]
        out.add(cur[0] - prev[0])                 # v(i, 0)
        pi = p[i - 1]
        for j in range(1, len(t) + 1):
            d = min(prev[j - 1] + sig(pi, t[j - 1]), prev[j] + g, cur[j - 1] + g)
            cur.append(d)
            out.add(d - prev[j])                  # v(i, j)
            out.add(d - cur[j - 1])               # h(i, j)
        prev = cur
    return prev[-1]


def dp_cost(s, p, t, search=False):
    """Plain DP cost. search=True: D(0,j)=0 and the minimum over end columns."""
    g, sig = s.gap, s.sigma
    col = [i * g for i in range(len(p) + 1)]
    best = col[-1]
    for c in t:
        diag, col[0] = col[0], (0 if search else col[0] + g)
        for i in range(1, len(p) + 1):
            cell = min(diag + sig(p[i - 1], c), col[i] + g, col[i - 1] + g)
            diag, col[i] = col[i], cell
        best = min(best, col[-1])
    return best if search else col[-1]


def all_strings(alphabet, max_len):
    out = [""]
    frontier = [""]
    for _ in range(max_len):
        frontier = [w + a for w in frontier for a in alphabet]
        out += frontier
    return out


def related_pair(rng, alphabet, m, n):
    """A long pattern and a text derived from it: optionally periodic, with edits and a shift.

    Uniform random pairs rarely produce the large-gap differences (they need
    near-identical sequences offset by gaps), so the sampler also draws these.
    """
    a = rng.randint(2, min(4, len(alphabet)))
    base = [rng.choice(alphabet[:a]) for _ in range(max(m, n) + 20)]
    if rng.random() < 0.5:
        per = rng.randint(1, 6)
        base = [base[k % per] for k in range(len(base))]
    t = list(base)
    for _ in range(rng.randint(0, 12)):
        k, op = rng.randrange(len(t)), rng.random()
        if op < 0.4:
            t[k] = rng.choice(alphabet)
        elif op < 0.7:
            t.insert(k, rng.choice(alphabet))
        else:
            del t[k]
    shift = rng.randint(0, 10)
    t = (t[shift:] + [rng.choice(alphabet) for _ in range(n)])[:n]
    return "".join(base[:m]), "".join(t)


def batch_differences(s, pairs, out):
    """dp_differences for many equal-length pairs at once (numpy over pairs)."""
    import numpy as np
    sym = {c: k for k, c in enumerate(s.alphabet)}
    P = np.array([[sym[c] for c in p] for p, _ in pairs])
    T = np.array([[sym[c] for c in t] for _, t in pairs])
    cost = np.array([[s.sigma(a, b) for b in s.alphabet] for a in s.alphabet])
    g, (B, m), n = s.gap, P.shape, T.shape[1]
    prev = np.tile(np.arange(n + 1) * g, (B, 1))
    out.update(np.unique(prev[:, 1:] - prev[:, :-1]).tolist())
    for i in range(m):
        cur = np.empty_like(prev)
        cur[:, 0] = prev[:, 0] + g
        c = cost[P[:, i:i + 1], T]
        for j in range(1, n + 1):
            cur[:, j] = np.minimum(np.minimum(prev[:, j - 1] + c[:, j - 1], prev[:, j] + g), cur[:, j - 1] + g)
        out.update(np.unique(cur - prev).tolist())
        out.update(np.unique(cur[:, 1:] - cur[:, :-1]).tolist())
        prev = cur


def empirical_delta(s, exhaustive_len, random_pairs, random_len, rng, long_batches=4):
    """Values of v and h seen in DPs (a lower bound on Δ): exhaustive short pairs,
    uniform random pairs, and batches of long related pairs."""
    seen = set()
    words = all_strings(s.alphabet, exhaustive_len)
    for p in words:
        for t in words:
            dp_differences(s, p, t, seen)
    for _ in range(random_pairs):
        m, n = rng.randint(1, random_len), rng.randint(0, random_len)
        p = "".join(rng.choice(s.alphabet) for _ in range(m))
        t = "".join(rng.choice(s.alphabet) for _ in range(n))
        dp_differences(s, p, t, seen)
    for _ in range(long_batches):
        m, n = rng.choice([20, 40, 60, 90]), rng.choice([20, 40, 60, 90])
        batch_differences(s, [related_pair(rng, s.alphabet, m, n) for _ in range(300)], seen)
    return seen


def cell(s_val, v, c, g):
    """Paper §3.1: d = min(c, s+g, v+g); h' = d - v; v' = d - s."""
    d = min(c, s_val + g, v + g)
    return d - v, d - s_val


def fixpoint_delta(s):
    """Close the boundary values {g} under the cell map, as §3.1 describes."""
    delta = {s.gap}
    while True:
        new = set(delta)
        for sv in delta:
            for v in delta:
                for c in s.costs:
                    new.update(cell(sv, v, c, s.gap))
        if new == delta:
            return delta
        delta = new


def cap_of(s, delta):
    """Theorem 2: cap = c_{C-1} - g. It is at least min Δ on the closure (§3.3, checked by theorems.py D12)."""
    assert s.costs[-1] - s.gap >= min(delta), (s, sorted(delta))
    return s.costs[-1] - s.gap


def chains_of(delta, cap):
    """§3.3: #{τ ∈ Δ : min Δ < τ <= cap}."""
    lo = min(delta)
    return sum(1 for t in delta if lo < t <= cap)


def level_of(x, delta):
    """[s >= x] as a level of Δ: None if constant (x <= min Δ or x > max Δ)."""
    if x <= min(delta) or x > max(delta):
        return None
    return min(y for y in delta if y >= x)


def has_cycle(edges, nodes):
    """Cycle among distinct nodes (self-loops ignored), by Kahn's algorithm."""
    indeg = {n: 0 for n in nodes}
    adj = {n: set() for n in nodes}
    for a, b in edges:
        if a != b and b not in adj[a]:
            adj[a].add(b)
            indeg[b] += 1
    queue = [n for n in nodes if indeg[n] == 0]
    seen = 0
    while queue:
        n = queue.pop()
        seen += 1
        for b in adj[n]:
            indeg[b] -= 1
            if indeg[b] == 0:
                queue.append(b)
    return seen != len(nodes)
