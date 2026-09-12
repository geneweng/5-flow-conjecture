"""Modulo-(2p+1) orientation tools.

A multigraph is a list of edges (u, v) with u != v (loops are irrelevant to
modulo orientations and are ignored).  Vertices are 0..n-1.

Main functions
  beta_orientation(n, edges, beta, k)  -> orientation (list of (tail, head)) or None
  has_mod_orientation(n, edges, k)     -> bool   (beta == 0)
  strongly_connected(n, edges, k)      -> (bool, bad_beta or None)
  nz_flow_Z5(n, edges)                 -> bool   (nowhere-zero Z5-flow of the graph)
  triple(edges)                        -> 3G

The encoding: one boolean x_e per edge (True = oriented u->v for the stored
pair (u,v)).  For each vertex v the constraint  sum_{e at v} s_e(x_e) == beta(v)
(mod k), where s_e = +1 if v is the tail, -1 if v is the head, is encoded with a
small mod-k counter automaton (unary state per residue) in CNF.
"""
from itertools import product
from pysat.solvers import Solver


class _CNF:
    def __init__(self):
        self.nv = 0
        self.clauses = []

    def var(self):
        self.nv += 1
        return self.nv

    def add(self, cl):
        self.clauses.append(list(cl))


def _encode(n, edges, beta, k):
    """Return (cnf, edge_vars)."""
    cnf = _CNF()
    x = [cnf.var() for _ in edges]
    inc = [[] for _ in range(n)]           # per vertex: list of (edge index, sign if x true)
    for i, (u, v) in enumerate(edges):
        if u == v:
            continue
        inc[u].append((i, +1))             # x true: u -> v, u is tail: +1 at u
        inc[v].append((i, -1))
    for v in range(n):
        # residue automaton: state r after processing edges 0..j
        # states[j][r] : var meaning residue == r after j edges (one-hot)
        prev = None
        for j, (i, s) in enumerate(inc[v]):
            cur = [cnf.var() for _ in range(k)]
            cnf.add(cur)                                   # at least one
            for r in range(k):
                for r2 in range(r + 1, k):
                    cnf.add([-cur[r], -cur[r2]])           # at most one
            for r in range(k):
                # x true: contribution s ; x false: contribution -s
                rt = (r + s) % k
                rf = (r - s) % k
                if prev is None:
                    if r != 0:
                        continue
                    # from residue 0
                    cnf.add([-x[i], cur[rt]])
                    cnf.add([x[i], cur[rf]])
                else:
                    cnf.add([-prev[r], -x[i], cur[rt]])
                    cnf.add([-prev[r], x[i], cur[rf]])
            prev = cur
        target = beta[v] % k
        if prev is None:
            if target != 0:
                cnf.add([])                                # unsatisfiable
        else:
            cnf.add([prev[target]])
    return cnf, x


def beta_orientation(n, edges, beta, k=5, solver="cadical153"):
    """Orientation D with d+(v) - d-(v) == beta(v) (mod k) for all v, or None."""
    if sum(beta) % k != 0:
        return None
    cnf, x = _encode(n, edges, beta, k)
    if any(len(c) == 0 for c in cnf.clauses):
        return None
    with Solver(name=solver, bootstrap_with=cnf.clauses) as s:
        if not s.solve():
            return None
        model = set(l for l in s.get_model() if l > 0)
    return [(u, v) if x[i] in model else (v, u) for i, (u, v) in enumerate(edges)]


def has_mod_orientation(n, edges, k=5):
    return beta_orientation(n, edges, [0] * n, k) is not None


def strongly_connected(n, edges, k=5):
    """Is G strongly Z_k-connected (every zero-sum beta achievable)?  Returns (bool, witness)."""
    for tail in product(range(k), repeat=n - 1):
        beta = list(tail) + [(-sum(tail)) % k]
        if beta_orientation(n, edges, beta, k) is None:
            return False, beta
    return True, None


def triple(edges, t=3):
    return [e for e in edges for _ in range(t)]


def nz_flow_Z5(n, edges):
    """Nowhere-zero Z5-flow exists iff 3G has a modulo 5-orientation."""
    return has_mod_orientation(n, triple(edges), 5)


def petersen():
    outer = [(i, (i + 1) % 5) for i in range(5)]
    spokes = [(i, i + 5) for i in range(5)]
    inner = [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
    return 10, outer + spokes + inner


if __name__ == "__main__":
    n, E = petersen()
    print("Petersen: nowhere-zero Z5-flow (via 3G mod-5 orientation):", nz_flow_Z5(n, E))
    # mod-3 orientation of G itself <=> nowhere-zero Z3-flow (values +-1); Petersen: False, K33: True
    print("Petersen: mod-3 orientation (= nowhere-zero Z3-flow):", has_mod_orientation(n, E, 3))
    print("K_{3,3}: mod-3 orientation (= nowhere-zero Z3-flow):", has_mod_orientation(6, [(i, 3 + j) for i in range(3) for j in range(3)], 3))
    print("Petersen: mod-5 orientation of G itself (impossible, cubic):", has_mod_orientation(n, E, 5))
    # K2 with a parallel edges: strongly Z5-connected iff a >= 4
    for a in range(2, 6):
        print(f"{a}K2 strongly Z5-connected:", strongly_connected(2, [(0, 1)] * a, 5)[0])
    # T_{2,2,3} and T_{1,3,3} should fail, T_{2,3,3} should succeed (a+b+c>=8, delta>=4)
    def T(a, b, c):
        return 3, [(0, 1)] * a + [(1, 2)] * b + [(0, 2)] * c
    for abc in [(2, 2, 3), (1, 3, 3), (2, 3, 3), (4, 4, 0), (3, 3, 2)]:
        n3, E3 = T(*abc)
        print(f"T{abc} strongly Z5-connected:", strongly_connected(n3, E3, 5))


# ---- fast exact strong Z_k-connectivity for small n via sumsets over Z_k^V ----
def achievable_boundaries(n, edges, k=5):
    """Boolean array over Z_k^n: which (d+ - d-) mod k vectors are realised by orientations."""
    import numpy as np
    R = np.zeros((k,) * n, dtype=bool)
    R[(0,) * n] = True
    for (u, v) in edges:
        if u == v:
            continue
        R = np.roll(np.roll(R, 1, axis=u), -1, axis=v) | np.roll(np.roll(R, -1, axis=u), 1, axis=v)
    return R


def strongly_connected_fast(n, edges, k=5):
    import numpy as np
    R = achievable_boundaries(n, edges, k)
    idx = np.indices((k,) * n).reshape(n, -1).T
    zero_sum = (idx.sum(axis=1) % k == 0)
    reach = R.reshape(-1)
    missing = np.where(zero_sum & ~reach)[0]
    return len(missing) == 0, (list(idx[missing[0]]) if len(missing) else None)
