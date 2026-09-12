"""Exact Z_k-connectivity test for small graphs by sumsets over Z_k^V.

G is Z_k-connected iff every zero-sum boundary b: V -> Z_k equals the boundary of
some nowhere-zero f: E -> Z_k.  The set of achievable boundaries is the Minkowski
sum over edges e=(u,v) of {a*(chi_u - chi_v) : a in Z_k \\ {0}}.  We compute it as a
boolean array over Z_k^V using np.roll.  Feasible for k=5 and |V| <= 12.
"""
import sys, itertools
import numpy as np
import networkx as nx


def achievable(n, edges, k=5):
    shape = (k,) * n
    R = np.zeros(shape, dtype=bool)
    R[(0,) * n] = True
    for (u, v) in edges:
        new = np.zeros_like(R)
        for a in range(1, k):
            new |= np.roll(np.roll(R, a, axis=u), -a, axis=v)
        R = new
    return R


def is_Zk_connected(n, edges, k=5):
    R = achievable(n, edges, k)
    idx = np.indices((k,) * n).reshape(n, -1).T
    zero_sum = (idx.sum(axis=1) % k == 0)
    reach = R.reshape(-1)
    missing = np.where(zero_sum & ~reach)[0]
    frac = reach[zero_sum].mean()
    return len(missing) == 0, frac, (idx[missing[0]] if len(missing) else None)


def petersen():
    outer = [(i, (i + 1) % 5) for i in range(5)]
    spokes = [(i, i + 5) for i in range(5)]
    inner = [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
    return 10, outer + spokes + inner


if __name__ == "__main__":
    tests = {
        "K4": (4, list(itertools.combinations(range(4), 2))),
        "K3,3": (6, [(i, 3 + j) for i in range(3) for j in range(3)]),
        "prism": (6, [(0,1),(1,2),(2,0),(3,4),(4,5),(5,3),(0,3),(1,4),(2,5)]),
        "cube": (8, [(u, v) for u in range(8) for v in range(u+1, 8) if bin(u ^ v).count('1') == 1]),
        "Petersen": petersen(),
    }
    for name, (n, E) in tests.items():
        for k in (5, 6):
            ok, frac, wit = is_Zk_connected(n, E, k)
            print(f"{name:9s} Z{k}-connected: {ok}  (fraction of zero-sum boundaries achievable: {frac:.4f})" + (f"  missing e.g. {wit}" if wit is not None else ""))
