"""Test (S1): min degree >= 9 and every proper set A with |A|,|A^c| >= 2 having d(A) >= 12
=> strongly Z5-connected?   Exhaustive over small multigraphs."""
import sys, itertools, time
import numpy as np
from mod5 import strongly_connected_fast, strongly_connected

n = int(sys.argv[1]); dmin = 9; dmax = int(sys.argv[2]) if len(sys.argv) > 2 else 9
pairs = list(itertools.combinations(range(n), 2))
seen = set(); stats = {}

def cuts_ok(m):
    for r in range(2, n - 1):
        for A in itertools.combinations(range(n), r):
            As = set(A)
            d = sum(k for (u, v), k in m.items() if (u in As) != (v in As))
            if d < 12:
                return False
    return True

def canon(m):
    rows = []
    for v in range(n):
        rows.append(tuple(sorted(m.get((min(u, v), max(u, v)), 0) for u in range(n) if u != v)))
    return tuple(sorted(rows))

deg = [0] * n; mult = {}
def rec(i):
    if i == len(pairs):
        if all(dmin <= d <= dmax for d in deg):
            key = canon(mult)
            if key in seen: return
            seen.add(key)
            if not cuts_ok(mult): return
            E = [p for p, k in mult.items() for _ in range(k)]
            ok, wit = strongly_connected_fast(n, E, 5)
            stats[ok] = stats.get(ok, 0) + 1
            if not ok:
                print("NOT strongly Z5-connected:", {p: k for p, k in mult.items() if k}, "degrees", deg, "bad beta", wit, flush=True)
        return
    u, v = pairs[i]
    for k in range(0, dmax + 1):
        if deg[u] + k > dmax or deg[v] + k > dmax: break
        deg[u] += k; deg[v] += k
        if k: mult[(u, v)] = k
        if not (v == n - 1 and deg[u] < dmin):
            rec(i + 1)
        deg[u] -= k; deg[v] -= k
        mult.pop((u, v), None)

t = time.time(); rec(0)
print(f"n={n}, degrees in [{dmin},{dmax}]: graphs passing the cut condition: {sum(stats.values())}, strongly Z5-connected: {stats.get(True,0)}, NOT: {stats.get(False,0)}  ({time.time()-t:.0f}s)")
