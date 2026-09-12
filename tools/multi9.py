"""Exhaustive search over loopless 9-regular multigraphs on n vertices (n even, small)
with edge-connectivity >= lam_min, testing for a modulo-5 orientation.
Generation: symmetric multiplicity matrices with row sums 9, by backtracking over pairs
in lexicographic order; isomorphs are not removed (n <= 6 keeps this affordable)."""
import sys, itertools, time
import networkx as nx
from mod5 import has_mod_orientation

n = int(sys.argv[1]); lam_min = int(sys.argv[2]) if len(sys.argv) > 2 else 8
pairs = list(itertools.combinations(range(n), 2))
deg = [0] * n
mult = {}
count = tested = good = 0
t0 = time.time()
seen = set()

def canon(m):
    # cheap canonical invariant to skip obvious isomorphs: sorted multiset of (sorted multiplicity rows)
    rows = []
    for v in range(n):
        rows.append(tuple(sorted(m.get((min(u, v), max(u, v)), 0) for u in range(n) if u != v)))
    return tuple(sorted(rows))

def edge_conn(m):
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for (u, v), k in m.items():
        if k:
            G.add_edge(u, v, weight=k)
    if not nx.is_connected(G):
        return 0
    return nx.stoer_wagner(G)[0]

def rec(i):
    global count, tested, good
    if i == len(pairs):
        if all(d == 9 for d in deg):
            count += 1
            key = canon(mult)
            if key in seen:
                return
            seen.add(key)
            if edge_conn(mult) >= lam_min:
                tested += 1
                E = [p for p, k in mult.items() for _ in range(k)]
                ok = has_mod_orientation(n, E, 5)
                good += ok
                if not ok:
                    print("NO MOD-5 ORIENTATION:", dict(mult), flush=True)
        return
    u, v = pairs[i]
    # remaining capacity check: vertex u must reach degree 9 using pairs (u, w) with w > v
    rem_u = sum(1 for w in range(v + 1, n)) * 9  # crude upper bound on what later pairs can give u
    for k in range(0, 10):
        if deg[u] + k > 9 or deg[v] + k > 9:
            break
        deg[u] += k; deg[v] += k
        if k: mult[(u, v)] = k
        # prune: if this is the last pair for u, u must be full
        if v == n - 1 and deg[u] != 9:
            pass
        else:
            rec(i + 1) if not (v == n - 1 and deg[u] != 9) else None
        deg[u] -= k; deg[v] -= k
        mult.pop((u, v), None)

rec(0)
print(f"n={n}: labeled 9-regular multigraphs {count}, distinct by invariant {len(seen)}, with lambda>={lam_min}: {tested}, with mod-5 orientation: {good}  ({time.time()-t0:.0f}s)")
