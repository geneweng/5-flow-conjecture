"""Estimate the TRUE averaging total  T = sum_S live(S) k(S)  where live(S) = number of 0-edge choices c for
which S is bad for some colouring of H(c), over all structured candidate sets S (whole circuits + <= 2 arcs
with the right matching count), by sampling M choices c.   usage: avg_true.py DUMP M"""
import sys, json, itertools, random, time
import networkx as nx
from oddness import canonical_colouring, H_components
d = json.load(open(sys.argv[1])); M_ = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]; n = G.number_of_nodes()
cedges = {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
Mset = {frozenset(e) for e in G.edges()} - cedges; mate = {}
for e in Mset: u, v = tuple(e); mate[u] = v; mate[v] = u
odd = [C for C in circuits if len(C) % 2]; N = 1
for C in odd: N *= len(C)
arcs = [(ci, frozenset(C[(p + k) % len(C)] for k in range(l))) for ci, C in enumerate(circuits) for p in range(len(C)) for l in range(1, len(C))]
cands = []   # (frozenset S, d(S), k)
ncirc = len(circuits)
for r in range(ncirc + 1):
    for W in itertools.combinations(range(ncirc), r):
        base = set(v for ci in W for v in circuits[ci]); free = [ci for ci in range(ncirc) if ci not in W]
        for na in (1, 2):
            for sel in itertools.combinations([a for a in arcs if a[0] in free], na):
                if len({a[0] for a in sel}) < na: continue
                S = base | set().union(*[a[1] for a in sel])
                if not S or 2 * len(S) > n: continue
                c1 = sum(1 for v in S if mate[v] not in S)
                if (na, c1) == (2, 7): cands.append((frozenset(S), 11, 1))
                elif (na, c1) == (1, 5): cands.append((frozenset(S), 7, 1))
                elif (na, c1) == (1, 4): cands.append((frozenset(S), 6, 2))
print(f"{len(cands)} candidate sets; sampling {M_} of {N} 0-edge choices", flush=True)
rng = random.Random(7); t0 = time.time(); total = 0.0; per_c_hist = {}
idx = {v: i for i, v in enumerate(G.nodes())}
for m in range(M_):
    zc = tuple(rng.randrange(len(C)) for C in odd)
    col = canonical_colouring(G, Mset, circuits, zc); comps = H_components(G, col)
    # per component, signed colour vector; imbalance of S within a component = |sum over v in S of sign|
    comp_sign = [{v: (1 if c else -1) for v, c in colouring.items()} for colouring, _ in comps]
    live_here = 0
    for S, dS, k in cands:
        best = 0
        for cs in comp_sign:
            s = 0
            for v in S:
                s += cs.get(v, 0)
            best += abs(s)
        if 5 * best > 3 * dS: total += k; live_here += 1
    per_c_hist[live_here] = per_c_hist.get(live_here, 0) + 1
T = total / M_ * N
print(f"estimated T = sum_S live(S) k(S) = {T:.0f}  vs 4N = {4*N}   (T/4N = {T/(4*N):.3f})  [{time.time()-t0:.0f}s]")
print("live candidate sets per sampled choice:", dict(sorted(per_c_hist.items())))
