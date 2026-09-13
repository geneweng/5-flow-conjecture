"""Empirical test of (L-circ): for every candidate tight set S (whole circuits + <= 2 arcs, right matching count)
of a 2-factor, compute the exact set of live 0-edge choices (S bad for some colouring of H(c)) over all
choices, and for each circuit the number of 0-edge positions occurring among the live choices.
Reports, over candidates with live > 0: distribution of live counts, and for circuits disjoint from S the
distribution of the marginal size (7 = unconstrained)."""
import sys, json, itertools, collections, time
import networkx as nx
from oddness import canonical_colouring, H_components
d = json.load(open(sys.argv[1])); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]; n = G.number_of_nodes()
cedges = {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
Mset = {frozenset(e) for e in G.edges()} - cedges; mate = {}
for e in Mset: u, v = tuple(e); mate[u] = v; mate[v] = u
odd = [C for C in circuits if len(C) % 2]
arcs = [(ci, frozenset(C[(p + k) % len(C)] for k in range(l))) for ci, C in enumerate(circuits) for p in range(len(C)) for l in range(1, len(C))]
cands = []; ncirc = len(circuits)
for r in range(ncirc + 1):
    for W in itertools.combinations(range(ncirc), r):
        base = set(v for ci in W for v in circuits[ci]); free = [ci for ci in range(ncirc) if ci not in W]
        for na in (1, 2):
            for sel in itertools.combinations([a for a in arcs if a[0] in free], na):
                if len({a[0] for a in sel}) < na: continue
                S = base | set().union(*[a[1] for a in sel])
                if not S or 2 * len(S) > n: continue
                c1 = sum(1 for v in S if mate[v] not in S)
                if (na, c1) in ((2, 7), (1, 5), (1, 4)):
                    # arc-end circuit edges (must be colour 2 for 7/11-cuts; for pair cuts allow colour 0/3 on one) and admissible z counts
                    bnd = [frozenset((u, v)) for u in S for v in G[u] if v not in S and frozenset((u, v)) in cedges]
                    cands.append((frozenset(S), c1 + 2 * na, set(W) | {a[0] for a in sel}, bnd, (2, 3, 4) if c1 == 4 else (3,), c1 == 4))
print(len(cands), "candidates", flush=True)
live = {i: [] for i in range(len(cands))}; t0 = time.time()
for m, zc in enumerate(itertools.product(*[range(len(C)) for C in odd])):
    col = canonical_colouring(G, Mset, circuits, zc); comps = H_components(G, col)
    comp_sign = None
    Z = {C[j] for C, j in zip(odd, zc)}
    for i, (S, dS, touched, bnd, zok, pair) in enumerate(cands):
        cols = [col[e] for e in bnd]
        if pair:
            if not (all(c == 2 for c in cols) or (sum(1 for c in cols if c == 2) == len(cols) - 1)): continue
        elif any(c != 2 for c in cols): continue
        if len(S & Z) not in zok: continue
        if comp_sign is None: comp_sign = [{v: (1 if c else -1) for v, c in colouring.items()} for colouring, _ in comps]
        best = 0
        for cs in comp_sign:
            s = 0
            for v in S: s += cs.get(v, 0)
            best += abs(s)
        if 5 * best > 3 * dS: live[i].append(zc)
    if m % 5000 == 0: print(f"  {m} choices [{time.time()-t0:.0f}s]", flush=True)
lc = collections.Counter(); marg_free = collections.Counter(); marg_touched = collections.Counter(); worst = []
for i, (S, dS, touched, bnd, zok, pair) in enumerate(cands):
    L = live[i]
    if not L: continue
    lc[len(L)] += 1
    oddidx = {ci: k for k, ci in enumerate(c for c in range(ncirc) if len(circuits[c]) % 2)}
    for ci in range(ncirc):
        if ci not in oddidx: continue                     # even circuits carry no 0-edge
        k = len({zc[oddidx[ci]] for zc in L})
        (marg_touched if ci in touched else marg_free)[k] += 1
        if ci not in touched and k == len(circuits[ci]): worst.append((len(S), dS, len(L), ci))
print(f"candidates with live > 0: {sum(lc.values())}; live-count distribution: {dict(sorted(lc.items()))}")
print("marginal sizes over circuits DISJOINT from S:", dict(sorted(marg_free.items())))
print("marginal sizes over circuits meeting S (whole or arc):", dict(sorted(marg_touched.items())))
print("examples of disjoint circuits with all 7 positions live (|S|, d, live, circuit):", worst[:10], "... total", len(worst))
out = sys.argv[1].rsplit("/", 1)[-1].replace(".json", "") + "_live.json"
json.dump([{"S": sorted(S), "d": dS, "k": 2 if dS == 6 else 1, "live": live[i]} for i, (S, dS, touched, bnd, zok, pair) in enumerate(cands) if live[i]], open(out, "w"))
T = sum((2 if dS == 6 else 1) * len(live[i]) for i, (S, dS, touched, bnd, zok, pair) in enumerate(cands))
N = 1
for C in odd: N *= len(C)
print(f"EXACT averaging total T = sum_S live(S) k(S) = {T}  vs 4N = {4*N}  (T/4N = {T/(4*N):.4f}); live sets dumped to {out}")
