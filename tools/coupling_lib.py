import itertools
def candidates(G, circuits):
    n = G.number_of_nodes(); cedges = {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    mate = {}
    for e in {frozenset(e) for e in G.edges()} - cedges: u, v = tuple(e); mate[u] = v; mate[v] = u
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
                        bnd = [frozenset((u, v)) for u in S for v in G[u] if v not in S and frozenset((u, v)) in cedges]
                        desc = tuple(sorted([(ci, "W") for ci in W] + [(a[0], len(a[1])) for a in sel]))
                        cands.append((frozenset(S), c1 + 2 * na, bnd, desc))
    return cands
