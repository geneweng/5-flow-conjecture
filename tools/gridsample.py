"""Sample random 0-edge choices c of an instance's 2-factor; for each re-pairing two-position shift c' (H three paths
at both), compute the badness (min-cut margin < 0) of the six grid colourings (rows a0 / E / E1, columns F / F') and
tally the pattern of fully-killed rows.  usage: gridsample.py DUMP NCHOICES SEED"""
import sys, json, random, itertools, collections
import networkx as nx
from oddness import canonical_colouring, H_components

def paths_of(G, col):
    H = nx.Graph([tuple(e) for e, c in col.items() if c in (1, 2)]); H.add_nodes_from(G.nodes()); out = []
    for cc in nx.connected_components(H):
        sub = H.subgraph(cc); ends = [v for v in cc if sub.degree(v) == 1]
        out.append(nx.shortest_path(sub, ends[0], ends[1]) if ends else None)
    return out
def margin(G, black):
    D = nx.DiGraph()
    for u, v in G.edges(): D.add_edge(u, v, capacity=3); D.add_edge(v, u, capacity=3)
    for v in G.nodes():
        if v in black: D.add_edge('s', v, capacity=5)
        else: D.add_edge(v, 't', capacity=5)
    cut, _ = nx.minimum_cut(D, 's', 't'); return cut - 5 * len(black)
fn, N, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); rng = random.Random(seed)
d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
odd = [C for C in circuits if len(C) % 2]
tally = collections.Counter(); nshift = 0
for _ in range(N):
    base = [rng.randrange(len(C)) for C in odd]
    col0 = canonical_colouring(G, M, circuits, base); comps0 = paths_of(G, col0)
    if any(c is None for c in comps0) or len(comps0) != 3: tally["base not three paths"] += 1; continue
    for i, C in enumerate(odd):
        L = len(C); j = base[i]
        for sgn in (1, -1):
            a0, a1, a2 = C[j], C[(j + sgn) % L], C[(j + 2 * sgn) % L]
            P = next(c for c in comps0 if a0 in c); Pp = next(c for c in comps0 if a1 in c)
            if Pp is P: continue
            zc = list(base); zc[i] = (j + 2 * sgn) % L
            comps1 = paths_of(G, canonical_colouring(G, M, circuits, zc))
            if any(c is None for c in comps1) or len(comps1) != 3: continue
            R = next(c for c in comps0 if c is not P and c is not Pp)
            ia1, ia2 = Pp.index(a1), Pp.index(a2); pi1 = Pp[:ia1 + 1] if ia1 < ia2 else Pp[ia1:]
            eps = {}
            for comp in comps0:
                for k, v in enumerate(comp): eps[v] = k % 2 == 0
            if not eps[a0]:
                for v in P: eps[v] = not eps[v]
            if not eps[a1]:
                for v in Pp: eps[v] = not eps[v]
            nshift += 1; rows = {}
            for row, flip in (("a0", ()), ("E", (P,)), ("E1", (pi1,))):
                bad = 0
                for colflip in ((), (R,)):
                    e = dict(eps)
                    for part in flip + colflip:
                        for v in part: e[v] = not e[v]
                    black = {v for v in G if e[v]}
                    if margin(G, black) < 0: bad += 1
                rows[row] = bad
            key = tuple(sorted((r, b) for r, b in rows.items())); tally[key] += 1
            full = [r for r, b in rows.items() if b == 2]
            if set(full) >= {"a0", "E1"}: print("rows a0 and E1 both fully killed:", fn, base, "circuit", i, "shift", 2 * sgn, rows, flush=True)
print(fn.split("/")[-1], f"{nshift} re-pairing two-shifts from {N} choices; tally of (row, #bad of 2):")
for k, v in sorted(tally.items(), key=lambda kv: -kv[1]): print("  ", v, k)
