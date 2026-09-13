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

def rows_for(G, M, circuits, odd, base, i, sgn):
    col0 = canonical_colouring(G, M, circuits, base); comps0 = paths_of(G, col0)
    if any(c is None for c in comps0) or len(comps0) != 3: return None
    C = odd[i]; L = len(C); j = base[i]
    a0, a1, a2 = C[j], C[(j + sgn) % L], C[(j + 2 * sgn) % L]
    P = next(c for c in comps0 if a0 in c); Pp = next(c for c in comps0 if a1 in c)
    if Pp is P: return None
    zc = list(base); zc[i] = (j + 2 * sgn) % L
    comps1 = paths_of(G, canonical_colouring(G, M, circuits, zc))
    if any(c is None for c in comps1) or len(comps1) != 3: return None
    R = next(c for c in comps0 if c is not P and c is not Pp)
    ia1, ia2 = Pp.index(a1), Pp.index(a2); pi1 = Pp[:ia1 + 1] if ia1 < ia2 else Pp[ia1:]
    eps = {}
    for comp in comps0:
        for k, v in enumerate(comp): eps[v] = k % 2 == 0
    if not eps[a0]:
        for v in P: eps[v] = not eps[v]
    if not eps[a1]:
        for v in Pp: eps[v] = not eps[v]
    rows = {}
    for row, flip in (("a0", ()), ("E", (P,)), ("E1", (pi1,))):
        bad = 0
        for colflip in ((), (R,)):
            e = dict(eps)
            for part in flip + colflip:
                for v in part: e[v] = not e[v]
            if margin(G, {v for v in G if e[v]}) < 0: bad += 1
        rows[row] = bad
    return rows, zc
