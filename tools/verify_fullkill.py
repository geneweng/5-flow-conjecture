"""Independent verification of a claimed template failure (all 2^t path colourings unbalanced)."""
import sys, json, itertools, collections
import networkx as nx
from oddness import canonical_colouring, H_components, partition
from badcuts import violating_set, classify
from gen_cyc6 import cyclically_6_connected, has_small_cyclic_cut
d = json.load(open(sys.argv[1])); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
n = G.number_of_nodes(); print("n =", n, "edges =", G.number_of_edges())
assert all(deg == 3 for _, deg in G.degree()), "not cubic"; print("cubic: ok")
assert nx.is_connected(G); print("connected: ok")
g = nx.girth(G); print("girth =", g); assert g >= 6
print("edge connectivity =", nx.edge_connectivity(G))
X = has_small_cyclic_cut(G); print("small cyclic cut:", X); assert X is None; print("cyclically 6-edge-connected: ok")
# 2-factor
cov = collections.Counter(v for C in circuits for v in C); assert set(cov) == set(G.nodes()) and all(c == 1 for c in cov.values())
for C in circuits:
    for i in range(len(C)): assert G.has_edge(C[i], C[(i + 1) % len(C)])
odd = [C for C in circuits if len(C) % 2]; print("2-factor: ok; circuit lengths", [len(C) for C in circuits], "; odd circuits:", len(odd))
assert len(odd) == 6
M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
assert len(M) == n // 2 and all(len(set(itertools.chain.from_iterable(M))) == n for _ in [0]); print("complementary perfect matching: ok")
col = canonical_colouring(G, M, circuits, d["zero"]); comps = H_components(G, col)
pidx = [i for i, (_, e) in enumerate(comps) if e]; t = len(pidx); print("H components:", len(comps), "; paths:", t, "; even H-circuits:", len(comps) - t)
ref = d["ref"]; f0 = [ref.get(str(min(c)), 0) for c, _ in comps]
allkilled = True; types = []
for bits in itertools.product((0, 1), repeat=t):
    f = list(f0)
    for i, b in zip(pidx, bits): f[i] = b
    black = partition(comps, f); assert 2 * len(black) == n, "unbalanced partition sizes"
    # every H-edge is bichromatic
    assert all((u in black) != (v in black) for e, c in col.items() if c in (1, 2) for u, v in [tuple(e)])
    S = violating_set(G, black)
    if S is None: allkilled = False; print(bits, "SURVIVES"); continue
    b = sum(1 for v in S if v in black); a = len(S) - b; dS = nx.cut_size(G, S)
    assert 5 * abs(b - a) > 3 * dS, "witness not bad"
    ty = classify(G, S, col, comps, black); types.append(ty)
    print(bits, f"killed: |S|={len(S)}, d(S)={dS}, |b-a|={abs(b-a)}, 5|b-a|-3d={5*abs(b-a)-3*dS}, type {ty[:5]} sep={ty[5]}")
print("ALL PATH COLOURINGS KILLED:", allkilled)
