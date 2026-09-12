"""Check on dumped instances: a killed colouring whose witness cut is crossed by no even circuit
(s = 0) must have its path-complement killed too (by V \\ S); with s > 0 the complement may survive."""
import sys, json, itertools
import networkx as nx
from oddness import canonical_colouring, H_components, partition
from badcuts import violating_set
for fn in sys.argv[1:]:
    d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
    M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    col = canonical_colouring(G, M, circuits, d["zero"]); comps = H_components(G, col)
    path_idx = [i for i, (_, ends) in enumerate(comps) if ends]; t = len(path_idx); ref = d["ref"]
    even_vertices = {v for i, (c, ends) in enumerate(comps) if not ends for v in c}
    killed = {}
    for bits in itertools.product((0, 1), repeat=t):
        f = list(ref)
        for i, b in zip(path_idx, bits): f[i] = b
        killed[bits] = violating_set(G, partition(comps, f))
    rows = []
    for bits, S in killed.items():
        if S is None: continue
        comp = tuple(1 - b for b in bits)
        s = sum(1 for u, v in G.edges() if (u in S) != (v in S) and col[frozenset((u, v))] == 1 and u in even_vertices)
        rows.append((bits, s, killed[comp] is not None))
    bad = [r for r in rows if r[1] == 0 and not r[2]]
    print(fn, "killed:", len(rows), "of", 2 ** t, "| (bits, s, complement killed):", rows, "| violations of the s=0 rule:", len(bad))
