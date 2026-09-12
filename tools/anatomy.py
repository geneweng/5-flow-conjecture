"""Anatomy of the witness cuts of a dumped instance: composition of S by circuit (arcs / whole circuits),
boundary edges by colour, matching edges inside vs leaving, z's and their colours."""
import sys, json, itertools, collections
import networkx as nx
from oddness import canonical_colouring, H_components, partition
d = json.load(open(sys.argv[1])); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
col = canonical_colouring(G, M, circuits, d["zero"]); comps = H_components(G, col)
odd = [C for C in circuits if len(C) % 2]; zs = {C[z]: ci for ci, (C, z) in enumerate(zip(odd, d["zero"]))}
circ_of = {v: ci for ci, C in enumerate(circuits) for v in C}
print("circuits:", [(ci, len(C), "odd" if len(C) % 2 else "even") for ci, C in enumerate(circuits)])
print("z's (vertex: circuit):", zs, " paths:", [tuple(e) for _, e in comps if e])
seen = set()
for k, Sl in d["S"].items():
    S = set(Sl)
    if len(S) > G.number_of_nodes() // 2: S = set(G.nodes()) - S
    if frozenset(S) in seen: continue
    seen.add(frozenset(S))
    parts = collections.defaultdict(list)
    for v in S: parts[circ_of[v]].append(v)
    desc = []
    for ci, vs in sorted(parts.items()):
        C = circuits[ci]; L = len(C); inside = [v in S for v in C]
        if all(inside): desc.append(f"C{ci}(whole,{L})")
        else:
            arcs = sum(1 for i in range(L) if inside[i] and not inside[i - 1]); desc.append(f"C{ci}({len(vs)} of {L}, {arcs} arc)")
    bnd = collections.Counter(col[frozenset((u, v))] for u, v in G.edges() if (u in S) != (v in S))
    minside = sum(1 for e in M if all(v in S for v in e)); mleave = sum(1 for e in M if sum(v in S for v in e) == 1)
    print(f"cut for colouring {k}: |S|={len(S)} boundary colours {dict(bnd)}; matching inside {minside}, leaving {mleave}; parts: {', '.join(desc)}; z in S: {[z for z in zs if z in S]}")
