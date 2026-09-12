"""Analyse a dumped worst instance of killsets_t.py: for every path-colouring class report the
witness bad cut (type, z's inside with colours, blocks touched) and for the surviving classes
the balance margin min_S [3 d(S) - 5(b_S - a_S)]."""
import sys, json, itertools, collections
import networkx as nx
from oddness import canonical_colouring, H_components, partition
from badcuts import violating_set, classify
d = json.load(open(sys.argv[1])); per_block = int(sys.argv[2]) if len(sys.argv) > 2 else 3
G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
block = {v: ci // per_block for ci, C in enumerate(circuits) for v in C}
col = canonical_colouring(G, M, circuits, d["zero"]); comps = H_components(G, col)
path_idx = [i for i, (_, ends) in enumerate(comps) if ends]; t = len(path_idx); ref = d["ref"]
paths = [tuple(comps[i][1]) for i in path_idx]
print(f"n={G.number_of_nodes()}, t={t}, paths (z ends -> blocks):", [(p, (block[p[0]], block[p[1]])) for p in paths])
print("path lengths (vertices):", [len(comps[i][0]) for i in path_idx])
def margin(black):
    D = nx.DiGraph()
    for u, v in G.edges(): D.add_edge(u, v, capacity=3); D.add_edge(v, u, capacity=3)
    for v in G.nodes():
        if v in black: D.add_edge('s', v, capacity=5)
        else: D.add_edge(v, 't', capacity=5)
    cut, _ = nx.minimum_cut(D, 's', 't'); return cut - 5 * len(black)
summary = "-s" in sys.argv; cuts = {}; nsurv = 0
for bits in itertools.product((0, 1), repeat=t):
    f = list(ref)
    for i, b in zip(path_idx, bits): f[i] = b
    black = partition(comps, f); S = violating_set(G, black)
    zcol = ["".join("B" if z in black else "w" for z in p) for p in paths]
    if S is None:
        nsurv += 1; print(f"class {bits}: SURVIVES  path-end colours {zcol}"); continue
    ty = classify(G, S, col, comps, black)
    if summary:
        key = frozenset(S); cuts.setdefault(key, [ty[:5], ty[5], [(z, block[z]) for p in paths for z in p if z in S], 0]); cuts[key][3] += 1; continue
    zin = [z for p in paths for z in p if z in S]
    blocks_in = collections.Counter(block[v] for v in S)
    print(f"class {bits}: killed by {ty[:5]} sep={ty[5]}  |S|={len(S)}  z in S: {[(z, 'B' if z in black else 'w', block[z]) for z in zin]}  blocks of S: {dict(sorted(blocks_in.items()))}  path-end colours {zcol}")

if summary:
    print(f"{nsurv} survivors of {2**t} path colourings; {len(cuts)} distinct witness cuts:")
    for key, (ty, sep, zs, cnt) in sorted(cuts.items(), key=lambda kv: -kv[1][3]):
        print(f"   {ty} sep={sep} kills {cnt:3d}  |S|={len(key):3d}  blocks={dict(sorted(collections.Counter(block[v] for v in key).items()))}  z in S: {zs}")
