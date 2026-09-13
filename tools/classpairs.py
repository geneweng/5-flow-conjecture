"""Pairwise class-kill correlations at a fixed choice from an exhaustive live-set dump (lcirc_test.py): for every
0-edge choice with H three paths, the set of classes killed by live 11/7-cuts (pattern-matched as in live_analysis.py);
tabulate P(i and j) against P(i)P(j) for adjacent classes (differ by one path flip) and for the full cover."""
import sys, json, itertools, collections
import networkx as nx
from oddness import canonical_colouring, H_components
d = json.load(open(sys.argv[1])); L = json.load(open(sys.argv[2]))
G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]; odd = [C for C in circuits if len(C) % 2]
M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
by_choice = collections.defaultdict(list)
for rec in L:
    if rec["d"] in (7, 11):
        for z in rec["live"]: by_choice[tuple(z)].append(set(rec["S"]))
N = 1
for C in odd: N *= len(C)
single = collections.Counter(); pair = collections.Counter(); nthree = 0; full = 0
for zc, sets in by_choice.items():
    col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
    paths = [tuple(ends) for _, ends in comps if ends]
    if len(paths) != 3 or len(comps) != 3: continue
    Z = {z: (i, e) for i, p in enumerate(paths) for e, z in enumerate(p)}
    killed = set()
    for S in sets:
        zin = [Z[z] for z in Z if z in S]
        if len(zin) == 3 and len({i for i, e in zin}) == 3:
            pat = tuple(e for i, e in sorted(zin))
            if pat[0] == 1: pat = tuple(1 - e for e in pat)      # normalize: path 0 end 0 black
            killed.add(pat)
    for k in killed: single[k] += 1
    for k1, k2 in itertools.combinations(sorted(killed), 2): pair[(k1, k2)] += 1
    if len(killed) == 4: full += 1
classes = [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1)]
print("choices (all with three paths assumed):", N, "; per-class kill counts:", {k: single[k] for k in classes}, "; full covers:", full)
for k1, k2 in itertools.combinations(classes, 2):
    ham = sum(a != b for a, b in zip(k1, k2))
    exp = single[k1] * single[k2] / N
    print(f"  classes {k1} {k2} (differ in {ham} path{'s' if ham > 1 else ''}): joint {pair[(k1, k2)]}, expected under independence {exp:.1f}, ratio {pair[(k1, k2)] / exp if exp else float('nan'):.2f}")
