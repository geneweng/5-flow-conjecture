"""Distribution of the number of killed colourings (all 2^{#components} colourings of H) over 0-edge
choices of a fixed graph and 2-factor.   usage: zero_killdist.py DUMP [N|all]"""
import sys, json, itertools, random, collections, time
import networkx as nx
from oddness import canonical_colouring, H_components, partition, balanced
d = json.load(open(sys.argv[1])); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
odd = [C for C in circuits if len(C) % 2]
space = list(itertools.product(*[range(len(C)) for C in odd])); mode = sys.argv[2] if len(sys.argv) > 2 else "5000"
if mode != "all": space = random.Random(2).sample(space, int(mode))
hist = collections.Counter(); comps_hist = collections.Counter(); t0 = time.time()
for zc in space:
    col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col); comps_hist[len(comps)] += 1
    k = sum(0 if balanced(G, partition(comps, f)) else 1 for f in itertools.product((0, 1), repeat=len(comps)))
    hist[(len(comps), k)] += 1
print(f"{len(space)} 0-edge choices [{time.time()-t0:.0f}s]; #H-components histogram {dict(sorted(comps_hist.items()))}")
print("(components, killed colourings) histogram:", dict(sorted(hist.items())))
