"""For the counterexample's graph and 2-factor: over 0-edge choices, does some 2-colouring of H (all
components) give a balanced valuation?  usage: zero_sweep.py DUMP [all|N]"""
import sys, json, itertools, random, time
import networkx as nx
from oddness import canonical_colouring, H_components, partition, balanced
d = json.load(open(sys.argv[1])); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
odd = [C for C in circuits if len(C) % 2]
space = list(itertools.product(*[range(len(C)) for C in odd]))
mode = sys.argv[2] if len(sys.argv) > 2 else "1000"
if mode != "all": space = random.Random(1).sample(space, int(mode))
ok = fail = 0; failures = []; t0 = time.time()
for k, zc in enumerate(space):
    col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
    found = any(balanced(G, partition(comps, f)) for f in itertools.product((0, 1), repeat=len(comps)))
    if found: ok += 1
    else: fail += 1; failures.append(zc)
    if (k + 1) % 2000 == 0: print(f"  {k+1}/{len(space)}: ok {ok}, fail {fail} [{time.time()-t0:.0f}s]", flush=True)
print(f"0-edge choices tested {len(space)}: some balanced colouring exists for {ok}; none for {fail}")
print("failing 0-edge choices:", failures[:50])
