"""Killed-colouring tail on RANDOM instances: generate cyclically 6-connected graphs of a profile, sample
0-edge choices, and tabulate the number of killed colourings (all components) per choice.
usage: tail_random.py PROFILE NGRAPHS NSAMPLES SEED"""
import sys, random, itertools, collections
from gen_cyc6 import generate
from oddness import canonical_colouring, H_components, partition, balanced
prof = {"six7": [7] * 6, "six7p8": [7] * 6 + [8]}[sys.argv[1]]; ng = int(sys.argv[2]); ns = int(sys.argv[3]); rng = random.Random(int(sys.argv[4]))
graphs = generate(prof, rng, ng)
for gi, (G, circuits) in enumerate(graphs):
    M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    odd = [C for C in circuits if len(C) % 2]; hist = collections.Counter()
    seen = set()
    while len(seen) < ns: seen.add(tuple(rng.randrange(len(C)) for C in odd))
    for zc in seen:
        col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
        k = sum(0 if balanced(G, partition(comps, f)) else 1 for f in itertools.product((0, 1), repeat=len(comps)))
        hist[(len(comps), k)] += 1
    print(f"graph {gi} n={G.number_of_nodes()}: (components, killed) histogram {dict(sorted(hist.items()))}", flush=True)
