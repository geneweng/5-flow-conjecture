"""For each killing cut S of a counterexample: over all 0-edge choices c, is S a bad cut for SOME
2-colouring of H(c)?  (max over colourings of |b_S - a_S| = sum over components K of H(c) of
|#K∩S black - #K∩S white| under K's canonical colouring.)  Reports the density of each killer and of
their union over the 0-edge space, and how many choices have all four cuts simultaneously 'live'."""
import sys, json, itertools, collections, time
import networkx as nx
from oddness import canonical_colouring, H_components
d = json.load(open(sys.argv[1])); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
odd = [C for C in circuits if len(C) % 2]
cuts = {}
for k, Sl in d["S"].items():
    S = frozenset(Sl); cuts.setdefault(S, []).append(k)
cuts = list(cuts); dS = [nx.cut_size(G, S) for S in cuts]; print(len(cuts), "distinct witness sets, sizes", [len(S) for S in cuts], "d =", dS)
live = collections.Counter(); per = [0] * len(cuts); t0 = time.time(); N = 0
for zc in itertools.product(*[range(len(C)) for C in odd]):
    N += 1; col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
    mask = 0
    for i, S in enumerate(cuts):
        best = 0
        for colouring, ends in comps:
            imb = sum(1 if c else -1 for v, c in colouring.items() if v in S); best += abs(imb)
        if 5 * best > 3 * dS[i]: per[i] += 1; mask |= 1 << i
    live[bin(mask).count("1")] += 1
print(f"{N} choices [{time.time()-t0:.0f}s]; per-cut density of 'bad for some colouring': {[round(p / N, 4) for p in per]}")
print("number of the witness cuts live simultaneously:", dict(sorted(live.items())))
