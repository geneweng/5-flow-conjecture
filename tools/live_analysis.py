"""Analyse the live-set dump of lcirc_test.py (<instance>_live.json) together with the instance:
(a) is each candidate's live set a product of its per-circuit marginals?  (b) for each 0-edge choice,
which classes are killed by live candidates whose pattern matches (an 11/7-cut with z's T kills the
class of T; a pair cut kills the two classes containing its two ends), and how many choices have all
four classes covered ('potential failures') vs the actual failing choice.   usage: live_analysis.py DUMP LIVEJSON"""
import sys, json, itertools, collections
import networkx as nx
d = json.load(open(sys.argv[1])); L = json.load(open(sys.argv[2]))
G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]; odd = [C for C in circuits if len(C) % 2]
ncirc = len(odd)
exact = 0; ratio = []
for rec in L:
    live = [tuple(z) for z in rec["live"]]; marg = [sorted({z[i] for z in live}) for i in range(ncirc)]
    prod = 1
    for m in marg: prod *= len(m)
    if prod == len(live): exact += 1
    ratio.append(len(live) / prod)
print(f"{len(L)} live candidates; live set equals the product of its marginals for {exact}; mean live/product = {sum(ratio)/len(ratio):.3f}, min {min(ratio):.3f}")
# classes: the paths are determined by the choice c (they depend on c!), so 'class of a triple' must be computed per c.
# For each choice c, compute the paths of H(c), then for each live candidate S at c, the set of z's inside S and the class(es) it can kill.
from oddness import canonical_colouring, H_components
M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
by_choice = collections.defaultdict(list)
for rec in L:
    for z in rec["live"]: by_choice[tuple(z)].append(rec)
cover = collections.Counter(); potential = []
for zc, recs in by_choice.items():
    col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
    paths = [tuple(ends) for _, ends in comps if ends]
    if len(paths) != 3: cover[("t", len(paths))] += 1; continue
    Z = {z: (i, e) for i, p in enumerate(paths) for e, z in enumerate(p)}
    killed = set()
    for rec in recs:
        S = set(rec["S"]); zin = [Z[z] for z in Z if z in S]
        if rec["d"] in (7, 11) and len(zin) == 3 and len({i for i, e in zin}) == 3:
            pat = tuple(sorted(zin)); comp = tuple(sorted((i, 1 - e) for i, e in zin))
            killed.add(min(pat, comp))
        elif rec["d"] == 6 and len(zin) in (2, 4):
            sep = [i for i in range(3) if sum(1 for z in Z if z in S and Z[z][0] == i) == 1]
            if len(sep) == 2:
                ends = {i: e for i, e in zin if i in sep}
                for e3 in (0, 1):   # both classes containing this pair pattern
                    i3 = [i for i in range(3) if i not in sep][0]
                    pat = tuple(sorted(list(ends.items()) + [(i3, e3)])); comp = tuple(sorted((i, 1 - e) for i, e in pat))
                    killed.add(min(pat, comp))
    cover[len(killed)] += 1
    if len(killed) == 4: potential.append(zc)
print("choices with >=1 live candidate:", len(by_choice), "; number of classes covered by pattern-matching live candidates:", dict(sorted(cover.items(), key=str)))
print("choices with all four classes covered ('potential failures'):", len(potential), potential[:10], "; actual failing choice:", d["zero"])
# (c) pins: for live candidates, circuits whose marginal is a single position (typically the arc circuits);
# per circuit, the set of pinned positions used by any live candidate; a position used by no pin, if chosen,
# kills every candidate that pins that circuit.
pins = collections.defaultdict(set); pin_count = collections.Counter()
for rec in L:
    live = [tuple(z) for z in rec["live"]]
    for i in range(ncirc):
        vals = {z[i] for z in live}
        if len(vals) == 1: pins[i].add(next(iter(vals))); pin_count[i] += 1
print("pinned positions per odd circuit (positions used as a pin by some live candidate):", {i: sorted(pins[i]) for i in range(ncirc)}, "; candidates pinning each circuit:", dict(pin_count))
free_pos = {i: [p for p in range(len(odd[i])) if p not in pins[i]] for i in range(ncirc)}
print("positions used by no pin, per circuit:", free_pos)
# choose, per circuit, a pin-free position if one exists; count live candidates at such choices
choices = list(itertools.product(*[free_pos[i] if free_pos[i] else range(len(odd[i])) for i in range(ncirc)]))
alive = collections.Counter(len(by_choice.get(zc, [])) for zc in choices)
print(f"{len(choices)} pin-avoiding choices; live candidates at them:", dict(sorted(alive.items())))
