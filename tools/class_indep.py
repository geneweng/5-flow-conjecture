"""Per-class kill probabilities and pairwise joints over the 0-edge space, from a live-set dump:
for each choice with live candidates, the set of classes covered by pattern-matching live candidates
(same rule as live_analysis.py); classes are labelled by the pattern (which end of P2, which end of P3
is inside together with an end of P1) -- but paths depend on the choice, so classes are compared via the
z-vertices: a class = the set of three z-vertices inside (canonical: the triple containing the smallest
z of P1's pair... here simply the frozenset of z's on the side containing the globally smallest z)."""
import sys, json, itertools, collections
import networkx as nx
from oddness import canonical_colouring, H_components
d = json.load(open(sys.argv[1])); L = json.load(open(sys.argv[2]))
G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]; odd = [C for C in circuits if len(C) % 2]
M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
N = 1
for C in odd: N *= len(C)
by_choice = collections.defaultdict(list)
for rec in L:
    for z in rec["live"]: by_choice[tuple(z)].append(rec)
per_choice_classes = {}
for zc, recs in by_choice.items():
    col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
    paths = [tuple(ends) for _, ends in comps if ends]
    if len(paths) != 3: continue
    Z = {z: (i, e) for i, p in enumerate(paths) for e, z in enumerate(p)}
    killed = set()
    for rec in recs:
        S = set(rec["S"]); zin = [Z[z] for z in Z if z in S]
        if rec["d"] in (7, 11) and len(zin) == 3 and len({i for i, e in zin}) == 3:
            pat = tuple(sorted(zin)); comp = tuple(sorted((i, 1 - e) for i, e in zin)); killed.add(min(pat, comp))
        elif rec["d"] == 6 and len(zin) in (2, 4):
            sep = [i for i in range(3) if sum(1 for z in Z if z in S and Z[z][0] == i) == 1]
            if len(sep) == 2:
                ends = {i: e for i, e in zin if i in sep}; i3 = [i for i in range(3) if i not in sep][0]
                for e3 in (0, 1):
                    pat = tuple(sorted(list(ends.items()) + [(i3, e3)])); comp = tuple(sorted((i, 1 - e) for i, e in pat)); killed.add(min(pat, comp))
    per_choice_classes[zc] = killed
# class labels (up to complement) as tuples of (path, end); the paths are labelled by the choice-dependent
# component order, so the label is only meaningful per choice.  Use the count structure instead:
cnt = collections.Counter(len(k) for k in per_choice_classes.values())
n_any = sum(cnt.values()); p_any = n_any / N
# per-class marginals: total class-kills / (4N)
p_class = sum(len(k) for k in per_choice_classes.values()) / (4 * N)
print(f"choices with >=1 killed class: {n_any} (p_any = {p_any:.4f}); average per-class kill probability p = {p_class:.4f}")
print("independent-classes model with this p: P[>=1] =", round(1 - (1 - p_class) ** 4, 4), " P[=2] =", round(6 * p_class**2 * (1 - p_class)**2, 5), " P[=3] =", round(4 * p_class**3 * (1 - p_class), 6), " P[=4] =", f"{p_class**4:.2e}")
print("observed: P[=1] =", round(cnt[1] / N, 4), " P[=2] =", round(cnt[2] / N, 5), " P[=3] =", round(cnt[3] / N, 6), " P[=4] =", f"{cnt[4] / N:.2e}")
