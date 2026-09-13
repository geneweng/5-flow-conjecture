"""(L-repair) test: for each failing 0-edge choice of an instance (given, or found among the instance's
'zero' and optional 'fails' list), classify every one-step move of a single 0-edge by its effect on H
(circuit-level pairing of the six ends; number of components) and whether the method succeeds at the
new choice.   usage: repair_test.py DUMP [DUMP ...]"""
import sys, json, itertools, collections
import networkx as nx
from oddness import canonical_colouring, H_components, partition, balanced
tot = collections.Counter()
for fn in sys.argv[1:]:
    d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
    M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    odd = [C for C in circuits if len(C) % 2]; circ_of = {v: k for k, C in enumerate(odd) for v in C}
    def analyse(zc):
        col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
        pairing = frozenset(frozenset(circ_of[v] for v in e) for _, e in comps if e); n = len(comps)
        fail = not any(balanced(G, partition(comps, f)) for f in itertools.product((0, 1), repeat=n))
        return pairing, n, fail
    fails = [list(z) for z in d.get("fails", [d["zero"]])]; cnt = collections.Counter(); nre = []
    for base in fails:
        p0, n0, f0 = analyse(base)
        if not f0: print(fn, "base", base, "is not a failure; skipped"); continue
        k = 0
        for i in range(len(odd)):
            for j in range(len(odd[i])):
                if j == base[i]: continue
                zc = list(base); zc[i] = j; p, n, f = analyse(zc)
                kind = "same pairing" if (p == p0 and n == n0) else ("re-paired" if n == n0 else "components changed")
                cnt[(kind, "fails" if f else "repaired")] += 1
                if kind == "re-paired": k += 1
        nre.append(k)
    print(fn.split("/")[-1], dict(cnt), "; re-pairing moves per failing choice:", nre, flush=True); tot.update(cnt)
print("TOTAL:", dict(tot))
