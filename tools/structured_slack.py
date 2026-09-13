"""For a dumped 6-of-8 instance and each surviving colouring: minimum of 3 d(S) - 5 k(S) over all
11-cut-shaped candidate sets S = (one whole z-circuit) + (arcs around the other two required ends) +
(any subset of the even circuits), with the three required (black) ends inside and the others outside."""
import sys, json, itertools
import networkx as nx
from oddness import canonical_colouring, H_components, partition
from badcuts import violating_set
d = json.load(open(sys.argv[1])); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
col = canonical_colouring(G, M, circuits, d["zero"]); comps = H_components(G, col)
pidx = [i for i, (_, e) in enumerate(comps) if e]; ref = d["ref"]
f0 = [ref.get(str(min(c)), 0) for c, _ in comps] if isinstance(ref, dict) else None
if f0 is None:   # recover list ref as in fullkill_sat
    want = {tuple(int(ch) for ch in k if ch in "01") for k in d["S"]}
    for cand in itertools.product((0, 1), repeat=len(comps)):
        if any(cand[i] for i in pidx): continue
        got = set()
        for bits in itertools.product((0, 1), repeat=len(pidx)):
            ff = list(cand)
            for i, b in zip(pidx, bits): ff[i] = b
            if violating_set(G, partition(comps, ff)) is not None: got.add(bits)
        if got == want: f0 = list(cand); break
odd = [ci for ci, C in enumerate(circuits) if len(C) % 2]; even = [ci for ci, C in enumerate(circuits) if len(C) % 2 == 0]
zof = {circuits[ci][zi]: ci for ci, zi in zip(odd, d["zero"])}; zpos = {ci: zi for ci, zi in zip(odd, d["zero"])}
paths = [tuple(comps[i][1]) for i in pidx]
def arcs(ci):
    C = circuits[ci]; L = len(C); z0 = zpos[ci]; out = []
    for back in range(0, L // 2 + 1):
        for fwd in range(0, L // 2 + 1):
            if back + fwd + 1 >= L: continue
            out.append([C[(z0 + k) % L] for k in range(-back, fwd + 1)])
    return out
for bits in itertools.product((0, 1), repeat=len(pidx)):
    f = list(f0)
    for i, b in zip(pidx, bits): f[i] = b
    black = partition(comps, f)
    if violating_set(G, black) is not None: continue
    req = [z for p in paths for z in p if z in black]; best = None
    for whole in req:
        others = [z for z in req if z != whole]
        for a1 in arcs(zof[others[0]]):
            for a2 in arcs(zof[others[1]]):
                base = set(circuits[zof[whole]]) | set(a1) | set(a2)
                for r in range(len(even) + 1):
                    for ev in itertools.combinations(even, r):
                        S = base | {v for ci in ev for v in circuits[ci]}
                        dS = nx.cut_size(G, S); k = sum(1 for v in S if v in black) - sum(1 for v in S if v not in black)
                        val = 3 * dS - 5 * k
                        if best is None or val < best[0]: best = (val, dS, k, whole, len(a1), len(a2), ev)
    print(f"survivor {bits}: best structured candidate slack {best[0]} (d={best[1]}, k={best[2]}, whole circuit of z{best[3]}, arcs {best[4]},{best[5]}, evens {best[6]})")
