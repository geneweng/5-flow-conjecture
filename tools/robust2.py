"""Two-point robustness: from a failing instance (graph, 2-factor, 0-edges c*), anneal over matching
2-swaps keeping c* failing (hard constraint) and pushing a second 0-edge choice c' (c* with one
0-edge moved by +1 on circuit CIRC) to fail as well: objective = (killed colourings under c',
-structured slack over survivors under c').   usage: robust2.py DUMP CIRC SEED STEPS"""
import sys, json, random, itertools, time
import networkx as nx
from oddness import canonical_colouring, H_components, partition
from badcuts import violating_set
from gen_cyc6 import cyclically_6_connected
fn, circ, seed, steps = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); rng = random.Random(seed)
d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
cedges = {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
odd = [ci for ci, C in enumerate(circuits) if len(C) % 2]; even = [ci for ci, C in enumerate(circuits) if len(C) % 2 == 0]
cstar = list(d["zero"]); cprime = list(cstar); cprime[circ] = (cprime[circ] + 1) % len(circuits[odd[circ]])
print("c* =", cstar, " c' =", cprime, flush=True)
EVSUB = [frozenset(v for ci in ev for v in circuits[ci]) for r in range(len(even) + 1) for ev in itertools.combinations(even, r)]
def arcs_around(ci, zi):
    C = circuits[ci]; L = len(C); out = []
    for back in range(0, L // 2 + 1):
        for fwd in range(0, L // 2 + 1):
            if back + fwd + 1 < L: out.append(frozenset(C[(zi + k) % L] for k in range(-back, fwd + 1)))
    return out
def all_killed(G, zc):
    M = {frozenset(e) for e in G.edges()} - cedges
    col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
    for f in itertools.product((0, 1), repeat=len(comps)):
        if violating_set(G, partition(comps, f)) is None: return False, comps, col
    return True, comps, col
def structured(G, zc, comps, black, paths):
    zof = {circuits[odd[k]][zc[k]]: (odd[k], zc[k]) for k in range(len(odd))}
    req = [z for p in paths for z in p if z in black]; best = 10 ** 9
    for whole in req:
        others = [z for z in req if z != whole]
        for a1 in arcs_around(*zof[others[0]]):
            for a2 in arcs_around(*zof[others[1]]):
                base = set(circuits[zof[whole][0]]) | a1 | a2
                for ev in EVSUB:
                    S = base | ev; dS = sum(1 for u, v in G.edges(S) if v not in S); k = 2 * sum(1 for v in S if v in black) - len(S)
                    best = min(best, 3 * dS - 5 * k)
    return best
def evaluate(G):
    ok, _, _ = all_killed(G, cstar)
    if not ok: return None
    M = {frozenset(e) for e in G.edges()} - cedges
    col = canonical_colouring(G, M, circuits, cprime); comps = H_components(G, col)
    pidx = [i for i, (_, e) in enumerate(comps) if e]; paths = [tuple(comps[i][1]) for i in pidx]
    killed = 0; worst = 0; total = 0
    for f in itertools.product((0, 1), repeat=len(comps)):
        total += 1; black = partition(comps, f)
        if violating_set(G, black) is not None: killed += 1
        else: worst = max(worst, structured(G, cprime, comps, black, paths))
    return (killed - total, -worst)      # 0 in the first coordinate = c' fails too
def swap(G):
    Ml = [tuple(e) for e in G.edges() if frozenset(e) not in cedges]
    (a, b), (c, dd) = rng.sample(Ml, 2)
    if rng.random() < 0.5: c, dd = dd, c
    H = G.copy(); H.remove_edge(a, b); H.remove_edge(c, dd)
    if nx.has_path(H, a, c) and nx.shortest_path_length(H, a, c) < 5: return None
    H.add_edge(a, c)
    if nx.has_path(H, b, dd) and nx.shortest_path_length(H, b, dd) < 5: return None
    H.add_edge(b, dd); return H
score = evaluate(G); assert score is not None, "c* does not fail on the input"; best = score; t0 = time.time(); T0 = 1.0
print("start", score, flush=True)
for step in range(1, steps + 1):
    T = T0 * (1 - step / steps); H = swap(G)
    if H is None: continue
    s2 = evaluate(H)
    if s2 is None: continue
    delta = (s2[0] - score[0]) * 20 + (s2[1] - score[1])
    if delta >= 0 or (T > 0 and rng.random() < 2.718 ** (delta / T)):
        if not cyclically_6_connected(H): continue
        G, score = H, s2
        if score > best:
            best = score; print(f"  step {step} [{time.time()-t0:.0f}s] best {best}", flush=True)
            if score[0] == 0:
                json.dump({"edges": [list(e) for e in G.edges()], "circuits": circuits, "zero": cstar, "zero2": cprime}, open(f"robust2_{circ}_{seed}.json", "w"))
                print("TWO FAILING 0-EDGE CHOICES ->", f"robust2_{circ}_{seed}.json", flush=True); break
print(f"done: best {best} [{time.time()-t0:.0f}s]", flush=True)
