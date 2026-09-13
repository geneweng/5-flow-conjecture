"""Targeted completion search from a near-miss instance: keep the 2-factor, 0-edges and even-circuit
colouring; anneal over matching 2-swaps (girth >= 6, cyclic 6-edge-connectivity rechecked on accepted
moves); objective = (killed colourings, -min structured 11-cut slack over survivors), where the structured
slack is min 3d(S)-5k(S) over S = whole circuit of one required end + arcs around the other two + any
subset of even circuits.  Slack < 0 for a survivor = a full kill.   usage: complete.py DUMP SEED STEPS"""
import sys, json, random, itertools, time
import networkx as nx
from oddness import canonical_colouring, H_components, partition
from badcuts import violating_set
from gen_cyc6 import cyclically_6_connected
fn, seed, steps = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); rng = random.Random(seed)
d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
cedges = {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
odd = [ci for ci, C in enumerate(circuits) if len(C) % 2]; even = [ci for ci, C in enumerate(circuits) if len(C) % 2 == 0]
zpos = {ci: zi for ci, zi in zip(odd, d["zero"])}; zof = {circuits[ci][zi]: ci for ci, zi in zpos.items()}
ARCS = {}
for ci in odd:
    C = circuits[ci]; L = len(C); z0 = zpos[ci]; out = []
    for back in range(0, L // 2 + 1):
        for fwd in range(0, L // 2 + 1):
            if back + fwd + 1 < L: out.append(frozenset(C[(z0 + k) % L] for k in range(-back, fwd + 1)))
    ARCS[ci] = out
EVSUB = [frozenset(v for ci in ev for v in circuits[ci]) for r in range(len(even) + 1) for ev in itertools.combinations(even, r)]
ref = d["ref"]
def base_flips(comps):
    if isinstance(ref, dict): return [ref.get(str(min(c)), 0) for c, _ in comps]
    return [0] * len(comps)     # list refs are not reproducible; start from the canonical colouring
def structured_slack(G, black, paths):
    req = [z for p in paths for z in p if z in black]; best = 10 ** 9
    bl = black
    for whole in req:
        others = [z for z in req if z != whole]
        for a1 in ARCS[zof[others[0]]]:
            for a2 in ARCS[zof[others[1]]]:
                base = set(circuits[zof[whole]]) | a1 | a2
                for ev in EVSUB:
                    S = base | ev
                    dS = sum(1 for u, v in G.edges(S) if v not in S); k = 2 * sum(1 for v in S if v in bl) - len(S)
                    best = min(best, 3 * dS - 5 * k)
    return best
def evaluate(G):
    M = {frozenset(e) for e in G.edges()} - cedges
    col = canonical_colouring(G, M, circuits, d["zero"]); comps = H_components(G, col)
    pidx = [i for i, (_, e) in enumerate(comps) if e]; paths = [tuple(comps[i][1]) for i in pidx]
    f0 = base_flips(comps); killed = 0; worst = 0
    for bits in itertools.product((0, 1), repeat=len(pidx)):
        f = list(f0)
        for i, b in zip(pidx, bits): f[i] = b
        black = partition(comps, f)
        if violating_set(G, black) is not None: killed += 1; continue
        worst = max(worst, structured_slack(G, black, paths))
    return (killed, -worst)
def swap(G):
    Ml = [tuple(e) for e in G.edges() if frozenset(e) not in cedges]
    (a, b), (c, dd) = rng.sample(Ml, 2)
    if rng.random() < 0.5: c, dd = dd, c
    H = G.copy(); H.remove_edge(a, b); H.remove_edge(c, dd)
    if nx.has_path(H, a, c) and nx.shortest_path_length(H, a, c) < 5: return None
    H.add_edge(a, c)
    if nx.has_path(H, b, dd) and nx.shortest_path_length(H, b, dd) < 5: return None
    H.add_edge(b, dd); return H
score = evaluate(G); best = score; t0 = time.time(); T0 = 1.0
print(f"{fn} seed {seed}: start {score}", flush=True)
for step in range(1, steps + 1):
    T = T0 * (1 - step / steps)
    H = swap(G)
    if H is None: continue
    s2 = evaluate(H); delta = (s2[0] - score[0]) * 20 + (s2[1] - score[1])
    if delta >= 0 or (T > 0 and rng.random() < 2.718 ** (delta / T)):
        if not cyclically_6_connected(H): continue
        G, score = H, s2
        if score > best:
            best = score; print(f"  step {step} [{time.time()-t0:.0f}s] best {best}", flush=True)
            if score[0] == 8 or score[1] > 0:
                json.dump({"edges": [list(e) for e in G.edges()], "circuits": circuits, "zero": d["zero"], "ref": ref, "score": score}, open(f"fullkill_complete_{seed}.json", "w"))
                print("FULL KILL ->", f"fullkill_complete_{seed}.json", flush=True); break
print(f"done: best {best} [{time.time()-t0:.0f}s]", flush=True)
