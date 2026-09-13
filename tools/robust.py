"""0-edge-robust failure search: from an instance (graph + 2-factor), anneal over matching 2-swaps
(girth >= 6, cyclic 6-edge-connectivity on accepted moves) maximizing the number of 0-edge choices,
in a fixed sample, for which NO 2-colouring of H (all components) is balanced; secondary objective:
total killed colourings over the sample.   usage: robust.py DUMP SEED STEPS [SAMPLE]"""
import sys, os, json, random, itertools, time
import networkx as nx
from oddness import canonical_colouring, H_components, partition, balanced
from gen_cyc6 import cyclically_6_connected
fn, seed, steps = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); K = int(sys.argv[4]) if len(sys.argv) > 4 else 40
rng = random.Random(seed); d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
cedges = {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
odd = [C for C in circuits if len(C) % 2]
def neigh(zc): return {tuple(list(zc[:i]) + [j] + list(zc[i + 1:])) for i in range(len(odd)) for j in range(len(odd[i]))}
if os.environ.get("NEIGH") == "1":   # the one-step neighbourhood of the base choice (all positions of each single 0-edge)
    sample = sorted(neigh(d["zero"]))
elif os.environ.get("NEIGH", "").startswith("axis:"):   # all positions of one circuit's 0-edge, the others fixed at the base
    ax = int(os.environ["NEIGH"].split(":")[1]); sample = [tuple(list(d["zero"][:ax]) + [j] + list(d["zero"][ax + 1:])) for j in range(len(odd[ax]))]
elif os.environ.get("NEIGH") == "2":   # union of one-step neighbourhoods of every currently failing choice (from the dump's sample, else the base)
    M0 = {frozenset(e) for e in G.edges()} - cedges; seeds = []
    for zc in d.get("sample", [d["zero"]]):
        col = canonical_colouring(G, M0, circuits, zc); comps = H_components(G, col)
        if not any(balanced(G, partition(comps, f)) for f in itertools.product((0, 1), repeat=len(comps))): seeds.append(tuple(zc))
    if not seeds: seeds = [tuple(d["zero"])]
    sample = sorted(set().union(*[neigh(z) for z in seeds])); print("failing seeds:", seeds, "sample size", len(sample), flush=True)
else:
    sample = [tuple(d["zero"])] + [tuple(rng.randrange(len(C)) for C in odd) for _ in range(K - 1)]
def evaluate(G):
    M = {frozenset(e) for e in G.edges()} - cedges; fails = 0; killed = 0
    for zc in sample:
        col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
        space = list(itertools.product((0, 1), repeat=len(comps)))
        if len(space) > 64: space = rng.sample(space, 64)
        k = sum(0 if balanced(G, partition(comps, f)) else 1 for f in space)
        killed += k / len(space)
        if k == len(space): fails += 1
    return (fails, round(killed, 3))
def swap(G):
    Ml = [tuple(e) for e in G.edges() if frozenset(e) not in cedges]
    (a, b), (c, dd) = rng.sample(Ml, 2)
    if rng.random() < 0.5: c, dd = dd, c
    H = G.copy(); H.remove_edge(a, b); H.remove_edge(c, dd)
    if nx.has_path(H, a, c) and nx.shortest_path_length(H, a, c) < 5: return None
    H.add_edge(a, c)
    if nx.has_path(H, b, dd) and nx.shortest_path_length(H, b, dd) < 5: return None
    H.add_edge(b, dd); return H
score = evaluate(G); best = score; t0 = time.time(); T0 = 0.3
print(f"{fn} seed {seed}: sample {K}, start {score}", flush=True)
for step in range(1, steps + 1):
    T = T0 * (1 - step / steps); H = swap(G)
    if H is None: continue
    s2 = evaluate(H); delta = (s2[0] - score[0]) + 0.05 * (s2[1] - score[1])
    if delta >= 0 or (T > 0 and rng.random() < 2.718 ** (delta / T)):
        if not cyclically_6_connected(H): continue
        G, score = H, s2
        if score > best:
            best = score; print(f"  step {step} [{time.time()-t0:.0f}s] best {best}", flush=True)
            json.dump({"edges": [list(e) for e in G.edges()], "circuits": circuits, "zero": d["zero"], "sample": sample, "score": score}, open(f"robust_best_{seed}.json", "w"))
print(f"done: best {best} [{time.time()-t0:.0f}s]", flush=True)
