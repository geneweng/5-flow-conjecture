"""Direct attack on (L-repair) for two-position shifts: fix a 2-factor, a 0-edge choice c and a +-2 shift c' on one
circuit; anneal over matching 2-swaps (girth >= 6, cyclic 6-edge-connectivity rechecked) to make ALL vertex
colourings admissible at c or at c' bad (12 when both H's are three paths).  Objective = (#bad classes, -sum of
balance margins of the surviving classes).  A full kill with H(c), H(c') three paths and different pairings would
refute (L-repair).   usage: climb12.py DUMP CIRCUIT SHIFT(+2|-2) SEED STEPS   (CIRCUIT = odd-circuit index, or 'auto')"""
import sys, json, random, itertools, time
import networkx as nx
from oddness import canonical_colouring, H_components, partition
from gen_cyc6 import cyclically_6_connected
fn, circ, shift, seed, steps = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]); rng = random.Random(seed)
d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
cedges = {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
odd = [C for C in circuits if len(C) % 2]
base = list(d["zero"])
def pairing(G, zc):
    M = {frozenset(e) for e in G.edges()} - cedges
    comps = H_components(G, canonical_colouring(G, M, circuits, zc))
    circ_of = {v: k for k, C in enumerate(odd) for v in C}
    return frozenset(frozenset(circ_of[v] for v in e) for _, e in comps if e), len(comps), sum(1 for _, e in comps if e)
if circ == "auto":
    p0 = pairing(G, base); opts = []
    for i in range(len(odd)):
        for s in (2, -2):
            zc = list(base); zc[i] = (zc[i] + s) % len(odd[i]); p = pairing(G, zc)
            if p[1] == p0[1] and p[0] != p0[0] and p[2] == 3: opts.append((i, s))
    if not opts: print("no re-pairing two-position shift"); sys.exit(0)
    circ, shift = rng.choice(opts)
circ = int(circ); alt = list(base); alt[circ] = (alt[circ] + shift) % len(odd[circ])
print(f"{fn}: c={base} c'={alt} (circuit {circ}, shift {shift:+d})", flush=True)
def margin(G, black):
    D = nx.DiGraph()
    for u, v in G.edges(): D.add_edge(u, v, capacity=3); D.add_edge(v, u, capacity=3)
    for v in G.nodes():
        if v in black: D.add_edge('s', v, capacity=5)
        else: D.add_edge(v, 't', capacity=5)
    cut, _ = nx.minimum_cut(D, 's', 't'); return cut - 5 * len(black)
def colourings(G, zc):
    M = {frozenset(e) for e in G.edges()} - cedges
    comps = H_components(G, canonical_colouring(G, M, circuits, zc)); out = set(); v0 = min(G.nodes())
    for f in itertools.product((0, 1), repeat=len(comps)):
        black = partition(comps, f)
        if v0 not in black: black = set(G.nodes()) - black      # representative up to complement
        out.add(frozenset(black))
    return out
def zdata(zc):
    zof = {}; ARCS = {}
    for k, C in enumerate(odd):
        L = len(C); z0 = zc[k]; zof[C[z0]] = k; out = []
        for back in range(0, L // 2 + 1):
            for fwd in range(0, L // 2 + 1):
                if back + fwd + 1 < L: out.append(frozenset(C[(z0 + j) % L] for j in range(-back, fwd + 1)))
        ARCS[k] = out
    return zof, ARCS
even = [C for C in circuits if len(C) % 2 == 0]
EVSUB = [frozenset(v for C in ev for v in C) for r in range(len(even) + 1) for ev in itertools.combinations(even, r)]
ZD = {tuple(base): zdata(base), tuple(alt): zdata(alt)}
def structured_slack(G, black, zc):
    zof, ARCS = ZD[tuple(zc)]; req = [z for z in zof if z in black]
    if len(req) > 3: req = [z for z in zof if z not in black]
    if len(req) != 3: return 10 ** 6
    best = 10 ** 9
    for whole in req:
        others = [z for z in req if z != whole]
        for a1 in ARCS[zof[others[0]]]:
            for a2 in ARCS[zof[others[1]]]:
                bs = set(odd[zof[whole]]) | a1 | a2
                for ev in EVSUB:
                    S = bs | ev
                    dS = sum(1 for u, v in G.edges(S) if v not in S); k = 2 * sum(1 for v in S if v in black) - len(S)
                    best = min(best, 3 * dS - 5 * abs(k))
    return best
def evaluate(G):
    A, B = colourings(G, base), colourings(G, alt); cols = A | B; bad = 0; surv = 0
    for black in cols:
        m = margin(G, black)
        if m < 0: bad += 1; continue
        s = min(structured_slack(G, black, zc) for zc, fam in ((base, A), (alt, B)) if black in fam)
        surv += min(s, 30)
    return (bad, -surv, len(cols))
def swap(G):
    Ml = [tuple(e) for e in G.edges() if frozenset(e) not in cedges]
    (a, b), (c, dd) = rng.sample(Ml, 2)
    if rng.random() < 0.5: c, dd = dd, c
    H = G.copy(); H.remove_edge(a, b); H.remove_edge(c, dd)
    if nx.has_path(H, a, c) and nx.shortest_path_length(H, a, c) < 5: return None
    H.add_edge(a, c)
    if nx.has_path(H, b, dd) and nx.shortest_path_length(H, b, dd) < 5: return None
    H.add_edge(b, dd); return H
score = evaluate(G); best = score; t0 = time.time(); T0 = 3.0
print(f"seed {seed}: start {score}", flush=True)
for step in range(1, steps + 1):
    T = T0 * (1 - step / steps); H = swap(G)
    if H is None: continue
    s2 = evaluate(H)
    if s2[2] != score[2]: continue                       # keep both H's three paths (same number of classes)
    delta = (s2[0] - score[0]) * 30 + (s2[1] - score[1])
    if delta >= 0 or (T > 0 and rng.random() < 2.718 ** (delta / T)):
        if not cyclically_6_connected(H): continue
        G, score = H, s2
        if score[:2] > best[:2]:
            best = score; print(f"  step {step} [{time.time()-t0:.0f}s] best {best}", flush=True)
            if score[0] == score[2]:
                p0, p1 = pairing(G, base), pairing(G, alt)
                out = f"double_{seed}.json"
                json.dump({"edges": [list(e) for e in G.edges()], "circuits": circuits, "zero": base, "zero2": alt, "circuit": circ, "shift": shift,
                           "pairings": [sorted(map(sorted, p0[0])), sorted(map(sorted, p1[0]))], "ncomp": [p0[1], p1[1]]}, open(out, "w"))
                print("ALL COLOURINGS BAD ->", out, "pairings", sorted(map(sorted, p0[0])), sorted(map(sorted, p1[0])), "components", p0[1], p1[1],
                      "*** (L-repair) VIOLATION ***" if (p0[1] == p1[1] == 3 and p0[0] != p1[0]) else "(not a re-pairing pair)", flush=True); break
print(f"done: best {best} [{time.time()-t0:.0f}s]", flush=True)
