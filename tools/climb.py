"""Directed search for a failure of the path-recolouring template at oddness 6 (or any t):
hill-climb / annealing over (matching M, 0-edge choice, even-circuit colouring) with the planted
2-factor fixed, maximizing the number of killed path colourings; tiebreak by the survivors' forced
slack.  Moves: matching 2-swap keeping girth >= 6 and cyclic 6-edge-connectivity (SAT-checked),
0-edge move, even-component flip.  A state with all 2^t colourings killed is dumped as fullkill_*.json.
usage: climb.py PROFILE SEED STEPS   (env: TEMP=0.5 initial temperature, RESTART=1500 stale steps)"""
import sys, os, json, random, itertools, time, collections
import networkx as nx
from gen_cyc6 import generate_straddle, generate, generate_rings, cyclically_6_connected
from oddness import canonical_colouring, H_components, partition
from badcuts import violating_set, classify, profiles
profiles.update({'t4ring': [[7, 7, 6]] * 4, 'ring3': [[7, 7, 6]] * 3, 'rand7': [7] * 6 + [8], 'rand9': [7] * 4 + [9, 9, 6]})
INF = 10 ** 6
def mincut(G, black, inside=(), outside=()):
    D = nx.DiGraph()
    for u, v in G.edges(): D.add_edge(u, v, capacity=3); D.add_edge(v, u, capacity=3)
    for v in G.nodes():
        if v in black: D.add_edge('s', v, capacity=5)
        else: D.add_edge(v, 't', capacity=5)
    for v in inside: D.add_edge('s', v, capacity=INF)
    for v in outside: D.add_edge(v, 't', capacity=INF)
    cut, (X, Y) = nx.minimum_cut(D, 's', 't'); return cut - 5 * len(black), set(X) - {'s'}
def evaluate(G, circuits, zc, refbits):
    """returns (n_killed, -sum survivor slack, details)"""
    M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
    keys = [min(c) for c, _ in comps]
    path_idx = [i for i, (_, ends) in enumerate(comps) if ends]; t = len(path_idx)
    paths = [tuple(comps[i][1]) for i in path_idx]
    killed = 0; slack_sum = 0; det = []
    for bits in itertools.product((0, 1), repeat=t):
        f = [refbits.get(k, 0) for k in keys]
        for i, b in zip(path_idx, bits): f[i] = b
        black = partition(comps, f); m, S = mincut(G, black)
        if m < 0: killed += 1; det.append((bits, "killed", classify(G, S, col, comps, black)[:5])); continue
        best = None
        for sub in [tuple(range(t))] + list(itertools.combinations(range(t), 2)):
            ins = [z for i in sub for z in paths[i] if z in black]; outs = [z for i in sub for z in paths[i] if z not in black]
            sl, S2 = mincut(G, black, ins, outs)
            if best is None or sl < best: best = sl
        slack_sum += best; det.append((bits, "alive", best))
    return (killed, -slack_sum), det, t, keys
def swap_move(G, circuits, rng, circuit_edges):
    """2-swap of two matching edges keeping girth >= 6; returns new graph or None"""
    M = [tuple(e) for e in G.edges() if frozenset(e) not in circuit_edges]
    (a, b), (c, d) = rng.sample(M, 2)
    if rng.random() < 0.5: c, d = d, c
    H = G.copy(); H.remove_edge(a, b); H.remove_edge(c, d)
    if nx.has_path(H, a, c) and nx.shortest_path_length(H, a, c) < 5: return None
    H.add_edge(a, c)
    if nx.has_path(H, b, d) and nx.shortest_path_length(H, b, d) < 5: return None
    H.add_edge(b, d)
    return H
name, seed, steps = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
rng = random.Random(seed); T0 = float(os.environ.get("TEMP", "0.5")); RESTART = int(os.environ.get("RESTART", "1500"))
pr = profiles[name]
def fresh():
    gs = generate_straddle(pr[0], pr[1], rng, 1) if isinstance(pr, tuple) else (generate_rings(pr, rng, 1) if isinstance(pr[0], list) else generate(pr, rng, 1))
    G, circuits = gs[0]; odd = [C for C in circuits if len(C) % 2]
    return G, circuits, odd, [rng.randrange(len(C)) for C in odd], {}
G, circuits, odd, zc, ref = fresh()
circuit_edges = {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
score, det, t, keys = evaluate(G, circuits, zc, ref); best = score; best_state = None; stale = 0; t0 = time.time()
print(f"{name} seed {seed}: n={G.number_of_nodes()} t={t} start score {score}", flush=True)
for step in range(1, steps + 1):
    T = T0 * (1 - step / steps)
    r = rng.random(); G2, zc2, ref2 = G, list(zc), dict(ref)
    if r < 0.5:
        G2 = swap_move(G, circuits, rng, circuit_edges)
        if G2 is None: continue
    elif r < 0.8:
        j = rng.randrange(len(odd)); zc2[j] = (zc2[j] + rng.choice((1, -1, rng.randrange(len(odd[j])))) ) % len(odd[j])
    else:
        k = rng.choice(keys); ref2[k] = 1 - ref2.get(k, 0)
    s2, det2, t2, keys2 = evaluate(G2, circuits, zc2, ref2)
    if t2 != t: continue
    delta = (s2[0] - score[0]) + 0.02 * (s2[1] - score[1])
    if delta >= 0 or (T > 0 and rng.random() < 2.718 ** (delta / T)):
        if G2 is not G and not cyclically_6_connected(G2): continue
        G, zc, ref, score, det, keys = G2, zc2, ref2, s2, det2, keys2
        if score > best:
            best = score; stale = 0
            print(f"  step {step} [{time.time()-t0:.0f}s] new best {best}: " + ", ".join(f"{b}:{w if w=='killed' else 'slack'}{'' if w=='killed' else x}" for b, w, x in det), flush=True)
            if score[0] >= 2 ** t - 2:                      # dump near-misses as solver hints
                json.dump({"profile": name, "seed": seed, "edges": [list(e) for e in G.edges()], "circuits": circuits, "zero": zc, "ref": ref, "score": score}, open(f"best_climb_{name}_{seed}.json", "w"))
            if score[0] == 2 ** t:
                fn = f"fullkill_climb_{name}_{seed}_{step}.json"
                json.dump({"profile": name, "seed": seed, "edges": [list(e) for e in G.edges()], "circuits": circuits, "zero": zc, "ref": ref, "details": [list(map(str, d)) for d in det]}, open(fn, "w"))
                print("FULL KILL FOUND ->", fn, flush=True); break
        else: stale += 1
    else: stale += 1
    if stale >= RESTART:
        G, circuits, odd, zc, ref = fresh(); circuit_edges = {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
        score, det, t, keys = evaluate(G, circuits, zc, ref); stale = 0
        print(f"  step {step}: restart, score {score}", flush=True)
print(f"done: best {best} in {time.time()-t0:.0f}s", flush=True)
