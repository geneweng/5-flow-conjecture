"""Validation of the count relaxation count_c6.py against real data: for every archived instance with H = three
paths (+ even circuits), take the min-cut killer of every unbalanced path colouring (black-heavy side), read
off its configuration (memberships of the path ends, (c1, c2), crossing mode and beta of the non-separated
paths, non-H boundary edge) and ask the model about every single killer, every pair and the whole set.
Every answer must be FEASIBLE: an INFEASIBLE here is a bug in the model (a constraint that is not necessary).
usage: validate_c6.py FILE...   [env EVENFLIPS=1: also all colourings of the even circuits of H]"""
import sys, json, itertools, os, collections
import networkx as nx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oddness import canonical_colouring, H_components, partition
from count_c6 import solve

def black_heavy_violator(G, black):
    D = nx.DiGraph()
    for u, v in G.edges(): D.add_edge(u, v, capacity=3); D.add_edge(v, u, capacity=3)
    for v in G.nodes():
        if v in black: D.add_edge('s', v, capacity=5)
        else: D.add_edge(v, 't', capacity=5)
    cut, (X, Y) = nx.minimum_cut(D, 's', 't')
    return None if cut - 5 * len(black) >= 0 else set(X) - {'s'}

def path_order(G, col, comp_ends):
    """vertex sequence of the H-path from comp_ends[0]"""
    H = nx.Graph([tuple(e) for e, c in col.items() if c in (1, 2)])
    return nx.shortest_path(H, comp_ends[0], comp_ends[1])

tally = collections.Counter(); failures = []
for fn in sys.argv[1:]:
    d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
    M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    col = canonical_colouring(G, M, circuits, d["zero"]); comps = H_components(G, col)
    pidx = [i for i, (_, e) in enumerate(comps) if e]; cidx = [i for i, (_, e) in enumerate(comps) if not e]
    if len(pidx) != 3: continue
    paths = [path_order(G, col, sorted(comps[i][1])) for i in pidx]
    evens = list(itertools.product((0, 1), repeat=len(cidx))) if os.environ.get("EVENFLIPS") == "1" else [tuple([0] * len(cidx))]
    for ev in evens:
        killers = {}
        for bits in itertools.product((0, 1), repeat=3):
            f = [0] * len(comps)
            for i, b in zip(pidx, bits): f[i] = b
            for i, b in zip(cidx, ev): f[i] = b
            black = partition(comps, f); S = black_heavy_violator(G, black)
            if S is None: continue
            bnd = [e for e in G.edges() if (e[0] in S) != (e[1] in S)]
            cc = collections.Counter(col[frozenset(e)] for e in bnd); k = sum(1 if v in black else -1 for v in S)
            nonH = len(bnd) - cc[1] - cc[2]
            typ = (len(bnd), k, cc[1], cc[2], nonH)
            if typ not in {(6, 4, 4, 2, 0), (6, 4, 4, 1, 1), (7, 5, 5, 2, 0), (11, 7, 7, 4, 0)}: tally["untyped " + str(typ)] += 1; continue
            cross = {}
            for i, P in enumerate(paths):
                ncr = sum(1 for a, b in zip(P, P[1:]) if (a in S) != (b in S))
                if (P[0] in S) == (P[-1] in S) and ncr: cross[i] = 1 if P[0] in black else 0
            key = frozenset(S)
            killers.setdefault(key, (dict(c1=cc[1], c2=cc[2], nonH=nonH, cross=cross), typ, []))[2].append(bits)
        ks = list(killers.items())
        if not ks: continue
        def cfg(sub):
            ends = {i: (tuple(int(P[0] in S) for S, _ in sub), tuple(int(P[-1] in S) for S, _ in sub)) for i, P in enumerate(paths)}
            return ends, [v[0] for _, v in sub]
        subsets = [s for r in (1, 2) for s in itertools.combinations(ks, r)] + ([tuple(ks)] if len(ks) > 2 else [])
        for sub in subsets:
            ends, cuts = cfg(sub); st, sol = solve(ends, cuts, 2, 120)
            desc = "+".join(f"{v[1][0]}{'x' if v[0]['cross'] else ''}{'n' if v[0]['nonH'] else ''}" for _, v in sub)
            tally[(len(sub), st)] += 1
            if st != "FEASIBLE":
                failures.append((fn, ev, desc, st)); print("!!", os.path.basename(fn), ev, desc, st, ends, cuts, flush=True)
        print(os.path.basename(fn), ev, "killers:", [(v[1], v[0]["cross"], len(v[2])) for _, v in ks], flush=True)
print("TALLY", dict(tally)); print("FAILURES", len(failures))
for f in failures: print(f)
