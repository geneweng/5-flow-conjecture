"""For a dumped instance and each SURVIVING path colouring: force the pattern (the S-end of every
path inside S, the other end outside) and compute min_S [3 d(S) - 5(b_S - a_S)] over such S by one
min cut with infinite source/sink capacities at the forced z's.  The optimum is the 'best attempt'
at a bad cut with that pattern; slack 0 would be a tight-but-not-bad cut, slack > 0 says how far."""
import sys, json, itertools
import networkx as nx
from oddness import canonical_colouring, H_components, partition
from badcuts import violating_set, classify
INF = 10 ** 6
def best_forced(G, black, inside, outside):
    D = nx.DiGraph()
    for u, v in G.edges(): D.add_edge(u, v, capacity=3); D.add_edge(v, u, capacity=3)
    for v in G.nodes():
        if v in black: D.add_edge('s', v, capacity=5)
        else: D.add_edge(v, 't', capacity=5)
    for v in inside: D.add_edge('s', v, capacity=INF)
    for v in outside: D.add_edge(v, 't', capacity=INF)
    cut, (X, Y) = nx.minimum_cut(D, 's', 't'); S = set(X) - {'s'}
    return cut - 5 * len(black), S
for fn in sys.argv[1:]:
    d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
    M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    col = canonical_colouring(G, M, circuits, d["zero"]); comps = H_components(G, col)
    path_idx = [i for i, (_, ends) in enumerate(comps) if ends]; t = len(path_idx); ref = d["ref"]
    paths = [tuple(comps[i][1]) for i in path_idx]
    print(fn)
    for bits in itertools.product((0, 1), repeat=t):
        f = list(ref)
        for i, b in zip(path_idx, bits): f[i] = b
        black = partition(comps, f)
        if violating_set(G, black) is not None: continue
        # pattern: for each path put its black end inside (all-separated pattern, q = t)
        inside = [z for p in paths for z in p if z in black]; outside = [z for p in paths for z in p if z not in black]
        slack, S = best_forced(G, black, inside, outside)
        ty = classify(G, S, col, comps, black)
        print(f"  survivor {bits}: best all-separated cut has slack {slack} (3d-5k), type {ty[:5]} |S|={len(S)}")
        # also the best pair patterns (two paths separated, third not)
        for a, b in itertools.combinations(range(t), 2):
            ins = [z for i in (a, b) for z in paths[i] if z in black]; outs = [z for i in (a, b) for z in paths[i] if z not in black]
            slack2, S2 = best_forced(G, black, ins, outs); ty2 = classify(G, S2, col, comps, black)
            print(f"      pair {(a, b)} forced: slack {slack2}, type {ty2[:5]} sep={ty2[5]} |S|={len(S2)}")
