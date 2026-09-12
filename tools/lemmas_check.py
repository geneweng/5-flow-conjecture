"""Sanity-check structural lemmas about tight bad cuts on generated instances:
 L1 (alternation): along every H-component, consecutive crossing edges of a bad cut alternate colours 1,2.
 L2: a path with both ends outside S never crosses a bad cut S; a path with both ends inside crosses evenly.
 L3: inside segments of H-components have odd numbers of vertices (both ends black).
Also report the codimension of each bad cut = number of paths that cross it (J*)."""
import sys, random, itertools, collections
import networkx as nx
from gen_cyc6 import generate_straddle, generate
from oddness import canonical_colouring, H_components, partition
from badcuts import violating_set, classify, profiles

def check(G, S, col, comps, black):
    H = nx.Graph([tuple(e) for e, c in col.items() if c in (1, 2)])
    viol = []
    codim = 0
    prof_paths = []; s_circ = 0
    for colouring, ends in comps:
        sub = H.subgraph(colouring.keys())
        # order the component's vertices along the path/circuit
        if ends:
            order = nx.shortest_path(sub, ends[0], ends[1])
            edges = list(zip(order, order[1:]))
        else:
            order = [v for v in nx.cycle_basis(sub)[0]]
            edges = list(zip(order, order[1:] + order[:1]))
        cross = [(i, e) for i, e in enumerate(edges) if (e[0] in S) != (e[1] in S)]
        if not cross:
            continue
        cols = [col[frozenset(e)] for _, e in cross]
        if ends: codim += 1; prof_paths.append(len(cross))
        else: s_circ += len(cross) // 2
        for a, b in zip(cols, cols[1:]):
            if a == b: viol.append("L1")
        if ends:
            inside = [v in S for v in ends]
            if not any(inside) and cross: viol.append("L2a")
            if all(inside) and len(cross) % 2: viol.append("L2b")
        # L3: segments between consecutive crossings
        for (i, e), (j, e2) in zip(cross, cross[1:]):
            seg = order[i + 1: j + 1]
            if seg and (seg[0] in S) and len(seg) % 2 == 0: viol.append("L3")
    return viol, codim, (tuple(sorted(prof_paths, reverse=True)), s_circ)

rng = random.Random(9); name = sys.argv[1] if len(sys.argv) > 1 else "straddleA"; pr = profiles[name]
graphs = generate_straddle(pr[0], pr[1], rng, 4) if isinstance(pr, tuple) else generate(pr, rng, 4)
viols = collections.Counter(); codims = collections.Counter(); n = 0
for G, circuits in graphs:
    M = {frozenset((u, v)) for u, v in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    odd = [C for C in circuits if len(C) % 2]
    for zc in rng.sample(list(itertools.product(*[range(len(C)) for C in odd])), 40):
        col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
        path_idx = [i for i, (_, ends) in enumerate(comps) if ends]
        ref = [rng.randint(0, 1) for _ in comps]
        for bits in itertools.product((0, 1), repeat=len(path_idx)):
            f = list(ref)
            for i, b in zip(path_idx, bits): f[i] = b
            black = partition(comps, f); S = violating_set(G, black)
            if S is None: continue
            n += 1
            v, cd, prof = check(G, S, col, comps, black)
            for x in v: viols[x] += 1
            t = classify(G, S, col, comps, black)
            codims[(t[0], t[1], t[2], prof)] += 1
print(f"{name}: {n} bad cuts checked; lemma violations: {dict(viols) or 'none'}")
print("(|dS|, c1, c2, (path crossing counts, circuit crossings s)) :")
for k_, v_ in sorted(codims.items()): print("   ", k_, v_)
