"""Statistics of bad cuts: for cyclically 6-connected graphs with a planted 6-odd 2-factor,
for each 0-edge choice look at all 8 path colourings (even components fixed), extract the
violating set S of every unbalanced one, and classify (|dS|, c1, c2, |S∩Z|, paths separated, sign)."""
import sys, random, itertools, collections
import networkx as nx
from gen_cyc6 import generate_straddle, generate
from oddness import canonical_colouring, H_components, partition

def violating_set(G, black):
    D = nx.DiGraph()
    for u, v in G.edges():
        D.add_edge(u, v, capacity=3); D.add_edge(v, u, capacity=3)
    for v in G.nodes():
        if v in black: D.add_edge('s', v, capacity=5)
        else: D.add_edge(v, 't', capacity=5)
    cut, (X, Y) = nx.minimum_cut(D, 's', 't')
    if cut - 5 * len(black) >= 0: return None
    S = set(X) - {'s'}
    if len(S) > G.number_of_nodes() // 2: S = set(G.nodes()) - S
    return S

def classify(G, S, col, comps, black):
    boundary = [(u, v) for u, v in G.edges() if (u in S) != (v in S)]
    c = collections.Counter(col[frozenset(e)] for e in boundary)
    zs = [v for (_, ends) in comps for v in ends]
    paths = [tuple(ends) for (_, ends) in comps if ends]
    sep = [i for i, (x, y) in enumerate(paths) if (x in S) != (y in S)]
    zin = [z for z in zs if z in S]
    q = abs(sum(1 if z in black else -1 for z in zin))
    return (len(boundary), c[1], c[2], len(zin), q, tuple(sep))

profiles = {
    "straddleA": ([[7, 6], [7], [7, 7], [6]], [(7, 3), (7, 3), None, None]),
    "straddle6A": ([[7, 6], [6], [7, 6], [6], [7, 6], [6]], [(7, 3), None, (7, 3), None, (7, 3), None]),
    "random6x7+8": [7] * 6 + [8],
    "straddleB": ([[7, 8], [7], [7, 7], [8]], [(7, 3), (7, 3), None, None]),
    "straddleC": ([[7, 6], [7, 6], [7, 7], [6]], [(7, 3), (7, 3), None, None]),
    "straddleD": ([[9, 6], [7], [7, 7], [6]], [(7, 3), (7, 3), None, None]),
    "straddleE": ([[7, 6], [7], [7, 7], [6]], [(7, 2), (7, 4), None, None]),
    "straddleF": ([[7, 6], [9], [7, 7], [6]], [(9, 5), (7, 3), None, None]),
    "straddleG": ([[7, 6], [7], [7, 7], [6]], [(7, 3), (7, 3), (6, 3), None]),
    "straddleH": ([[7, 6], [7], [7, 7], [8]], [(7, 3), (7, 3), (8, 4), (6, 3)]),
    "six12": ([[7, 6], [12], [7, 6], [12], [7, 6], [12]], [(7, 3), None, (7, 3), None, (7, 3), None]),
    "six10": ([[7, 8], [10], [7, 8], [10], [7, 8], [10]], [(7, 3), None, (7, 3), None, (7, 3), None]),
}
rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
name = sys.argv[1] if len(sys.argv) > 1 else "straddleA"
pr = profiles[name]
graphs = generate_straddle(pr[0], pr[1], rng, 6) if isinstance(pr, tuple) else generate(pr, rng, 6)
hist = collections.Counter(); types = collections.Counter(); n_inst = 0
for G, circuits in graphs:
    M = {frozenset((u, v)) for u, v in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    odd = [C for C in circuits if len(C) % 2]
    zero_space = list(itertools.product(*[range(len(C)) for C in odd])); zero_space = rng.sample(zero_space, 60)
    for zc in zero_space:
        col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
        path_idx = [i for i, (_, ends) in enumerate(comps) if ends]
        ref = [rng.randint(0, 1) for _ in comps]; nbal = 0
        for bits in itertools.product((0, 1), repeat=len(path_idx)):
            f = list(ref)
            for i, b in zip(path_idx, bits): f[i] = b
            black = partition(comps, f); S = violating_set(G, black)
            if S is None: nbal += 1
            else: types[classify(G, S, col, comps, black)] += 1
        hist[nbal] += 1; n_inst += 1
print(f"{name}: {len(graphs)} graphs, {n_inst} (2-factor, 0-edge) instances; #balanced among the 8 path colourings:", dict(sorted(hist.items())))
print("bad cut types (|dS|, c1, c2, |S∩Z|, q, separated paths): count")
for k, v in types.most_common(15): print("  ", k, v)
