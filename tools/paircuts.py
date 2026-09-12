"""Systematic census of tight PAIR cuts at oddness 6.  For every instance and every path colouring,
for every pair {a,b} and each side for the third path, force the pattern (black ends of P_a, P_b
inside, white ends outside, both ends of P_c on the chosen side) and take the min cut; a negative
value is a bad cut with that pair pattern.  Record (pair, sign s_ab, s = even-circuit crossings,
third path crosses?, type) and count how many distinct pairs coexist per instance."""
import sys, random, itertools, collections
import networkx as nx
from gen_cyc6 import generate_straddle, generate
from oddness import canonical_colouring, H_components, partition
from badcuts import classify, profiles
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
def sample_zero_choices(odd, k, rng):
    seen = set()
    while len(seen) < k: seen.add(tuple(rng.randrange(len(C)) for C in odd))
    return sorted(seen)
name = sys.argv[1]; rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 5); ng = int(sys.argv[3]) if len(sys.argv) > 3 else 6
pr = profiles[name]
graphs = generate_straddle(pr[0], pr[1], rng, ng) if isinstance(pr, tuple) else generate(pr, rng, ng)
coexist = collections.Counter(); kinds = collections.Counter(); n_inst = 0; examples = []
for gi, (G, circuits) in enumerate(graphs):
    M = {frozenset((u, v)) for u, v in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    odd = [C for C in circuits if len(C) % 2]
    for zc in sample_zero_choices(odd, 40, rng):
        col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
        path_idx = [i for i, (_, ends) in enumerate(comps) if ends]
        if len(path_idx) != 3: continue
        paths = [tuple(comps[i][1]) for i in path_idx]; ref = [rng.randint(0, 1) for _ in comps]
        even_vertices = {v for (c, ends) in comps for v in c if not ends}
        path_vertices = [set(comps[i][0]) for i in path_idx]
        found = {}   # frozenset(S) -> record
        for bits in itertools.product((0, 1), repeat=3):
            f = list(ref)
            for i, b in zip(path_idx, bits): f[i] = b
            black = partition(comps, f)
            for a, b in itertools.combinations(range(3), 2):
                c = 3 - a - b
                ins = [z for i in (a, b) for z in paths[i] if z in black]; outs = [z for i in (a, b) for z in paths[i] if z not in black]
                for side in (0, 1):
                    val, S = mincut(G, black, ins + (list(paths[c]) if side else []), outs + ([] if side else list(paths[c])))
                    if val >= 0: continue
                    if len(S) > G.number_of_nodes() // 2: S = set(G.nodes()) - S
                    key = frozenset(S)
                    if key in found: continue
                    ty = classify(G, S, col, comps, black)
                    if ty[0] != 6: continue
                    s_even = sum(1 for u, v in G.edges() if (u in S) != (v in S) and col[frozenset((u, v))] == 1 and u in even_vertices)
                    crosses3 = any((u in S) != (v in S) for u, v in G.edges() if u in path_vertices[c] and v in path_vertices[c])
                    # sign: which end of each separated path is inside, relative to the path's canonical colouring
                    sign = tuple(int(paths[i][0] in S) for i in (a, b))
                    found[key] = ((a, b), sign, s_even, crosses3, ty[:5])
        n_inst += 1
        pairs = {rec[0] for rec in found.values()}
        codim2 = {rec[0] for rec in found.values() if rec[2] == 0 and not rec[3]}
        coexist[(len(pairs), len(codim2))] += 1
        for rec in found.values(): kinds[(rec[0], rec[2], rec[3])] += 1
        if len(pairs) >= 2 and len(examples) < 5: examples.append((gi, zc, sorted(found.values())))
print(f"{name}: {n_inst} instances; (distinct pairs with a tight pair cut, distinct pairs with a codim-2 s=0 pair cut) histogram: {dict(sorted(coexist.items()))}")
print("  pair cuts by (pair, even crossings s, third path crosses):", dict(sorted(kinds.items())))
for ex in examples: print("  example with >=2 pairs: graph", ex[0], "zero", ex[1], ex[2])
