"""For pairs of bad cuts (for different colourings) in the same instance whose four uncrossing parts
each contain a z, record the part-cut sizes |E(U_i,U_j)| and the part boundaries, by cut-type pair.
Validates the MS 3/0 structure for pairs of (4,2)-6-cuts and shows what holds for 7/11-cuts."""
import sys, random, itertools, collections
import networkx as nx
from gen_cyc6 import generate_straddle, generate
from oddness import canonical_colouring, H_components, partition
from badcuts import violating_set, classify, profiles
rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
name = sys.argv[1]; pr = profiles[name]
graphs = generate_straddle(pr[0], pr[1], rng, 6) if isinstance(pr, tuple) else generate(pr, rng, 6)
stats = collections.Counter(); npairs = 0; nz4 = 0
for G, circuits in graphs:
    M = {frozenset((u, v)) for u, v in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    odd = [C for C in circuits if len(C) % 2]
    for zc in rng.sample(list(itertools.product(*[range(len(C)) for C in odd])), 60):
        col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
        path_idx = [i for i, (_, ends) in enumerate(comps) if ends]
        if len(path_idx) != 3: continue
        zs = [v for (_, e) in comps for v in e]
        ref = [rng.randint(0, 1) for _ in comps]; cuts = []
        for bits in itertools.product((0, 1), repeat=3):
            f = list(ref)
            for i, b in zip(path_idx, bits): f[i] = b
            black = partition(comps, f); S = violating_set(G, black)
            if S is None: continue
            t = classify(G, S, col, comps, black); cuts.append((bits, frozenset(S), (t[0], t[1], t[2])))
        for (b1, S1, t1), (b2, S2, t2) in itertools.combinations(cuts, 2):
            if (b1[0]^b1[1], b1[0]^b1[2]) == (b2[0]^b2[1], b2[0]^b2[2]): continue
            npairs += 1
            parts = {'A': S1 & S2, 'C': S1 - S2, 'D': S2 - S1, 'B': set(G.nodes()) - S1 - S2}
            if any(not (p & set(zs)) for p in parts.values()): continue
            nz4 += 1
            names = ['A', 'B', 'C', 'D']
            sizes = tuple(sum(1 for u, v in G.edges() if (u in parts[x] and v in parts[y]) or (u in parts[y] and v in parts[x])) for x, y in [('A','B'),('C','D'),('A','C'),('A','D'),('B','C'),('B','D')])
            stats[(tuple(sorted([t1, t2])), sizes)] += 1
print("different-class pairs:", npairs, "with four z-parts:", nz4)
print(name, "pairs of bad cuts (different colourings) with four z-parts, by (types, |AB|,|CD|,|AC|,|AD|,|BC|,|BD|):")
for k, v in sorted(stats.items(), key=lambda kv: (kv[0][0], -kv[1])): print("  ", k, v)
