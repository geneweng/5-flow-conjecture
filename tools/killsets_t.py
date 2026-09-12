"""killsets for general t: number of killed classes (of 2^(t-1)) per (2-factor, 0-edge) instance."""
import sys, random, itertools, collections, os
from gen_cyc6 import generate_straddle, generate, generate_rings
from oddness import canonical_colouring, H_components, partition
from badcuts import violating_set, classify
profiles = {
    "t4str": ([[7, 6], [7], [7, 7], [6], [7, 7]], [(7, 3), (7, 3), None, None, None]),
    "t4str2": ([[7, 6], [7], [7, 7], [6], [7, 7], [8]], [(7, 3), (7, 3), None, None, None, None]),
    "t4rand": [7] * 8 + [8, 6],
    "t4ring": [[7, 7, 6]] * 4,
    "t5str": ([[7, 6], [7], [7, 7], [6], [7, 7], [7, 7]], [(7, 3), (7, 3), None, None, None, None]),
}
name = sys.argv[1]; rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 5); ng = int(sys.argv[3]) if len(sys.argv) > 3 else 4
pr = profiles[name]
graphs = generate_straddle(pr[0], pr[1], rng, ng) if isinstance(pr, tuple) else (generate_rings(pr, rng, ng) if isinstance(pr[0], list) else generate(pr, rng, ng))
hist = collections.Counter(); types = collections.Counter(); n_inst = 0; worst = 0
for G, circuits in graphs:
    M = {frozenset((u, v)) for u, v in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    odd = [C for C in circuits if len(C) % 2]
    zero_space = list(itertools.product(*[range(len(C)) for C in odd])); zero_space = rng.sample(zero_space, min(40, len(zero_space)))
    for zc in zero_space:
        col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
        path_idx = [i for i, (_, ends) in enumerate(comps) if ends]; t = len(path_idx)
        ref = [rng.randint(0, 1) for _ in comps]; killed = set()
        for bits in itertools.product((0, 1), repeat=t):
            if bits[0] == 1: continue                      # classes up to complement
            f = list(ref)
            for i, b in zip(path_idx, bits): f[i] = b
            black = partition(comps, f); S = violating_set(G, black)
            if S is not None:
                killed.add(bits); ty = classify(G, S, col, comps, black); types[(ty[0], ty[1], ty[2], ty[4], len(ty[5]))] += 1
        n_inst += 1; hist[len(killed)] += 1; worst = max(worst, len(killed))
print(f"{name} (t={t}, {2**(t-1)} classes): {n_inst} instances; killed classes histogram: {dict(sorted(hist.items()))}; worst {worst}")
print("bad cut types (|dS|, c1, c2, q, #separated paths):", dict(types.most_common(12)))
