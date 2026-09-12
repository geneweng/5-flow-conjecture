"""killsets for general t: number of killed classes (of 2^(t-1)) per (2-factor, 0-edge) instance."""
import sys, random, itertools, collections, os
from gen_cyc6 import generate_straddle, generate, generate_rings
from oddness import canonical_colouring, H_components, partition
from badcuts import violating_set, classify
def sample_zero_choices(odd, k, rng):
    """k distinct 0-edge choices (one position per odd circuit) without materializing the product
    (7^12 tuples for twelve 7-circuits exhausts memory)."""
    total = 1
    for C in odd: total *= len(C)
    if total <= 4 * k:
        space = list(itertools.product(*[range(len(C)) for C in odd])); return rng.sample(space, min(k, total))
    seen = set()
    while len(seen) < k: seen.add(tuple(rng.randrange(len(C)) for C in odd))
    return sorted(seen)

profiles = {
    "t4str": ([[7, 6], [7], [7, 7], [6], [7, 7]], [(7, 3), (7, 3), None, None, None]),
    "t4str2": ([[7, 6], [7], [7, 7], [6], [7, 7], [8]], [(7, 3), (7, 3), None, None, None, None]),
    "t4rand": [7] * 8 + [8, 6],
    "t4ring": [[7, 7, 6]] * 4,
    "t6rand": [7] * 12 + [8, 6],
    "t7rand": [7] * 14 + [8],
    "t6ring": [[7, 7, 6]] * 6,
    "t7ring": [[7, 7, 6]] * 7,
    "t8ring": [[7, 7, 6]] * 8,
    "t6ring8": [[7, 7, 8]] * 6,
    "t6ring2": [[7, 7]] * 6,
    "t8ring2": [[7, 7]] * 8,
    "t5ring": [[7, 7, 6]] * 5,
    "t6ring4": [[7, 7, 7, 7, 6]] * 3,
    "t8ring4": [[7, 7, 7, 7, 6]] * 4,
    "t5str": ([[7, 6], [7], [7, 7], [6], [7, 7], [7, 7]], [(7, 3), (7, 3), None, None, None, None]),
}
name = sys.argv[1]; seed = int(sys.argv[2]) if len(sys.argv) > 2 else 5; rng = random.Random(seed); ng = int(sys.argv[3]) if len(sys.argv) > 3 else 4
nzero = int(os.environ.get("NZERO", "40"))
pr = profiles[name]
graphs = generate_straddle(pr[0], pr[1], rng, ng) if isinstance(pr, tuple) else (generate_rings(pr, rng, ng) if isinstance(pr[0], list) else generate(pr, rng, ng))
hist = collections.Counter(); types = collections.Counter(); n_inst = 0; worst = 0
worst_sig = None
for gi, (G, circuits) in enumerate(graphs):
    M = {frozenset((u, v)) for u, v in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    odd = [C for C in circuits if len(C) % 2]
    zero_space = sample_zero_choices(odd, nzero, rng)
    for zc in zero_space:
        col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
        path_idx = [i for i, (_, ends) in enumerate(comps) if ends]; t = len(path_idx)
        ref = [rng.randint(0, 1) for _ in comps]; killed = set(); sig = {}
        for bits in itertools.product((0, 1), repeat=t):
            if bits[0] == 1: continue                      # classes up to complement
            f = list(ref)
            for i, b in zip(path_idx, bits): f[i] = b
            black = partition(comps, f); S = violating_set(G, black)
            if S is not None:
                killed.add(bits); ty = classify(G, S, col, comps, black); types[(ty[0], ty[1], ty[2], ty[4], len(ty[5]))] += 1
                sig[bits] = (ty[0], ty[1], ty[2], ty[4], len(ty[5]))
        n_inst += 1; hist[len(killed)] += 1
        if len(killed) > worst:
            worst = len(killed); worst_sig = collections.Counter(sig[b] for b in killed)
            worst_inst = {"profile": name, "seed": seed, "graph": gi, "edges": [list(e) for e in G.edges()],
                          "circuits": circuits, "zero": list(zc), "ref": ref, "killed": sorted(killed), "sig": {str(b): sig[b] for b in killed}}
        if len(killed) == 2 ** (t - 1):                    # every class killed: dump the instance
            import json
            fn = f"fullkill_{name}_{seed}_{gi}_{n_inst}.json"
            json.dump({"profile": name, "seed": seed, "graph": gi, "edges": [list(e) for e in G.edges()],
                       "circuits": circuits, "zero": list(zc), "ref": ref, "killed": sorted(killed)}, open(fn, "w"))
            print("FULL KILL ->", fn, flush=True)
if n_inst == 0: print(f"{name}: generator produced no cyclically 6-connected graph"); sys.exit(0)
print(f"{name} (t={t}, {2**(t-1)} classes): {n_inst} instances; killed classes histogram: {dict(sorted(hist.items()))}; worst {worst}")
print("bad cut types (|dS|, c1, c2, q, #separated paths):", dict(types.most_common(12)))
print("worst instance: killing cut types per killed class:", dict(worst_sig) if worst_sig else None, flush=True)
if worst_sig:
    import json; json.dump(worst_inst, open(f"worst_{name}_{seed}.json", "w"))
