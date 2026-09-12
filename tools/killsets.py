"""Per (2-factor, 0-edge) instance: which of the 4 colouring classes (up to complement) are killed,
and by which cut types.  Class label of a path colouring (b1,b2,b3) is (b1^b2, b1^b3)."""
import sys, random, itertools, collections
from gen_cyc6 import generate_straddle, generate, cyclically_5_connected
import os, json
CYC5 = os.environ.get('CYC5') == '1'
from oddness import canonical_colouring, H_components, partition
from badcuts import violating_set, classify, profiles
rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 5)
name = sys.argv[1]; pr = profiles[name]
graphs = (generate_straddle(pr[0], pr[1], rng, int(sys.argv[3]) if len(sys.argv) > 3 else 6, min_dist=(4 if CYC5 else 5), check=(cyclically_5_connected if CYC5 else None)) if isinstance(pr, tuple) else generate(pr, rng, 6))
killed_hist = collections.Counter(); combo = collections.Counter(); n_inst = 0; col_hist = collections.Counter()
seed = int(sys.argv[2]) if len(sys.argv) > 2 else 5
for gi, (G, circuits) in enumerate(graphs):
    M = {frozenset((u, v)) for u, v in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    odd = [C for C in circuits if len(C) % 2]
    zero_space = rng.sample(list(itertools.product(*[range(len(C)) for C in odd])), 60)
    for zc in zero_space:
        col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
        path_idx = [i for i, (_, ends) in enumerate(comps) if ends]
        if len(path_idx) != 3: continue
        ref = [rng.randint(0, 1) for _ in comps]; killed = {}; sets = {}
        for bits in itertools.product((0, 1), repeat=3):
            f = list(ref)
            for i, b in zip(path_idx, bits): f[i] = b
            black = partition(comps, f); S = violating_set(G, black)
            if S is not None:
                lab = (bits[0] ^ bits[1], bits[0] ^ bits[2])
                t = classify(G, S, col, comps, black)
                killed.setdefault(lab, set()).add((t[0], t[1], t[2], t[5])); sets[str(bits)] = sorted(S)
        n_inst += 1; killed_hist[len(killed)] += 1; col_hist[len(sets)] += 1
        if len(killed) >= 3:                                   # dump for uncrossing analysis
            json.dump({"profile": name, "seed": seed, "graph": gi, "edges": [list(e) for e in G.edges()], "circuits": circuits,
                       "zero": list(zc), "ref": ref, "killed": {str(k): sorted(v) for k, v in killed.items()}, "S": sets},
                      open(f"threekill_{name}_{seed}_{gi}_{n_inst}.json", "w"))
        sig = tuple(sorted((lab, tuple(sorted(ts))) for lab, ts in killed.items()))
        combo[sig] += 1
print(f"{name}: {n_inst} instances; number of killed classes (of 4):", dict(sorted(killed_hist.items())))
print(f"{name}: number of killed path colourings (of 8, even circuits fixed):", dict(sorted(col_hist.items())))
print("most common (killed class -> cut types) signatures:")
for sig, cnt in combo.most_common(12):
    print(f"  {cnt:4d}  " + ("; ".join(f"class{lab}: {list(ts)}" for lab, ts in sig) if sig else "(nothing killed)"))
