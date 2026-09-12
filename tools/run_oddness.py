"""Decisive experiment for the oddness-6 route: on cyclically 6-edge-connected cubic graphs
with a planted 2-factor F2 having 6 (or 4) odd circuits, does the flow-partition method
(Steffen / Mazzuoccolo-Steffen) find a balanced +-5/3 valuation?
For each graph and F2 we measure, over 0-edge choices c:
   paths   : some 2-colouring differing from a reference only on the t paths is balanced
   allcomp : some 2-colouring of all components of H is balanced
and whether the graph has a nowhere-zero 5-flow at all (SAT)."""
import sys, random, time, itertools
import networkx as nx
from gen_cyc6 import generate, generate_rings, generate_straddle
from oddness import canonical_colouring, H_components, partition, balanced
from mod5 import nz_flow_Z5, triple

profiles = {
    "6x7+8":      [7] * 6 + [8],
    "6x7+2x6":    [7] * 6 + [6, 6],
    "4x7+2x9+6":  [7] * 4 + [9, 9, 6],
    "6x7+14":     [7] * 6 + [14],
    "6x9":        [9] * 6,
    "4x7+8+6":    [7] * 4 + [8, 6],
    "4x9+6+6":    [9] * 4 + [6, 6],
}
ring_profiles = {
    "ring3[7,7,6]":   [[7, 7, 6]] * 3,
    "ring3[7,9,8]":   [[7, 9, 8]] * 3,
    "ring3[7,7,8]":   [[7, 7, 8]] * 3,
    "ring3[7,7,6,6]": [[7, 7, 6, 6]] * 3,
    "ring3[9,9,6]":   [[9, 9, 6]] * 3,
    "ring3[7,7]":     [[7, 7]] * 3,
}
profiles.update(ring_profiles)
straddle_profiles = {   # MS crossing structure: 4 blocks, odd circuits straddling boundaries 0 and 1
    "straddleA": ([[7, 6], [7], [7, 7], [6]], [(7, 3), (7, 3), None, None]),
    "straddleB": ([[7, 8], [7], [7, 7], [8]], [(7, 3), (7, 3), None, None]),
    "straddleC": ([[7, 6], [7, 6], [7, 7], [6]], [(7, 3), (7, 3), None, None]),
    "straddleD": ([[9, 6], [7], [7, 7], [6]], [(7, 3), (7, 3), None, None]),
    "straddleE": ([[7, 6], [7], [7, 7], [6]], [(7, 2), (7, 4), None, None]),
    "straddleF": ([[7, 6], [9], [7, 7], [6]], [(9, 5), (7, 3), None, None]),
    # 6-ring with straddles on alternate boundaries: three distinct crossing bad-type 6-cuts
    "straddle6A": ([[7, 6], [6], [7, 6], [6], [7, 6], [6]], [(7, 3), None, (7, 3), None, (7, 3), None]),
    "straddle6B": ([[7, 6], [8], [7, 6], [6], [7, 8], [6]], [(7, 3), None, (7, 3), None, (7, 3), None]),
    "straddle6C": ([[7], [6], [7], [8], [7], [6]], [(7, 3), None, (7, 3), None, (7, 3), None]),
    "six12": ([[7, 6], [12], [7, 6], [12], [7, 6], [12]], [(7, 3), None, (7, 3), None, (7, 3), None]),
    "six10": ([[7, 8], [10], [7, 8], [10], [7, 8], [10]], [(7, 3), None, (7, 3), None, (7, 3), None]),
    "straddle6D": ([[7, 6], [6], [7, 6], [6], [7, 6], [6]], [(7, 2), None, (7, 4), None, (7, 3), None]),
}
profiles.update(straddle_profiles)
which = sys.argv[1:] or list(profiles)
rng = random.Random(11)
for name in which:
    lengths = profiles[name]
    if isinstance(lengths, tuple):
        graphs = generate_straddle(lengths[0], lengths[1], rng, int(__import__("os").environ.get("NGRAPHS", "8")))
    elif name.startswith("ring"):
        graphs = generate_rings(lengths, rng, 5)
    else:
        graphs = generate(lengths, rng, 6)
    for gi, (G, circuits) in enumerate(graphs):
        n = G.number_of_nodes(); t0 = time.time()
        M = {frozenset((u, v)) for u, v in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
        odd = [C for C in circuits if len(C) % 2]
        E = list(G.edges())
        has5 = nz_flow_Z5(n, E)
        zero_space = list(itertools.product(*[range(len(C)) for C in odd]))
        zero_space = rng.sample(zero_space, min(int(__import__("os").environ.get("NZERO", "120")), len(zero_space)))
        paths_ok = all_ok = 0
        for zc in zero_space:
            col = canonical_colouring(G, M, circuits, zc)
            comps = H_components(G, col)
            path_idx = [i for i, (_, ends) in enumerate(comps) if ends]
            # paths only (other components at a fixed random reference colouring)
            ref = [rng.randint(0, 1) for _ in comps]
            pok = False
            for bits in itertools.product((0, 1), repeat=len(path_idx)):
                f = list(ref)
                for i, b in zip(path_idx, bits): f[i] = b
                if balanced(G, partition(comps, f)): pok = True; break
            paths_ok += pok
            if pok: all_ok += 1
            else:
                space = list(itertools.product((0, 1), repeat=len(comps)))
                if len(space) > 512: space = rng.sample(space, 512)
                all_ok += any(balanced(G, partition(comps, f)) for f in space)
        print(f"{name} #{gi} n={n} |odd|={len(odd)} #H-comps~{len(comps)} 5-flow={has5} | 0-edge choices tried {len(zero_space)}: paths-only OK {paths_ok}, all-components OK {all_ok}  [{time.time()-t0:.0f}s]", flush=True)
