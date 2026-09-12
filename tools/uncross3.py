"""Uncrossing analysis of a dumped three-kill instance (killsets.py -> threekill_*.json):
for the witness sets S of the killed classes report the cut type, the z's inside, the
pairwise part-cut sizes d(U_i, U_j) of the four parts of S, S', and the atom sizes of all cuts."""
import sys, json, itertools, collections
import networkx as nx
from oddness import canonical_colouring, H_components, partition
from badcuts import classify
d = json.load(open(sys.argv[1]))
G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
col = canonical_colouring(G, M, circuits, d["zero"]); comps = H_components(G, col)
path_idx = [i for i, (_, ends) in enumerate(comps) if ends]; ref = d["ref"]
paths = [tuple(comps[i][1]) for i in path_idx]; Z = [z for p in paths for z in p]
print(f"{d['profile']} graph {d['graph']}: n={G.number_of_nodes()}, paths {paths}, path lengths {[len(comps[i][0]) for i in path_idx]}")
def edges_between(A, B): return sum(1 for u, v in G.edges() if (u in A and v in B) or (u in B and v in A))
cuts = []
for bits_s, Sl in d["S"].items():
    bits = eval(bits_s); S = set(Sl); f = list(ref)
    for i, b in zip(path_idx, bits): f[i] = b
    black = partition(comps, f); ty = classify(G, S, col, comps, black)
    lab = (bits[0] ^ bits[1], bits[0] ^ bits[2])
    zin = [(z, "B" if z in black else "w") for z in Z if z in S]
    print(f"  class {lab} (bits {bits}): cut {ty[:5]} sep={ty[5]} |S|={len(S)} z in S: {zin}  path-end colours {[''.join('B' if z in black else 'w' for z in p) for p in paths]}")
    cuts.append((lab, S))
V = set(G.nodes())
for (la, A), (lb, B) in itertools.combinations(cuts, 2):
    parts = {"A&B": A & B, "A-B": A - B, "B-A": B - A, "rest": V - A - B}
    print(f"  pair {la},{lb}: part sizes {[len(p) for p in parts.values()]}; part cuts:",
          {f"{x}|{y}": edges_between(parts[x], parts[y]) for x, y in itertools.combinations(parts, 2)},
          "z's:", {k: [z for z in Z if z in p] for k, p in parts.items()})
if len(cuts) >= 3:
    atoms = collections.Counter()
    for v in V: atoms[tuple(int(v in S) for _, S in cuts)] += 1
    print("  atoms (membership vector -> size):", dict(sorted(atoms.items())))
    print("  z atoms:", {z: tuple(int(z in S) for _, S in cuts) for z in Z})
