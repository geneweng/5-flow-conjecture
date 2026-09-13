import sys, json, random
import networkx as nx
from gridsample_lib import rows_for
fn, N, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); rng = random.Random(seed)
d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
odd = [C for C in circuits if len(C) % 2]
ok = mism = 0
for _ in range(N):
    base = [rng.randrange(len(C)) for C in odd]
    for i in range(len(odd)):
        for sgn in (1, -1):
            r = rows_for(G, M, circuits, odd, base, i, sgn)
            if r is None: continue
            rows, zc = r
            r2 = rows_for(G, M, circuits, odd, zc, i, -sgn)
            if r2 is None: print("reverse not re-pairing?!", base, i, sgn); continue
            rows2, back = r2
            assert back == base
            if (rows["a0"], rows["E"], rows["E1"]) == (rows2["E1"], rows2["E"], rows2["a0"]): ok += 1
            else: mism += 1; print("MISMATCH", base, i, sgn, rows, rows2)
print("ok", ok, "mismatch", mism)
