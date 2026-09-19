"""Scan random 0-edge choices of archived instances for killers of type 7-cut or (6,4,1,3) (rare types), and
dump every choice that has such a killer together with at least one other killer, as an instance file for
validate_c6.py.  usage: scan_types_c6.py OUTDIR NCHOICES FILE..."""
import sys, json, itertools, os, random, collections
import networkx as nx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oddness import canonical_colouring, H_components, partition
def violator(G, black):
    D = nx.DiGraph()
    for u, v in G.edges(): D.add_edge(u, v, capacity=3); D.add_edge(v, u, capacity=3)
    for v in G.nodes():
        if v in black: D.add_edge('s', v, capacity=5)
        else: D.add_edge(v, 't', capacity=5)
    cut, (X, Y) = nx.minimum_cut(D, 's', 't')
    return None if cut - 5 * len(black) >= 0 else set(X) - {'s'}
out, N = sys.argv[1], int(sys.argv[2]); rng = random.Random(7); tally = collections.Counter(); dumped = 0
for fn in sys.argv[3:]:
    d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
    M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    odd = [C for C in circuits if len(C) % 2]
    for it in range(N):
        zero = [rng.randrange(len(C)) for C in odd]
        col = canonical_colouring(G, M, circuits, zero); comps = H_components(G, col)
        pidx = [i for i, (_, e) in enumerate(comps) if e]
        if len(pidx) != 3: continue
        types = []
        for bits in itertools.product((0,), (0, 1), (0, 1)):
            f = [0] * len(comps)
            for i, b in zip(pidx, bits): f[i] = b
            black = partition(comps, f); S = violator(G, black)
            if S is None: continue
            bnd = [e for e in G.edges() if (e[0] in S) != (e[1] in S)]
            cc = collections.Counter(col[frozenset(e)] for e in bnd); types.append((len(bnd), cc[1], cc[2]))
        for t in types: tally[t] += 1
        if (len(types) >= 2 and any(t in {(7, 5, 2), (6, 4, 1)} for t in types)) or (os.environ.get('SMALL2') == '1' and sum(1 for t in types if t[0] != 11) >= 2):
            dd = dict(d); dd["zero"] = zero; dumped += 1
            json.dump(dd, open(os.path.join(out, f"{os.path.basename(fn)[:-5]}_{it}.json"), "w"))
    print(os.path.basename(fn), dict(tally), "dumped", dumped, flush=True)
