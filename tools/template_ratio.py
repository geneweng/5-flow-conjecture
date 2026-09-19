"""For archived template failures: the worst ratio |b-a|/d(S) per path colouring, and the circular flow
number r = 2/(1-rho) that the template achieves at the failing 0-edge choice."""
import sys, json, itertools, glob
from fractions import Fraction
import networkx as nx
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from oddness import canonical_colouring, H_components, partition

def worst(G, black):
    rho = Fraction(0); best = None
    while True:
        p, q = rho.numerator, rho.denominator          # test q(b-a) <= p d
        D = nx.DiGraph()
        for u, v in G.edges():
            D.add_edge(u, v, capacity=p); D.add_edge(v, u, capacity=p)
        for v in G.nodes():
            if v in black: D.add_edge('s', v, capacity=q)
            else: D.add_edge(v, 't', capacity=q)
        cut, (Sx, _) = nx.minimum_cut(D, 's', 't')
        if cut - q * len(black) >= 0: return rho, best
        Sx = set(Sx) - {'s'}
        b = len(Sx & black); a = len(Sx) - b; d = nx.cut_size(G, Sx)
        new = Fraction(b - a, d)
        assert new > rho
        rho, best = new, (b - a, d, len(Sx))

for fn in sys.argv[1:]:
    d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
    M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    col = canonical_colouring(G, M, circuits, d["zero"]); comps = H_components(G, col)
    pidx = [i for i, (_, e) in enumerate(comps) if e]
    out = []
    for bits in itertools.product((0, 1), repeat=len(comps) - 1):
        f = [0] + list(bits)
        rho, w = worst(G, partition(comps, f)); out.append((rho, w))
    rho, w = min(out)
    print(f"{fn.split('/')[-1]:50s} n={G.number_of_nodes()} comps={len(comps)} best rho={rho} (k,d,|S|)={w} r={float(2/(1-rho)):.4f}  all: {sorted(set(str(o[0]) for o in out))}")
