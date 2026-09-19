"""Petersen fragments in the minimisers of F(G;5) (experiment D).  Fragments = induced copies of
P-v (9 vertices, 3-pole), P-uv (uv an edge; 8 vertices, 4-pole = dot-product block), P-P3 (7 vertices, 5-pole)."""
import sys, os, json, glob, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from networkx.algorithms import isomorphism as iso
from oddness import petersen
from flowdata import g6_edges

def fragments():
    P = petersen(); out = {}
    for name, rem in (('P-v', [0]), ('P-uv', [0, 1]), ('P-P3', [0, 1, 2])):
        H = P.copy(); H.remove_nodes_from(rem); out[name] = H
    return out
FR = fragments()

def count(G):
    res = {}
    for name, H in FR.items():
        s = set()
        for m in iso.GraphMatcher(G, H).subgraph_isomorphisms_iter():
            s.add(frozenset(m))
            if len(s) >= 500: break
        res[name] = len(s)
    c5 = [c for c in nx.simple_cycles(G, length_bound=5) if len(c) == 5]
    res['C5'] = len(c5); res['v_on_C5'] = len({v for c in c5 for v in c})
    return res

if __name__ == '__main__':
    print('reference: Petersen', count(petersen()))
    for fn in ('flowdata_enum_18-ctf.json', 'flowdata_enum_20-ctf.json', 'flowdata_enum_22-ctf.json'):
        d = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), fn)))
        for c in d['classified']:
            n, E = g6_edges(c['g6']); print(f"exhaustive n={n} F={c['F5']} cyc={c['cyc']} F4={c['F4']}:", count(nx.Graph(E)))
    for fn in sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'flowclimb_best_*.json'))):
        d = json.load(open(fn)); G = nx.Graph([tuple(e) for e in d['edges']])
        if len(G) >= 24: print(f"climber {os.path.basename(fn)} n={len(G)} F={d['F5']}:", count(G))
    rng = random.Random(5)
    for n in (22, 30):
        tot = {}; k = 0
        while k < 40:
            G = nx.random_regular_graph(3, n, seed=rng.randrange(10 ** 9))
            if nx.girth(G) < 5 or nx.edge_connectivity(G) < 3: continue
            k += 1
            for a, b in count(G).items(): tot[a] = tot.get(a, 0) + b / 40
        print(f'random girth>=5 cubic, n={n}, mean over 40:', {a: round(b, 2) for a, b in tot.items()})
