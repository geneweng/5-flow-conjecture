"""Random 9-regular multigraphs (pairing model), filtered to 8-edge-connected,
tested for a modulo-5 orientation (Jaeger's conjecture, p=2)."""
import random, sys, time
import networkx as nx
from mod5 import has_mod_orientation


def pairing(n, d, rng):
    pts = [v for v in range(n) for _ in range(d)]
    rng.shuffle(pts)
    return [(pts[i], pts[i + 1]) for i in range(0, len(pts), 2)]


def edge_connectivity(n, edges):
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for u, v in edges:
        if u == v:
            continue
        if G.has_edge(u, v):
            G[u][v]['weight'] += 1
        else:
            G.add_edge(u, v, weight=1)
    if not nx.is_connected(G):
        return 0
    return nx.stoer_wagner(G)[0]


if __name__ == "__main__":
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    budget = float(sys.argv[2]) if len(sys.argv) > 2 else 600
    rng = random.Random(seed)
    t0 = time.time(); stats = {}
    while time.time() - t0 < budget:
        n = rng.choice([10, 12, 14, 16, 18, 20, 24])
        E = [(u, v) for u, v in pairing(n, 9, rng) if u != v]
        if len(E) != 9 * n // 2:
            continue                      # had loops; skip
        lam = edge_connectivity(n, E)
        if lam < 8:
            continue
        ok = has_mod_orientation(n, E, 5)
        key = (n, lam)
        stats.setdefault(key, [0, 0])
        stats[key][0] += 1
        stats[key][1] += ok
        if not ok:
            print("NO MOD-5 ORIENTATION:", n, lam, E, flush=True)
    for key in sorted(stats):
        print(f"n={key[0]:2d} lambda={key[1]:2d}: tested {stats[key][0]:5d}, with mod-5 orientation {stats[key][1]:5d}")
