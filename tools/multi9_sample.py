"""Random loopless 9-regular multigraphs on n vertices with heavy multiplicities
(random symmetric multiplicity matrices via a Markov chain of 'switches'),
filtered to edge-connectivity >= 8, tested for modulo-5 orientations."""
import sys, random, time, itertools
import networkx as nx
from mod5 import has_mod_orientation

n = int(sys.argv[1]); seed = int(sys.argv[2]); budget = float(sys.argv[3])
rng = random.Random(seed)
pairs = list(itertools.combinations(range(n), 2))

def start():
    # 3G for a random cubic multigraph pairing (may have multiplicity up to 3 per pair), then perturb
    while True:
        pts = [v for v in range(n) for _ in range(3)]
        rng.shuffle(pts)
        m = {}
        ok = True
        for i in range(0, len(pts), 2):
            u, v = pts[i], pts[i + 1]
            if u == v:
                ok = False; break
            key = (min(u, v), max(u, v))
            m[key] = m.get(key, 0) + 3
        if ok:
            return m

def switch(m):
    # pick two disjoint pairs (a,b),(c,d) with m>0 and move one unit to (a,c),(b,d)  (degree-preserving)
    es = [p for p, k in m.items() if k > 0]
    for _ in range(50):
        (a, b), (c, d) = rng.sample(es, 2)
        if len({a, b, c, d}) < 4:
            continue
        if rng.random() < 0.5:
            c, d = d, c
        for x, y in ((a, b), (c, d)):
            p = (min(x, y), max(x, y)); m[p] -= 1
        for x, y in ((a, c), (b, d)):
            p = (min(x, y), max(x, y)); m[p] = m.get(p, 0) + 1
        return

def edge_conn(m):
    G = nx.Graph(); G.add_nodes_from(range(n))
    for (u, v), k in m.items():
        if k: G.add_edge(u, v, weight=k)
    return nx.stoer_wagner(G)[0] if nx.is_connected(G) else 0

t0 = time.time(); tested = good = 0; maxmult = 0
m = start()
while time.time() - t0 < budget:
    for _ in range(rng.randint(1, 8)):
        switch(m)
    if edge_conn(m) < 8:
        continue
    E = [p for p, k in m.items() for _ in range(k)]
    tested += 1
    mm = max(m.values()); maxmult = max(maxmult, mm)
    if has_mod_orientation(n, E, 5):
        good += 1
    else:
        print("NO MOD-5 ORIENTATION:", {p: k for p, k in m.items() if k}, flush=True)
print(f"n={n} seed={seed}: tested {tested} 8-edge-connected 9-regular multigraphs, with mod-5 orientation {good}; max multiplicity seen {maxmult}")
