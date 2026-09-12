"""Flow partitions (Steffen 2010, Mazzuoccolo-Steffen 2017) and the balanced-valuation
test for nowhere-zero 5-flows on cubic graphs.

Given a cubic graph G, a 2-factor F2 with odd circuits C_1..C_{2t}, and a choice of one
"0-edge" on each odd circuit, the canonical colouring c gives: colour 1 = the perfect
matching M = E - F2; on each circuit colours 2,3 alternate (odd circuits: starting after
the 0-edge).  H = c^{-1}(1) u c^{-1}(2) is a disjoint union of even circuits and t paths
whose ends are the 2t vertices missing colour 2.  A *flow partition* is a proper
2-colouring (A = white, B = black) of H (every component 2-coloured independently);
Jaeger's theorem then says: G has a nowhere-zero 5-flow with this black/white pattern
iff  5*(b_S - a_S) <= 3*|d(S)| for all S, which we test by one min-cut computation.
"""
import random, itertools, sys, time
import networkx as nx


# ---------- graphs ----------
def petersen():
    E = [(i, (i + 1) % 5) for i in range(5)] + [(i, i + 5) for i in range(5)] + [(5 + i, 5 + (i + 2) % 5) for i in range(5)]
    return nx.Graph(E)

def flower(k):
    a = lambda i: ('a', i % k); b = lambda i: ('b', i % k); c = lambda i: ('c', i % k); d = lambda i: ('d', i % k)
    E = []
    for i in range(k):
        E += [(a(i), b(i)), (a(i), c(i)), (a(i), d(i)), (b(i), b(i + 1))]
    for i in range(k - 1):
        E += [(c(i), c(i + 1)), (d(i), d(i + 1))]
    E += [(c(k - 1), d(0)), (d(k - 1), c(0))]
    return nx.convert_node_labels_to_integers(nx.Graph(E))

def gp(n, k):
    E = [(('u', i), ('u', (i + 1) % n)) for i in range(n)] + [(('u', i), ('v', i)) for i in range(n)] + [(('v', i), ('v', (i + k) % n)) for i in range(n)]
    return nx.convert_node_labels_to_integers(nx.Graph(E))

def insert_P2(G, e, tag):
    """Insert the 2-pole P2 (Petersen with a subdivided edge split off) into edge e=(x,y)."""
    G = G.copy(); x, y = e; G.remove_edge(x, y)
    P = petersen(); m = {v: (tag, v) for v in P}
    for u, v in P.edges():
        G.add_edge(m[u], m[v])
    # subdivide edge (0,1) of the Petersen copy by two terminals s (joined to x) and t (joined to y)
    G.remove_edge(m[0], m[1])
    G.add_edge(m[0], x); G.add_edge(m[1], y)
    return G

def R2():
    """Petersen with P2 inserted into the three edges at vertex 0: order 40, oddness 6 (Lukot'ka et al.)."""
    G = petersen()
    for j, (x, y) in enumerate(list(G.edges(0))):
        G = insert_P2(G, (x, y), f'p{j}')
    return nx.convert_node_labels_to_integers(G)


# ---------- 2-factors and colourings ----------
def random_2factor(G, rng):
    """Random perfect matching (max-weight matching with random weights); returns (M, circuits of G - M)."""
    while True:
        for u, v in G.edges():
            G[u][v]['w'] = rng.random()
        M = nx.max_weight_matching(G, maxcardinality=True, weight='w')
        M = {frozenset(e) for e in M}
        if len(M) * 2 != G.number_of_nodes():
            continue
        nbr = {v: [] for v in G.nodes()}
        for u, v in G.edges():
            if frozenset((u, v)) not in M:
                nbr[u].append(v); nbr[v].append(u)
        if any(len(l) != 2 for l in nbr.values()):
            continue
        unseen = set(G.nodes()); circuits = []
        while unseen:
            start = next(iter(unseen)); C = [start]; prev, cur = None, start
            while True:
                nxt = nbr[cur][0] if nbr[cur][0] != prev else nbr[cur][1]
                if nxt == start: break
                C.append(nxt); prev, cur = cur, nxt
            unseen -= set(C); circuits.append(C)
        return M, circuits

def circuit_edges(C):
    return [(C[i], C[(i + 1) % len(C)]) for i in range(len(C))]

def canonical_colouring(G, M, circuits, zero_choice):
    """zero_choice[j] = index of the 0-edge on the j-th odd circuit. Returns colour dict on frozenset edges."""
    col = {e: 1 for e in M}
    odd = [C for C in circuits if len(C) % 2]
    for C in circuits:
        es = circuit_edges(C)
        if len(C) % 2 == 0:
            for i, e in enumerate(es):
                col[frozenset(e)] = 2 if i % 2 == 0 else 3
        else:
            j = odd.index(C); z = zero_choice[j]
            es = es[z:] + es[:z]                     # es[0] is the 0-edge
            col[frozenset(es[0])] = 0
            for i, e in enumerate(es[1:]):
                col[frozenset(e)] = 2 if i % 2 == 0 else 3
    return col

def H_components(G, col):
    H = nx.Graph([tuple(e) for e, c in col.items() if c in (1, 2)])
    H.add_nodes_from(G.nodes())
    comps = []
    for cc in nx.connected_components(H):
        sub = H.subgraph(cc)
        colouring = nx.bipartite.color(sub)          # proper 2-colouring, one canonical
        ends = [v for v in cc if sub.degree(v) == 1]
        comps.append((colouring, ends))              # ends empty for even circuits, two z's for paths
    return comps

def partition(comps, flips):
    black = set()
    for (colouring, ends), f in zip(comps, flips):
        for v, c in colouring.items():
            if c ^ f: black.add(v)
    return black

def balanced(G, black):
    """Is w = +5/3 on black, -5/3 on white a balanced valuation?  Equivalent to
    min_S [3 d(S) - 5(b_S - a_S)] >= 0, computed by one s-t min cut."""
    n = G.number_of_nodes()
    if 2 * len(black) != n:
        return False
    D = nx.DiGraph()
    for u, v in G.edges():
        D.add_edge(u, v, capacity=3); D.add_edge(v, u, capacity=3)
    for v in G.nodes():
        if v in black: D.add_edge('s', v, capacity=5)        # pay 5 if a black vertex is left out of S
        else:          D.add_edge(v, 't', capacity=5)        # pay 5 if a white vertex is put into S
    cut, _ = nx.minimum_cut(D, 's', 't')
    # cut = min_S [ 5|B\S| + 5|A∩S| + 3 d(S) ] = min_S [3 d(S) - 5(b_S - a_S)] + 5|B|
    return cut - 5 * len(black) >= 0


def test_graph(G, rng, want_odd=(4, 6), samples=200, max_zero=400, max_flips=64, verbose=True):
    """Sample 2-factors; for those with the wanted number of odd circuits, try 0-edge choices and
    2-colourings of H; report whether some flow partition is balanced."""
    stats = {}
    seen = set()
    for _ in range(samples):
        M, circuits = random_2factor(G, rng)
        key = frozenset(M)
        if key in seen: continue
        seen.add(key)
        odd = [C for C in circuits if len(C) % 2]
        if len(odd) not in want_odd: continue
        zero_space = list(itertools.product(*[range(len(C)) for C in odd]))
        if len(zero_space) > max_zero: zero_space = rng.sample(zero_space, max_zero)
        f2_ok = False; ok_per_c = 0; tot_c = 0
        for zc in zero_space:
            col = canonical_colouring(G, M, circuits, zc)
            comps = H_components(G, col)
            flip_space = list(itertools.product((0, 1), repeat=len(comps)))
            if len(flip_space) > max_flips: flip_space = rng.sample(flip_space, max_flips)
            c_ok = any(balanced(G, partition(comps, f)) for f in flip_space)
            tot_c += 1; ok_per_c += c_ok; f2_ok |= c_ok
        k = len(odd)
        st = stats.setdefault(k, [0, 0, 0, 0]); st[0] += 1; st[1] += f2_ok; st[2] += ok_per_c; st[3] += tot_c
        if verbose and not f2_ok:
            print(f"  2-factor with {k} odd circuits (lengths {sorted(len(C) for C in odd)}): NO balanced flow partition", flush=True)
    return stats


if __name__ == "__main__":
    rng = random.Random(1)
    which = sys.argv[1:] or ["Petersen", "GP(12,3)", "GP(14,3)", "GP(16,3)", "R2", "J7", "J9"]
    builders = {"Petersen": petersen, "J7": lambda: flower(7), "J9": lambda: flower(9), "J11": lambda: flower(11),
                "GP(12,3)": lambda: gp(12, 3), "GP(14,3)": lambda: gp(14, 3), "GP(16,3)": lambda: gp(16, 3), "GP(18,3)": lambda: gp(18, 3), "GP(13,3)": lambda: gp(13, 3), "GP(17,3)": lambda: gp(17, 3), "GP(15,4)": lambda: gp(15, 4), "R2": R2}
    for name in which:
        G = builders[name](); t = time.time()
        import os
        want = tuple(int(x) for x in os.environ.get('ODD', '2,4,6').split(','))
        st = test_graph(G, rng, want_odd=want, samples=int(os.environ.get('SAMPLES', '300')))
        print(f"{name} (n={G.number_of_nodes()}): " + "; ".join(f"{k} odd: {v[1]}/{v[0]} 2-factors OK, {v[2]}/{v[3]} (2-factor,0-edge) choices OK" for k, v in sorted(st.items())) + f"  [{time.time()-t:.0f}s]", flush=True)
