"""Generate cubic graphs with a planted 2-factor of prescribed circuit lengths, girth >= 6,
and certified cyclic 6-edge-connectivity (SAT check: no vertex set X with 6 <= |X| <= n-6
and d(X) <= 5; with girth >= 6 that is exactly cyclic 6-edge-connectivity)."""
import random, itertools
import networkx as nx
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver


def planted_graph(lengths, rng, min_dist=5, tries=200):
    """Circuits of the given lengths + a random perfect matching joining vertices at distance >= min_dist."""
    n = sum(lengths)
    for _ in range(tries):
        G = nx.Graph(); circuits = []; idx = 0
        for L in lengths:
            C = list(range(idx, idx + L)); idx += L
            G.add_edges_from((C[i], C[(i + 1) % L]) for i in range(L)); circuits.append(C)
        free = set(range(n)); ok = True
        while free:
            u = rng.choice(sorted(free))
            dist = nx.single_source_shortest_path_length(G, u, cutoff=min_dist - 1)
            cand = [v for v in free if v != u and v not in dist]
            if not cand: ok = False; break
            v = rng.choice(cand); G.add_edge(u, v); free -= {u, v}
        if ok:
            return G, circuits
    return None


def has_small_cyclic_cut(G, kmax=5, side_min=6):
    """SAT: is there X with side_min <= |X| <= n - side_min and d(X) <= kmax?"""
    nodes = list(G.nodes()); n = len(nodes); vid = {v: i + 1 for i, v in enumerate(nodes)}
    top = n; clauses = []; cutvars = []
    for u, v in G.edges():
        top += 1; y = top; cutvars.append(y)          # y <-> x_u xor x_v
        a, b = vid[u], vid[v]
        clauses += [[-a, b, y], [a, -b, y], [-a, -b, -y], [a, b, -y]]
    xs = [vid[v] for v in nodes]
    enc = CardEnc.atmost(lits=cutvars, bound=kmax, top_id=top, encoding=EncType.seqcounter); clauses += enc.clauses; top = enc.nv
    enc = CardEnc.atleast(lits=xs, bound=side_min, top_id=top, encoding=EncType.seqcounter); clauses += enc.clauses; top = enc.nv
    enc = CardEnc.atmost(lits=xs, bound=n - side_min, top_id=top, encoding=EncType.seqcounter); clauses += enc.clauses; top = enc.nv
    with Solver(name='cadical153', bootstrap_with=clauses) as s:
        if s.solve():
            m = set(l for l in s.get_model() if l > 0)
            return [v for v in nodes if vid[v] in m]
    return None


def cyclically_6_connected(G):
    return nx.girth(G) >= 6 and nx.edge_connectivity(G) >= 3 and has_small_cyclic_cut(G) is None


def cyclically_5_connected(G):
    return nx.girth(G) >= 5 and nx.edge_connectivity(G) >= 3 and has_small_cyclic_cut(G, kmax=4, side_min=5) is None


def generate(lengths, rng, count, verbose=False):
    out = []; attempts = 0
    while len(out) < count and attempts < 100 * count:
        attempts += 1
        r = planted_graph(lengths, rng)
        if r is None: continue
        G, circuits = r
        if cyclically_6_connected(G):
            out.append((G, circuits))
            if verbose: print(f"  found {len(out)} (attempt {attempts})", flush=True)
    return out


if __name__ == "__main__":
    import time
    rng = random.Random(1); t = time.time()
    gs = generate([7] * 6 + [8], rng, 3, verbose=True)
    print(len(gs), "graphs, n =", [g.number_of_nodes() for g, _ in gs], f"[{time.time()-t:.0f}s]")


def _match_free(G, free, rng, min_dist):
    """Backtracking: perfect matching of `free` inside G keeping all distances >= min_dist (girth >= min_dist+1)."""
    if not free:
        return True
    u = min(free)
    dist = nx.single_source_shortest_path_length(G, u, cutoff=min_dist - 1)
    cand = [v for v in free if v != u and v not in dist]; rng.shuffle(cand)
    for v in cand:
        G.add_edge(u, v)
        if _match_free(G, free - {u, v}, rng, min_dist):
            return True
        G.remove_edge(u, v)
    return False


def ring_graph(blocks, rng, min_dist=5, tries=60):
    """blocks: list of lists of circuit lengths.  Blocks are arranged in a ring; consecutive blocks
    are joined by exactly 3 matching edges; the remaining vertices are matched inside their block,
    keeping girth >= 6.  Every block must have an even number of vertices."""
    for _ in range(tries):
        G = nx.Graph(); circuits = []; blockverts = []; idx = 0
        for bl in blocks:
            bv = []
            for L in bl:
                C = list(range(idx, idx + L)); idx += L
                G.add_edges_from((C[i], C[(i + 1) % L]) for i in range(L)); circuits.append(C); bv += C
            blockverts.append(bv)
        m = len(blocks); free = [set(bv) for bv in blockverts]; ok = True
        for j in range(m):
            k = (j + 1) % m
            for _ in range(3):
                us = sorted(free[j]); rng.shuffle(us); placed = False
                for u in us:
                    dist = nx.single_source_shortest_path_length(G, u, cutoff=min_dist - 1)
                    cand = [v for v in free[k] if v not in dist]
                    if cand:
                        v = rng.choice(cand); G.add_edge(u, v); free[j].discard(u); free[k].discard(v); placed = True; break
                if not placed: ok = False; break
            if not ok: break
        if not ok: continue
        for j in range(m):
            if not _match_free(G, frozenset(free[j]), rng, min_dist): ok = False; break
        if ok and all(d == 3 for _, d in G.degree()):
            return G, circuits
    return None


def generate_rings(blocks, rng, count):
    out = []; attempts = 0
    while len(out) < count and attempts < 100 * count:
        attempts += 1
        r = ring_graph(blocks, rng)
        if r is None: continue
        G, circuits = r
        if cyclically_6_connected(G): out.append((G, circuits))
    return out


def ring_graph_straddle(blocks, straddles, rng, min_dist=5, tries=60):
    """Ring of blocks; straddles[j] is None or (L, a): a circuit of length L with a vertices in
    block j and L-a in block j+1 (2 circuit edges cross boundary j).  Boundary j additionally gets
    3 - 2*[straddled] matching links, so every block boundary carries exactly 3 edges and the cut
    around two adjacent blocks is a 6-cut with 4 or 2 matching edges."""
    m = len(blocks)
    for _ in range(tries):
        G = nx.Graph(); circuits = []; free = [set() for _ in range(m)]; idx = 0
        for j, bl in enumerate(blocks):
            for L in bl:
                C = list(range(idx, idx + L)); idx += L
                G.add_edges_from((C[i], C[(i + 1) % L]) for i in range(L)); circuits.append(C); free[j] |= set(C)
        for j, st in enumerate(straddles):
            if st is None: continue
            L, a = st; C = list(range(idx, idx + L)); idx += L
            G.add_edges_from((C[i], C[(i + 1) % L]) for i in range(L)); circuits.append(C)
            free[j] |= set(C[:a]); free[(j + 1) % m] |= set(C[a:])
        ok = True
        for j in range(m):
            k = (j + 1) % m
            for _ in range(3 - (2 if straddles[j] is not None else 0)):
                us = sorted(free[j]); rng.shuffle(us); placed = False
                for u in us:
                    dist = nx.single_source_shortest_path_length(G, u, cutoff=min_dist - 1)
                    cand = [v for v in free[k] if v not in dist]
                    if cand:
                        v = rng.choice(cand); G.add_edge(u, v); free[j].discard(u); free[k].discard(v); placed = True; break
                if not placed: ok = False; break
            if not ok: break
        if not ok: continue
        for j in range(m):
            if len(free[j]) % 2: raise ValueError(f"block {j} has an odd number of free vertices")
            if not _match_free(G, frozenset(free[j]), rng, min_dist): ok = False; break
        if ok and all(d == 3 for _, d in G.degree()):
            return G, circuits
    return None


def generate_straddle(blocks, straddles, rng, count, min_dist=5, check=None):
    check = check or cyclically_6_connected
    out = []; attempts = 0
    while len(out) < count and attempts < 100 * count:
        attempts += 1
        r = ring_graph_straddle(blocks, straddles, rng, min_dist=min_dist)
        if r is None: continue
        if check(r[0]): out.append(r)
    return out
