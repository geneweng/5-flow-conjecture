"""Simulated annealing over cubic graphs of fixed order n, MINIMISING F(G;5) (experiment D).

Moves: edge 2-swaps  ab, cd -> ac, bd  (keeping the graph simple and cubic).  A move is first
evaluated (exact F by tensor contraction, reusing the current contraction plan), then put to the
Metropolis test on log F, and only then checked against the class constraints:
   --cyc 3   3-edge-connected (= 3-connected for cubic graphs)
   --cyc 4|5|6   cyclically c-edge-connected (SAT test of tools/gen_cyc6.py, exact, see flowdata.cyc_conn)
   --girth g     girth >= g
Usage:  python3 flowclimb.py N [N ...] --cyc 5 --girth 5 --minutes 4 --seed 1
Writes tools/flowclimb_best_c{cyc}g{girth}_n{N}.json (best graph found) and prints one summary line per N.
"""
import os, sys, json, math, random, time, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from flowcount import best_plan, plan_cost, _contract, flow_count
from gen_cyc6 import has_small_cyclic_cut

HERE = os.path.dirname(os.path.abspath(__file__))


def short_cycle_through(adj, a, c, g):
    """after adding edge ac: is there a cycle of length < g through it? (BFS from a to c avoiding the edge ac)"""
    if g <= 3: return False
    dist = {a: 0}; frontier = [a]; d = 0
    while frontier and d < g - 2:
        d += 1; nxt = []
        for x in frontier:
            for y in adj[x]:
                if x == a and y == c: continue
                if y == c: return True          # path of length d, cycle length d+1 <= g-1
                if y not in dist: dist[y] = d; nxt.append(y)
        frontier = nxt
    return False


def class_ok(adj, cyc, girth, full_girth=False):
    G = nx.Graph((u, v) for u in adj for v in adj[u] if u < v)
    if full_girth and girth > 3 and nx.girth(G) < girth: return False
    if not nx.is_connected(G) or has_small_cyclic_cut(G, kmax=2, side_min=1) is not None: return False   # 3-edge-connected
    g = max(girth, 3)
    for c in range(3, cyc):
        side = max(c, g)
        if len(G) >= 2 * side and has_small_cyclic_cut(G, kmax=c, side_min=side) is not None: return False
    return True


def random_start(n, cyc, girth, rng):
    for t in range(20000):
        G = nx.random_regular_graph(3, n, seed=rng.randrange(10 ** 9))
        if girth > 3 and nx.girth(G) < girth: continue
        adj = {v: set(G[v]) for v in G}
        if class_ok(adj, cyc, girth): return adj
    return None


def structure(E, n):
    """local structures of a graph: short cycles and Petersen fragments"""
    from networkx.algorithms import isomorphism as iso
    from oddness import petersen
    G = nx.Graph(E)
    cyc = {k: 0 for k in (3, 4, 5, 6)}
    for c in nx.simple_cycles(G, length_bound=6): cyc[len(c)] += 1
    out = dict(girth=nx.girth(G), cycles=cyc, F4=flow_count(E, 4))
    P = petersen()
    Pv = P.copy(); Pv.remove_node(0)                       # Petersen minus a vertex (9 vertices)
    Pe = P.copy(); Pe.remove_nodes_from([0, 1])            # Petersen minus two adjacent vertices (8 vertices, the dot-product block)
    for name, H in (('P-v', Pv), ('P-uv', Pe)):
        sets = set()
        if len(G) <= 60:
            gm = iso.GraphMatcher(G, H)
            for m in gm.subgraph_isomorphisms_iter():     # induced copies
                sets.add(frozenset(m));
                if len(sets) > 200: break
        out[name] = len(sets)
    from flowdata import cyc_conn
    out['cyc'] = cyc_conn(G, 7)
    return out


def climb(n, cyc, girth, minutes, seed, T0=0.25, T1=0.01, start=None, log=print):
    rng = random.Random(seed)
    adj = start or random_start(n, cyc, girth, rng)
    if adj is None:
        log(f'n={n}: no start graph in the class found'); return None
    verts = sorted(adj)
    def edges_of(a): return [(u, v) for u in verts for v in a[u] if u < v]
    E = edges_of(adj); plan, (w, _) = best_plan(E, verts, restarts=20, rng=rng)
    f = _contract(E, verts, plan, 5); lf = math.log(f)
    best = (f, list(E)); t0 = time.time(); budget = minutes * 60.0
    steps = acc = rej_class = 0; since_replan = 0
    while True:
        el = time.time() - t0
        if el > budget: break
        T = T0 * (T1 / T0) ** (el / budget)
        steps += 1
        (a, b), (c, d) = rng.sample(E, 2)
        if rng.random() < 0.5: c, d = d, c
        if len({a, b, c, d}) < 4 or c in adj[a] or d in adj[b]: continue
        # new edges ac, bd
        adj[a].remove(b); adj[b].remove(a); adj[c].remove(d); adj[d].remove(c)
        adj[a].add(c); adj[c].add(a); adj[b].add(d); adj[d].add(b)
        ok = not (short_cycle_through(adj, a, c, girth) or short_cycle_through(adj, b, d, girth))
        if ok:
            E2 = edges_of(adj)
            w2, _ = plan_cost(E2, verts, plan, 5)
            p2 = plan
            if w2 > max(w + 1, 8) or w2 > 10:
                p2, (w2, _) = best_plan(E2, verts, restarts=6, rng=rng)
            if w2 > 12:                         # too wide to evaluate quickly: skip the move
                f2 = 0
            else:
                f2 = _contract(E2, verts, p2, 5)
            ok = f2 > 0
            if ok:
                lf2 = math.log(f2)
                ok = lf2 <= lf or rng.random() < math.exp(-(lf2 - lf) / T)
            if ok and not class_ok(adj, cyc, girth):
                ok = False; rej_class += 1
        if ok:
            acc += 1; E, plan, w, f, lf = E2, p2, w2, f2, lf2; since_replan += 1
            if since_replan >= 40:
                p3, (w3, c3) = best_plan(E, verts, restarts=4, rng=rng)
                if (w3, c3) < plan_cost(E, verts, plan, 5): plan, w = p3, w3
                else: w = plan_cost(E, verts, plan, 5)[0]
                since_replan = 0
            if f < best[0]: best = (f, list(E))
        else:
            adj[a].remove(c); adj[c].remove(a); adj[b].remove(d); adj[d].remove(b)
            adj[a].add(b); adj[b].add(a); adj[c].add(d); adj[d].add(c)
    fb, Eb = best
    assert flow_count(Eb, 5, restarts=10) == fb                      # recount with a fresh plan
    adjb = {v: set() for v in verts}
    for u, v in Eb: adjb[u].add(v); adjb[v].add(u)
    assert class_ok(adjb, cyc, girth, full_girth=True)
    st = structure(Eb, n)
    res = dict(n=n, cyc_req=cyc, girth_req=girth, seed=seed, minutes=minutes, F5=fb, g=fb ** (1.0 / n), steps=steps, accepted=acc,
               rejected_by_class=rej_class, structure=st, edges=Eb)
    log(f"n={n} class(cyc>={cyc},girth>={girth}) seed={seed}: min F5 = {fb}  F^(1/n)={fb ** (1 / n):.4f}  F/1.5^n={fb / 1.5 ** n:.3f}  "
        f"F mod 5={fb % 5}  [steps {steps}, accepted {acc}, class-rejected {rej_class}]  girth={st['girth']} cyc={st['cyc']} "
        f"cycles(3,4,5,6)={tuple(st['cycles'][k] for k in (3, 4, 5, 6))} F4={st['F4']} P-v copies={st['P-v']} P-uv copies={st['P-uv']}")
    return res


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('n', type=int, nargs='+'); ap.add_argument('--cyc', type=int, default=3); ap.add_argument('--girth', type=int, default=3)
    ap.add_argument('--minutes', type=float, default=4); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--start', default=None, help='json file with field "edges" (start graph)')
    ap.add_argument('--T0', type=float, default=0.25); ap.add_argument('--T1', type=float, default=0.01)
    a = ap.parse_args()
    for n in a.n:
        start = None
        if a.start:
            G = nx.convert_node_labels_to_integers(nx.Graph([tuple(e) for e in json.load(open(a.start))['edges']]))
            assert len(G) == n; start = {v: set(G[v]) for v in G}
        r = climb(n, a.cyc, a.girth, a.minutes, a.seed, a.T0, a.T1, start=start, log=lambda s: print(s, flush=True))
        if r is None: continue
        fn = os.path.join(HERE, f'flowclimb_best_c{a.cyc}g{a.girth}_n{n}.json')
        old = json.load(open(fn)) if os.path.exists(fn) else None
        if old is None or old['F5'] > r['F5']: json.dump(r, open(fn, 'w'))
