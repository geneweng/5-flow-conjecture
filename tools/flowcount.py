"""Exact number F(G;q) of nowhere-zero Z_q-flows of a (multi)graph, by dynamic programming
over a vertex order with small edge frontier (experiment D of notes/oddness.md, 18.4).

Method.  F(G;q) is a tensor-network contraction: fix an orientation; every vertex v carries the
0/1 tensor K_v[values on its edges] = [inflow = outflow mod q] with one axis of length q-1 (the
values 1..q-1) per incident edge, and F is the full contraction.  Clusters of vertices are merged
pairwise along a contraction tree ("plan"); the tensor of a cluster S has one axis per edge of the
cut d(S), so the cost is governed by max d(S) over the tree (carving-width style; a vertex order
with small frontier is the special case of a caterpillar tree).  Plans are found by a randomised
greedy heuristic with restarts; a plan stays valid when the graph is changed on the same vertex set
(used by the climber).

Arithmetic.  float64 with BLAS is exact as long as every entry stays below 2^52 (all numbers are
non-negative integers, so partial sums are bounded by the entry); this is checked after every
contraction, and otherwise the computation is redone in int64 modulo primes near 2^20 and
recombined by CRT (rigorous: F <= (q-1)^(m-n+1)).

Loops contribute a factor q-1; parallel edges are allowed (edges are given as a list).

CLI:   python3 flowcount.py --test          validation suite
       python3 flowcount.py --json FILE     F(G;5) of the graph in FILE (field "edges")
"""
import os
for _v in ('OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sys, random, itertools, math, json, time
import numpy as np

_KERNELS = {}


def kernel(k, f, q):
    key = (k, f, q)
    if key not in _KERNELS:
        d = k + f
        K = np.zeros((q - 1,) * d, dtype=np.float64)
        for t in itertools.product(range(1, q), repeat=d):
            if (sum(t[:k]) - sum(t[k:])) % q == 0:
                K[tuple(x - 1 for x in t)] = 1.0
        _KERNELS[key] = K
    return _KERNELS[key]


# ---------------------------------------------------------------- contraction plans
# A plan is a list of pairs (i, j) of cluster ids; clusters 0..n-1 are the vertices (in the order of
# `verts`), the cluster created by the t-th contraction has id n+t.  The tensor of a cluster S has one
# axis per edge of the cut d(S), so the width of a plan is max d(S) (a carving-width style parameter).
def _verts_adj(edges):
    verts = []; seen = set()
    for u, v in edges:
        for x in (u, v):
            if x not in seen: seen.add(x); verts.append(x)
    return verts


def plan_cost(edges, verts, plan, q=5):
    """(width, flops) of a plan on this graph (the plan may come from another graph on the same vertices)."""
    idx = {v: i for i, v in enumerate(verts)}
    nb = [dict() for _ in verts]
    for u, v in edges:
        if u == v: continue
        a, b = idx[u], idx[v]
        nb[a][b] = nb[a].get(b, 0) + 1; nb[b][a] = nb[b].get(a, 0) + 1
    deg = [sum(d.values()) for d in nb]
    width = max(deg) if deg else 0; flops = 0.0
    for (a, b) in plan:
        sh = nb[a].get(b, 0)
        flops += float(q - 1) ** (deg[a] + deg[b] - sh)
        new = {}
        for x in (a, b):
            for y, c in nb[x].items():
                if y != a and y != b: new[y] = new.get(y, 0) + c
        c = len(nb); nb.append(new); deg.append(deg[a] + deg[b] - 2 * sh)
        for y, cnt in new.items():
            nb[y].pop(a, None); nb[y].pop(b, None); nb[y][c] = cnt
        nb[a] = nb[b] = None
        width = max(width, deg[c])
    return width, flops


def greedy_plan(edges, verts, rng, alpha=1.0, tau=0.3, q=5):
    idx = {v: i for i, v in enumerate(verts)}
    nb = {i: {} for i in range(len(verts))}
    for u, v in edges:
        if u == v: continue
        a, b = idx[u], idx[v]
        nb[a][b] = nb[a].get(b, 0) + 1; nb[b][a] = nb[b].get(a, 0) + 1
    deg = {i: sum(d.values()) for i, d in nb.items()}
    plan = []; nxt = len(verts); lq = math.log(q - 1)
    while len(nb) > 1:
        best = None
        for a, d in nb.items():
            for b, sh in d.items():
                if a < b:
                    dab = deg[a] + deg[b] - 2 * sh
                    sc = dab - alpha * (max(deg[a], deg[b]) + math.log1p(math.exp(-abs(deg[a] - deg[b]) * lq)) / lq)
                    if tau: sc += tau * rng.gauss(0, 1)
                    if best is None or sc < best[0]: best = (sc, a, b)
        if best is None:                       # disconnected: join two arbitrary clusters
            ks = sorted(nb, key=lambda x: deg[x]); a, b = ks[0], ks[1]
        else:
            _, a, b = best
        sh = nb[a].get(b, 0); new = {}
        for x in (a, b):
            for y, c in nb[x].items():
                if y != a and y != b: new[y] = new.get(y, 0) + c
        del nb[a], nb[b]
        for y, cnt in new.items():
            nb[y].pop(a, None); nb[y].pop(b, None); nb[y][nxt] = cnt
        nb[nxt] = new; deg[nxt] = deg[a] + deg[b] - 2 * sh
        plan.append((a, b)); nxt += 1
    return plan


def best_plan(edges, verts=None, restarts=40, rng=None, q=5):
    rng = rng or random.Random(0)
    verts = verts or _verts_adj(edges)
    best = None; bp = None
    for r in range(restarts):
        alpha = 1.0 if r == 0 else rng.choice([0.0, 0.5, 1.0, 1.0, 1.5])
        tau = 0.0 if r == 0 else rng.choice([0.1, 0.3, 0.6, 1.0])
        p = greedy_plan(edges, verts, rng, alpha, tau, q)
        c = plan_cost(edges, verts, p, q)
        if best is None or c < best: best, bp = c, p
    return bp, best


# ---------------------------------------------------------------- the contraction
class Overflow(Exception):
    pass


def _small_primes(lo=2 ** 20 - 3000, hi=2 ** 20):
    return [p for p in range(hi - 1, lo, -1) if all(p % d for d in range(2, int(p ** 0.5) + 1))]


def _contract(edges, verts, plan, q, mod=None):
    idx = {v: i for i, v in enumerate(verts)}
    inc = [([], []) for _ in verts]            # per vertex: (incoming edge ids, outgoing edge ids)
    loops = 0
    for i, (u, v) in enumerate(edges):
        if u == v: loops += 1; continue
        inc[idx[u]][1].append(i); inc[idx[v]][0].append(i)
    cl = []
    dt = np.float64 if mod is None else np.int64
    for ins, outs in inc:
        cl.append((kernel(len(ins), len(outs), q).astype(dt), ins + outs))
    for (a, b) in plan:
        TA, xa = cl[a]; TB, xb = cl[b]
        sh = [e for e in xa if e in xb]
        if mod is not None: assert len(sh) <= 11
        T = np.tensordot(TA, TB, axes=([xa.index(e) for e in sh], [xb.index(e) for e in sh]))
        if mod is None:
            if T.size and T.max() >= 2.0 ** 52: raise Overflow()
        else:
            T %= mod
        cl.append((T, [e for e in xa if e not in sh] + [e for e in xb if e not in sh]))
        cl[a] = cl[b] = None
    T, x = cl[-1]
    assert not x
    val = int(round(float(T))) if mod is None else int(T)
    if mod is None: return val * (q - 1) ** loops
    return (val * pow(q - 1, loops, mod)) % mod


def flow_count(edges, q=5, plan=None, verts=None, restarts=40, rng=None, return_info=False, max_width=14):
    """Exact F(G;q).  edges: list of pairs (parallel edges and loops allowed)."""
    edges = [tuple(e) for e in edges]
    verts = verts or _verts_adj(edges)
    if not [e for e in edges if e[0] != e[1]]:
        val = (q - 1) ** len(edges)
        return (val, 0, []) if return_info else val
    if plan is None:
        plan, _ = best_plan(edges, verts, restarts, rng, q)
    width, flops = plan_cost(edges, verts, plan, q)
    if width > max_width: raise MemoryError(f'plan width {width} > {max_width}')
    try:
        val = _contract(edges, verts, plan, q)
    except Overflow:
        n = len(verts); m = len(edges)
        bound = (q - 1) ** (max(m - n + 1, 1) + 2)
        res = []; M = 1
        for p in _small_primes():
            res.append((_contract(edges, verts, plan, q, mod=p), p)); M *= p
            if M > bound: break
        assert M > bound, 'not enough primes'
        val = 0
        for r, p in res:
            Mi = M // p; val += r * Mi * pow(Mi, -1, p)
        val %= M
    if return_info:
        return val, width, plan
    return val


def F(G, q=5, **kw):
    """networkx graph -> F(G;q)"""
    return flow_count(list(G.edges()), q, **kw)


# ---------------------------------------------------------------- independent checks
def flow_count_dc(edges, q=5):
    """Deletion-contraction F(G) = F(G/e) - F(G-e) (e not a loop), F(loop+G) = (q-1) F(G). Small graphs only."""
    edges = [tuple(e) for e in edges]
    if not edges: return 1
    for i, (u, v) in enumerate(edges):
        if u == v:
            return (q - 1) * flow_count_dc(edges[:i] + edges[i + 1:], q)
    # bridge/pendant shortcut: a vertex of degree 1 gives 0
    deg = {}
    for u, v in edges: deg[u] = deg.get(u, 0) + 1; deg[v] = deg.get(v, 0) + 1
    if min(deg.values()) == 1: return 0
    (u, v) = edges[0]; rest = edges[1:]
    contracted = [(u if a == v else a, u if b == v else b) for a, b in rest]
    return flow_count_dc(contracted, q) - flow_count_dc(rest, q)


def flow_count_brute(edges, q=5):
    """Enumerate values on the co-tree edges of a spanning forest, solve the tree edges."""
    import networkx as nx
    edges = [tuple(e) for e in edges]
    G = nx.MultiGraph(); G.add_edges_from((u, v, i) for i, (u, v) in enumerate(edges))
    tree = set()
    Tg = nx.Graph()
    for u, v, i in G.edges(keys=True):
        if u != v and not (Tg.has_node(u) and Tg.has_node(v) and nx.has_path(Tg, u, v)):
            Tg.add_edge(u, v, key=i); tree.add(i)
    cot = [i for i in range(len(edges)) if i not in tree]
    # leaf-elimination order of the forest
    Tc = Tg.copy(); elim = []
    while Tc.number_of_edges():
        for x in [x for x in Tc if Tc.degree(x) == 1]:
            if Tc.degree(x) != 1: continue
            y = next(iter(Tc[x])); elim.append((x, Tc[x][y]['key'])); Tc.remove_node(x)
    count = 0
    for vals in itertools.product(range(1, q), repeat=len(cot)):
        exc = {}
        for i, a in zip(cot, vals):
            u, v = edges[i]
            if u == v: continue
            exc[u] = exc.get(u, 0) - a; exc[v] = exc.get(v, 0) + a      # a flows u -> v
        ok = True
        for x, i in elim:
            u, v = edges[i]; y = v if u == x else u
            e = exc.get(x, 0) % q                                      # surplus at x must leave along the tree edge
            if e == 0: ok = False; break
            exc[y] = exc.get(y, 0) + e; exc[x] = 0
        if ok and all(e % q == 0 for e in exc.values()): count += 1
    return count


def heuristic(n, m, q=5):
    return (q - 1) ** m / q ** (n - 1)


def selftest():
    import networkx as nx
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from oddness import petersen, flower, gp
    P = petersen()
    for q, want in ((3, 0), (4, 0), (5, 240), (6, None)):
        got = F(P, q); dc = flow_count_dc(list(P.edges()), q)
        print(f'Petersen q={q}: DP {got}  DC {dc}  want {want}'); assert got == dc and (want is None or got == want)
    assert flow_count_brute(list(P.edges()), 5) == 240
    K4 = nx.complete_graph(4)
    for q in (3, 4, 5, 6, 7):
        assert F(K4, q) == (q - 1) * (q - 2) * (q - 3), q
    print('K4 ok: (q-1)(q-2)(q-3)')
    named = {'K33': nx.complete_bipartite_graph(3, 3), 'prism': nx.circular_ladder_graph(3), 'cube': nx.hypercube_graph(3),
             'J3': flower(3), 'gp(7,2)': gp(7, 2), 'heawood': nx.heawood_graph()}
    for name, G in named.items():
        E = list(G.edges())
        for q in (3, 4, 5, 6):
            a = F(G, q); b = flow_count_dc(E, q) if len(E) <= 18 else flow_count_brute(E, q) if q <= 4 else a
            assert a == b, (name, q, a, b)
        print(f'{name}: F3,F4,F5,F6 =', [F(G, q) for q in (3, 4, 5, 6)])
    # K33 flow polynomial (q-1)(q-2)(q^2-6q+10)
    for q in (3, 4, 5, 6): assert F(named['K33'], q) == (q - 1) * (q - 2) * (q * q - 6 * q + 10)
    print('K33 polynomial ok')
    rng = random.Random(1)
    for t in range(40):
        n = rng.choice([6, 8, 10]); G = nx.random_regular_graph(3, n, seed=rng.randrange(10 ** 6))
        E = list(G.edges()); q = rng.choice([3, 4, 5, 6])
        a = F(G, q, rng=rng); b = flow_count_dc(E, q); assert a == b, (E, q, a, b)
        o = list(G.nodes()); rng.shuffle(o); assert flow_count(E, q, verts=o, plan=[(0, 1)] + [(len(o) + i, i + 2) for i in range(len(o) - 2)]) == a
    print('40 random cubic graphs: DP == deletion-contraction, order-independent')
    # multigraphs with loops / parallel edges
    for t in range(30):
        n = rng.randint(2, 5); m = rng.randint(3, 9)
        E = [(rng.randrange(n), rng.randrange(n)) for _ in range(m)]; q = rng.choice([3, 4, 5])
        assert flow_count(E, q) == flow_count_dc(E, q) == flow_count_brute(E, q), E
    print('30 random multigraphs: DP == DC == brute force')
    # CRT path (int64 modulo primes near 2^20) against the float path
    G = nx.random_regular_graph(3, 40, seed=5); E = list(G.edges()); V = _verts_adj(E); pl, _ = best_plan(E, V)
    a = _contract(E, V, pl, 5); M = 1; res = []
    for p in _small_primes()[:3]: res.append((_contract(E, V, pl, 5, mod=p), p)); M *= p
    v = sum(r * (M // p) * pow(M // p, -1, p) for r, p in res) % M
    print('CRT check n=40:', a, v); assert a == v
    for n in (60, 80):
        t = time.time(); G = nx.random_regular_graph(3, n, seed=7)
        val, w, _ = flow_count(list(G.edges()), 5, return_info=True, restarts=100)
        print(f'random cubic n={n}: F5 = {val} (width {w}, {time.time() - t:.1f}s), heuristic {heuristic(n, 3 * n // 2):.4e}')
    print('ALL TESTS PASSED')


if __name__ == '__main__':
    if '--test' in sys.argv:
        selftest()
    elif '--json' in sys.argv:
        d = json.load(open(sys.argv[sys.argv.index('--json') + 1])); E = [tuple(e) for e in d['edges']]
        t = time.time(); val, w, _ = flow_count(E, 5, return_info=True, restarts=200)
        n = len({x for e in E for x in e})
        print(f'n={n} m={len(E)} width={w} F5={val} F^(1/n)={val ** (1 / n):.5f} ratio={val / heuristic(n, len(E)):.4f} ({time.time() - t:.1f}s)')
