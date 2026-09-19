"""Data for experiment D (notes/exp-D-counting.md): F(G;5) on named families, random cubic graphs, the
project's cyclically 6-connected graphs, and exhaustive enumeration with geng.

  python3 flowdata.py families
  python3 flowdata.py random
  python3 flowdata.py enum N [geng flags, default -c]     e.g.  enum 16 ;  enum 22 -ctf
"""
import os, sys, json, glob, random, subprocess, time, math, heapq
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from flowcount import flow_count, best_plan, plan_cost, _contract, _verts_adj, heuristic
from oddness import petersen, flower, gp, R2

HERE = os.path.dirname(os.path.abspath(__file__))


# ------------------------------------------------------------ constructions
def dot_product(G1, G2, e1=None, e2=None, f=None):
    """Isaacs dot product: remove two independent edges ab, cd of G1, an edge xy of G2 with its ends;
    join a,b to the other neighbours of x and c,d to those of y."""
    G1 = nx.relabel_nodes(G1, {v: (1, v) for v in G1}); G2 = nx.relabel_nodes(G2, {v: (2, v) for v in G2})
    E1 = list(G1.edges())
    if e1 is None:
        e1 = E1[0]; e2 = next(e for e in E1 if not set(e) & set(e1))
    else:
        e1 = ((1, e1[0]), (1, e1[1])); e2 = ((1, e2[0]), (1, e2[1]))
    f = list(G2.edges())[0] if f is None else ((2, f[0]), (2, f[1]))
    x, y = f
    nx_ = [u for u in G2[x] if u != y]; ny = [u for u in G2[y] if u != x]
    G = nx.union(G1, G2); G.remove_edge(*e1); G.remove_edge(*e2); G.remove_nodes_from([x, y])
    G.add_edges_from([(e1[0], nx_[0]), (e1[1], nx_[1]), (e2[0], ny[0]), (e2[1], ny[1])])
    assert all(d == 3 for _, d in G.degree())
    return nx.convert_node_labels_to_integers(G)


def blanusa_chain(k):
    """P . P . ... . P (k factors): a Blanusa-type snark of order 8k+2."""
    G = petersen()
    for _ in range(k - 1): G = dot_product(G, petersen())
    return G


def girth(G):
    return nx.girth(G)


def cyc_conn(G, upto=6):
    """cyclic edge-connectivity of a cubic graph, capped at `upto` (SAT; both sides must contain a cycle:
    for a cut of size c in a cubic graph a connected side has a cycle iff it has >= c vertices; we use the
    sufficient formulation 'side with >= max(c, girth) vertices' on minimal cuts)."""
    from gen_cyc6 import has_small_cyclic_cut
    if nx.edge_connectivity(G) < 3: return nx.edge_connectivity(G)
    g = nx.girth(G)
    for c in range(3, upto):
        # a cycle-separating c-cut: both sides contain a cycle => both sides have >= g vertices, and a side
        # with s vertices and cut c has (3s-c)/2 edges >= s iff s >= c.  Minimal cuts have connected sides.
        side = max(c, g)
        if len(G) >= 2 * side and has_small_cyclic_cut(G, kmax=c, side_min=side) is not None:
            return c
    return upto


def row(name, G, q4=True, restarts=60):
    E = list(G.edges()); n = len(G); m = len(E)
    t = time.time(); f5, w, _ = flow_count(E, 5, return_info=True, restarts=restarts, rng=random.Random(1))
    f4 = flow_count(E, 4, restarts=restarts, rng=random.Random(1)) if q4 else None
    return dict(name=name, n=n, F5=f5, F4=f4, g=f5 ** (1.0 / n) if f5 else 0, r_naive=f5 / heuristic(n, m),
                r_bethe=f5 / 1.5 ** n, mod5=f5 % 5, F_over_4_mod=[(f5 // 4) % k for k in (2, 3, 5)], girth=nx.girth(G),
                width=w, secs=round(time.time() - t, 2))


def show(r):
    print(f"{r['name']:<28} n={r['n']:<3} F5={r['F5']:<22} F4={str(r['F4']):<14} F^(1/n)={r['g']:.4f} "
          f"F/1.5^n={r['r_bethe']:<9.4g} F/(4^m/5^(n-1))={r['r_naive']:<9.3g} F mod 5={r['mod5']} girth={r['girth']} w={r['width']} "
          + (f"cyc={r['cyc']}" if 'cyc' in r else ''), flush=True)


def families():
    out = []
    def add(name, G, cyc=True):
        try: r = row(name, G)
        except MemoryError as ex:
            print(f'{name}: skipped ({ex})', flush=True); return
        if cyc and len(G) <= 80:
            try: r['cyc'] = cyc_conn(G, 7 if nx.girth(G) >= 6 else 6)
            except Exception as ex: r['cyc'] = '?'
        out.append(r); show(r)
    add('K4', nx.complete_graph(4)); add('K33', nx.complete_bipartite_graph(3, 3)); add('prism', nx.circular_ladder_graph(3))
    add('cube', nx.hypercube_graph(3)); add('Petersen', petersen()); add('Heawood', nx.heawood_graph())
    add('dodecahedron', nx.dodecahedral_graph()); add('Tutte-Coxeter', nx.LCF_graph(30, [-13, -9, 7, -7, 9, 13], 5))
    add('McGee (7-cage)', nx.LCF_graph(24, [12, 7, -7], 8)); add('Desargues', nx.desargues_graph())
    add('Coxeter', nx.Graph([(f'a{i}', f'z{i}') for i in range(7)] + [(f'b{i}', f'z{i}') for i in range(7)] + [(f'c{i}', f'z{i}') for i in range(7)] + [(f'a{i}', f'a{(i + 1) % 7}') for i in range(7)] + [(f'b{i}', f'b{(i + 2) % 7}') for i in range(7)] + [(f'c{i}', f'c{(i + 3) % 7}') for i in range(7)]))
    for k in range(3, 27, 2): add(f'flower J{k}', flower(k))
    for k in range(1, 8): add(f'Blanusa chain P^.{k}', blanusa_chain(k))
    add('R2 (oddness 6, n=40)', R2())
    for (n, k) in [(5, 2), (6, 2), (7, 2), (8, 3), (9, 2), (10, 3), (11, 2), (12, 5), (13, 5), (14, 5), (15, 4), (16, 6), (17, 4),
                   (18, 6), (20, 6), (21, 7), (24, 6), (28, 7), (30, 6), (35, 7), (36, 6), (42, 6), (42, 7)]:
        add(f'gp({n},{k})', gp(n, k))
    for k in (4, 6, 8, 10, 15, 20, 30): add(f'prism C{k}xK2', nx.circular_ladder_graph(k))
    for k in (4, 6, 8, 10, 15, 20, 30): add(f'Moebius ladder M{2 * k}', nx.circulant_graph(2 * k, [1, k]))
    for f in sorted(glob.glob(os.path.join(HERE, '..', 'data', '*.json'))):
        d = json.load(open(f)); G = nx.Graph([tuple(e) for e in d['edges']])
        add('data/' + os.path.basename(f).replace('oddness6_', '').replace('.json', ''), G)
    fk = sorted(glob.glob(os.path.join(HERE, 'fullkill_*.json')))
    for f in fk:
        try:
            d = json.load(open(f)); G = nx.Graph([tuple(e) for e in d['edges']])
            if len(G) <= 80: add('tools/' + os.path.basename(f).replace('.json', ''), G, cyc=False)
        except Exception as ex:
            print('skip', f, ex)
    json.dump(out, open(os.path.join(HERE, 'flowdata_families.json'), 'w'), indent=0)


def random_graphs():
    out = []
    for n in (20, 24, 30, 36, 42, 50, 56, 60):
        rows = []; seed = 0
        while len(rows) < (30 if n <= 50 else 12):
            seed += 1
            G = nx.random_regular_graph(3, n, seed=1000 * n + seed)
            if nx.edge_connectivity(G) < 3: continue
            r = row(f'random n={n} seed={seed}', G, q4=False); r['tri'] = sum(nx.triangles(G).values()) // 3
            rows.append(r); out.append(r)
        gs = [r['g'] for r in rows]; rb = [r['r_bethe'] for r in rows]
        lo = min(rows, key=lambda r: r['F5'])
        print(f"n={n}: {len(rows)} 3-connected random cubic graphs: F^(1/n) min {min(gs):.4f} mean {sum(gs) / len(gs):.4f} max {max(gs):.4f};"
              f" F/1.5^n min {min(rb):.3f} mean {sum(rb) / len(rb):.3f} max {max(rb):.3f}; min F={lo['F5']} (girth {lo['girth']}, triangles {lo['tri']})"
              f"; mean log(F)/n by girth: " + ', '.join(f"g{k}:{sum(math.log(r['F5']) / n for r in rows if r['girth'] == k) / max(1, sum(1 for r in rows if r['girth'] == k)):.4f}" for k in sorted({r['girth'] for r in rows})), flush=True)
    json.dump(out, open(os.path.join(HERE, 'flowdata_random.json'), 'w'), indent=0)


# ------------------------------------------------------------ exhaustive enumeration
def g6_edges(s):
    s = s.strip(); n = ord(s[0]) - 63; bits = []
    for ch in s[1:]:
        x = ord(ch) - 63
        bits += [(x >> k) & 1 for k in (5, 4, 3, 2, 1, 0)]
    E = []; k = 0
    for j in range(1, n):
        for i in range(j):
            if bits[k]: E.append((i, j))
            k += 1
    return n, E


def enum(n, flags='-c', resmod=None, keep=400):
    cmd = ['geng', '-q', flags, '-d3', '-D3', str(n)] + ([resmod] if resmod else [])
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    rng = random.Random(1); verts = list(range(n))
    cnt = 0; zero = 0; mod5 = Counter(); q4mod = Counter(); heap = []; hist = Counter(); gcd = 0; t0 = time.time()
    for line in proc.stdout:
        _, E = g6_edges(line)
        plan, _ = best_plan(E, verts, restarts=2, rng=rng)
        f = _contract(E, verts, plan, 5)
        cnt += 1
        if f == 0: zero += 1; continue
        gcd = math.gcd(gcd, f); mod5[f % 5] += 1; q4mod[(f // 4) % 6] += 1
        if len(heap) < keep: heapq.heappush(heap, (-f, line.strip()))
        elif -heap[0][0] > f: heapq.heapreplace(heap, (-f, line.strip()))
    low = sorted((-a, b) for a, b in heap)
    res = dict(n=n, flags=flags, resmod=resmod, count=cnt, bridged=zero, gcd=gcd, mod5=dict(mod5), F4mod6=dict(q4mod),
               lowest=low, secs=round(time.time() - t0, 1))
    return res


def classify(g6):
    n, E = g6_edges(g6); G = nx.Graph(E)
    lam = nx.edge_connectivity(G); gi = nx.girth(G)
    tri = sum(nx.triangles(G).values()) // 3
    c4 = sum(1 for c in nx.simple_cycles(G, length_bound=4) if len(c) == 4)
    cyc = cyc_conn(G, 6) if lam >= 3 else lam
    f4 = flow_count(E, 4)
    return dict(lam=lam, girth=gi, tri=tri, c4=c4, cyc=cyc, F4=f4)


def enum_report(n, flags='-c', resmod=None):
    res = enum(n, flags, resmod)
    tag = f"enum_{n}{flags}" + (('_' + resmod.replace('/', 'of')) if resmod else '')
    print(f"n={n} geng {flags} {resmod or ''}: {res['count']} graphs, {res['bridged']} with F=0 (bridge), gcd of nonzero F = {res['gcd']}, "
          f"F mod 5: {sorted(res['mod5'].items())}, (F/4) mod 6: {sorted(res['F4mod6'].items())}  [{res['secs']}s]", flush=True)
    # classify the lowest ones until we have seen a few of every class
    seen = Counter(); rows = []
    for f, g6 in res['lowest']:
        need = [k for k in ('c3', 'c4', 'c5g5') if seen[k] < 5]
        if not need: break
        c = classify(g6)
        keys = []
        if c['lam'] >= 3: keys.append('c3')
        if c['cyc'] >= 4: keys.append('c4')
        if c['cyc'] >= 5 and c['girth'] >= 5: keys.append('c5g5')
        if not any(seen[k] < 5 for k in keys): continue
        for k in keys: seen[k] += 1
        c.update(F5=f, g6=g6, classes=keys); rows.append(c)
        print(f"   F5={f:<10} F^(1/n)={f ** (1 / n):.4f} lam={c['lam']} cyc={c['cyc']} girth={c['girth']} tri={c['tri']} C4={c['c4']} F4={c['F4']} {g6}", flush=True)
    res['classified'] = rows; res['lowest'] = res['lowest'][:100]
    json.dump(res, open(os.path.join(HERE, f'flowdata_{tag}.json'), 'w'))
    for k, lab in (('c3', '3-connected'), ('c4', 'cyclically 4-connected'), ('c5g5', 'cyclically 5-connected, girth>=5')):
        xs = [r['F5'] for r in rows if k in r['classes']]
        print(f"   MIN over {lab}: {min(xs) if xs else 'none among the ' + str(len(res['lowest'])) + ' lowest'}", flush=True)


def enum_below(n, flags, bound):
    """exact minimum of F over cyclically 4-connected graphs: collect all graphs with F <= bound, test them all"""
    proc = subprocess.Popen(['geng', '-q', flags, '-d3', '-D3', str(n)], stdout=subprocess.PIPE, text=True)
    rng = random.Random(1); verts = list(range(n)); rows = []; cnt = 0
    for line in proc.stdout:
        _, E = g6_edges(line); cnt += 1
        plan, _ = best_plan(E, verts, restarts=2, rng=rng)
        f = _contract(E, verts, plan, 5)
        if 0 < f <= bound: rows.append((f, line.strip()))
    rows.sort(); out = []
    for f, g in rows:
        G = nx.Graph(g6_edges(g)[1]); c = cyc_conn(G, 6)
        if c >= 4: out.append((f, c, nx.girth(G), g))
    print(f"n={n} geng {flags}: {cnt} graphs, {len(rows)} with 0 < F <= {bound}, of which {len(out)} cyclically 4-connected", flush=True)
    for f, c, gi, g in out[:15]: print(f"   F5={f} cyc={c} girth={gi} {g}", flush=True)
    for cc in (4, 5):
        xs = [f for f, c, gi, g in out if c >= cc]
        print(f"   MIN over cyclically {cc}-connected: {min(xs) if xs else None}", flush=True)


def snark_census(n, flags='-ctf'):
    """all graphs from geng with F(G;4)=0 and F(G;5)>0 (snarks, including those with 3-cuts): F(G;5) statistics"""
    proc = subprocess.Popen(['geng', '-q', flags, '-d3', '-D3', str(n)], stdout=subprocess.PIPE, text=True)
    rng = random.Random(1); verts = list(range(n)); rows = []; cnt = 0
    for line in proc.stdout:
        _, E = g6_edges(line); cnt += 1
        plan, _ = best_plan(E, verts, restarts=2, rng=rng)
        if _contract(E, verts, plan, 4) != 0: continue
        f = _contract(E, verts, plan, 5)
        if f == 0: continue
        rows.append((f, line.strip()))
    rows.sort()
    print(f"n={n} geng {flags}: {cnt} graphs, {len(rows)} snarks (F4=0, bridgeless); F5 mod 5: {sorted(Counter(f % 5 for f, _ in rows).items())}; "
          f"F5 min {rows[0][0] if rows else None} max {rows[-1][0] if rows else None}", flush=True)
    for f, g in rows:
        G = nx.Graph(g6_edges(g)[1]); c = cyc_conn(G, 6)
        if c >= 4: print(f"   cyc>={c} snark: F5={f} F/1.5^n={f / 1.5 ** n:.3f} {g}", flush=True)
    json.dump(rows, open(os.path.join(HERE, f'flowdata_snarks_{n}.json'), 'w'))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'families': families()
    elif cmd == 'random': random_graphs()
    elif cmd == 'below':
        enum_below(int(sys.argv[2]), sys.argv[3], int(sys.argv[4]))
    elif cmd == 'snarks':
        for n in sys.argv[2:]: snark_census(int(n))
    elif cmd == 'enum':
        n = int(sys.argv[2]); flags = sys.argv[3] if len(sys.argv) > 3 else '-c'
        enum_report(n, flags, sys.argv[4] if len(sys.argv) > 4 else None)
