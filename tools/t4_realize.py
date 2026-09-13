"""Realize a (T4) blueprint: atoms with prescribed sizes (membership pattern in S_1..S_4), one z in
each of the six z-atoms, paths z1z2, z3z4, z5z6.  Everything else is free: the cubic graph with its
edge colours (1 = matching, 2/3 = 2-factor circuits, 0 = the 0-edge at each z), the base 2-colouring
x of H, path labels with monotone potentials.  Edges between atoms must be H-edges on paths (s = 0).
For each j, S_j (a fixed vertex set) must be a bad cut for the pattern colouring (S_j-ends black);
by s = 0 the complement colourings are killed by V \\ S_j, so all 8 colourings die.  Girth >= 6 and
cyclic 6-edge-connectivity are enforced lazily (short cycles / small cyclic cuts excluded, resolve).
usage: t4_realize.py "<blueprint>" [time] [workers]   blueprint like 0000:23,0011:3,0100:2,...
"""
import sys, os, time, itertools, json
import networkx as nx
from ortools.sat.python import cp_model
from gen_cyc6 import has_small_cyclic_cut
bp = {k: int(v) for k, v in (kv.split(":") for kv in sys.argv[1].split(","))}
TIME = float(sys.argv[2]) if len(sys.argv) > 2 else 600; WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else 10
zatoms = {"1111": 1, "0000": 1, "1100": 2, "0011": 2, "1010": 3, "0101": 3}     # atom -> path label of its z
for a in zatoms: assert bp.get(a, 0) >= 3, f"z-atom {a} needs >= 3 vertices"
verts = []; atom_of = {}; zs = {}
for a, k in bp.items():
    for i in range(k):
        v = len(verts); verts.append(v); atom_of[v] = a
        if i == 0 and a in zatoms: zs[v] = zatoms[a]
n = len(verts); t = 3
inS = {j: {v for v in verts if atom_of[v][j] == "1"} for j in range(4)}
print(f"n={n}, atoms {bp}, |S_j| = {[len(inS[j]) for j in range(4)]}, z's {zs}", flush=True)
m = cp_model.CpModel()
pairs = list(itertools.combinations(verts, 2))
y = {}   # (u,v) -> dict colour -> Bool ; colours 1,2,3,0
pres = {}
for (u, v) in pairs:
    same = atom_of[u] == atom_of[v]
    y[u, v] = {c: m.NewBoolVar(f"y{c}_{u}_{v}") for c in (1, 2, 3, 0)}
    pres[u, v] = m.NewBoolVar(f"p_{u}_{v}"); m.Add(sum(y[u, v].values()) == pres[u, v])
    if not same: m.Add(y[u, v][3] == 0); m.Add(y[u, v][0] == 0)        # inter-atom edges are H-edges
def inc(v, c): return [y[min(u, v), max(u, v)][c] for u in verts if u != v]
for v in verts:
    m.Add(sum(inc(v, 1)) == 1)
    if v in zs: m.Add(sum(inc(v, 2)) == 0); m.Add(sum(inc(v, 3)) == 1); m.Add(sum(inc(v, 0)) == 1)
    else:       m.Add(sum(inc(v, 2)) == 1); m.Add(sum(inc(v, 3)) + sum(inc(v, 0)) == 1)
for (u, v) in pairs:     # every 0-edge joins a z to a non-z vertex
    if (u in zs) == (v in zs): m.Add(y[u, v][0] == 0)
# no triangles (girth >= 6 completed lazily)
for a, b, c in (itertools.combinations(verts, 3) if os.environ.get('NOTRI') != '1' else []):
    m.AddBoolOr([pres[a, b].Not(), pres[a, c].Not(), pres[b, c].Not()])
# H-edge indicator
h = {p: m.NewBoolVar("") for p in pairs}
for p in pairs: m.Add(h[p] == y[p][1] + y[p][2])
# base colouring proper on H
x = {v: m.NewBoolVar(f"x{v}") for v in verts}
for (u, v) in pairs: m.Add(x[u] != x[v]).OnlyEnforceIf(h[u, v])
# labels: one-hot 0..3, propagated along H-edges; z labels fixed; potentials monotone on nonzero labels
lab = {v: [m.NewBoolVar(f"L{v}_{l}") for l in range(t + 1)] for v in verts}
for v in verts:
    m.AddExactlyOne(lab[v])
    if v in zs: m.Add(lab[v][zs[v]] == 1)
for (u, v) in pairs:
    for l in range(t + 1): m.Add(lab[u][l] == lab[v][l]).OnlyEnforceIf(h[u, v])
    if atom_of[u] != atom_of[v]: m.AddImplication(pres[u, v], lab[u][0].Not())   # inter-atom edges lie on paths
p = {v: m.NewIntVar(0, n, f"pot{v}") for v in verts}
lo = {v: [] for v in verts}; hi = {v: [] for v in verts}
for (u, v) in pairs:
    for a, b in ((u, v), (v, u)):
        L_ = m.NewBoolVar(""); H_ = m.NewBoolVar("")
        m.AddImplication(L_, h[u, v]); m.AddImplication(H_, h[u, v])
        m.Add(p[b] == p[a] - 1).OnlyEnforceIf(L_); m.Add(p[b] == p[a] + 1).OnlyEnforceIf(H_)
        m.AddBoolOr([h[u, v].Not(), lab[a][0], L_, H_]); lo[a].append(L_); hi[a].append(H_)
for v in verts:
    nz = lab[v][0].Not()
    if os.environ.get('NOPOT') == '1': break
    if v in zs: m.Add(sum(lo[v]) + sum(hi[v]) == 1).OnlyEnforceIf(nz)
    else: m.Add(sum(lo[v]) == 1).OnlyEnforceIf(nz); m.Add(sum(hi[v]) == 1).OnlyEnforceIf(nz)
# the four bad cuts under their pattern colourings
zend = {}   # (j, label) -> the z of that path inside S_j
for v, l in zs.items():
    for j in range(4):
        if v in inS[j]: zend[j, l] = v
cut_expr = {}
for j in range(4):
    black = {}
    for v in verts:
        # flip path l iff its S_j-end is white under x:  black = x_v xor (lab_v = l and not x_{zend})
        flips = []
        for l in range(1, t + 1):
            f = m.NewBoolVar(""); m.AddBoolAnd([lab[v][l], x[zend[j, l]].Not()]).OnlyEnforceIf(f)
            m.AddBoolOr([lab[v][l].Not(), x[zend[j, l]], f]); flips.append(f)
        fl = m.NewBoolVar(""); m.Add(fl == sum(flips))
        b = m.NewBoolVar(""); m.AddBoolXOr([x[v], fl, b.Not()]); black[v] = b
    d = sum(pres[u, v] for (u, v) in pairs if (u in inS[j]) != (v in inS[j]))
    imbalance = sum(2 * black[v] - 1 for v in inS[j])
    if os.environ.get('NOCUT') != '1': m.Add(5 * imbalance - 3 * d >= 1)
    cut_expr[j] = d
print("model built; solving...", flush=True)
solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = TIME; solver.parameters.num_workers = WORKERS
t0 = time.time(); rounds = 0
while True:
    st = solver.Solve(m); rounds += 1
    print(f"round {rounds}: {solver.StatusName(st)} [{time.time()-t0:.0f}s]", flush=True)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE): break
    G = nx.Graph(); G.add_nodes_from(verts); col = {}
    for (u, v) in pairs:
        for c in (1, 2, 3, 0):
            if solver.Value(y[u, v][c]): G.add_edge(u, v); col[frozenset((u, v))] = c
    g = nx.girth(G); added = 0
    if g < 6:   # exclude every short cycle found
        for v in verts:
            for cyc in nx.cycle_basis(G):
                pass
        Gd = G.copy()
        for cyc in nx.simple_cycles(G, length_bound=5):
            es = [(min(a, b), max(a, b)) for a, b in zip(cyc, cyc[1:] + cyc[:1])]
            m.AddBoolOr([pres[e].Not() for e in es]); added += 1
            if added >= 200: break
        print(f"   girth {g}: excluded {added} short cycles", flush=True); continue
    X = has_small_cyclic_cut(G)
    if X is not None:
        Xs = set(X); m.Add(sum(pres[u, v] for (u, v) in pairs if (u in Xs) != (v in Xs)) >= 6)
        print(f"   small cyclic cut |X|={len(Xs)} d={nx.cut_size(G, Xs)}: excluded", flush=True); continue
    out = {"edges": [list(e) for e in G.edges()], "colours": {f"{min(e)},{max(e)}": c for e, c in ((tuple(k), c) for k, c in col.items())},
           "zs": zs, "x": {v: solver.Value(x[v]) for v in verts}, "labels": {v: [solver.Value(b) for b in lab[v]].index(1) for v in verts},
           "atoms": atom_of, "cuts": [solver.Value(cut_expr[j]) for j in range(4)]}
    json.dump(out, open("t4_realized.json", "w")); print("REALIZED: cyclically 6-connected, girth >= 6 -> t4_realized.json", flush=True); break
