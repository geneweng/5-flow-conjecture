"""Count relaxation of the (T4) atom model (no walk order, no bound on segments): a complete proof
of the pure case if infeasible.  Key fact: at a crossing of S_j by path i the colour of the crossing
edge is forced by the direction and sigma_ij (leaving S_j: colour 1 iff sigma=0; entering: colour 2 iff
sigma=0), so every directed atom transition (a->b) of path i has a forced colour or is forbidden
(conflicting requirements from two cuts).  Variables x[i,a,b] = number of transitions a->b of path i.
Constraints: flow conservation of a trail from the start atom to the end atom; 7 colour-1 and 4
colour-2 crossings per cut, per-path crossing numbers odd; atom parities; n(a) >= visits(a) with
n(a) = z(a) mod 2; for atom unions U (<= UMAX atoms): d(U) <= 3 n(U), n(U)<=5 => d(U) >= n(U)+2,
n(U)>=1 => d(U) >= 3, d(U) <= 5 => n(U) <= 3; total transitions <= 44.
usage: t4_count.py [UMAX] [time]"""
import sys, itertools, time
from ortools.sat.python import cp_model
UMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3; TIME = float(sys.argv[2]) if len(sys.argv) > 2 else 600
cells = list(itertools.product((0, 1), repeat=4)); A = range(16)
ends = {0: ((1, 1, 1, 1), (0, 0, 0, 0)), 1: ((1, 1, 0, 0), (0, 0, 1, 1)), 2: ((1, 0, 1, 0), (0, 1, 0, 1))}
zc = {a: sum(1 for i in ends for c in ends[i] if c == cells[a]) for a in A}
def forced_colour(i, a, b):
    """colour (1 or 2) forced for a transition a->b of path i, or None if forbidden; a==b never used"""
    col = None
    for j in range(4):
        if cells[a][j] == cells[b][j]: continue
        sigma = 0 if ends[i][0][j] == 1 else 1
        leaving = cells[a][j] == 1
        c = (1 if sigma == 0 else 2) if leaving else (2 if sigma == 0 else 1)
        if col is None: col = c
        elif col != c: return None
    return col
m = cp_model.CpModel(); x = {}
for i in range(3):
    for a in A:
        for b in A:
            if a == b: continue
            c = forced_colour(i, a, b)
            if c is None: continue
            x[i, a, b] = (m.NewIntVar(0, 9, f"x{i}_{a}_{b}"), c)
def out(i, a): return [x[i, a, b][0] for b in A if (i, a, b) in x]
def inn(i, a): return [x[i, b, a][0] for b in A if (i, b, a) in x]
for i in range(3):
    sa = cells.index(ends[i][0]); ea = cells.index(ends[i][1])
    for a in A:
        m.Add(sum(out(i, a)) - sum(inn(i, a)) == (1 if a == sa else -1 if a == ea else 0))
# cut counts and per-path crossing numbers
for j in range(4):
    c1 = []; c2 = []
    for i in range(3):
        cr = [(v, c) for (ii, a, b), (v, c) in x.items() if ii == i and cells[a][j] != cells[b][j]]
        r = m.NewIntVar(1, 9, f"r{i}_{j}"); h = m.NewIntVar(0, 4, ""); m.Add(r == 2 * h + 1); m.Add(r == sum(v for v, c in cr))
        c1 += [v for v, c in cr if c == 1]; c2 += [v for v, c in cr if c == 2]
    m.Add(sum(c1) == 7); m.Add(sum(c2) == 4)
# atoms: visits, n, parities
n = {}; visits = {}
for a in A:
    visits[a] = sum(inn(i, a) for i in range(3)) if False else sum(v for i in range(3) for v in inn(i, a)) + sum(1 for i in range(3) if cells.index(ends[i][0]) == a)
    n[a] = m.NewIntVar(0, 500, f"n{a}"); m.Add(n[a] >= visits[a])
    h = m.NewIntVar(0, 250, ""); m.Add(n[a] == 2 * h + zc[a] % 2)
    inc1 = [v for (i, p, q), (v, c) in x.items() if (p == a or q == a) and c == 1]
    inc2 = [v for (i, p, q), (v, c) in x.items() if (p == a or q == a) and c == 2]
    h1 = m.NewIntVar(0, 100, ""); m.Add(sum(inc1) == 2 * h1 + zc[a] % 2)
    h2 = m.NewIntVar(0, 100, ""); m.Add(sum(inc2) == 2 * h2)
    # a nonempty atom that no path enters consists of even H-circuit vertices only; then n(a) >= 6 (girth) or 0
    if all(cells.index(ends[i][0]) != a for i in range(3)):
        nov = m.NewBoolVar(""); m.Add(sum(v for i in range(3) for v in inn(i, a)) == 0).OnlyEnforceIf(nov)
        m.Add(sum(v for i in range(3) for v in inn(i, a)) >= 1).OnlyEnforceIf(nov.Not())
        big = m.NewBoolVar(""); m.Add(n[a] >= 6).OnlyEnforceIf([nov, big]); m.Add(n[a] == 0).OnlyEnforceIf([nov, big.Not()])
nU = 0
for r in range(1, UMAX + 1):
    for U in itertools.combinations(A, r):
        Us = set(U); dU = sum(v for (i, p, q), (v, c) in x.items() if (p in Us) != (q in Us)); nUv = sum(n[a] for a in U); nU += 1
        d = m.NewIntVar(0, 200, ""); m.Add(d == dU); nn = m.NewIntVar(0, 8000, ""); m.Add(nn == nUv)
        m.Add(d <= 3 * nn)
        empty = m.NewBoolVar(""); m.Add(nn == 0).OnlyEnforceIf(empty); m.Add(nn >= 1).OnlyEnforceIf(empty.Not())
        le5 = m.NewBoolVar(""); m.Add(nn <= 5).OnlyEnforceIf(le5); m.Add(nn >= 6).OnlyEnforceIf(le5.Not()); m.Add(d >= nn + 2).OnlyEnforceIf([le5, empty.Not()])
        m.Add(d >= 3).OnlyEnforceIf(empty.Not()); m.Add(d == 0).OnlyEnforceIf(empty)
        if r <= 2:   # complement holds >= 4 z's, so the small side (if any) is U
            small = m.NewBoolVar(""); m.Add(nn <= 3).OnlyEnforceIf(small); m.Add(d >= 6).OnlyEnforceIf(small.Not())
        else:        # either side may be the small (<= 3 vertices, acyclic) one
            small = m.NewBoolVar(""); csmall = m.NewBoolVar("")
            m.Add(nn <= 3).OnlyEnforceIf(small); m.Add(sum(n[a] for a in A if a not in Us) <= 3).OnlyEnforceIf(csmall)
            m.Add(d >= 6).OnlyEnforceIf([small.Not(), csmall.Not()])
m.Add(sum(v for v, c in x.values()) <= 44)
# ---- 2-factor circuits as closed walks through atoms via colour-2 crossings.
# Classes k = 0..5: the odd circuit of z_{k+1} (z-atoms in order 1111,0000,1100,0011,1010,0101); class 6: all even circuits.
zatom_of = {0: cells.index((1,1,1,1)), 1: cells.index((0,0,0,0)), 2: cells.index((1,1,0,0)), 3: cells.index((0,0,1,1)), 4: cells.index((1,0,1,0)), 5: cells.index((0,1,0,1))}
upairs = [(a, b) for a in A for b in A if a < b]
t2 = {(a, b): sum(v for (i, p, q), (v, c) in x.items() if c == 2 and {p, q} == {a, b}) for (a, b) in upairs}
w = {(k, a, b): m.NewIntVar(0, 9, f"w{k}_{a}_{b}") for k in range(7) for (a, b) in upairs}
for (a, b) in upairs: m.Add(sum(w[k, a, b] for k in range(7)) == t2[a, b])
deg = {}
for k in range(7):
    for a in A:
        deg[k, a] = sum(w[k, p, q] for (p, q) in upairs if a in (p, q))
        hh = m.NewIntVar(0, 20, ""); m.Add(deg[k, a] == 2 * hh)            # closed walk: even degree at every atom
for a in A:
    m.Add(sum(deg[k, a] for k in range(7)) + zc[a] <= n[a])                # every colour-2 crossing uses its own vertex; z has none
for k in range(6):
    za = zatom_of[k]
    small = m.NewBoolVar(""); m.Add(n[za] <= 6).OnlyEnforceIf(small); m.Add(n[za] >= 7).OnlyEnforceIf(small.Not())
    m.Add(deg[k, za] >= 2).OnlyEnforceIf(small)                              # a circuit of length >= 7 must leave a small atom
    # a 3-vertex z-atom {u,z,w} is traversed only by z's own circuit
    three = m.NewBoolVar(""); m.Add(n[za] == 3).OnlyEnforceIf(three); m.Add(n[za] != 3).OnlyEnforceIf(three.Not())
    for kk in range(7):
        if kk != k: m.Add(deg[kk, za] == 0).OnlyEnforceIf(three)
    m.Add(deg[k, za] == 2).OnlyEnforceIf(three)
# room for circuit visits: each visit of a circuit to an atom uses >= 2 vertices (the visit through z uses 3)
for a in A:
    m.Add(sum(deg[k, a] for k in range(7)) + zc[a] <= n[a])
m.Add(sum(n[a] for a in A) >= 42)                     # six odd circuits of length >= 7
for a in A:
    if zc[a]: m.Add(n[a] >= 3)                        # a z has two non-H edges inside its atom
NMAX = int(__import__("os").environ.get("NMAX", "60")); m.Add(sum(n[a] for a in A) <= NMAX)
if __import__("os").environ.get("MINATOM"):
    MA = int(__import__("os").environ["MINATOM"])
    for a in A:
        e = m.NewBoolVar(""); m.Add(n[a] == 0).OnlyEnforceIf(e); m.Add(n[a] >= MA).OnlyEnforceIf(e.Not())
def add_union(U):
    Us = set(U); dU = sum(v for (i, p, q), (v, c) in x.items() if (p in Us) != (q in Us)); nUv = sum(n[a] for a in U)
    d = m.NewIntVar(0, 200, ""); m.Add(d == dU); nn = m.NewIntVar(0, 8000, ""); m.Add(nn == nUv); cn = m.NewIntVar(0, 8000, ""); m.Add(cn == sum(n[a] for a in A if a not in Us))
    m.Add(d <= 3 * nn); m.Add(d <= 3 * cn)
    for side in (nn, cn):
        s15 = m.NewBoolVar(""); m.AddLinearConstraint(side, 1, 5).OnlyEnforceIf(s15)
        m.AddBoolOr([s15, m.NewBoolVar("")]) if False else None
        out15 = m.NewBoolVar(""); m.Add(side == 0).OnlyEnforceIf([s15.Not(), out15]); m.Add(side >= 6).OnlyEnforceIf([s15.Not(), out15.Not()])
        m.Add(d >= side + 2).OnlyEnforceIf(s15)
    e1 = m.NewBoolVar(""); m.Add(nn == 0).OnlyEnforceIf(e1); m.Add(nn >= 1).OnlyEnforceIf(e1.Not())
    e2 = m.NewBoolVar(""); m.Add(cn == 0).OnlyEnforceIf(e2); m.Add(cn >= 1).OnlyEnforceIf(e2.Not())
    m.Add(d >= 3).OnlyEnforceIf([e1.Not(), e2.Not()]); m.Add(d == 0).OnlyEnforceIf(e1); m.Add(d == 0).OnlyEnforceIf(e2)
    small = m.NewBoolVar(""); csmall = m.NewBoolVar("")
    m.Add(nn <= 3).OnlyEnforceIf(small); m.Add(cn <= 3).OnlyEnforceIf(csmall); m.Add(d >= 6).OnlyEnforceIf([small.Not(), csmall.Not()])
def violated(sol_x, sol_n):
    """all unions (bitmasks) violating cubic / girth / 3-edge-conn / cyclic-6 conditions"""
    N = sum(sol_n); bad = []
    tr = [(1 << p | 0, 1 << q, v) for (i, p, q), v in sol_x.items() if v]
    for mask in range(1, 1 << 16):
        nU = sum(sol_n[a] for a in A if mask >> a & 1); cn = N - nU
        d = sum(v for pm, qm, v in tr if bool(pm & mask) != bool(qm & mask))
        ok = d <= 3 * nU and d <= 3 * cn
        for side in (nU, cn):
            if 1 <= side <= 5 and d < side + 2: ok = False
        if nU >= 1 and cn >= 1 and d < 3: ok = False
        if (nU == 0 or cn == 0) and d != 0: ok = False
        if d <= 5 and nU > 3 and cn > 3: ok = False
        if not ok: bad.append((mask, nU, d))
    return bad
print(f"transition variables {len(x)}, unions {nU}; solving...", flush=True)
solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = TIME; solver.parameters.num_workers = 10
t0 = time.time(); rounds = 0
while True:
    st = solver.Solve(m); rounds += 1; print(f"round {rounds}: status {solver.StatusName(st)} [{time.time()-t0:.0f}s]", flush=True)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE): break
    sol_x = {k: solver.Value(v) for k, (v, c) in x.items()}; sol_n = [solver.Value(n[a]) for a in A]
    bad = violated(sol_x, sol_n)
    if not bad: print("no violated union: relaxation genuinely feasible", flush=True); break
    bad.sort(key=lambda t: (bin(t[0]).count("1"), t[2]))
    for mask, nUv, d in bad[:40]: add_union([a for a in A if mask >> a & 1])
    print(f"   {len(bad)} violated unions (smallest: {bin(bad[0][0]).count('1')} atoms, n={bad[0][1]}, d={bad[0][2]}); added {min(40, len(bad))}", flush=True)
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    for i in range(3):
        print(f"P{i+1}:", {(''.join(map(str, cells[a])), ''.join(map(str, cells[b])), c): solver.Value(v) for (ii, a, b), (v, c) in x.items() if ii == i and solver.Value(v)})
    print("n:", {''.join(map(str, cells[a])): solver.Value(n[a]) for a in A if solver.Value(n[a])})
