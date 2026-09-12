"""Abstract feasibility of (T4): four tight 11-cuts S_1..S_4 (type (11,7,4,3), s = 0) whose z-triples
realize the four classes, in a cyclically 6-edge-connected cubic graph with the canonical colouring.
Atoms = 16 cells a in {0,1}^4 (membership in S_1..S_4).  Every edge between different cells lies in
some dS_j, hence is an H-edge, hence a path edge (s = 0).  So each path P_i is a walk through cells:
segments (cell, length).  Colours along a path alternate 1,2,1,... from z; under the pattern colouring
of S_j the vertices are black iff position = sigma_ij (mod 2), sigma_ij = 0 iff the start z is in S_j.
Constraints: per cut 7 colour-1 and 4 colour-2 crossings, S_j-end of every crossing edge black;
per cell parities from the three perfect matchings (colour 1, colour 2 U {z}, colours 0/3);
d(U) >= 6 or n(U) <= 3 for unions U of at most UMAX cells.
usage: t4_abstract.py [K segments per path] [UMAX] [time]"""
import sys, itertools, time
from ortools.sat.python import cp_model
K = int(sys.argv[1]) if len(sys.argv) > 1 else 45; UMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 2; TIME = float(sys.argv[3]) if len(sys.argv) > 3 else 600
cells = list(itertools.product((0, 1), repeat=4)); cid = {c: i for i, c in enumerate(cells)}
# end cells: z1 in all four, z2 in none; z3 in S1,S2; z4 in S3,S4; z5 in S1,S3; z6 in S2,S4
ends = {0: ((1, 1, 1, 1), (0, 0, 0, 0)), 1: ((1, 1, 0, 0), (0, 0, 1, 1)), 2: ((1, 0, 1, 0), (0, 1, 0, 1))}
zcount = {c: 0 for c in cells}
for i in ends:
    for c in ends[i]: zcount[c] += 1
m = cp_model.CpModel(); LMAX = 16   # lengths >= 15 are recorded as 15/16 by parity: valid for every constraint used (thresholds 3 and 5; d <= 3n vacuous since d <= 44)
used = {}; cell = {}; length = {}; par = {}; pos = {}
for i in range(3):
    for s in range(K):
        used[i, s] = m.NewBoolVar(f"u{i}_{s}"); length[i, s] = m.NewIntVar(0, LMAX, f"l{i}_{s}")
        par[i, s] = m.NewBoolVar(f"p{i}_{s}"); pos[i, s] = m.NewBoolVar(f"q{i}_{s}")   # pos = parity of positions before segment s
        cell[i, s] = [m.NewBoolVar(f"c{i}_{s}_{a}") for a in range(16)]
for i in range(3):
    for s in range(K):
        m.Add(sum(cell[i, s]) == 1)
        if s > 0: m.AddImplication(used[i, s], used[i, s - 1])
        m.Add(length[i, s] >= 1).OnlyEnforceIf(used[i, s]); m.Add(length[i, s] == 0).OnlyEnforceIf(used[i, s].Not())
        q = m.NewIntVar(0, LMAX, ""); m.Add(length[i, s] == 2 * q + par[i, s]); m.AddImplication(used[i, s].Not(), par[i, s].Not())
        if s == 0: m.Add(pos[i, 0] == 0); m.Add(used[i, 0] == 1); m.Add(cell[i, 0][cid[ends[i][0]]] == 1)
        else: m.AddBoolXOr([pos[i, s - 1], par[i, s - 1], pos[i, s].Not()])   # pos_s = pos_{s-1} + par_{s-1}
        # last used segment ends in the end cell; total length even (last vertex index odd)
        last = m.NewBoolVar("")
        if s + 1 < K: m.AddBoolAnd([used[i, s], used[i, s + 1].Not()]).OnlyEnforceIf(last); m.AddBoolOr([used[i, s].Not(), used[i, s + 1], last])
        else: m.Add(last == used[i, s])
        m.Add(cell[i, s][cid[ends[i][1]]] == 1).OnlyEnforceIf(last)
        tot_par = m.NewBoolVar(""); m.AddBoolXOr([pos[i, s], par[i, s], tot_par.Not()]); m.Add(tot_par == 0).OnlyEnforceIf(last)
        # consecutive used segments lie in different cells
        if s > 0:
            for a in range(16): m.AddBoolOr([used[i, s].Not(), cell[i, s - 1][a].Not(), cell[i, s][a].Not()])
def coord(i, s, j): return sum(cell[i, s][a] for a in range(16) if cells[a][j] == 1)
def inU(i, s, U): return sum(cell[i, s][a] for a in U)
# transitions: after segment s (exists iff used[i,s+1]); edge parity pm = pos_{s+1} + 1
trans = []
for i in range(3):
    for s in range(K - 1):
        ex = used[i, s + 1]; pm = pos[i, s + 1].Not()          # pm = 1 - pos_{s+1}  (pm even <=> colour 1)
        col1 = pos[i, s + 1]                                    # colour 1 iff pm even iff pos_{s+1} = 1
        trans.append((i, s, ex, col1))
cnt1 = {j: [] for j in range(4)}; cnt2 = {j: [] for j in range(4)}
for (i, s, ex, col1) in trans:
    for j in range(4):
        sigma = 0 if ends[i][0][j] == 1 else 1
        c0 = m.NewBoolVar(""); m.Add(c0 == coord(i, s, j)); c1 = m.NewBoolVar(""); m.Add(c1 == coord(i, s + 1, j))
        leave = m.NewBoolVar(""); m.AddBoolAnd([ex, c0, c1.Not()]).OnlyEnforceIf(leave); m.AddBoolOr([ex.Not(), c0.Not(), c1, leave])
        enter = m.NewBoolVar(""); m.AddBoolAnd([ex, c0.Not(), c1]).OnlyEnforceIf(enter); m.AddBoolOr([ex.Not(), c0, c1.Not(), enter])
        # black S_j-end: leaving -> pm == sigma ; entering -> pm == 1 - sigma.   pm = 1 - pos_{s+1}
        pm_val = 1 - sigma  # required pos_{s+1} value when leaving
        m.Add(pos[i, s + 1] == pm_val).OnlyEnforceIf(leave); m.Add(pos[i, s + 1] == 1 - pm_val).OnlyEnforceIf(enter)
        x1 = m.NewBoolVar(""); m.AddMultiplicationEquality(x1, [leave, col1]); y1 = m.NewBoolVar(""); m.AddMultiplicationEquality(y1, [enter, col1])
        x2 = m.NewBoolVar(""); m.AddMultiplicationEquality(x2, [leave, col1.Not()]); y2 = m.NewBoolVar(""); m.AddMultiplicationEquality(y2, [enter, col1.Not()])
        cnt1[j] += [x1, y1]; cnt2[j] += [x2, y2]
for j in range(4): m.Add(sum(cnt1[j]) == 7); m.Add(sum(cnt2[j]) == 4)
# per cut and path: an odd number r in {1,3,5,7,9} of crossings, summing to 11 over the three paths
for j in range(4):
    rs = []
    for i in range(3):
        cr = []
        for (ii, s, ex, col1) in trans:
            if ii != i: continue
            c0 = m.NewBoolVar(""); m.Add(c0 == coord(i, s, j)); c1 = m.NewBoolVar(""); m.Add(c1 == coord(i, s + 1, j))
            x = m.NewBoolVar(""); m.AddBoolXOr([c0, c1, x.Not()]); y = m.NewBoolVar(""); m.AddMultiplicationEquality(y, [x, ex]); cr.append(y)
        r = m.NewIntVar(1, 9, f"r{i}_{j}"); h = m.NewIntVar(0, 4, ""); m.Add(r == 2 * h + 1); m.Add(r == sum(cr)); rs.append(r)
    m.Add(sum(rs) == 11)
# segments per path: at most 1 + total crossings of that path over the four cuts (each transition crosses a cut)
for i in range(3):
    m.Add(sum(used[i, s] for s in range(K)) <= 1 + sum(m.GetIntVarFromProtoIndex(v.Index()) for v in [] ) + 36)
# per-cell parities: n_a = path vertices + extra;  n_a = z_a (mod 2), colour-1 crossings at a = z_a (mod 2), colour-2 crossings at a even
n = {}; extra = {}
for a in range(16):
    extra[a] = m.NewIntVar(0, 16, f"e{a}")
    pv = []
    for i in range(3):
        for s in range(K):
            v = m.NewIntVar(0, LMAX, ""); m.Add(v == length[i, s]).OnlyEnforceIf(cell[i, s][a]); m.Add(v == 0).OnlyEnforceIf(cell[i, s][a].Not()); pv.append(v)
    n[a] = m.NewIntVar(0, 600, f"n{a}"); m.Add(n[a] == sum(pv) + extra[a])
    h = m.NewIntVar(0, 1000, ""); m.Add(n[a] == 2 * h + zcount[cells[a]] % 2)
    inc1 = []; inc2 = []
    for (i, s, ex, col1) in trans:
        touch = m.NewBoolVar(""); m.AddBoolOr([cell[i, s][a], cell[i, s + 1][a]]).OnlyEnforceIf(touch); m.AddImplication(touch, ex)
        m.AddBoolOr([touch, ex.Not(), cell[i, s][a].Not()]); m.AddBoolOr([touch, ex.Not(), cell[i, s + 1][a].Not()])
        t1 = m.NewBoolVar(""); m.AddMultiplicationEquality(t1, [touch, col1]); t2 = m.NewBoolVar(""); m.AddMultiplicationEquality(t2, [touch, col1.Not()])
        inc1.append(t1); inc2.append(t2)
    h1 = m.NewIntVar(0, 100, ""); m.Add(sum(inc1) == 2 * h1 + zcount[cells[a]] % 2)
    h2 = m.NewIntVar(0, 100, ""); m.Add(sum(inc2) == 2 * h2)
# cyclic 6-edge-connectivity for unions of at most UMAX cells
nU = 0
for r in range(1, UMAX + 1):
    for U in itertools.combinations(range(16), r):
        Uset = set(U); cross = []
        for (i, s, ex, col1) in trans:
            a_in = m.NewBoolVar(""); m.Add(a_in == inU(i, s, U)); b_in = m.NewBoolVar(""); m.Add(b_in == inU(i, s + 1, U))
            x = m.NewBoolVar(""); m.AddBoolXOr([a_in, b_in, x.Not()]); c = m.NewBoolVar(""); m.AddMultiplicationEquality(c, [x, ex]); cross.append(c)
        nUv = m.NewIntVar(0, 1200, ""); m.Add(nUv == sum(n[a] for a in U)); dU = m.NewIntVar(0, 200, ""); m.Add(dU == sum(cross)); nU += 1
        m.Add(dU <= 3 * nUv)                                                    # cubic
        le5 = m.NewBoolVar(""); m.Add(nUv <= 5).OnlyEnforceIf(le5); m.Add(nUv >= 6).OnlyEnforceIf(le5.Not())
        m.Add(dU >= nUv + 2).OnlyEnforceIf(le5)                                 # girth >= 6: at most 5 vertices span a forest
        empty = m.NewBoolVar(""); m.Add(nUv == 0).OnlyEnforceIf(empty); m.Add(nUv >= 1).OnlyEnforceIf(empty.Not())
        m.Add(dU >= 3).OnlyEnforceIf(empty.Not())                               # 3-edge-connected (complement holds >= 4 z's)
        m.Add(dU == 0).OnlyEnforceIf(empty)
        # cyclic 6-edge-connectivity: d(U) <= 5 => U has at most 3 vertices (complement contains >= 6-|U| >= 4 z's)
        assert r <= 2
        small = m.NewBoolVar(""); m.Add(nUv <= 3).OnlyEnforceIf(small); m.Add(dU >= 6).OnlyEnforceIf(small.Not())
m.Add(sum(ex for (_, _, ex, _) in trans) <= 44)
print(f"K={K}, unions checked {nU}; solving...", flush=True)
solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = TIME; solver.parameters.num_workers = 10
t0 = time.time(); st = solver.Solve(m); print("status", solver.StatusName(st), f"[{time.time()-t0:.0f}s]", flush=True)
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    for i in range(3):
        walk = [(cells[[solver.Value(b) for b in cell[i, s]].index(1)], solver.Value(length[i, s])) for s in range(K) if solver.Value(used[i, s])]
        print(f"P{i+1}: " + " -> ".join(f"{''.join(map(str, c))}x{l}" for c, l in walk))
    print("cells n:", {''.join(map(str, cells[a])): solver.Value(n[a]) for a in range(16) if solver.Value(n[a])})
