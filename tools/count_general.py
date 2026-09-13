"""Count relaxation of the atom model for a general configuration of tight cuts (all s = 0, boundaries in
E(H)).  Config: J cuts; for each path i the membership vectors of its two ends in the J cuts; for each
cut its (c1, c2) and which paths cross it (others have zero crossings).  Infeasible => the configuration
cannot occur in a cyclically 6-edge-connected cubic graph with the canonical colouring.
usage: count_general.py CONFIG [UMAX] [time]      CONFIG in {sit2, t4, sit3a, sit4}"""
import sys, os, itertools, time
from ortools.sat.python import cp_model

def solve(ends, cuts, UMAX=2, TIME=300, NOCIRC=False, NOUNION=False, NOPAR=False, NO42=False, NOCUTS=False, workers=8):
    """count relaxation for a configuration; returns (status, atom sizes or None)"""
    J = len(cuts)
    os.environ["NOCUTS"] = "1" if NOCUTS else "0"; os.environ["NOPAR"] = "1" if NOPAR else "0"; os.environ["NO42"] = "1" if NO42 else "0"
    cells = list(itertools.product((0, 1), repeat=J)); A = range(len(cells)); nc = len(cells)
    zc = {a: sum(1 for i in ends for c in ends[i] if c == cells[a]) for a in A}
    sep = {j: [i for i in ends if ends[i][0][j] != ends[i][1][j]] for j in range(J)}
    def forced_colour(i, a, b):
        col = None
        for j in range(J):
            if cells[a][j] == cells[b][j]: continue
            if i not in sep[j]: return None                 # a path not separated by S_j does not cross it (codim assumption)
            sigma = 0 if ends[i][0][j] == 1 else 1
            c = (1 if sigma == 0 else 2) if cells[a][j] == 1 else (2 if sigma == 0 else 1)
            if col is None: col = c
            elif col != c: return None
        return col
    m = cp_model.CpModel(); x = {}
    for i in ends:
        for a in A:
            for b in A:
                if a == b: continue
                c = forced_colour(i, a, b)
                if c is not None: x[i, a, b] = (m.NewIntVar(0, 11, f"x{i}_{a}_{b}"), c)
    def out(i, a): return [x[i, a, b][0] for b in A if (i, a, b) in x]
    def inn(i, a): return [x[i, b, a][0] for b in A if (i, b, a) in x]
    for i in ends:
        sa = cells.index(ends[i][0]); ea = cells.index(ends[i][1])
        for a in A: m.Add(sum(out(i, a)) - sum(inn(i, a)) == (1 if a == sa else 0) - (1 if a == ea else 0))
    for j in range(J):
        c1 = []; c2 = []
        for i in sep[j]:
            cr = [(v, c) for (ii, a, b), (v, c) in x.items() if ii == i and cells[a][j] != cells[b][j]]
            r = m.NewIntVar(1, 11, ""); h = m.NewIntVar(0, 5, ""); m.Add(r == 2 * h + 1); m.Add(r == sum(v for v, c in cr))
            c1 += [v for v, c in cr if c == 1]; c2 += [v for v, c in cr if c == 2]
        if os.environ.get('NOCUTS') != '1': m.Add(sum(c1) == cuts[j][0]); m.Add(sum(c2) == cuts[j][1])
    n = {}
    for a in A:
        visits = sum(v for i in ends for v in inn(i, a)) + sum(1 for i in ends if cells.index(ends[i][0]) == a)
        n[a] = m.NewIntVar(0, 500, f"n{a}"); m.Add(n[a] >= visits)
        h = m.NewIntVar(0, 250, ""); m.Add(n[a] == 2 * h + zc[a] % 2)
        inc1 = [v for (i, p, q), (v, c) in x.items() if (p == a or q == a) and c == 1]
        inc2 = [v for (i, p, q), (v, c) in x.items() if (p == a or q == a) and c == 2]
        if os.environ.get('NOPAR') != '1':
            h1 = m.NewIntVar(0, 100, ""); m.Add(sum(inc1) == 2 * h1 + zc[a] % 2)
            h2 = m.NewIntVar(0, 100, ""); m.Add(sum(inc2) == 2 * h2)
        m.Add(sum(inc2) + zc[a] <= n[a])                           # colour-2 crossings use distinct non-z vertices
        if zc[a]: m.Add(n[a] >= 3)
        if all(cells.index(ends[i][0]) != a for i in ends):
            nov = m.NewBoolVar(""); m.Add(sum(v for i in ends for v in inn(i, a)) == 0).OnlyEnforceIf(nov); m.Add(sum(v for i in ends for v in inn(i, a)) >= 1).OnlyEnforceIf(nov.Not())
            big = m.NewBoolVar(""); m.Add(n[a] >= 6).OnlyEnforceIf([nov, big]); m.Add(n[a] == 0).OnlyEnforceIf([nov, big.Not()])
    if os.environ.get('NO42') != '1': m.Add(sum(n[a] for a in A) >= 42)
    # circuits: odd circuit of each z as a closed walk (colour-2 crossings partitioned into classes)
    zlist = [(i, e) for i in ends for e in (0, 1)]
    zatom = {k: cells.index(ends[i][e]) for k, (i, e) in enumerate(zlist)}
    upairs = [(a, b) for a in A for b in A if a < b]
    t2 = {(a, b): sum(v for (i, p, q), (v, c) in x.items() if c == 2 and {p, q} == {a, b}) for (a, b) in upairs}
    K = len(zlist) + 1
    w = {(k, a, b): m.NewIntVar(0, 11, "") for k in range(K) for (a, b) in upairs}
    if not NOCIRC:
        for (a, b) in upairs: m.Add(sum(w[k, a, b] for k in range(K)) == t2[a, b])
    deg = {(k, a): sum(w[k, p, q] for (p, q) in upairs if a in (p, q)) for k in range(K) for a in A}
    for k in range(K):
        for a in A: hh = m.NewIntVar(0, 20, ""); m.Add(deg[k, a] == 2 * hh)
    for k, za in (zatom.items() if not NOCIRC else []):
        small = m.NewBoolVar(""); m.Add(n[za] <= 6).OnlyEnforceIf(small); m.Add(n[za] >= 7).OnlyEnforceIf(small.Not()); m.Add(deg[k, za] >= 2).OnlyEnforceIf(small)
        three = m.NewBoolVar(""); m.Add(n[za] == 3).OnlyEnforceIf(three); m.Add(n[za] != 3).OnlyEnforceIf(three.Not())
        for kk in range(K):
            if kk != k: m.Add(deg[kk, za] == 0).OnlyEnforceIf(three)
        m.Add(deg[k, za] == 2).OnlyEnforceIf(three)
    def add_union(U):
        Us = set(U); dU = sum(v for (i, p, q), (v, c) in x.items() if (p in Us) != (q in Us))
        d = m.NewIntVar(0, 200, ""); m.Add(d == dU); nn = m.NewIntVar(0, 8000, ""); m.Add(nn == sum(n[a] for a in U)); cn = m.NewIntVar(0, 8000, ""); m.Add(cn == sum(n[a] for a in A if a not in Us))
        m.Add(d <= 3 * nn); m.Add(d <= 3 * cn)
        for side in (nn, cn):
            s15 = m.NewBoolVar(""); m.AddLinearConstraint(side, 1, 5).OnlyEnforceIf(s15)
            o = m.NewBoolVar(""); m.Add(side == 0).OnlyEnforceIf([s15.Not(), o]); m.Add(side >= 6).OnlyEnforceIf([s15.Not(), o.Not()])
            m.Add(d >= side + 2).OnlyEnforceIf(s15)
        e1 = m.NewBoolVar(""); m.Add(nn == 0).OnlyEnforceIf(e1); m.Add(nn >= 1).OnlyEnforceIf(e1.Not())
        e2 = m.NewBoolVar(""); m.Add(cn == 0).OnlyEnforceIf(e2); m.Add(cn >= 1).OnlyEnforceIf(e2.Not())
        m.Add(d >= 3).OnlyEnforceIf([e1.Not(), e2.Not()]); m.Add(d == 0).OnlyEnforceIf(e1); m.Add(d == 0).OnlyEnforceIf(e2)
        small = m.NewBoolVar(""); csmall = m.NewBoolVar("")
        m.Add(nn <= 3).OnlyEnforceIf(small); m.Add(cn <= 3).OnlyEnforceIf(csmall); m.Add(d >= 6).OnlyEnforceIf([small.Not(), csmall.Not()])
    for r in (range(1, UMAX + 1) if not NOUNION else []):
        for U in itertools.combinations(A, r): add_union(U)
    def violated(sol_x, sol_n):
        N = sum(sol_n); bad = []
        tr = [(1 << p, 1 << q, v) for (i, p, q), v in sol_x.items() if v]
        for mask in range(1, 1 << nc):
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

    solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = TIME; solver.parameters.num_workers = workers
    t0 = time.time(); rounds = 0
    while rounds < 60:
        st = solver.Solve(m); rounds += 1; 
        if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE): return solver.StatusName(st), None
        sol_x = {k: solver.Value(v) for k, (v, c) in x.items()}; sol_n = [solver.Value(n[a]) for a in A]
        bad = violated(sol_x, sol_n)
        # walk connectivity (subtour elimination): a path's transitions must form one trail from its start atom;
        # an odd circuit's colour-2 crossings must form one closed walk through its z-atom
        import networkx as nx
        sub_added = 0
        for i in ends:
            Gp = nx.Graph(); Gp.add_edges_from((p, q) for (ii, p, q), v in sol_x.items() if ii == i and v)
            sa = cells.index(ends[i][0])
            for comp in nx.connected_components(Gp):
                if sa in comp: continue
                W = set(comp); inside = [v for (ii, p, q), (v, c) in x.items() if ii == i and p in W and q in W]
                bnd = [v for (ii, p, q), (v, c) in x.items() if ii == i and (p in W) != (q in W)]
                if not inside: continue
                use = m.NewBoolVar(""); m.Add(sum(inside) == 0).OnlyEnforceIf(use.Not()); m.Add(sum(bnd) >= 1).OnlyEnforceIf(use); sub_added += 1
        for k, za in zatom.items():
            Gc = nx.Graph(); Gc.add_edges_from((p, q) for (p, q) in upairs if solver.Value(w[k, p, q]))
            for comp in nx.connected_components(Gc):
                if za in comp: continue
                W = set(comp); inside = [w[k, p, q] for (p, q) in upairs if p in W and q in W]
                bnd = [w[k, p, q] for (p, q) in upairs if (p in W) != (q in W)]
                use = m.NewBoolVar(""); m.Add(sum(inside) == 0).OnlyEnforceIf(use.Not()); m.Add(sum(bnd) >= 1).OnlyEnforceIf(use); sub_added += 1
        if sub_added and not bad: continue
        if not bad:
            return "FEASIBLE", {''.join(map(str, cells[a])): sol_n[a] for a in A if sol_n[a]}

        bad.sort(key=lambda t: (bin(t[0]).count("1"), t[2]))
        for mask, nUv, d in bad[:40]: add_union([a for a in A if mask >> a & 1])

    return solver.StatusName(st), None

if __name__ == "__main__":
    CFG = sys.argv[1]; UMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 2; TIME = float(sys.argv[3]) if len(sys.argv) > 3 else 600
    if CFG == "t4":      # four 11-cuts, classes z1z3z5, z1z3z6, z1z4z5, z1z4z6
        J = 4; ends = {0: ((1,1,1,1), (0,0,0,0)), 1: ((1,1,0,0), (0,0,1,1)), 2: ((1,0,1,0), (0,1,0,1))}
        cuts = [(7, 4)] * 4
    elif CFG == "sit2":  # three pair 6-cuts: S12 ∋ z1,z3 ; S13 ∋ z1,z5 ; S23 ∋ z3,z5 ; third path avoids each
        J = 3; ends = {0: ((1,1,0), (0,0,0)), 1: ((1,0,1), (0,0,0)), 2: ((0,1,1), (0,0,0))}
        cuts = [(4, 2)] * 3
    elif CFG == "sit3a": # two pair cuts on different pairs (S12 ∋ z1,z3 ; S13 ∋ z1,z5) + an 11-cut killing a remaining class (z2,z4,z6 -> complement: z1,z3,z5 pattern? use z1,z4,z6)
        J = 3; ends = {0: ((1,1,1), (0,0,0)), 1: ((1,0,0), (0,0,1)), 2: ((0,1,0), (0,0,1))}
        cuts = [(4, 2), (4, 2), (7, 4)]
    elif CFG == "sit4":  # one pair cut S12 ∋ z1,z3 + two 11-cuts z1,z4,z5 and z1,z4,z6
        J = 3; ends = {0: ((1,1,1), (0,0,0)), 1: ((1,0,0), (0,1,1)), 2: ((0,1,0), (0,0,1))}
        cuts = [(4, 2), (7, 4), (7, 4)]
    elif CFG == "one":   # a single pair cut S12 ∋ z1,z3 (must be feasible)
        J = 1; ends = {0: ((1,), (0,)), 1: ((1,), (0,)), 2: ((0,), (0,))}; cuts = [(4, 2)]
    elif CFG == "two":   # two pair cuts S12 ∋ z1,z3 and S13 ∋ z1,z5 (observed in data; must be feasible)
        J = 2; ends = {0: ((1,1), (0,0)), 1: ((1,0), (0,0)), 2: ((0,1), (0,0))}; cuts = [(4, 2), (4, 2)]
    elif CFG == "twosame":   # two pair cuts on the same pair with the same sign: S ∋ z1,z3 ; S' ∋ z1,z3
        J = 2; ends = {0: ((1,1), (0,0)), 1: ((1,1), (0,0)), 2: ((0,0), (0,0))}; cuts = [(4, 2), (4, 2)]
    elif CFG == "ms":    # MS situation 1: same pair, opposite signs: S ∋ z1,z3 ; S' ∋ z1,z4  (must be infeasible: Lemma 6.4)
        J = 2; ends = {0: ((1,1), (0,0)), 1: ((1,0), (0,1)), 2: ((0,0), (0,0))}; cuts = [(4, 2), (4, 2)]
    elif CFG == "two_in":   # observed: S12 ∋ z1,z3 ; S13 ∋ z1,z5 with P2 entirely INSIDE S13 (must be feasible)
        J = 2; ends = {0: ((1,1), (0,0)), 1: ((1,1), (0,1)), 2: ((0,1), (0,0))}; cuts = [(4, 2), (4, 2)]
    elif CFG.startswith("sit2_"):   # three pair cuts, triangle z1z3z5; placement bits: P1 wrt S23, P2 wrt S13, P3 wrt S12
        p1, p2, p3 = (int(ch) for ch in CFG[5:])
        J = 3; ends = {0: ((1,1,p1), (0,0,p1)), 1: ((1,p2,1), (0,p2,0)), 2: ((p3,1,1), (p3,0,0))}; cuts = [(4, 2)] * 3
    else: raise SystemExit("unknown config")
    NOCIRC = os.environ.get("NOCIRC") == "1"; NOUNION = os.environ.get("NOUNION") == "1"
    print(CFG, solve(ends, cuts, UMAX, TIME, NOCIRC=NOCIRC, NOUNION=NOUNION))
