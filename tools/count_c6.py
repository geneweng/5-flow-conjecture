"""Count relaxation of the atom model, generalized for (C6) (notes §18.1): configurations of tight 6- and
7-cuts (and 11-cuts) INCLUDING the impure cases
  * even circuits of H crossing the cuts (s > 0): variables u[a,b,c] = number of colour-c edges of even
    H-circuits from atom a (black end) to atom b (white end); by Lemma 6.1 the S_j-end of every crossing
    edge is black for every cut j it crosses, by Lemma 6.2 each cut has as many colour-1 as colour-2 such
    edges, and in every atom the segments pair up so that B1 - B2 = W1 - W2;
  * non-separated paths crossing a pair cut (codimension-3 pair cuts): mode 'cross' with beta = 1 iff the
    start of the path is black under the colouring killed by the cut; an even number >= 2 of crossings;
  * cuts of type (6,4,1,3): one boundary edge of colour 0 or 3 (an F2-circuit edge outside H), variable
    y_j[k,a,b] (circuit class k, atom pair a,b), shared between two such cuts if it crosses both.
Every constraint is a necessary condition, so INFEASIBLE is a proof that the configuration cannot occur in a
cyclically 6-edge-connected cubic graph (n >= 42, canonical colouring of a 2-factor with six odd circuits).

cut = dict(c1=, c2=, nonH=0/1, cross={path: beta})   (paths not separated and not in 'cross' avoid the cut)
PURE=True forbids even-circuit crossings (reproduces count_general.py)."""
import sys, os, itertools, time
import networkx as nx
from ortools.sat.python import cp_model

def solve(ends, cuts, UMAX=2, TIME=300, PURE=False, NOCIRC=False, NOUNION=False, workers=4, rounds_max=400, verbose=False):
    J = len(cuts)
    cells = list(itertools.product((0, 1), repeat=J)); A = range(len(cells)); nc = len(cells)
    zc = {a: sum(1 for i in ends for c in ends[i] if c == cells[a]) for a in A}
    sep = {j: [i for i in ends if ends[i][0][j] != ends[i][1][j]] for j in range(J)}
    crs = {j: dict(cuts[j].get("cross", {})) for j in range(J)}
    for j in range(J): assert not set(crs[j]) & set(sep[j])
    def diff(a, b): return [j for j in range(J) if cells[a][j] != cells[b][j]]
    def forced_colour(i, a, b):
        col = None
        for j in diff(a, b):
            if i in sep[j]: beta = 1 if ends[i][0][j] == 1 else 0
            elif i in crs[j]: beta = crs[j][i]
            else: return None
            c = 1 if (cells[a][j] == 1) == (beta == 1) else 2
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
    # even circuits of H: a = black end, b = white end; needs a_j = 1, b_j = 0 on every cut crossed
    u = {}
    if not PURE:
        for a in A:
            for b in A:
                if a != b and all(cells[a][j] == 1 for j in diff(a, b)):
                    for c in (1, 2): u[a, b, c] = m.NewIntVar(0, 5, f"u{a}_{b}_{c}")
    # non-H edges of the (6,4,1,3) cuts
    NH = [j for j in range(J) if cuts[j].get("nonH")]
    K = 2 * len(ends) + 1
    upairs = [(a, b) for a in A for b in A if a < b]
    y = {}; e = {}
    for (a, b) in upairs:
        D = diff(a, b)
        if D and set(D) <= set(NH):
            for k in range(K):
                e[k, a, b] = m.NewBoolVar(f"e{k}_{a}_{b}")
    for j in NH: m.Add(sum(v for (k, a, b), v in e.items() if j in diff(a, b)) == 1)
    def out(i, a): return [x[i, a, b][0] for b in A if (i, a, b) in x]
    def inn(i, a): return [x[i, b, a][0] for b in A if (i, b, a) in x]
    for i in ends:
        sa = cells.index(ends[i][0]); ea = cells.index(ends[i][1])
        for a in A: m.Add(sum(out(i, a)) - sum(inn(i, a)) == (1 if a == sa else 0) - (1 if a == ea else 0))
    for j in range(J):
        c1 = []; c2 = []
        for i in list(sep[j]) + list(crs[j]):
            cr = [(v, c) for (ii, a, b), (v, c) in x.items() if ii == i and cells[a][j] != cells[b][j]]
            r = m.NewIntVar(1, 11, ""); h = m.NewIntVar(0, 5, "")
            if i in sep[j]: m.Add(r == 2 * h + 1)
            else: m.Add(r == 2 * h); m.Add(h >= 1)
            m.Add(r == sum(v for v, c in cr))
            c1 += [v for v, c in cr if c == 1]; c2 += [v for v, c in cr if c == 2]
        s1 = [v for (a, b, c), v in u.items() if c == 1 and j in diff(a, b)]
        s2 = [v for (a, b, c), v in u.items() if c == 2 and j in diff(a, b)]
        m.Add(sum(s1) == sum(s2))
        m.Add(sum(c1) + sum(s1) == cuts[j]["c1"]); m.Add(sum(c2) + sum(s2) == cuts[j]["c2"])
    n = {}
    for a in A:
        B = {c: sum(v for (p, q, cc), v in u.items() if p == a and cc == c) for c in (1, 2)}
        W = {c: sum(v for (p, q, cc), v in u.items() if q == a and cc == c) for c in (1, 2)}
        if u: m.Add(B[1] - B[2] == W[1] - W[2])
        ecinc = B[1] + B[2] + W[1] + W[2]
        nhinc = sum(v for (k, p, q), v in e.items() if a in (p, q))
        pvis = sum(v for i in ends for v in inn(i, a)) + sum(1 for i in ends if cells.index(ends[i][0]) == a)
        n[a] = m.NewIntVar(0, 500, f"n{a}")
        m.Add(2 * n[a] >= 2 * pvis + ecinc)
        inc1 = sum(v for (i, p, q), (v, c) in x.items() if (p == a or q == a) and c == 1) + B[1] + W[1]
        inc2 = sum(v for (i, p, q), (v, c) in x.items() if (p == a or q == a) and c == 2) + B[2] + W[2]
        h1 = m.NewIntVar(0, 300, ""); m.Add(n[a] - inc1 == 2 * h1)                   # matching inside the atom
        h2 = m.NewIntVar(0, 300, ""); m.Add(n[a] - inc2 - zc[a] == 2 * h2)           # colour-2 edges inside the atom
        h3 = m.NewIntVar(0, 300, ""); m.Add(inc2 + nhinc == 2 * h3)                  # F2-circuits cross evenly
        m.Add(inc2 + zc[a] <= n[a]); m.Add(inc1 <= n[a])
        if zc[a]: m.Add(n[a] + nhinc >= 3)
        allinc = inc1 + inc2 + nhinc
        if not any(cells.index(ends[i][0]) == a for i in ends):
            nov = m.NewBoolVar(""); m.Add(allinc == 0).OnlyEnforceIf(nov); m.Add(allinc >= 1).OnlyEnforceIf(nov.Not())
            big = m.NewBoolVar(""); m.Add(n[a] >= 6).OnlyEnforceIf([nov, big]); m.Add(n[a] == 0).OnlyEnforceIf([nov, big.Not()])
    m.Add(sum(n[a] for a in A) >= 42)
    # F2-circuits as closed walks: colour-2 crossings and non-H edges partitioned into classes
    zlist = [(i, ee) for i in ends for ee in (0, 1)]
    zatom = {k: cells.index(ends[i][ee]) for k, (i, ee) in enumerate(zlist)}
    t2 = {(a, b): sum(v for (i, p, q), (v, c) in x.items() if c == 2 and {p, q} == {a, b})
                  + sum(v for (p, q, c), v in u.items() if c == 2 and {p, q} == {a, b}) for (a, b) in upairs}
    w = {(k, a, b): m.NewIntVar(0, 11, "") for k in range(K) for (a, b) in upairs}
    if not NOCIRC:
        for (a, b) in upairs: m.Add(sum(w[k, a, b] for k in range(K)) == t2[a, b])
    deg = {(k, a): sum(w[k, p, q] for (p, q) in upairs if a in (p, q)) + sum(v for (kk, p, q), v in e.items() if kk == k and a in (p, q))
           for k in range(K) for a in A}
    if not NOCIRC:
        for k in range(K):
            for a in A: hh = m.NewIntVar(0, 30, ""); m.Add(deg[k, a] == 2 * hh)
    for k, za in (zatom.items() if not NOCIRC else []):
        small = m.NewBoolVar(""); m.Add(n[za] <= 6).OnlyEnforceIf(small); m.Add(n[za] >= 7).OnlyEnforceIf(small.Not()); m.Add(deg[k, za] >= 2).OnlyEnforceIf(small)
        nhinc = sum(v for (kk, p, q), v in e.items() if za in (p, q))
        three = m.NewBoolVar(""); m.Add(n[za] == 3).OnlyEnforceIf(three); m.Add(n[za] != 3).OnlyEnforceIf(three.Not())
        clean = m.NewBoolVar(""); m.Add(nhinc == 0).OnlyEnforceIf(clean); m.Add(nhinc >= 1).OnlyEnforceIf(clean.Not())
        for kk in range(K):
            if kk != k: m.Add(deg[kk, za] == 0).OnlyEnforceIf([three, clean])
        m.Add(deg[k, za] == 2).OnlyEnforceIf([three, clean])
    def trans():
        """all (a, b, var) inter-atom edge variables"""
        return [(p, q, v) for (i, p, q), (v, c) in x.items()] + [(p, q, v) for (p, q, c), v in u.items()] + [(p, q, v) for (k, p, q), v in e.items()]
    TR = trans()
    def add_union(U):
        Us = set(U); dU = sum(v for p, q, v in TR if (p in Us) != (q in Us))
        d = m.NewIntVar(0, 300, ""); m.Add(d == dU); nn = m.NewIntVar(0, 8000, ""); m.Add(nn == sum(n[a] for a in U)); cn = m.NewIntVar(0, 8000, ""); m.Add(cn == sum(n[a] for a in A if a not in Us))
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
    def violated(tr, sol_n):
        N = sum(sol_n); bad = []
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
    rounds = 0; t0 = time.time()
    while rounds < rounds_max and time.time() - t0 < 6 * TIME:
        st = solver.Solve(m); rounds += 1
        if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE): return solver.StatusName(st), None
        sol_x = {k: solver.Value(v) for k, (v, c) in x.items()}; sol_n = [solver.Value(n[a]) for a in A]
        tr = [(1 << p, 1 << q, solver.Value(v)) for p, q, v in TR]; tr = [t for t in tr if t[2]]
        bad = [] if NOUNION else violated(tr, sol_n)
        sub_added = 0
        for i in ends:
            Gp = nx.Graph(); Gp.add_edges_from((p, q) for (ii, p, q), v in sol_x.items() if ii == i and v)
            sa = cells.index(ends[i][0])
            for comp in nx.connected_components(Gp):
                if sa in comp: continue
                Wc = set(comp); inside = [v for (ii, p, q), (v, c) in x.items() if ii == i and p in Wc and q in Wc]
                bnd = [v for (ii, p, q), (v, c) in x.items() if ii == i and (p in Wc) != (q in Wc)]
                if not inside: continue
                use = m.NewBoolVar(""); m.Add(sum(inside) == 0).OnlyEnforceIf(use.Not()); m.Add(sum(bnd) >= 1).OnlyEnforceIf(use); sub_added += 1
        for k, za in (zatom.items() if not NOCIRC else []):
            Gc = nx.Graph(); Gc.add_edges_from((p, q) for (p, q) in upairs if solver.Value(w[k, p, q]) or ((k, p, q) in e and solver.Value(e[k, p, q])))
            for comp in nx.connected_components(Gc):
                if za in comp: continue
                Wc = set(comp); inside = [w[k, p, q] for (p, q) in upairs if p in Wc and q in Wc] + [v for (kk, p, q), v in e.items() if kk == k and p in Wc and q in Wc]
                bnd = [w[k, p, q] for (p, q) in upairs if (p in Wc) != (q in Wc)] + [v for (kk, p, q), v in e.items() if kk == k and (p in Wc) != (q in Wc)]
                use = m.NewBoolVar(""); m.Add(sum(inside) == 0).OnlyEnforceIf(use.Not()); m.Add(sum(bnd) >= 1).OnlyEnforceIf(use); sub_added += 1
        if sub_added and not bad: continue
        if not bad:
            sol = {''.join(map(str, cells[a])): sol_n[a] for a in A if sol_n[a]}
            sol["s"] = sum(solver.Value(v) for v in u.values()); return "FEASIBLE", sol
        bad.sort(key=lambda t: (bin(t[0]).count("1"), t[2]))
        if verbose: print(f"  round {rounds}: {len(bad)} violated unions, {sub_added} walk cuts, n={sum(sol_n)} [{time.time()-t0:.0f}s]", flush=True)
        for mask, nUv, d in bad[:300]: add_union([a for a in A if mask >> a & 1])
    return "UNKNOWN", None
