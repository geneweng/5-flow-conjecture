"""Atom-count relaxation for the 3x2 grid of a two-position shift (notes S16).  Cuts: tight 11-cuts S_{x,y} for the
grid colourings, x in rows (a0 | E | E1), y in columns (F | F').  Pieces P (E->a0), pi1 (a1->E1), pi2 (E2->a2), R (F->F')
are walks through the 2^J atoms; crossing colours are forced by the piece's start colour under the cut's colouring
(black start: exits are matching edges (1), entries visits (2); white start: the reverse); glue edges a1a2 (in H(c))
and a0a1 (in H(c')) are added explicitly.  Every constraint is a necessary condition, so INFEASIBLE proves that the
chosen rows cannot all be fully killed by 11-cuts.   usage: grid_count.py ROWS [UMAX] [TIME]   e.g. ROWS=a0,E1 or a0,E,E1"""
import sys, os, itertools, time, random
import networkx as nx
from ortools.sat.python import cp_model
def solve(rows, UMAX=2, TIME=600, workers=8, verbose=True):
    cuts = [(r, c) for r in rows for c in ("F", "Fp")]; J = len(cuts)
    cells = list(itertools.product((0, 1), repeat=J)); A = range(len(cells)); nc = len(cells); idx = {c: a for a, c in enumerate(cells)}
    def vec(f): return tuple(f(r, c) for r, c in cuts)
    FREE = None
    # fixed memberships of the special vertices
    mem = {"E": vec(lambda r, c: int(r == "E")), "E1": vec(lambda r, c: int(r == "E1")), "E2": vec(lambda r, c: 1),
           "F": vec(lambda r, c: int(c == "F")), "Fp": vec(lambda r, c: int(c == "Fp")),
           "a1": vec(lambda r, c: int(r == "a0")),
           "a0": vec(lambda r, c: 1 if r == "a0" else 0 if r == "E" else FREE),
           "a2": vec(lambda r, c: FREE if r == "a0" else 0)}
    # pieces: (name, start, end, black_start per cut)
    pieces = {"P": ("E", "a0", vec(lambda r, c: r == "E")), "pi1": ("a1", "E1", vec(lambda r, c: r in ("a0", "E"))),
              "pi2": ("E2", "a2", vec(lambda r, c: True)), "R": ("F", "Fp", vec(lambda r, c: c == "F"))}
    m = cp_model.CpModel()
    # atom selection for vertices with free coordinates
    def atoms_matching(v): return [a for a in A if all(f is None or f == cells[a][j] for j, f in enumerate(mem[v]))]
    sel = {}
    for v in ("a0", "a2"):
        cand = atoms_matching(v); sel[v] = {a: m.NewBoolVar(f"sel_{v}_{a}") for a in cand}; m.AddExactlyOne(sel[v].values())
    def in_cut(v, j):
        """0/1 or a BoolVar-like linear expr: is special vertex v inside cut j"""
        f = mem[v][j]
        if f is not None: return f
        return sum(b for a, b in sel[v].items() if cells[a][j] == 1)
    def forced_colour(i, a, b):
        col = None
        for j in range(J):
            if cells[a][j] == cells[b][j]: continue
            bs = pieces[i][2][j]; leaving = cells[a][j] == 1
            c = (1 if bs else 2) if leaving else (2 if bs else 1)
            if col is None: col = c
            elif col != c: return None
        return col
    x = {}
    for i in pieces:
        for a in A:
            for b in A:
                if a == b: continue
                c = forced_colour(i, a, b)
                if c is not None: x[i, a, b] = (m.NewIntVar(0, 11, f"x_{i}_{a}_{b}"), c)
    def out(i, a): return [x[i, a, b][0] for b in A if (i, a, b) in x]
    def inn(i, a): return [x[i, b, a][0] for b in A if (i, b, a) in x]
    # flow conservation per piece: start atom fixed or selected, end atom fixed or selected
    def src_ind(v, a):
        if v in sel: return sel[v].get(a, 0)
        return 1 if cells[a] == mem[v] else 0
    for i, (s, e, bs) in pieces.items():
        for a in A: m.Add(sum(out(i, a)) - sum(inn(i, a)) == src_ind(s, a) - src_ind(e, a))
    # glue crossings: g12[j] = a1a2 crosses cut j (rows a0 only: a1 in, a2 out); g01[j] = a0a1 crosses cut j (rows E1 only: a0 in)
    g12 = {}; g01 = {}
    for j, (r, c) in enumerate(cuts):
        g12[j] = (1 - in_cut("a2", j)) if r == "a0" else 0
        g01[j] = in_cut("a0", j) if r == "E1" else 0
    # cut counts: 7 matching + 4 visit crossings, black S-end (colours already forced); per piece odd/even crossing parity
    for j in range(J):
        c1 = []; c2 = []
        for i, (s, e, bs) in pieces.items():
            cr = [(v, c) for (ii, a, b), (v, c) in x.items() if ii == i and cells[a][j] != cells[b][j]]
            c1 += [v for v, c in cr if c == 1]; c2 += [v for v, c in cr if c == 2]
            # crossings ≡ start_in + end_in (mod 2): encode with an auxiliary integer k: crossings = start_in + end_in + 2k - 2*start_in*end_in... simpler: crossings + start_in + end_in even
            h = m.NewIntVar(0, 8, ""); m.Add(sum(v for v, c in cr) + in_cut(s, j) + in_cut(e, j) == 2 * h)
        m.Add(sum(c1) == 7); m.Add(sum(c2) + g12[j] + g01[j] == 4)
    # atom vertex counts
    n = {a: m.NewIntVar(0, 500, f"n{a}") for a in A}
    zc_c = {a: sum(1 for v in ("E", "E1", "E2", "F", "Fp") if cells[a] == mem[v]) for a in A}     # fixed z's
    def a_in(v, a):          # indicator that special vertex v lies in atom a
        return src_ind(v, a)
    for a in A:
        visits = sum(v for i in pieces for v in inn(i, a)) + sum(src_ind(s, a) for i, (s, e, bs) in pieces.items())
        m.Add(n[a] >= visits)
        inc1 = [v for (i, p, q), (v, c) in x.items() if (p == a or q == a) and c == 1]
        inc2 = [v for (i, p, q), (v, c) in x.items() if (p == a or q == a) and c == 2]
        # glue crossings incident to atom a: a1a2 crossing leaves a1's atom and enters a2's atom
        # a1's atom is fixed; a1a2 crosses cut j iff a2 outside; it is incident to atom(a1) and atom(a2) when they differ
        a1_atom = idx[mem["a1"]]
        diff_a2 = 1 - sel["a2"].get(a1_atom, 0) if a == a1_atom else sel["a2"].get(a, 0) if a != a1_atom else 0
        # a0a1 crossing: incident to atom(a0) and atom(a1) when they differ
        a0_atom_here = sel["a0"].get(a, 0)
        diff_a0 = (1 - sel["a0"].get(a1_atom, 0)) if a == a1_atom else (a0_atom_here if a != a1_atom else 0)
        # matching parity: inc1 ≡ n[a]
        d1 = m.NewIntVar(-250, 250, ""); m.Add(n[a] - sum(inc1) == 2 * d1)
        # visit parity at c: n - zc_c - [a0 here] == 2*internal + inc2 + [a1a2 crossing incident here]
        zc_here_c = zc_c[a] + a_in("a0", a)
        d2 = m.NewIntVar(-250, 250, ""); m.Add(n[a] - zc_here_c - sum(inc2) - diff_a2 == 2 * d2)
        # visit parity at c': z's of c' are the fixed ones plus a2; glue a0a1
        zc_here_cp = zc_c[a] + a_in("a2", a)
        d3 = m.NewIntVar(-250, 250, ""); m.Add(n[a] - zc_here_cp - sum(inc2) - diff_a0 == 2 * d3)
        # room: colour-2 incidences (at c) + z's of c <= n
        m.Add(sum(inc2) + diff_a2 + zc_here_c <= n[a])
        # z-atoms: a fixed z has its two non-H circuit neighbours in its atom
        if zc_c[a]: m.Add(n[a] >= 3 * zc_c[a])
        # a0's atom contains a6 (edge a6a0 never crosses); a2's atom contains a3
        m.Add(n[a] >= 2).OnlyEnforceIf(sel["a0"][a]) if a in sel["a0"] else None
        m.Add(n[a] >= 2).OnlyEnforceIf(sel["a2"][a]) if a in sel["a2"] else None
    m.Add(sum(n[a] for a in A) >= 42)
    # circuits: colour-2 crossings partitioned into 6 circuit classes with even degree at every atom; C's class also
    # contains the glue crossings.  Class 0 = C, classes 1..5 = circuits of E, E1, E2, F, F' (each contains its z-atom).
    upairs = [(a, b) for a in A for b in A if a < b]
    t2 = {(a, b): sum(v for (i, p, q), (v, c) in x.items() if c == 2 and {p, q} == {a, b}) for (a, b) in upairs}
    K = 6; w = {(k, a, b): m.NewIntVar(0, 11, "") for k in range(K) for (a, b) in upairs}
    for (a, b) in upairs: m.Add(sum(w[k, a, b] for k in range(K)) == t2[a, b])
    a1_atom = idx[mem["a1"]]
    for k in range(K):
        for a in A:
            deg = sum(w[k, p, q] for (p, q) in upairs if a in (p, q))
            extra = 0
            if k == 0:
                # glue edges incident to atom a within circuit C's closed walk
                if a == a1_atom: extra = (1 - sel["a2"].get(a1_atom, 0)) + (1 - sel["a0"].get(a1_atom, 0))
                else: extra = sel["a2"].get(a, 0) + sel["a0"].get(a, 0)
            hh = m.NewIntVar(0, 30, ""); m.Add(deg + extra == 2 * hh)
    # union facts: cubic, girth >= 6 (a side with 1..5 vertices has d = n + 2 only if acyclic; d >= n+2 anyway), cyclic 6-edge-connectivity
    def add_union(U):
        Us = set(U)
        d_pieces = sum(v for (i, p, q), (v, c) in x.items() if (p in Us) != (q in Us))
        # glue crossings: a1a2 between atom(a1) and atom(a2); a0a1 between atom(a0) and atom(a1)
        a1in = a1_atom in Us
        g_a2 = sum(b for a, b in sel["a2"].items() if (a in Us) != a1in)
        g_a0 = sum(b for a, b in sel["a0"].items() if (a in Us) != a1in)
        d = m.NewIntVar(0, 300, ""); m.Add(d == d_pieces + g_a2 + g_a0)
        nn = m.NewIntVar(0, 8000, ""); m.Add(nn == sum(n[a] for a in U)); cn = m.NewIntVar(0, 8000, ""); m.Add(cn == sum(n[a] for a in A if a not in Us))
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
    for r in range(1, UMAX + 1):
        for U in itertools.combinations(A, r): add_union(U)
    # structured unions: each cut's inside, pairwise intersections / differences / unions of cuts
    sides = {j: [a for a in A if cells[a][j] == 1] for j in range(J)}
    seen = set()
    for j in range(J):
        for jj in range(J):
            if j == jj: continue
            for U in (set(sides[j]) & set(sides[jj]), set(sides[j]) - set(sides[jj]), set(sides[j]) | set(sides[jj])):
                if U and len(U) < nc and frozenset(U) not in seen: seen.add(frozenset(U)); add_union(sorted(U))
    solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = TIME; solver.parameters.num_workers = workers
    rng = random.Random(1); t0 = time.time()
    for rnd in range(40):
        st = solver.Solve(m)
        if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            return solver.StatusName(st), None
        sol_x = {k: solver.Value(v) for k, (v, c) in x.items()}; sol_n = [solver.Value(n[a]) for a in A]
        sa0 = next(a for a, b in sel["a0"].items() if solver.Value(b)); sa2 = next(a for a, b in sel["a2"].items() if solver.Value(b))
        tr = [(p, q, v) for (i, p, q), v in sol_x.items() if v]
        glue = [(a1_atom, sa2, 1 if sa2 != a1_atom else 0), (sa0, a1_atom, 1 if sa0 != a1_atom else 0)]
        N = sum(sol_n)
        def check(Us):
            nU = sum(sol_n[a] for a in Us); cn = N - nU
            d = sum(v for p, q, v in tr if (p in Us) != (q in Us)) + sum(v for p, q, v in glue if (p in Us) != (q in Us))
            ok = d <= 3 * nU and d <= 3 * cn
            for side in (nU, cn):
                if 1 <= side <= 5 and d < side + 2: ok = False
            if nU >= 1 and cn >= 1 and d < 3: ok = False
            if (nU == 0 or cn == 0) and d != 0: ok = False
            if d <= 5 and nU > 3 and cn > 3: ok = False
            return ok, nU, d
        bad = []
        occupied = [a for a in A if sol_n[a]]
        # exhaustive over unions of occupied atoms if few, else sampled + connected unions of the transition graph
        if len(occupied) <= 16:
            for r in range(1, len(occupied)):
                for U in itertools.combinations(occupied, r):
                    ok, nU, d = check(set(U))
                    if not ok: bad.append((set(U), nU, d))
        else:
            for _ in range(20000):
                r = rng.randint(1, len(occupied) - 1); U = set(rng.sample(occupied, r)); ok, nU, d = check(U)
                if not ok: bad.append((U, nU, d))
        # walk connectivity (subtour elimination) for pieces and circuit classes
        sub_added = 0
        for i, (s, e, bs) in pieces.items():
            Gp = nx.Graph(); Gp.add_edges_from((p, q) for (ii, p, q), v in sol_x.items() if ii == i and v)
            sa = idx[mem[s]] if s not in sel else next(a for a, b in sel[s].items() if solver.Value(b))
            for comp in nx.connected_components(Gp):
                if sa in comp: continue
                W = set(comp); inside = [v for (ii, p, q), (v, c) in x.items() if ii == i and p in W and q in W]
                bnd = [v for (ii, p, q), (v, c) in x.items() if ii == i and (p in W) != (q in W)]
                if not inside: continue
                use = m.NewBoolVar(""); m.Add(sum(inside) == 0).OnlyEnforceIf(use.Not()); m.Add(sum(bnd) >= 1).OnlyEnforceIf(use); sub_added += 1
        if verbose: print(f"  round {rnd}: status {solver.StatusName(st)}, N={N}, occupied atoms {len(occupied)}, violated unions {len(bad)}, subtours {sub_added} [{time.time()-t0:.0f}s]", flush=True)
        if not bad and not sub_added:
            return "FEASIBLE", {"n": {''.join(map(str, cells[a])): sol_n[a] for a in occupied}, "a0": ''.join(map(str, cells[sa0])), "a2": ''.join(map(str, cells[sa2]))}
        bad.sort(key=lambda t: (len(t[0]), t[2]))
        for U, nU, d in bad[:60]: add_union(sorted(U))
    return "UNDECIDED", None
if __name__ == "__main__":
    rows = sys.argv[1].split(","); UMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 2; TIME = float(sys.argv[3]) if len(sys.argv) > 3 else 600
    print(rows, solve(rows, UMAX, TIME))
