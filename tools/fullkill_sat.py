"""Exact search (CP-SAT) for a failure of the path-recolouring template on a fixed planted 2-factor.
Variables: perfect matching M (girth >= 6 by clauses), 0-edge per odd circuit, base 2-colouring x of
H, path labels L in {0,1,2,3} with monotone potentials (sound: nonzero labels are exactly the three
z-z paths), and for each of the 2^t path colourings a set S with 5(b_S - a_S) >= 3 d(S) + 1.
Objective: maximize the number of killed colourings; 2^t = template failure.  Cyclic 6-edge-
connectivity is added lazily (CEGAR) as linear cut constraints on M.
usage: fullkill_sat.py PROFILE [TIME_S] [WORKERS]     profiles: lengths lists, e.g. 6x7, 6x7+8"""
import sys, time, itertools, json
import networkx as nx
from ortools.sat.python import cp_model
from gen_cyc6 import has_small_cyclic_cut
prof = sys.argv[1]; TIME = float(sys.argv[2]) if len(sys.argv) > 2 else 600; WORKERS = int(sys.argv[3]) if len(sys.argv) > 3 else 8
HINT = sys.argv[4] if len(sys.argv) > 4 else None
FIX = __import__("os").environ.get("FIX") == "1"      # FIX=1: impose the hint as hard constraints (model validation)
def hint(var, val):
    if FIX: m.Add(var == val)
    else: m.AddHint(var, val)
lengths = []
for part in prof.split("+"):
    if "x" in part: k, L = part.split("x"); lengths += [int(L)] * int(k)
    else: lengths.append(int(part))
circuits = []; idx = 0
for L in lengths: circuits.append(list(range(idx, idx + L))); idx += L
n = idx; odd = [C for C in circuits if len(C) % 2]; t = len(odd) // 2
circ_of = {v: ci for ci, C in enumerate(circuits) for v in C}; pos_of = {v: i for C in circuits for i, v in enumerate(C)}
def cdist(u, v):
    if circ_of[u] != circ_of[v]: return 99
    L = len(circuits[circ_of[u]]); d = abs(pos_of[u] - pos_of[v]); return min(d, L - d)
cedges = [(C[i], C[(i + 1) % len(C)]) for C in circuits for i in range(len(C))]
m = cp_model.CpModel()
# ---- matching
pairs = [(u, v) for u in range(n) for v in range(u + 1, n) if cdist(u, v) >= 5]
M = {p: m.NewBoolVar(f"m{p}") for p in pairs}
inc = {v: [] for v in range(n)}
for (u, v), var in M.items(): inc[u].append(var); inc[v].append(var)
for v in range(n): m.AddExactlyOne(inc[v])
# girth >= 6: two matching edges (u,v),(x,y) with cdist(u,x)+cdist(v,y) <= 3 close a cycle of length <= 5
near = {v: [w for w in range(n) if w != v and cdist(v, w) <= 2] for v in range(n)}
cnt = 0
for (u, v), var in M.items():
    for x in near[u]:
        for y in near[v]:
            if x == y or cdist(u, x) + cdist(v, y) > 3: continue
            q = (min(x, y), max(x, y))
            if q in M and q > (u, v): m.AddBoolOr([var.Not(), M[q].Not()]); cnt += 1
# ---- 0-edges, colour-2 edges, z's
c2 = {}; zvar = {v: None for v in range(n)}
for C in circuits:
    L = len(C); es = [(C[i], C[(i + 1) % L]) for i in range(L)]
    if L % 2 == 0:
        for i, e in enumerate(es): c2[e] = m.NewConstant(1 if i % 2 == 0 else 0)
    else:
        o = [m.NewBoolVar(f"o{C[0]}_{i}") for i in range(L)]; m.AddExactlyOne(o)
        for i, e in enumerate(es):
            c2[e] = m.NewBoolVar(f"c2{e}"); m.Add(c2[e] == sum(o[j] for j in range(L) if (i - j) % L % 2 == 1))
        for i in range(L): zvar[C[i]] = o[i]
for v in range(n):
    if zvar[v] is None: zvar[v] = m.NewConstant(0)
# H-edges: list of (u, v, indicator)
hedges = [(u, v, c2[(u, v)]) for (u, v) in cedges] + [(u, v, var) for (u, v), var in M.items()]
# ---- base colouring x, proper on H
x = [m.NewBoolVar(f"x{v}") for v in range(n)]
for u, v, h in hedges: m.Add(x[u] != x[v]).OnlyEnforceIf(h)
# ---- labels (one-hot over 0..t) and potentials
lab = [[m.NewBoolVar(f"lab{v}_{l}") for l in range(t + 1)] for v in range(n)]
for v in range(n): m.AddExactlyOne(lab[v])
for u, v, h in hedges:
    for l in range(t + 1): m.Add(lab[u][l] == lab[v][l]).OnlyEnforceIf(h)
for v in range(n): m.AddImplication(zvar[v], lab[v][0].Not())
for l in range(1, t + 1):
    zl = []
    for v in range(n):
        if zvar[v].Name() == "": continue          # constant 0
        a = m.NewBoolVar(f"zl{v}_{l}"); m.AddMultiplicationEquality(a, [zvar[v], lab[v][l]]); zl.append(a)
    m.Add(sum(zl) == 2)
p = [m.NewIntVar(0, n, f"p{v}") for v in range(n)]
lo = {v: [] for v in range(n)}; hi = {v: [] for v in range(n)}
for u, v, h in hedges:
    for a, b in ((u, v), (v, u)):
        L_ = m.NewBoolVar(""); H_ = m.NewBoolVar("")
        m.AddImplication(L_, h); m.AddImplication(H_, h)
        m.Add(p[b] == p[a] - 1).OnlyEnforceIf(L_); m.Add(p[b] == p[a] + 1).OnlyEnforceIf(H_)
        # an H-edge at a nonzero-labelled vertex is a +-1 step
        m.AddBoolOr([h.Not(), lab[a][0], L_, H_])
        lo[a].append(L_); hi[a].append(H_)
for v in range(n):
    nzv = lab[v][0].Not()
    # internal path vertex: one lower, one higher; z: exactly one step in total
    m.Add(sum(lo[v]) == 1).OnlyEnforceIf([nzv, zvar[v].Not()]); m.Add(sum(hi[v]) == 1).OnlyEnforceIf([nzv, zvar[v].Not()])
    m.Add(sum(lo[v]) + sum(hi[v]) == 1).OnlyEnforceIf([nzv, zvar[v]])
# ---- the 2^t colourings and their bad cuts
kills = []; Svars = {}
for eps in itertools.product((0, 1), repeat=t):
    kill = m.NewBoolVar(f"kill{eps}"); kills.append(kill)
    black = []
    for v in range(n):
        sel = m.NewBoolVar(""); m.Add(sel == sum(lab[v][l + 1] for l in range(t) if eps[l]))
        b = m.NewBoolVar(""); m.AddBoolXOr([x[v], sel, b.Not()]); black.append(b)
    S = [m.NewBoolVar(f"S{eps}_{v}") for v in range(n)]; Svars[eps] = (S, kill)
    cut_terms = []
    for u, v, h in hedges + [(u, v, None) for (u, v) in cedges if (u, v) in c2 and c2[(u, v)].Name() == ""]:
        pass
    # circuit edges (always present) and matching pairs (present iff M)
    cut1 = []; cut2 = []; cut03 = []
    for (u, v) in cedges:
        xr = m.NewBoolVar(""); m.AddBoolXOr([S[u], S[v], xr.Not()]); cut_terms.append(xr)
        c2e = c2[(u, v)]
        # crossing colour-2 edge: the S-end is black;  crossing colour-0/3 edge counted separately
        is2 = m.NewBoolVar(""); m.AddMultiplicationEquality(is2, [xr, c2e]); cut2.append(is2)
        is03 = m.NewBoolVar(""); m.Add(is03 == xr - is2); cut03.append(is03)
        m.AddBoolOr([is2.Not(), S[u].Not(), black[u]]); m.AddBoolOr([is2.Not(), S[v].Not(), black[v]])
    for (u, v), var in M.items():
        xr = m.NewBoolVar(""); m.AddBoolXOr([S[u], S[v], xr.Not()])
        c = m.NewBoolVar(""); m.AddImplication(c, var); m.AddImplication(c, xr); m.AddBoolOr([var.Not(), xr.Not(), c]); cut_terms.append(c); cut1.append(c)
        m.AddBoolOr([c.Not(), S[u].Not(), black[u]]); m.AddBoolOr([c.Not(), S[v].Not(), black[v]])   # crossing matching edge: S-end black
    m.Add(sum(cut1) >= 4).OnlyEnforceIf(kill); m.Add(sum(cut1) <= 3 * t - 2).OnlyEnforceIf(kill)
    m.Add(sum(cut2) >= 1).OnlyEnforceIf(kill); m.Add(sum(cut2) <= 2 * t - 2).OnlyEnforceIf(kill)
    m.Add(sum(cut03) <= 1).OnlyEnforceIf(kill)
    bs = []
    for v in range(n):
        a = m.NewBoolVar(""); m.AddMultiplicationEquality(a, [S[v], black[v]]); bs.append(a)
    imbalance = sum(2 * a for a in bs) - sum(S)
    m.Add(5 * imbalance - 3 * sum(cut_terms) >= 1).OnlyEnforceIf(kill)
    # structure of tight cuts (Lemma 6.1): at most 5t-4 cut edges, imbalance at least 4, at least 6 vertices on each side
    m.Add(sum(cut_terms) <= 5 * t - 4).OnlyEnforceIf(kill); m.Add(imbalance >= 4).OnlyEnforceIf(kill)
    m.Add(sum(S) >= 6).OnlyEnforceIf(kill); m.Add(sum(S) <= n - 6).OnlyEnforceIf(kill)
m.Maximize(sum(kills))
if HINT:   # warm start from a dumped instance (killsets/climb json): map its circuits onto ours in order
    from oddness import canonical_colouring, H_components, partition
    d = json.load(open(HINT)); Gh = nx.Graph(); Gh.add_edges_from(map(tuple, d["edges"])); hc = d["circuits"]
    assert [len(C) for C in hc] == [len(C) for C in circuits], "hint profile mismatch"
    vmap = {hv: v for hC, C in zip(hc, circuits) for hv, v in zip(hC, C)}
    Mh = {frozenset(e) for e in Gh.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in hc for i in range(len(C))}
    col = canonical_colouring(Gh, Mh, hc, d["zero"]); comps = H_components(Gh, col)
    ref = d["ref"]; pidx = [i for i, (_, e) in enumerate(comps) if e]
    if isinstance(ref, dict): f = [ref.get(str(min(c)), 0) for (c, _) in comps]
    else:   # old list dumps: component order is not reproducible; recover the colouring from the recorded killed set
        from badcuts import violating_set
        want = {tuple(int(ch) for ch in k if ch in "01") for k in d["S"]}; f = None
        for cand in itertools.product((0, 1), repeat=len(comps)):
            if any(cand[i] for i in pidx): continue
            got = set()
            for bits in itertools.product((0, 1), repeat=len(pidx)):
                ff = list(cand)
                for i, b in zip(pidx, bits): ff[i] = b
                if violating_set(Gh, partition(comps, ff)) is not None: got.add(bits)
            if got == want: f = list(cand); break
        assert f is not None, "could not recover the hint colouring"
    black = partition(comps, f)
    labh = {v: (pidx.index(i) + 1 if i in pidx else 0) for i, (c, _) in enumerate(comps) for v in c}
    mset = {frozenset((vmap[a], vmap[b])) for e in Mh for a, b in [tuple(e)]}
    for (u, v), var in M.items(): hint(var, 1 if frozenset((u, v)) in mset else 0)
    for hv, v in vmap.items(): hint(x[v], 1 if hv in black else 0)
    for hv, v in vmap.items():
        for l in range(t + 1): hint(lab[v][l], 1 if labh[hv] == l else 0)
    hodd = [C for C in hc if len(C) % 2]
    for hC, zi in zip(hodd, d["zero"]): hint(zvar[vmap[hC[zi]]], 1)
    print("hint loaded from", HINT, flush=True)
    if __import__("os").environ.get("FIXS") == "1":     # debug: also fix the witness sets of the dump
        for k, Sl in d["S"].items():
            eps = tuple(int(ch) for ch in k if ch in "01"); Sm = {vmap[v] for v in Sl}
            S, kill = Svars[eps]; m.Add(kill == 1)
            blk = {vmap[hv] for hv in Gh.nodes() if (hv in black) ^ (labh[hv] > 0 and eps[labh[hv] - 1] == 1)}
            b = sum(1 for v in Sm if v in blk); a = len(Sm) - b
            if b < a: Sm = set(range(n)) - Sm          # the dump keeps the smaller side; we need the black-heavy side
            for v in range(n): m.Add(S[v] == (1 if v in Sm else 0))
        print("witness sets fixed", flush=True)
print(f"{prof}: n={n}, t={t}, matching vars {len(M)}, girth clauses {cnt}", flush=True)
class CB(cp_model.CpSolverSolutionCallback):
    def __init__(s): super().__init__(); s.t0 = time.time()
    def on_solution_callback(s): print(f"  [{time.time()-s.t0:.0f}s] solution with {int(s.ObjectiveValue())} of {2**t} colourings killed", flush=True)
extra_cuts = 0
while True:
    solver = cp_model.CpSolver(); solver.parameters.max_time_in_seconds = TIME; solver.parameters.num_workers = WORKERS
    st = solver.Solve(m, CB()); name = solver.StatusName(st)
    print(f"status {name}, objective {solver.ObjectiveValue() if st in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None}, bound {solver.BestObjectiveBound()}", flush=True)
    if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE): break
    G = nx.Graph(); G.add_edges_from(cedges); G.add_edges_from(pr for pr, var in M.items() if solver.Value(var))
    # external cross-check of the solution with the reference tools
    from badcuts import violating_set
    xs = [solver.Value(b) for b in x]; labs = [[solver.Value(l) for l in lab[v]].index(1) for v in range(n)]
    ext = 0
    for eps in itertools.product((0, 1), repeat=t):
        blk = {v for v in range(n) if xs[v] ^ (1 if labs[v] and eps[labs[v] - 1] else 0)}
        if violating_set(G, blk) is not None: ext += 1
    print(f"  external check: girth {nx.girth(G)}, labels used {sorted(set(labs))}, z count {sum(solver.Value(zvar[v]) for v in range(n))}, killed colourings by reference tools = {ext}", flush=True)
    if solver.ObjectiveValue() < 2 ** t: break
    X = has_small_cyclic_cut(G)
    if X is None:
        out = {"profile": prof, "edges": [list(e) for e in G.edges()], "circuits": circuits,
               "zero": [[solver.Value(zvar[v]) for v in C].index(1) for C in odd], "x": [solver.Value(b) for b in x],
               "labels": [[solver.Value(l) for l in lab[v]].index(1) for v in range(n)]}
        json.dump(out, open(f"fullkill_sat_{prof}.json", "w")); print("FULL KILL, cyclically 6-connected ->", f"fullkill_sat_{prof}.json", flush=True); break
    Xs = set(X); const = sum(1 for u, v in cedges if (u in Xs) != (v in Xs))
    m.Add(const + sum(var for (u, v), var in M.items() if (u in Xs) != (v in Xs)) >= 6); extra_cuts += 1
    print(f"  small cyclic cut of size {nx.cut_size(G, Xs)} found (|X|={len(Xs)}); added constraint #{extra_cuts}", flush=True)
