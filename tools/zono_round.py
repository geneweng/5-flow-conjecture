"""Experiment C (notes/oddness.md 18.3): rounding inside the zonotope.

Z(G) = { w in R^V : |w(X)| <= d(X) for all X } = boundaries of real flows with |f(e)| <= 1.
G has a circular nowhere-zero r-flow iff Z(G) contains a vector with all coordinates +-c,
c = r/(r-2).  We start at w = 0 and try to reach such a vector by

  walk : Lovett-Meka-style random walk.  Gaussian steps in the subspace orthogonal to the
         all-ones vector, to e_v for frozen v (|w_v| = c) and to 1_X for every cut X that has
         become tight.  A step that would leave the box or Z is shrunk to the boundary (the cut
         is found by max-flow, Newton/Dinkelbach iteration on the step length); the coordinate /
         cut that stopped the step is added to the constraint list (sticky, never released).
         delta = step length per coordinate (in units of c); delta = inf gives a "ray shooting"
         walk in which every step ends on a new constraint.
  lp   : iterative LP rounding in the flow formulation (f_e in [-1,1], w = boundary of f,
         |w_v| <= c): maximise a random linear objective, freeze the coordinates at +-c of the
         vertex solution, repeat.  Stuck = LP_TRIES consecutive random objectives freeze nothing.

For a stuck point the lattice of tight sets is read off the residual digraph of a feasible flow
(X is tight iff no residual arc leaves X); its strongly connected components are the links of a
maximal chain of tight sets.  CP-SAT (exact) is used for ground truth: does a +-c vector exist
in Z(G), does the frozen pattern of a stuck run extend to one, and how many frozen signs must be
flipped to reach one.

Usage:  python3 tools/zono_round.py [--restarts 50] [--out FILE.json] [--graphs 'a;b;...'] [--jobs 2]
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
import sys, json, time, argparse, math
from fractions import Fraction
from collections import Counter
import numpy as np
import networkx as nx
from ortools.graph.python import max_flow
from ortools.linear_solver import pywraplp
from ortools.sat.python import cp_model

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from oddness import petersen, flower, gp, R2

ROOT = os.path.dirname(HERE)
K = 1 << 40            # integer scaling for max-flow
FEAS_TOL = 1e-9        # w is "inside" Z if max_X w(X)-d(X) <= FEAS_TOL
TIGHT_TOL = 1e-6       # residual capacity below this = saturated (stuck-point analysis)
LP_TRIES = 10
R_GRID = [Fraction(4), Fraction(9, 2), Fraction(5), Fraction(11, 2), Fraction(6), Fraction(7), Fraction(8),
          Fraction(10), Fraction(16)]     # 10 and 16 added to see where the procedures start to work


# ---------------------------------------------------------------- graphs
def load_json_graph(path):
    d = json.load(open(path))
    return nx.Graph([tuple(e) for e in d["edges"]])

def test_graphs():
    gs = {
        "Petersen": petersen(),
        "J5": flower(5), "J7": flower(7), "J9": flower(9),
        "GP(13,3)": gp(13, 3), "GP(17,4)": gp(17, 4),
        "R2": R2(),
        "cx42": load_json_graph(os.path.join(ROOT, "data/oddness6_template_counterexample_n42.json")),
        "cx50": load_json_graph(os.path.join(ROOT, "data/oddness6_template_counterexample_n50.json")),
        "cx60": load_json_graph(os.path.join(ROOT, "data/oddness6_template_counterexample.json")),
    }
    out = {}
    for name, G in gs.items():
        G = nx.convert_node_labels_to_integers(G, ordering="sorted") if all(isinstance(v, int) for v in G) \
            else nx.convert_node_labels_to_integers(G)
        assert all(d == 3 for _, d in G.degree()), name
        out[name] = G
    return out


# ---------------------------------------------------------------- zonotope oracle
class Zono:
    def __init__(self, G):
        self.G = G
        self.n = n = G.number_of_nodes()
        self.E = np.array(sorted(tuple(sorted(e)) for e in G.edges()), dtype=np.int64)
        self.m = m = len(self.E)
        self.s, self.t = n, n + 1
        tails = np.concatenate([self.E[:, 0], self.E[:, 1], np.full(n, self.s), np.arange(n)])
        heads = np.concatenate([self.E[:, 1], self.E[:, 0], np.arange(n), np.full(n, self.t)])
        caps = np.concatenate([np.full(2 * m, K), np.zeros(2 * n)]).astype(np.int64)
        self.smf = max_flow.SimpleMaxFlow()
        self.arcs = self.smf.add_arcs_with_capacity(tails.astype(np.int32), heads.astype(np.int32), caps)
        self.src_arcs = self.arcs[2 * m:2 * m + n]
        self.snk_arcs = self.arcs[2 * m + n:]
        self.nflow = 0

    def d(self, mask):
        return int(np.sum(mask[self.E[:, 0]] != mask[self.E[:, 1]]))

    def mincut(self, w):
        """(viol, X): viol = max_X w(X) - d(X) >= 0 (0 iff w in Z, given sum w = 0), X a maximiser (bool mask)."""
        iw = np.rint(w * K).astype(np.int64)
        pos = np.maximum(iw, 0); neg = np.maximum(-iw, 0)
        self.smf.set_arcs_capacity(self.src_arcs, pos)
        self.smf.set_arcs_capacity(self.snk_arcs, neg)
        st = self.smf.solve(self.s, self.t)
        assert st == self.smf.OPTIMAL
        self.nflow += 1
        viol = (int(pos.sum()) - self.smf.optimal_flow()) / K
        X = np.zeros(self.n, dtype=bool)
        if viol > 0:
            side = [v for v in self.smf.get_source_side_min_cut() if v < self.n]
            X[side] = True
        return viol, X

    def edge_flow(self, w):
        """A flow f (|f|<=1) whose boundary is w up to the violation; f[e] > 0 means E[e,0] -> E[e,1]."""
        viol, _ = self.mincut(w)
        fl = np.array(self.smf.flows(self.arcs[:2 * self.m]), dtype=np.float64) / K
        return fl[:self.m] - fl[self.m:], viol

    def tight_lattice(self, f):
        """Residual digraph of f on V; tight sets = successor-closed sets.  Returns the SCCs in an order
        such that every suffix is successor-closed... (we return prefix order: X_i = union of first i SCCs is tight)."""
        D = nx.DiGraph(); D.add_nodes_from(range(self.n))
        for (u, v), fe in zip(self.E, f):
            if fe < 1 - TIGHT_TOL: D.add_edge(int(u), int(v))     # can push more u->v
            if fe > -1 + TIGHT_TOL: D.add_edge(int(v), int(u))
        C = nx.condensation(D)
        order = list(nx.topological_sort(C))[::-1]               # sinks first: prefixes are successor-closed
        return [sorted(C.nodes[i]["members"]) for i in order]


def c_of(r):
    r = Fraction(r)
    return r / (r - 2)


# ---------------------------------------------------------------- Lovett-Meka style walk
class Basis:
    """Orthonormal basis of the span of the active constraint normals."""
    def __init__(self, n):
        self.n = n; self.Q = np.zeros((n, 0))
    def add(self, a):
        a = np.asarray(a, dtype=float).copy()
        for _ in range(2):
            a -= self.Q @ (self.Q.T @ a)
        nrm = np.linalg.norm(a)
        if nrm > 1e-8:
            self.Q = np.hstack([self.Q, (a / nrm)[:, None]]); return True
        return False
    def rank(self): return self.Q.shape[1]
    def project(self, g):
        for _ in range(2):
            g = g - self.Q @ (self.Q.T @ g)
        return g


def walk(Z, c, rng, delta=0.1, max_steps=200000):
    n = Z.n; c = float(c)
    w = np.zeros(n); frozen = np.zeros(n, dtype=bool); cuts = []
    B = Basis(n); B.add(np.ones(n))
    steps = 0; null_steps = 0; status = None
    while not frozen.all() and B.rank() < n and steps < max_steps:
        steps += 1
        g = B.project(rng.standard_normal(n))
        g[frozen] = 0.0
        ng = np.linalg.norm(g)
        if ng < 1e-9: continue
        if math.isinf(delta):
            g /= ng; t = np.inf
        else:
            g *= delta * c; t = 1.0
        # box
        free = ~frozen
        with np.errstate(divide="ignore", invalid="ignore"):
            tb = np.where(g > 0, (c - w) / g, np.where(g < 0, (-c - w) / g, np.inf))
        tb[~free] = np.inf
        t = min(t, float(tb.min()))
        # cuts: shrink t until w + t g is in Z.  A step is only ever accepted if max-flow certifies it.
        hit = None; ok = False
        for _ in range(60):
            viol, X = Z.mincut(w + t * g)
            if viol <= FEAS_TOL: ok = True; break
            gX = g[X].sum(); slack = Z.d(X) - w[X].sum()
            hit = X
            if gX <= 1e-12: t = 0.0; break              # this cut cannot be relieved along g
            tn = slack / gX
            if tn >= t: tn = 0.5 * t                     # numerical: bisect instead of trusting the Newton step
            t = max(tn, 0.0)
            if t <= 1e-15: t = 0.0; break
        if not ok: t = 0.0                               # reject: never move to an uncertified point
        if t > 0.0:
            w = w + t * g
        newf = free & (np.abs(w) >= c - 1e-9)
        for v in np.flatnonzero(newf):
            w[v] = math.copysign(c, w[v]); frozen[v] = True
            e = np.zeros(n); e[v] = 1; B.add(e)
        progressed = t > 0.0 or newf.any()
        # register the cut only if it is really tight at the point we ended on
        if hit is not None and Z.d(hit) - w[hit].sum() <= 1e-7 and B.add(hit.astype(float)):
            cuts.append(hit); progressed = True
        if progressed:
            null_steps = 0
        else:
            null_steps += 1
            if null_steps >= 200: status = "numerical"; break
    viol, _ = Z.mincut(w)
    success = bool(frozen.all() and viol <= 1e-6)
    if status is None:
        status = "success" if success else ("stuck" if B.rank() >= n and viol <= 1e-6 else "numerical")
    return dict(w=w, frozen=frozen, cuts=cuts, steps=steps, viol=viol, success=success, rank=B.rank(), status=status)


# ---------------------------------------------------------------- iterative LP rounding
class FlowLP:
    def __init__(self, Z, c):
        self.Z = Z; self.c = float(c)
        self.S = pywraplp.Solver.CreateSolver("GLOP")
        self.S.SetNumThreads(1)
        self.f = [self.S.NumVar(-1.0, 1.0, f"f{e}") for e in range(Z.m)]
        self.inc = [[] for _ in range(Z.n)]
        for e, (u, v) in enumerate(Z.E):
            self.inc[u].append((e, 1.0)); self.inc[v].append((e, -1.0))
        self.con = []
        for v in range(Z.n):
            ct = self.S.Constraint(-self.c, self.c)
            for e, sg in self.inc[v]: ct.SetCoefficient(self.f[e], sg)
            self.con.append(ct)
    def fix(self, v, val): self.con[v].SetBounds(val, val)
    def solve(self, obj_w):
        coef = np.zeros(self.Z.m)
        for v in range(self.Z.n):
            for e, sg in self.inc[v]: coef[e] += sg * obj_w[v]
        o = self.S.Objective(); o.Clear()
        for e in range(self.Z.m): o.SetCoefficient(self.f[e], float(coef[e]))
        o.SetMaximization()
        st = self.S.Solve()
        if st != pywraplp.Solver.OPTIMAL: return None, None
        f = np.array([x.solution_value() for x in self.f])
        w = np.zeros(self.Z.n)
        np.add.at(w, self.Z.E[:, 0], f); np.subtract.at(w, self.Z.E[:, 1], f)
        return f, w


def lp_round(Z, c, rng, frozen0=None, sign0=None):
    lp = FlowLP(Z, c); n = Z.n; cf = float(c)
    frozen = np.zeros(n, dtype=bool); sign = np.zeros(n)
    if frozen0 is not None:
        frozen = frozen0.copy(); sign = np.where(frozen, np.sign(sign0), 0.0)
        for v in np.flatnonzero(frozen): lp.fix(int(v), sign[v] * cf)
    rounds = 0; w = np.zeros(n); f = np.zeros(Z.m); first_frozen = None; spread = 0.0
    while not frozen.all():
        tries = 0; ws = []
        while tries < LP_TRIES:
            tries += 1; rounds += 1
            f, w = lp.solve(rng.standard_normal(n))
            if f is None:
                return dict(w=w, frozen=frozen, success=False, rounds=rounds, error="LP not optimal")
            ws.append(w)
            newf = (~frozen) & (np.abs(w) >= cf - 1e-7)
            if newf.any(): break
        if not newf.any():
            spread = float(max(np.abs(a - ws[0]).max() for a in ws))
            break
        for v in np.flatnonzero(newf):
            frozen[v] = True; sign[v] = math.copysign(1, w[v]); lp.fix(v, sign[v] * cf)
        if first_frozen is None: first_frozen = int(frozen.sum())
    w = np.where(frozen, sign * cf, w)
    return dict(w=w, f=f, frozen=frozen, success=bool(frozen.all()), rounds=rounds,
                first_frozen=first_frozen, spread=spread)


# ---------------------------------------------------------------- exact (CP-SAT)
EXACT_R = {Fraction(5), Fraction(11, 2), Fraction(6), Fraction(8)}   # stuck runs are sent to CP-SAT only for these r

def exact(Z, r, frozen=None, sign=None, mode="exists", time_limit=5.0):
    """mode 'exists': is there a +-c vector in Z (if frozen given: agreeing with the frozen signs)?
       mode 'dist'  : min number of frozen signs to flip.  Returns (status, value)."""
    cfr = c_of(r); a, b = cfr.numerator, cfr.denominator
    M = cp_model.CpModel()
    F = [M.NewIntVar(-b, b, f"F{e}") for e in range(Z.m)]
    x = [M.NewBoolVar(f"x{v}") for v in range(Z.n)]
    inc = [[] for _ in range(Z.n)]
    for e, (u, v) in enumerate(Z.E):
        inc[u].append(F[e]); inc[v].append(-F[e])
    for v in range(Z.n):
        M.Add(sum(inc[v]) == a * (2 * x[v] - 1))
    if frozen is None: M.Add(x[0] == 1)       # symmetry w -> -w
    if frozen is not None:
        idx = [int(v) for v in np.flatnonzero(frozen)]
        if mode == "exists":
            for v in idx: M.Add(x[v] == (1 if sign[v] > 0 else 0))
        else:
            M.Minimize(sum((1 - x[v]) if sign[v] > 0 else x[v] for v in idx))
    S = cp_model.CpSolver(); S.parameters.num_workers = 1; S.parameters.max_time_in_seconds = time_limit
    st = S.Solve(M)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        val = int(round(S.ObjectiveValue())) if mode == "dist" else 1
        xs = np.array([1.0 if S.Value(x[v]) else -1.0 for v in range(Z.n)])
        viol, _ = Z.mincut(xs * float(cfr))                 # independent check by max-flow
        assert viol <= 1e-6 and abs(xs.sum()) < 1e-9
        return ("opt" if st == cp_model.OPTIMAL else "feas"), val, xs
    if st == cp_model.INFEASIBLE: return "infeasible", None, None
    return "unknown", None, None


# ---------------------------------------------------------------- stuck-point analysis
def analyse(Z, r, res, n_exact_left, rng=None):
    c = c_of(r); b = c.denominator; cf = float(c)
    w = res["w"]; frozen = res["frozen"]
    U = np.flatnonzero(~frozen)
    vals = w[U] * b
    lat_err = float(np.abs(vals - np.rint(vals)).max()) if len(U) else 0.0
    f = res.get("f")
    if f is None: f, _ = Z.edge_flow(w)
    links = Z.tight_lattice(f)
    Uset = set(int(u) for u in U)
    link_sizes = [len(L) for L in links]
    link_unfrozen = [len(Uset & set(L)) for L in links]
    chain = []
    mask = np.zeros(Z.n, dtype=bool)
    for L in links[:-1]:
        mask[L] = True; chain.append((int(mask.sum()), Z.d(mask)))
    # how spread out are the unfrozen vertices
    diam = 0
    if len(U) > 1:
        for i, u in enumerate(U):
            dist = nx.single_source_shortest_path_length(Z.G, int(u))
            diam = max(diam, max(dist[int(v)] for v in U))
    out = dict(n_unfrozen=int(len(U)), vals_b=[int(x) for x in np.rint(vals)] if lat_err < 1e-5 else [float(x) for x in vals],
               lat_err=lat_err, n_links=len(links), link_sizes=link_sizes, link_unfrozen=link_unfrozen,
               max_unfrozen_per_link=max(link_unfrozen) if link_unfrozen else 0,
               chain=chain, unfrozen_diam=int(diam), n_walk_cuts=len(res.get("cuts", [])),
               walk_cuts=[(int(X.sum()), Z.d(X)) for X in res.get("cuts", [])])
    if n_exact_left > 0:
        sign = np.sign(w)
        st2, dist, xs = exact(Z, r, frozen, sign, "dist")
        out["flip_dist"] = dist; out["flip_status"] = st2
        out["extends"] = (dist == 0) if dist is not None else None
        if dist:
            # where are the flipped vertices (one optimal solution; not unique) relative to the unfrozen ones
            flipped = [int(v) for v in np.flatnonzero(frozen & (xs != sign))]
            dd = []
            for v in flipped:
                dist_v = nx.single_source_shortest_path_length(Z.G, v)
                dd.append(min(dist_v[int(u)] for u in U))
            out["flip_to_unfrozen"] = dd
        # null model: uniformly random signs on the same frozen set
        if rng is not None and dist is not None:
            rs = np.where(frozen, rng.choice([-1.0, 1.0], size=Z.n), 0.0)
            nst, nd, _ = exact(Z, r, frozen, rs, "dist")
            out["null_flip_dist"] = nd; out["null_status"] = nst
    return out


# ---------------------------------------------------------------- driver
METHODS = ["walk0.1", "walk0.03", "walkray", "lp", "walk0.1+lp"]

def run_one(Z, r, method, rng):
    c = c_of(r)
    if method == "lp": return lp_round(Z, c, rng)
    if method == "walkray": return walk(Z, c, rng, delta=np.inf)
    if method.endswith("+lp"):
        r1 = walk(Z, c, rng, delta=float(method[4:-3]))
        if r1["success"] or r1["status"] == "numerical": return r1
        r2 = lp_round(Z, c, rng, r1["frozen"], r1["w"])
        r2["steps"] = r1["steps"]; r2["walk_frozen"] = int(r1["frozen"].sum())
        return r2
    return walk(Z, c, rng, delta=float(method[4:]))

def run_graph(args):
    name, restarts, rgrid, seed, n_exact, outdir, deadline = args
    G = test_graphs()[name]; Z = Zono(G)
    t0 = time.time()
    rec = dict(name=name, n=Z.n, m=Z.m, diam=nx.diameter(G), exists={}, runs={})
    for r in rgrid:
        st, _, _ = exact(Z, r, time_limit=120.0)
        rec["exists"][str(r)] = st
    rec["complete"] = True
    for method in METHODS:
        for r in rgrid:
            if deadline and time.process_time() > deadline:
                rec["complete"] = False; print(f"[{name}] deadline reached before {method} r={r}", flush=True); continue
            rng = np.random.default_rng([seed, METHODS.index(method), r.numerator, r.denominator])
            succ = 0; stuck = []; steps = []; left = n_exact if r in EXACT_R else 0; bad = 0; numerical = 0
            nrest = max(1, restarts // 4) if method == "walk0.03" else restarts   # slow step-size control: a quarter of the restarts
            for it in range(nrest):
                res = run_one(Z, r, method, rng)
                steps.append(res.get("steps", res.get("rounds")))
                if res["success"]:
                    viol, _ = Z.mincut(res["w"])          # independent check of the final +-c vector
                    assert viol <= 1e-6 and abs(res["w"].sum()) < 1e-6 and np.allclose(np.abs(res["w"]), float(c_of(r)))
                    succ += 1
                else:
                    if "error" in res: bad += 1; continue
                    if res.get("status") == "numerical": numerical += 1; continue     # not a stuck point: excluded
                    a = analyse(Z, r, res, left, rng); left -= 1
                    if "lp" in method: a["first_frozen"] = res["first_frozen"]; a["spread"] = res["spread"]
                    a["rank"] = res.get("rank"); a["viol"] = float(res.get("viol", 0.0))
                    stuck.append(a)
            rec["runs"][f"{method}|{r}"] = dict(method=method, r=str(r), restarts=nrest, success=succ,
                                                 errors=bad, numerical=numerical, mean_steps=float(np.mean(steps)), stuck=stuck)
            print(f"[{name}] {method} r={r}: {succ}/{nrest}  ({time.time()-t0:.0f}s)", flush=True)
    rec["seconds"] = time.time() - t0
    if outdir:
        json.dump(rec, open(os.path.join(outdir, "zono_" + name.replace("(", "").replace(")", "").replace(",", "_") + ".json"), "w"))
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--restarts", type=int, default=50)
    ap.add_argument("--out", default=os.path.join(HERE, "zono_round_results.json"))
    ap.add_argument("--graphs", default="")
    ap.add_argument("--jobs", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--rs", default="", help="';'-separated r values, e.g. '8;6;11/2;5' (default: full grid)")
    ap.add_argument("--snarks", default="", help="';'-separated graph names that additionally get --snark-rs")
    ap.add_argument("--snark-rs", default="9/2")
    ap.add_argument("--minutes", type=float, default=0.0, help="CPU-time cap in minutes; cells not started by then are skipped")
    ap.add_argument("--partial-dir", default="", help="also write one JSON per graph here as soon as it is done")
    ap.add_argument("--exact-per-cell", type=int, default=10, help="stuck runs per (graph, method, r) sent to CP-SAT")
    a = ap.parse_args()
    names = [g for g in a.graphs.split(";") if g] or list(test_graphs())     # ';'-separated: names contain commas
    unknown = [g for g in names if g not in test_graphs()]
    assert not unknown, unknown
    deadline = time.process_time() + 60 * a.minutes if a.minutes else 0   # CPU time: immune to machine sleep
    base = [Fraction(x) for x in a.rs.split(";") if x] or R_GRID
    extra = [Fraction(x) for x in a.snark_rs.split(";") if x]
    snarks = set(g for g in a.snarks.split(";") if g)
    tasks = [(nm, a.restarts, sorted(set(base + (extra if nm in snarks else []))), a.seed, a.exact_per_cell,
              a.partial_dir, deadline) for nm in names]
    if a.jobs > 1:
        import multiprocessing as mp
        with mp.Pool(min(a.jobs, 2)) as P:
            recs = P.map(run_graph, tasks, chunksize=1)
    else:
        recs = [run_graph(t) for t in tasks]
    json.dump(recs, open(a.out, "w"))
    print("wrote", a.out)


if __name__ == "__main__":
    main()
