"""Experiment B (notes/oddness.md 18.2): the "defect distance" from the MS template to a true 5-flow.

For a cubic graph G, 2-factor, 0-edge choice and canonical colouring, H = colours 1,2.  The template uses
bipartitions that are proper 2-colourings of H.  Here we minimise, over ALL balanced bipartitions
(= nowhere-zero 5-flows), the number of monochromatic H-edges ("defects").

Balancedness (exponentially many cut conditions) is modelled exactly by a flow: black is balanced iff there is
an integer flow f with |f(e)| <= 3 and net outflow +5 at black, -5 at white vertices (max-flow/min-cut on the
network of oddness.balanced; capacities integral => integral flow).  CP-SAT model:
    x_v in {0,1}, f_e in [-3,3], d_e in {0,1} for H-edges,  x_u xor x_v xor d_e = 1,
    sum_out f - sum_in f = 10 x_v - 5,   minimise sum d_e.
Every solution is re-checked with oddness.balanced.

usage:  python3 tools/defect_min.py            (all 26 archived failures + controls; writes tools/defect_min_results.json)
"""
import sys, os, json, glob, itertools, random, time, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from ortools.sat.python import cp_model
from oddness import canonical_colouring, H_components, partition, balanced
from badcuts import violating_set

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKERS = 2
ENUM_CAP = 200


def load(fn, zero=None):
    d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
    M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    zero = list(d["zero"]) if zero is None else list(zero)
    col = canonical_colouring(G, M, circuits, zero); comps = H_components(G, col)
    return d, G, circuits, M, zero, col, comps


def build(G, col, allowed=(1, 2)):
    m = cp_model.CpModel()
    x = {v: m.NewBoolVar(f"x{v}") for v in G.nodes()}
    f = {}; out = {v: [] for v in G.nodes()}
    for u, v in G.edges():
        f[(u, v)] = m.NewIntVar(-3, 3, f"f{u}_{v}")
        out[u].append(f[(u, v)]); out[v].append(-f[(u, v)])
    for v in G.nodes():
        m.Add(sum(out[v]) == 10 * x[v] - 5)
    dvar = {}
    for e, c in col.items():
        if c not in (1, 2): continue
        u, v = tuple(e)
        if c in allowed:
            dv = m.NewBoolVar(f"d{u}_{v}"); m.AddBoolXOr([x[u], x[v], dv]); dvar[e] = dv
        else:
            m.Add(x[u] + x[v] == 1)
    return m, x, dvar


def solver():
    s = cp_model.CpSolver(); s.parameters.num_workers = WORKERS; s.parameters.max_time_in_seconds = 600
    return s


def min_defects(G, col, allowed=(1, 2)):
    """Returns (opt, black, defect list) or (None, None, None) if infeasible."""
    m, x, dvar = build(G, col, allowed)
    m.Minimize(sum(dvar.values())) if dvar else None
    s = solver(); st = s.Solve(m)
    if st == cp_model.INFEASIBLE: return None, None, None
    assert st == cp_model.OPTIMAL, s.StatusName(st)
    black = {v for v in G.nodes() if s.Value(x[v])}
    D = [e for e, dv in dvar.items() if s.Value(dv)]
    assert balanced(G, black), "CP-SAT solution fails the min-cut check"
    assert len(D) == sum(1 for e, c in col.items() if c in (1, 2) and len(set(tuple(e)) & black) != 1)
    return len(D), black, D


def enum_defect_sets(G, col, opt, cap=ENUM_CAP):
    """All defect sets of size opt that admit a balanced bipartition (no-good cuts on the d-variables)."""
    m, x, dvar = build(G, col)
    m.Add(sum(dvar.values()) == opt)
    found = []
    while len(found) < cap:
        s = solver(); st = s.Solve(m)
        if st == cp_model.INFEASIBLE: return found, True
        assert st in (cp_model.OPTIMAL, cp_model.FEASIBLE)
        black = {v for v in G.nodes() if s.Value(x[v])}; assert balanced(G, black)
        D = frozenset(e for e, dv in dvar.items() if s.Value(dv)); assert len(D) == opt
        found.append(D)
        m.Add(sum(dvar[e] for e in D) <= opt - 1)
    return found, False


def colourings_with_defects(G, col, comps, D):
    """All bipartitions whose monochromatic H-edges are exactly D, up to global complement
    (component 0 not flipped).  Yields (flips, black)."""
    H = nx.Graph(); H.add_nodes_from(G.nodes())
    for e, c in col.items():
        if c in (1, 2): u, v = tuple(e); H.add_edge(u, v, same=(e in D))
    base = []
    for colouring, _ in comps:
        root = min(colouring); val = {root: 0}; stack = [root]; ok = True
        while stack:
            u = stack.pop()
            for w in H[u]:
                c = val[u] ^ (0 if H[u][w]["same"] else 1)
                if w in val: ok &= (val[w] == c)
                else: val[w] = c; stack.append(w)
        if not ok: return
        base.append(val)
    for bits in itertools.product((0, 1), repeat=len(comps) - 1):
        flips = (0,) + bits
        yield flips, {v for val, f in zip(base, flips) for v, c in val.items() if c ^ f}


def retemplate(G, M, circuits, zero, black):
    """Is black a TEMPLATE bipartition for another 0-edge choice on the same 2-factor?  Needs every matching edge
    bichromatic and, on each odd circuit, a 0-edge position p with edges p+1,p+3,.. bichromatic (canonical convention of
    oddness.canonical_colouring; 'mirrored' = p+2,p+4,..), on each even circuit one alternating class bichromatic.
    Returns (moves_canonical, moves_any): number of circuits whose 0-edge/colour class must change, or None."""
    bi = lambda u, v: (u in black) != (v in black)
    if not all(bi(*tuple(e)) for e in M): return None, None
    mc = ma = 0; j = 0
    for C in circuits:
        L = len(C); b = [bi(C[i], C[(i + 1) % L]) for i in range(L)]
        if L % 2 == 0:
            ok0 = all(b[i] for i in range(0, L, 2)); ok1 = all(b[i] for i in range(1, L, 2))
            can = 0 if ok0 else None; any_ = 0 if ok0 else (1 if ok1 else None)
        else:
            z = zero[j]; j += 1
            canon = [p for p in range(L) if all(b[(p + k) % L] for k in range(1, L, 2))]
            mirr = [p for p in range(L) if all(b[(p + k) % L] for k in range(2, L, 2))]
            can = None if not canon else (0 if z in canon else 1)
            any_ = None if not (canon or mirr) else (0 if z in canon else 1)
        mc = None if (mc is None or can is None) else mc + can
        ma = None if (ma is None or any_ is None) else ma + any_
    return mc, ma


def path_order(G, col, colouring, ends):
    Hn = {v: [] for v in colouring}
    for e, c in col.items():
        if c in (1, 2):
            u, v = tuple(e)
            if u in Hn: Hn[u].append(v); Hn[v].append(u)
    start = min(ends) if ends else min(colouring)
    P = [start]; prev = None
    while True:
        nxt = [w for w in Hn[P[-1]] if w != prev]
        if not nxt or nxt[0] == start: break
        prev = P[-1]; P.append(nxt[0])
    return P


def analyse(fn, zero=None, enum=True, variants=True):
    d, G, circuits, M, zero, col, comps = load(fn, zero)
    n = G.number_of_nodes(); t0 = time.time()
    npaths = sum(1 for _, e in comps if e); ncirc = len(comps) - npaths
    orders = [path_order(G, col, c, e) for c, e in comps]
    where = {}                                   # H-edge -> (component index, position, #edges in comp)
    for ci, P in enumerate(orders):
        L = len(P) - 1 if comps[ci][1] else len(P)
        for i in range(L): where[frozenset((P[i], P[(i + 1) % len(P)]))] = (ci, i, L)
    circ_of = {v: j for j, C in enumerate(circuits) for v in C}
    # template colourings and their killers (computed by min cut; smaller side)
    killers = {}; template_ok = 0
    if len(comps) <= 9:
        for bits in itertools.product((0, 1), repeat=len(comps) - 1):
            fl = (0,) + bits; S = violating_set(G, partition(comps, fl))
            if S is None: template_ok += 1
            else: killers[fl] = frozenset(S)
    stored = []
    if isinstance(d.get("S"), dict) and zero == list(d["zero"]):
        for S in d["S"].values():
            S = frozenset(S); S = S if 2 * len(S) <= n else frozenset(G.nodes()) - S
            stored.append(S)
    comp_kill = set(killers.values()); stored_set = set(stored)
    res = dict(file=os.path.relpath(fn, ROOT), n=n, zero=zero, circuits=[len(C) for C in circuits], paths=npaths, Hcircuits=ncirc,
               path_ends=[sorted(e) for _, e in comps if e], path_edges=[len(P) - 1 for P, (_, e) in zip(orders, comps) if e],
               template_balanced=template_ok, killers_distinct=len(comp_kill),
               killer_sizes=sorted((len(S), nx.cut_size(G, S)) for S in comp_kill),
               stored_S=len(stored_set), stored_equal_computed=(stored_set == comp_kill) if stored else None,
               graph_key=hash((frozenset(map(frozenset, G.edges())), tuple(zero), tuple(map(tuple, circuits)))))
    opt, black, D = min_defects(G, col)
    res["min_defects"] = opt
    if variants:
        res["min_colour2_only"] = min_defects(G, col, allowed=(2,))[0]
        res["min_colour1_only"] = min_defects(G, col, allowed=(1,))[0]
    if enum and opt:
        Ds, complete = enum_defect_sets(G, col, opt)
        res["n_defect_sets"] = len(Ds); res["enum_complete"] = complete
        allk = comp_kill | stored_set
        sols = []
        for Dset in Ds:
            flips_ok = []; retemps = []
            for fl, bl in colourings_with_defects(G, col, comps, Dset):
                if balanced(G, bl):
                    # nearest template colouring (per component majority) and its killer
                    near = []
                    for (colouring, _), f in zip(comps, fl):
                        agree = sum(1 for v, c in colouring.items() if (c ^ f) == (v in bl))
                        near.append(f if 2 * agree >= len(colouring) else 1 - f)
                    if near[0] == 1: near = [1 - b for b in near]      # killers are complement-invariant
                    flips_ok.append((fl, tuple(near), bl)); retemps.append(retemplate(G, M, circuits, zero, bl))
            assert flips_ok
            edges = []
            for e in sorted(Dset, key=lambda e: where[e]):
                u, v = sorted(e); ci, pos, L = where[e]
                rel = collections.Counter("boundary" if len(e & S) == 1 else ("inside" if len(e & S) == 2 else "outside") for S in allk)
                near_rel = collections.Counter()
                for _, near, _ in flips_ok:
                    S = killers.get(near)
                    if S is not None: near_rel["boundary" if len(e & S) == 1 else ("inside" if len(e & S) == 2 else "outside")] += 1
                edges.append(dict(edge=[u, v], colour=col[e], comp=ci, is_path=bool(comps[ci][1]), ends=sorted(comps[ci][1]), pos=pos, L=L,
                                  circuits=sorted({circ_of[u], circ_of[v]}),
                                  killers_boundary=rel["boundary"], killers_inside=rel["inside"], killers_outside=rel["outside"],
                                  nearest_killer=dict(near_rel)))
            sols.append(dict(defects=edges, n_bipartitions=len(flips_ok), retemplate=retemps))
        res["solutions"] = sols
        if opt == 1 and complete:                # baseline: every colour-2 H-edge, does it work as the single defect, and where is it
            works = {next(iter(Dset)) for Dset in Ds}; tab = []
            for e, c in sorted(col.items(), key=lambda t: sorted(t[0])):
                if c != 2: continue
                rel = "boundary" if any(len(e & S) == 1 for S in allk) else ("inside" if any(len(e & S) == 2 for S in allk) else "outside")
                tab.append(dict(edge=sorted(e), works=e in works, rel=rel, circuit_len=len(circuits[circ_of[min(e)]])))
            res["colour2_edges"] = tab
        res["n_bipartitions"] = sum(s["n_bipartitions"] for s in sols)
    res["seconds"] = round(time.time() - t0, 1)
    return res


def control(fn, k, rng):
    """k random non-failing 0-edge choices: min defects (expected 0)."""
    d, G, circuits, *_ = load(fn); odd = [C for C in circuits if len(C) % 2]
    out = []; tried = 0
    while len(out) < k:
        zc = [rng.randrange(len(C)) for C in odd]; tried += 1
        _, _, _, _, _, col, comps = load(fn, zc)
        ok = any(balanced(G, partition(comps, (0,) + b)) for b in itertools.product((0, 1), repeat=len(comps) - 1))
        if not ok: continue
        out.append((zc, min_defects(G, col)[0]))
    return dict(file=os.path.relpath(fn, ROOT), sampled=tried, nonfailing=len(out), min_defects=collections.Counter(o[1] for o in out))


if __name__ == "__main__":
    files = sorted(glob.glob(os.path.join(ROOT, "data", "oddness6_*.json"))) + sorted(glob.glob(os.path.join(ROOT, "tools", "fullkill_*.json")))
    results = dict(instances=[], controls=[], robust_fails=[])
    for fn in files:
        r = analyse(fn); results["instances"].append(r)
        print(f"{r['file']:55s} n={r['n']} H={r['paths']}P+{r['Hcircuits']}C template_ok={r['template_balanced']} min={r['min_defects']} "
              f"c2only={r['min_colour2_only']} c1only={r['min_colour1_only']} sets={r.get('n_defect_sets')} bip={r.get('n_bipartitions')} [{r['seconds']}s]", flush=True)
    rng = random.Random(2026)
    for fn in ["data/oddness6_template_counterexample_n42.json", "tools/fullkill_climb_six7_1110_1643.json", "tools/fullkill_complete_505.json"]:
        c = control(os.path.join(ROOT, fn), 30, rng); c["min_defects"] = dict(c["min_defects"]); results["controls"].append(c); print("control", c, flush=True)
    fn = os.path.join(ROOT, "data", "oddness6_robust_stage2b_n50.json")
    for zc in json.load(open(fn))["fails"]:
        r = analyse(fn, zero=zc); results["robust_fails"].append(r)
        print(f"robust n50 zero={zc} template_ok={r['template_balanced']} min={r['min_defects']} c2only={r['min_colour2_only']} c1only={r['min_colour1_only']} sets={r.get('n_defect_sets')}", flush=True)
    json.dump(results, open(os.path.join(ROOT, "tools", "defect_min_results.json"), "w"), indent=1, default=str)
