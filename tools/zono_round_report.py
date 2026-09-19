"""Tables for notes/exp-C-rounding.md from the JSON written by zono_round.py.
Usage: python3 tools/zono_round_report.py [results.json]  > tables.md"""
import sys, json, os
from fractions import Fraction
from collections import Counter
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "zono_round_results.json")
recs = json.load(open(path))
ORDER = ["Petersen", "J5", "J7", "J9", "GP(13,3)", "GP(17,4)", "R2", "cx42", "cx50", "cx60"]
recs.sort(key=lambda r: ORDER.index(r["name"]) if r["name"] in ORDER else 99)
def fr(x): return float(Fraction(x))
RS = sorted({r for rec in recs for r in rec["exists"]}, key=fr)      # union: some r (the 9/2 control) exist only for the snarks
METHODS = []
for rec in recs:
    for k in rec["runs"]:
        m = k.split("|")[0]
        if m not in METHODS: METHODS.append(m)
def ex_cell(rec, r):
    st = rec["exists"].get(r)
    return "n/a" if st is None else ("yes" if st in ("opt", "feas") else ("no" if st == "infeasible" else "?"))
def restarts_of(m):
    return sorted({R["restarts"] for rec in recs for k, R in rec["runs"].items() if k.split("|")[0] == m})

def med(xs): return f"{np.median(xs):g}" if len(xs) else "-"
def rng_(xs): return f"{min(xs)}-{max(xs)}" if len(xs) else "-"
def hist(c, keys=None):
    keys = sorted(c) if keys is None else keys
    return " ".join(f"{k}:{c[k]}" for k in keys if c.get(k))

print("### T1. Ground truth (CP-SAT): does Z(G) contain a +-r/(r-2) vector?\n")
print("| graph | n | diam | " + " | ".join(RS) + " |")
print("|---|---|---|" + "---|" * len(RS))
for rec in recs:
    print(f"| {rec['name']} | {rec['n']} | {rec['diam']} | " + " | ".join(ex_cell(rec, r) for r in RS) + " |")

for m in METHODS:
    print(f"\n### T2[{m}]. Successes out of {'/'.join(map(str, restarts_of(m)))} restarts ('x' = no +-c vector exists, 'n/a' = cell not run)\n")
    print("| graph | " + " | ".join(RS) + " |")
    print("|---|" + "---|" * len(RS))
    for rec in recs:
        cells = []
        for r in RS:
            R = rec["runs"].get(f"{m}|{r}")
            if R is None: cells.append("n/a"); continue
            cells.append("x" if rec["exists"].get(r) == "infeasible" and R["success"] == 0 else str(R["success"]))
        print(f"| {rec['name']} | " + " | ".join(cells) + " |")

print("\n### T3. Smallest r of the grid at which some restart succeeds\n")
print("| graph | exists (exact) | " + " | ".join(METHODS) + " |")
print("|---|---|" + "---|" * len(METHODS))
for rec in recs:
    ex = [r for r in RS if rec["exists"].get(r) in ("opt", "feas")]
    cells = []
    for m in METHODS:
        ok = [r for r in RS if rec["runs"].get(f"{m}|{r}", {"success": 0})["success"] > 0]
        cells.append(min(ok, key=fr) if ok else "none")
    print(f"| {rec['name']} | {min(ex, key=fr) if ex else 'none'} | " + " | ".join(cells) + " |")

print("\n### T4. Stuck runs.  u = number of unfrozen coordinates; links = strongly connected components of the residual digraph "
      "(= links of a maximal chain of tight sets); d = boundary sizes of the chain sets; "
      "flip = min number of frozen signs to change so that the frozen pattern extends to a +-c vector in Z (exact, first 10 stuck runs per cell); "
      "null = same for uniformly random signs on the same frozen set\n")
for m in METHODS:
    print(f"\n#### T4[{m}]\n")
    print("| graph | r | stuck | u med (range) | in (1/b)Z | max unfrozen per link | links med | links with >1 vertex med | sticky cuts med | chain d range | chain d<=7 share | unfrozen diam med / graph diam | extends | flip hist | null hist |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for rec in recs:
        for r in RS:
            if r not in ("5", "11/2", "6", "8"): continue
            R = rec["runs"].get(f"{m}|{r}")
            if R is None: continue
            S = R["stuck"]
            if not S: print(f"| {rec['name']} | {r} | 0 |" + " |" * 12); continue
            u = [s["n_unfrozen"] for s in S]
            lat = sum(1 for s in S if s["lat_err"] < 1e-6)
            mul = max(s["max_unfrozen_per_link"] for s in S)
            links = [s["n_links"] for s in S]
            big = [sum(1 for z in s["link_sizes"] if z > 1) for s in S]
            wc = [s["n_walk_cuts"] for s in S]
            ds = [d for s in S for (_, d) in s["chain"]]
            small = sum(1 for d in ds if d <= 7) / max(1, len(ds))
            ud = [s["unfrozen_diam"] for s in S]
            ex = [s for s in S if "flip_dist" in s]
            ext = sum(1 for s in ex if s.get("extends"))
            fh = Counter(s["flip_dist"] for s in ex); nh = Counter(s.get("null_flip_dist") for s in ex)
            print(f"| {rec['name']} | {r} | {len(S)} | {med(u)} ({rng_(u)}) | {lat}/{len(S)} | {mul} | {med(links)} | {med(big)} | "
                  f"{med(wc) if 'walk' in m and '+' not in m else '-'} | {rng_(ds)} | {small:.2f} | {med(ud)} / {rec['diam']} | {ext}/{len(ex)} | {hist(fh)} | {hist(nh)} |")

print("\n### T5. Values of the unfrozen coordinates at stuck points, in units of 1/b (c = a/b), pooled over graphs\n")
print("| method | r | c | histogram of b*|w_v| over unfrozen v |")
print("|---|---|---|---|")
for m in METHODS:
    for r in RS:
        if r not in ("5", "11/2", "6", "8"): continue
        c = Fraction(r) / (Fraction(r) - 2); h = Counter(); bad = 0
        for rec in recs:
            for s in rec["runs"].get(f"{m}|{r}", {"stuck": []})["stuck"]:
                if s["lat_err"] < 1e-6: h.update(abs(int(v)) for v in s["vals_b"])
                else: bad += 1
        print(f"| {m} | {r} | {c} | {hist(h)}" + (f" (+{bad} stuck points off the lattice)" if bad else "") + " |")

print("\n### T6. Sizes of the links that contain an unfrozen vertex, and (|X|, d(X)) of the sticky cuts of the walk, pooled over graphs with n >= 40, r = 5\n")
for m in METHODS:
    ls = Counter(); cs = Counter()
    for rec in recs:
        if rec["n"] < 40: continue
        for s in rec["runs"].get(f"{m}|5", {"stuck": []})["stuck"]:
            for z, k in zip(s["link_sizes"], s["link_unfrozen"]):
                if k: ls[z] += 1
            for (x, d) in s["walk_cuts"]: cs[d] += 1
    print(f"- **{m}**: link sizes {hist(ls)}" + (f"; sticky cut d(X): {hist(cs)}" if cs else ""))

print("\n### T7. Distance (in G) from each flipped vertex of one optimal repair to the nearest unfrozen vertex, pooled, r in {5, 6}\n")
for m in METHODS:
    h = Counter()
    for rec in recs:
        for r in ("5", "6"):
            for s in rec["runs"].get(f"{m}|{r}", {"stuck": []})["stuck"]:
                h.update(s.get("flip_to_unfrozen", []))
    print(f"- **{m}**: {hist(h)}")
