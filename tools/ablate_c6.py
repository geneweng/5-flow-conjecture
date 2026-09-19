"""Ablation of the infeasible PAIRS of (C6) items: which constraint groups does each exclusion need?
Re-solves every cached INFEASIBLE pair (one representative per type signature) without circuit facts, without
union facts, and in PURE mode (no even circuits of H crossing).  Output: a table by pair type."""
import sys, os, json, collections
sys.argv = ["x", "noop"]; src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "enum_c6.py")).read().split('if MODE == "pairs":')[0]
ns = {"__file__": os.path.abspath("enum_c6.py")}; exec(src, ns)
from count_c6 import solve
cache = ns["cache"]
pairs = [k.split() for k, v in cache.items() if v == "INFEASIBLE" and len(k.split()) == 2 and not k.startswith("PURE")]
tab = collections.Counter()
for a, b in pairs:
    ends, cuts = ns["config"]([a, b])
    typ = tuple(sorted(n.split("_")[0][:2] if n[0] == "P" else "L" for n in (a, b)))
    nc = solve(ends, cuts, 2, 60, NOCIRC=True)[0]; nu = solve(ends, cuts, 2, 60, NOUNION=True)[0]
    tab[(typ, "needs circuits" if nc == "FEASIBLE" else "no circuits needed", "needs unions" if nu == "FEASIBLE" else "no unions needed")] += 1
    print(a, b, "NOCIRC:", nc, "NOUNION:", nu, flush=True)
print()
for k in sorted(tab): print(k, tab[k])
