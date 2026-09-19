"""Re-solve the configurations left UNKNOWN in enum_c6_cache.json with the lazy model (UMAX = 1) and a long limit."""
import sys, os, json, time
sys.argv = ["x", "noop"]; src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "enum_c6.py")).read().split('if MODE == "pairs":')[0]
ns = {"__file__": os.path.abspath("enum_c6.py")}; exec(src, ns)
from count_c6 import solve
T = float(os.environ.get("T", "900")); cache = ns["cache"]
for key in [k for k, v in cache.items() if v == "UNKNOWN"]:
    ends, cuts = ns["config"](key.split()); t0 = time.time()
    st, sol = solve(ends, cuts, 1, T, workers=8, verbose=True)
    print(f"{key:60s} {st} [{time.time()-t0:.0f}s] {sol or ''}", flush=True)
    if st != "UNKNOWN":
        c = json.load(open(ns["CACHE"])); c[key] = st; c[key + " T"] = T
        json.dump(c, open(ns["CACHE"] + ".tmp", "w"), indent=0); os.replace(ns["CACHE"] + ".tmp", ns["CACHE"])
