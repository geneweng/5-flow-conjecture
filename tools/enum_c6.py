"""(C6), notes §18.1: can the four colouring classes be covered by tight 6- and 7-cuts?  Items:
  lines  L(pair,b,p)      pair cut (6,4,2,2), third path avoids the cut (placement p), kills two classes
  points P7(cls)          7-cut (7,5,2,3)
         P641(cls)        6-cut (6,4,1,3), one non-H boundary edge
         Pc(cls,i,p)      pair cut (6,4,2,2) crossed by its non-separated path i (placement p), kills one class
all with even circuits of H allowed to cross (s free).  A sub-configuration of a realizable configuration is
realizable, so the search is hierarchical: all pairs of items first, then every minimal cover that contains no
infeasible pair (then no infeasible triple).
usage: enum_c6.py pairs|covers [TIME] [PURE]      (results cached in tools/enum_c6_cache.json)"""
import itertools, sys, time, json, os
from count_c6 import solve
HERE = os.path.dirname(os.path.abspath(__file__)); CACHE = os.path.join(HERE, "enum_c6_cache.json")
MODE = sys.argv[1] if len(sys.argv) > 1 else "pairs"; TIME = float(sys.argv[2]) if len(sys.argv) > 2 else 60
PURE = len(sys.argv) > 3 and sys.argv[3] == "PURE"
ALL = [(0, 0), (0, 1), (1, 0), (1, 1)]; PAIRS = [(0, 1), (0, 2), (1, 2)]
def line(pair, b, p):
    if pair == (0, 1): mem = {1: 1, 2: 0, 3: 1 - b, 4: b, 5: p, 6: p}; killed = {(b, 0), (b, 1)}
    elif pair == (0, 2): mem = {1: 1, 2: 0, 5: 1 - b, 6: b, 3: p, 4: p}; killed = {(0, b), (1, b)}
    else: mem = {3: 1, 4: 0, 5: 1 - b, 6: b, 1: p, 2: p}; killed = {(0, 0), (1, 1)} if b == 0 else {(0, 1), (1, 0)}
    return mem, killed
def point_mem(cls):
    b2, b3 = cls; return {1: 1, 2: 0, 3: 1 - b2, 4: b2, 5: 1 - b3, 6: b3}
items = {}
for pr in PAIRS:
    for b in (0, 1):
        for p in (0, 1):
            mem, killed = line(pr, b, p); items[f"L{pr[0]}{pr[1]}b{b}p{p}"] = (mem, dict(c1=4, c2=2), frozenset(killed))
for cls in ALL:
    nm = f"{cls[0]}{cls[1]}"
    items[f"P7_{nm}"] = (point_mem(cls), dict(c1=5, c2=2), frozenset([cls]))
    items[f"P641_{nm}"] = (point_mem(cls), dict(c1=4, c2=1, nonH=1), frozenset([cls]))
    for i in range(3):
        for p in (0, 1):
            mem = point_mem(cls); mem[2 * i + 1] = mem[2 * i + 2] = p
            beta = [1, 1 - cls[0], 1 - cls[1]][i]            # start of path i (z1, z3, z5) black under cls?
            items[f"Pc{i}p{p}_{nm}"] = (mem, dict(c1=4, c2=2, cross={i: beta}), frozenset([cls]))
def config(names):
    J = len(names); ends = {}
    for i, (za, zb) in enumerate([(1, 2), (3, 4), (5, 6)]):
        ends[i] = (tuple(items[nm][0][za] for nm in names), tuple(items[nm][0][zb] for nm in names))
    return ends, [items[nm][1] for nm in names]
cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
def status(names, T=TIME):
    key = ("PURE " if PURE else "") + " ".join(sorted(names))
    if key in cache and (cache[key] != "UNKNOWN" or T <= cache.get(key + " T", 0)): return cache[key]
    ends, cuts = config(sorted(names)); t0 = time.time()
    st, sol = solve(ends, cuts, 2, T, PURE=PURE)
    cache[key] = st; cache[key + " T"] = T
    tmp = CACHE + ".tmp"; json.dump(cache, open(tmp, "w"), indent=0); os.replace(tmp, CACHE)
    print(f"{key:60s} {st:10s} [{time.time()-t0:.0f}s] {sol if sol else ''}", flush=True)
    return st
def covers():
    names = sorted(items); L = [n for n in names if n[0] == "L"]; P = {c: [n for n in names if n[0] == "P" and items[n][2] == frozenset([c])] for c in ALL}
    out = []
    for r in (2, 3):
        for ls in itertools.combinations(L, r):
            ks = [items[l][2] for l in ls]
            if frozenset().union(*ks) == set(ALL) and all(frozenset().union(*(ks[:q] + ks[q + 1:])) != set(ALL) for q in range(r)): out.append(ls)
    for l1, l2 in itertools.combinations(L, 2):
        miss = set(ALL) - items[l1][2] - items[l2][2]
        if len(miss) == 1 and items[l1][2] != items[l2][2]:
            for p in P[next(iter(miss))]: out.append((l1, l2, p))
    for l in L:
        miss = sorted(set(ALL) - items[l][2])
        for p, q in itertools.product(P[miss[0]], P[miss[1]]): out.append((l, p, q))
    for ps in itertools.product(*(P[c] for c in ALL)): out.append(ps)
    return out
if MODE == "pairs":
    names = sorted(items); print(len(names), "items", flush=True)
    for a, b in itertools.combinations(names, 2):
        if items[a][2] == items[b][2] and a[0] == b[0] == "P": continue      # two points on the same class: not needed in a minimal cover
        status((a, b))
else:
    cv = covers(); print(len(cv), "minimal covers", flush=True); left = []
    for c in cv:
        if any(cache.get(("PURE " if PURE else "") + " ".join(sorted(s))) == "INFEASIBLE" for r in (2, 3) for s in itertools.combinations(c, r)): continue
        left.append(c)
    print(len(left), "covers contain no known infeasible pair/triple", flush=True)
    tally = {}
    for c in left:
        if len(c) == 4:      # try triples first
            if any(status(s) == "INFEASIBLE" for s in itertools.combinations(c, 3)): tally["INFEASIBLE"] = tally.get("INFEASIBLE", 0) + 1; continue
        st = status(c); tally[st] = tally.get(st, 0) + 1
    print("tally over remaining covers:", tally)
