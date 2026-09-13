"""Enumerate every minimal cover of the four colouring classes at oddness 6 by tight cuts with s = 0:
lines = pair 6-cuts (6,4,2,2) of codimension 2 (kill two classes), points = 11-cuts (11,7,4,3) or 7-cuts
(7,5,2,3) (kill one class); all placements of the avoided paths; run the count relaxation on each.
Classes = triples containing z1, encoded (b2, b3): b2 = 0 means z3 (else z4), b3 = 0 means z5 (else z6).
Output: one line per configuration with the relaxation status; a second pass on the infeasible ones
without circuit / union constraints shows which ingredients the proof needs."""
import itertools, sys, time
from count_general import solve
UMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 2; TIME = float(sys.argv[2]) if len(sys.argv) > 2 else 120
import re; FILTER = re.compile(sys.argv[3]) if len(sys.argv) > 3 else None; SECOND = len(sys.argv) <= 4
PT = {"11": (7, 4), "7": (5, 2)}
def line(pair, b, p):
    """membership of z1..z6 for a pair cut; pair in {(0,1),(0,2),(1,2)}; b = which end of the second path (or sign for (1,2)); p = placement of the avoided path"""
    mem = {}
    if pair == (0, 1): mem = {1: 1, 2: 0, 3: 1 - b, 4: b, 5: p, 6: p}; killed = {(b, 0), (b, 1)}
    elif pair == (0, 2): mem = {1: 1, 2: 0, 5: 1 - b, 6: b, 3: p, 4: p}; killed = {(0, b), (1, b)}
    else:   # (1,2): S ∋ z3 and z5 (b = 0, kills b2 == b3) or z3 and z6 (b = 1, kills b2 != b3)
        mem = {3: 1, 4: 0, 5: 1 - b, 6: b, 1: p, 2: p}; killed = {(0, 0), (1, 1)} if b == 0 else {(0, 1), (1, 0)}
    return mem, killed
def point(cls):
    b2, b3 = cls; return {1: 1, 2: 0, 3: 1 - b2, 4: b2, 5: 1 - b3, 6: b3}, {cls}
def config(items):
    """items: list of (membership dict, (c1,c2)); returns ends, cuts"""
    J = len(items); ends = {}
    for i, (za, zb) in enumerate([(1, 2), (3, 4), (5, 6)]):
        ends[i] = (tuple(items[j][0][za] for j in range(J)), tuple(items[j][0][zb] for j in range(J)))
    return ends, [it[1] for it in items]
ALL = {(0, 0), (0, 1), (1, 0), (1, 1)}
configs = []
pairs = [(0, 1), (0, 2), (1, 2)]
# situation 1: two lines on the same pair with opposite signs
for pr in pairs:
    for p, q in itertools.product((0, 1), repeat=2):
        (m1, k1), (m2, k2) = line(pr, 0, p), line(pr, 1, q); configs.append((f"S1 {pr} pl{p}{q}", [(m1, (4, 2)), (m2, (4, 2))]))
# situation 2: three lines, one per pair, covering everything
for b, bb in itertools.product((0, 1), repeat=2):
    for sg in (0, 1):
        (m1, k1), (m2, k2), (m3, k3) = line((0, 1), b, 0), line((0, 2), bb, 0), line((1, 2), sg, 0)
        if k1 | k2 | k3 != ALL: continue
        for p1, p2, p3 in itertools.product((0, 1), repeat=3):
            items = [(line((0, 1), b, p1)[0], (4, 2)), (line((0, 2), bb, p2)[0], (4, 2)), (line((1, 2), sg, p3)[0], (4, 2))]
            configs.append((f"S2 b{b}{bb} s{sg} pl{p1}{p2}{p3}", items))
# situation 3: two lines on different pairs + one point
for pa, pb in itertools.combinations(pairs, 2):
    for b, bb in itertools.product((0, 1), repeat=2):
        (m1, k1), (m2, k2) = line(pa, b, 0), line(pb, bb, 0); miss = ALL - k1 - k2
        if len(miss) != 1: continue
        for pt in ("11", "7"):
            for p1, p2 in itertools.product((0, 1), repeat=2):
                items = [(line(pa, b, p1)[0], (4, 2)), (line(pb, bb, p2)[0], (4, 2)), (point(next(iter(miss)))[0], PT[pt])]
                configs.append((f"S3 {pa}{pb} b{b}{bb} P{pt} pl{p1}{p2}", items))
# situation 4: one line + two points
for pr in pairs:
    for b in (0, 1):
        m1, k1 = line(pr, b, 0); miss = sorted(ALL - k1)
        for pt1, pt2 in itertools.product(("11", "7"), repeat=2):
            for p in (0, 1):
                items = [(line(pr, b, p)[0], (4, 2)), (point(miss[0])[0], PT[pt1]), (point(miss[1])[0], PT[pt2])]
                configs.append((f"S4 {pr} b{b} P{pt1},{pt2} pl{p}", items))
# situation 5: four points
for pts in itertools.product(("11", "7"), repeat=4):
    items = [(point(c)[0], PT[t]) for c, t in zip(sorted(ALL), pts)]
    configs.append((f"S5 P{','.join(pts)}", items))
print(len(configs), "configurations", flush=True)
results = {}
configs = [c for c in configs if FILTER is None or FILTER.search(c[0])]
for name, items in configs:
    ends, cuts = config(items); t0 = time.time()
    st, n = solve(ends, cuts, UMAX, TIME)
    results[name] = st
    print(f"{name:34s} {st:11s} [{time.time()-t0:.0f}s]" + (f"  n={n}" if n else ""), flush=True)
if not SECOND: raise SystemExit
print("\nsecond pass on infeasible configurations without circuit / union constraints:", flush=True)
for name, items in configs:
    if results[name] != "INFEASIBLE": continue
    ends, cuts = config(items)
    a = solve(ends, cuts, UMAX, TIME, NOCIRC=True)[0]; b = solve(ends, cuts, UMAX, TIME, NOUNION=True)[0]; c = solve(ends, cuts, UMAX, TIME, NOCIRC=True, NOUNION=True)[0]
    print(f"{name:34s} no-circuits:{a:11s} no-unions:{b:11s} neither:{c}", flush=True)
