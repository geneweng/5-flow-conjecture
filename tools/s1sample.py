"""Random multigraphs with degrees in [dmin,dmax] via degree-preserving switches,
filtered to: every proper set A with |A|,|A^c|>=2 has d(A) >= 12 (the (S1) hypothesis),
tested for strong Z5-connectivity (exact, sumset).  Usage: n dmin dmax seed budget"""
import sys, random, time, itertools
from mod5 import strongly_connected_fast

n, dmin, dmax, seed, budget = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
rng = random.Random(seed)

def start():
    while True:
        degs = [rng.randint(dmin, dmax) for _ in range(n)]
        if sum(degs) % 2: continue
        pts = [v for v in range(n) for _ in range(degs[v])]
        rng.shuffle(pts)
        m = {}; ok = True
        for i in range(0, len(pts), 2):
            u, v = pts[i], pts[i + 1]
            if u == v: ok = False; break
            key = (min(u, v), max(u, v)); m[key] = m.get(key, 0) + 1
        if ok: return m

def switch(m):
    es = [p for p, k in m.items() if k > 0]
    for _ in range(50):
        (a, b), (c, d) = rng.sample(es, 2)
        if len({a, b, c, d}) < 4: continue
        if rng.random() < 0.5: c, d = d, c
        for x, y in ((a, b), (c, d)):
            m[(min(x, y), max(x, y))] -= 1
        for x, y in ((a, c), (b, d)):
            p = (min(x, y), max(x, y)); m[p] = m.get(p, 0) + 1
        return

def cuts_ok(m):
    for r in range(2, n - 1):
        for A in itertools.combinations(range(n), r):
            As = set(A)
            if sum(k for (u, v), k in m.items() if (u in As) != (v in As)) < 12:
                return False
    return True

t0 = time.time(); tested = good = 0; seen = set()
m = start()
while time.time() - t0 < budget:
    for _ in range(rng.randint(1, 6)): switch(m)
    if rng.random() < 0.02: m = start()
    if not cuts_ok(m): continue
    key = tuple(sorted((p, k) for p, k in m.items() if k))
    if key in seen: continue
    seen.add(key)
    E = [p for p, k in m.items() for _ in range(k)]
    ok, wit = strongly_connected_fast(n, E, 5)
    tested += 1; good += ok
    if not ok:
        print("NOT strongly Z5-connected:", {p: k for p, k in m.items() if k}, "bad beta", wit, flush=True)
print(f"n={n} degrees [{dmin},{dmax}] seed={seed}: distinct graphs meeting the cut condition {tested}, strongly Z5-connected {good}")
