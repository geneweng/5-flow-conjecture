"""Stress test of (L-repair): repeatedly (1) climb on a fresh random instance of a profile until a 6-of-8 state,
(2) run the completion search from it, (3) verify any full kill, (4) classify all one-step moves from the failing
choice (repair_test logic) and append to a running tally.   usage: repair_pipeline.py PROFILE HOURS SEED"""
import sys, json, time, random, subprocess, itertools, collections, os, glob
import networkx as nx
from oddness import canonical_colouring, H_components, partition, balanced
prof, hours, seed0 = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]); t_end = time.time() + 3600 * hours
tot = collections.Counter(); found = 0; tried = 0; seed = seed0
def classify(fn):
    d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
    M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    odd = [C for C in circuits if len(C) % 2]; circ_of = {v: k for k, C in enumerate(odd) for v in C}
    def analyse(zc):
        col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
        pairing = frozenset(frozenset(circ_of[v] for v in e) for _, e in comps if e); n = len(comps)
        fail = not any(balanced(G, partition(comps, f)) for f in itertools.product((0, 1), repeat=n))
        return pairing, n, fail
    base = list(d["zero"]); p0, n0, f0 = analyse(base)
    if not f0: return None
    cnt = collections.Counter()
    for i in range(len(odd)):
        for j in range(len(odd[i])):
            if j == base[i]: continue
            zc = list(base); zc[i] = j; p, n, f = analyse(zc)
            kind = "same pairing" if (p == p0 and n == n0) else ("re-paired" if n == n0 else "components changed")
            cnt[(kind, "fails" if f else "repaired")] += 1
    return cnt
while time.time() < t_end:
    seed += 1; tried += 1
    r = subprocess.run(["python3", "climb.py", prof, str(seed), "3000"], capture_output=True, text=True)
    dumps = glob.glob(f"best_climb_{prof}_{seed}.json") + glob.glob(f"fullkill_climb_{prof}_{seed}_*.json")
    kills = [f for f in dumps if f.startswith("fullkill")]
    if not kills and dumps:
        r2 = subprocess.run(["python3", "complete.py", dumps[0], str(seed), "1500"], capture_output=True, text=True)
        kills = glob.glob(f"fullkill_complete_{seed}.json")
    for fn in kills:
        cnt = classify(fn)
        if cnt is None: print(f"[{time.strftime('%H:%M')}] {fn}: not a failure at its own choice?!", flush=True); continue
        found += 1; tot.update(cnt)
        print(f"[{time.strftime('%H:%M')}] new failure #{found} ({fn}): {dict(cnt)}", flush=True)
        if cnt.get(("re-paired", "fails"), 0): print("*** (L-repair) VIOLATION ***", fn, flush=True)
    print(f"[{time.strftime('%H:%M')}] tried {tried} instances, failures {found}; running tally {dict(tot)}", flush=True)
