"""Per-shape summary of a live-set dump: for each shape (whole circuits, arc lengths): number of live
candidates, pin density A/N, live density, routing factor; plus the exact total T and T/4N.
usage: shape_table.py DUMP LIVEJSON"""
import sys, json, collections
d = json.load(open(sys.argv[1])); L = json.load(open(sys.argv[2])); circuits = d["circuits"]
odd = [C for C in circuits if len(C) % 2]; N = 1
for C in odd: N *= len(C)
def circuit_poly(C, S):
    L_ = len(C); inside = [v in S for v in C]; k = sum(inside)
    if k == 0: return [L_ if L_ % 2 else 1, 0]
    if k == L_: return [0, L_] if L_ % 2 else [1, 0]
    if L_ % 2 == 0:   # even circuit: fixed colours (edge i colour 2 iff i even)
        ok = all(i % 2 == 0 for i in range(L_) if inside[i] != inside[(i + 1) % L_]); return [1 if ok else 0, 0]
    out = [0, 0]
    for j in range(L_):
        def col(m): return 2 if ((m - j) % L_) % 2 == 1 else (0 if m == j else 3)
        if all(col(i) == 2 for i in range(L_) if inside[i] != inside[(i + 1) % L_]): out[1 if inside[j] else 0] += 1
    return out
def A(S, zreq):
    poly = [1]
    for C in circuits:
        q = circuit_poly(C, S); new = [0] * (len(poly) + 1)
        for a, x in enumerate(poly): new[a] += x * q[0]; new[a + 1] += x * q[1]
        poly = new
    return sum(poly[z] for z in zreq if z < len(poly))
def shape(S):
    whole = 0; arcs = []
    for C in circuits:
        k = sum(1 for v in C if v in S)
        if k == len(C): whole += 1
        elif k: arcs.append(k)
    return (whole, tuple(sorted(arcs)))
agg = collections.defaultdict(lambda: [0, 0, 0]); T = 0
for r in L:
    S = set(r["S"]); sh = shape(S); a = A(S, (2, 3, 4) if r["d"] == 6 else (3,)); agg[sh][0] += 1; agg[sh][1] += a; agg[sh][2] += len(r["live"]); T += len(r["live"]) * r["k"]
print(f"{len(L)} live candidates; T = {T}; T/4N = {T/(4*N):.4f}; live per candidate (mean) = {sum(v[2] for v in agg.values())/max(1,len(L)):.1f} = N/{N/max(1,sum(v[2] for v in agg.values())/max(1,len(L))):.0f}")
for sh, (c, a, l) in sorted(agg.items(), key=lambda kv: -kv[1][2]):
    print(f"  {sh}: {c:4d} sets  A/N={a/c/N:.5f}  live/N={l/c/N:.6f}  routing={l/a if a else float('nan'):.4f}")
