"""Averaging bound for a fixed graph and 2-factor.  Candidate tight cuts: S = (union of whole circuits)
U (at most two arcs), with the exact matching-edge count of its type: pair 6-cut (c1=4, 1 arc),
7-cut (c1=5, 1 arc), 11-cut (c1=7, 2 arcs).  For each S, A(S) = number of 0-edge choices c for which
the arc-end colours are right (both colour 2; for pair cuts also one colour 2 + one colour 0/3) and
|S ∩ Z(c)| is admissible (3 for 7- and 11-cuts; 2 or 4 for pair cuts) -- an upper bound on the number
of c at which S is tight.  Reports  B = sum_S A(S) k(S)  (k = classes killed: 2 for pair cuts, else 1)
against 4 * prod |C_i|.  B < 4N proves that some c has a surviving class (a balanced colouring).
usage: avg_bound.py DUMP"""
import sys, json, itertools, collections
import networkx as nx
d = json.load(open(sys.argv[1])); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
n = G.number_of_nodes(); circ_of = {v: ci for ci, C in enumerate(circuits) for v in C}; pos = {v: i for C in circuits for i, v in enumerate(C)}
cedges = {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
mate = {}
for u, v in G.edges():
    if frozenset((u, v)) not in cedges: mate[u] = v; mate[v] = u
odd = [ci for ci, C in enumerate(circuits) if len(C) % 2]; N = 1
for ci in odd: N *= len(circuits[ci])
# arcs: (circuit, start p, length l) -> vertex set; proper arcs 1 <= l <= L-1
arcs = []
for ci, C in enumerate(circuits):
    L = len(C)
    for p in range(L):
        for l in range(1, L): arcs.append((ci, p, l, frozenset(C[(p + k) % L] for k in range(l))))
def poly_for_circuit(ci, inS_whole, arc):
    """generating polynomial over 0-edge positions j of circuit ci: coefficient list [#j with z outside S, #j with z inside S]
    restricted to j satisfying the arc-end colour condition; arc = (p, l) or None; whole even circuits have no j"""
    C = circuits[ci]; L = len(C)
    if L % 2 == 0:   # even circuit: colours fixed (edge i colour 2 iff i even), no z
        if arc is None: return [1, 0]
        p, l = arc; ok = (((p - 1) % L) % 2 == 0) and ((p + l - 1) % L % 2 == 0)   # both ends colour 2
        ok2 = (((p - 1) % L) % 2 == 0) != ((p + l - 1) % L % 2 == 0)              # one end colour 2, one colour 3
        return [1 if (ok or (ALLOW03 and ok2)) else 0, 0]
    if arc is None: return [0, L] if inS_whole else [L, 0]
    p, l = arc; out = [0, 0]
    for j in range(L):
        def col(m): return 2 if ((m - j) % L) % 2 == 1 else (0 if m == j else 3)
        c_a, c_b = col((p - 1) % L), col((p + l - 1) % L)
        good = (c_a == 2 and c_b == 2) or (ALLOW03 and ((c_a == 2) != (c_b == 2)) and (c_a in (0, 3) or c_b in (0, 3)))
        if not good: continue
        zin = any(C[(p + k) % L] == C[j] for k in range(l)); out[1 if zin else 0] += 1
    return out
def count_choices(S_desc, zcounts):
    """number of 0-edge choices with admissible |S ∩ Z|; S_desc: dict circuit -> 'whole' | ('arc', p, l) | absent (outside)"""
    poly = [1]
    for ci in range(len(circuits)):
        spec = S_desc.get(ci); arc = None if spec in (None, "whole") else (spec[1], spec[2])
        q = poly_for_circuit(ci, spec == "whole", arc)
        new = [0] * (len(poly) + 1)
        for a, x in enumerate(poly):
            new[a] += x * q[0]; new[a + 1] += x * q[1]
        poly = new
    return sum(poly[z] for z in zcounts if z < len(poly))
results = collections.Counter(); B = 0; examples = collections.Counter()
ncirc = len(circuits)
for r_whole in range(0, ncirc + 1):
    for W in itertools.combinations(range(ncirc), r_whole):
        base = set(v for ci in W for v in circuits[ci])
        free = [ci for ci in range(ncirc) if ci not in W]
        for na in (1, 2):
            for arcsel in itertools.combinations([a for a in arcs if a[0] in free], na):
                if len({a[0] for a in arcsel}) < na: continue          # at most one arc per circuit here
                S = base | set().union(*[a[3] for a in arcsel])
                if not S or 2 * len(S) > n: continue                    # each cut once, via its smaller side
                c1 = sum(1 for v in S if mate[v] not in S)
                desc = {ci: "whole" for ci in W}; desc.update({a[0]: ("arc", a[1], a[2]) for a in arcsel})
                if na == 2 and c1 == 7:
                    ALLOW03 = False; A = count_choices(desc, (3,)); results["11-cut"] += 1; B += A; examples["11"] += (A > 0)
                elif na == 1 and c1 == 5:
                    ALLOW03 = False; A = count_choices(desc, (3,)); results["7-cut"] += 1; B += A; examples["7"] += (A > 0)
                elif na == 1 and c1 == 4:
                    ALLOW03 = True; A = count_choices(desc, (2, 3, 4)); results["6-cut"] += 1; B += 2 * A; examples["6"] += (A > 0)
print(f"n={n}, N={N}, 4N={4*N}; candidate sets by type: {dict(results)}; with A>0: {dict(examples)}")
print(f"B = sum A(S) k(S) = {B}  ->  {'B < 4N: averaging closes for this instance' if B < 4 * N else 'B >= 4N: routing factor needed'}  (B/4N = {B/(4*N):.3f})")
