"""Coupling analysis for (L-repair): at every re-pairing one-step move from a failing 0-edge choice (new H three
paths), list ALL killer sets (structured candidates: whole circuits + <=2 arcs) at the old and the new choice,
split the classes killed at the new choice into those killed by OLD killer sets and by NEWBORN sets, and type each
class as 'reach' (triple transversal to both pairings, Lemma 13.1) or 'split' (contains both ends of an old path).
usage: coupling.py DUMP [DUMP ...]   (env DETAIL=1 prints every move with total >= 3)"""
import sys, os, json, itertools, collections
import networkx as nx
from oddness import canonical_colouring, H_components
DETAIL = os.environ.get("DETAIL")
tally = collections.Counter(); typet = collections.Counter(); moves_out = []
def candidates(G, circuits):
    n = G.number_of_nodes(); cedges = {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    mate = {}
    for e in {frozenset(e) for e in G.edges()} - cedges: u, v = tuple(e); mate[u] = v; mate[v] = u
    arcs = [(ci, frozenset(C[(p + k) % len(C)] for k in range(l))) for ci, C in enumerate(circuits) for p in range(len(C)) for l in range(1, len(C))]
    cands = []; ncirc = len(circuits)
    for r in range(ncirc + 1):
        for W in itertools.combinations(range(ncirc), r):
            base = set(v for ci in W for v in circuits[ci]); free = [ci for ci in range(ncirc) if ci not in W]
            for na in (1, 2):
                for sel in itertools.combinations([a for a in arcs if a[0] in free], na):
                    if len({a[0] for a in sel}) < na: continue
                    S = base | set().union(*[a[1] for a in sel])
                    if not S or 2 * len(S) > n: continue
                    c1 = sum(1 for v in S if mate[v] not in S)
                    if (na, c1) in ((2, 7), (1, 5), (1, 4)):
                        bnd = [frozenset((u, v)) for u in S for v in G[u] if v not in S and frozenset((u, v)) in cedges]
                        desc = tuple(sorted([(ci, "W") for ci in W] + [(a[0], len(a[1])) for a in sel]))
                        cands.append((frozenset(S), c1 + 2 * na, bnd, desc))
    return cands
def killers(G, cands, col, comps, odd, zc):
    """{S: set of classes killed}; class = frozenset of odd-circuit indices whose z is black, canonical up to complement."""
    paths = [(colouring, ends) for colouring, ends in comps if ends]
    evens = [(colouring, ends) for colouring, ends in comps if not ends]
    zcirc = {C[j]: k for k, (C, j) in enumerate(zip(odd, zc))}
    allc = frozenset(range(len(odd)))
    out = collections.defaultdict(set)
    for S, dS, bnd, desc in cands:
        cols = [col[e] for e in bnd]
        if sum(1 for c in cols if c == 2) < len(cols) - (1 if dS == 6 else 0): continue
        zin = [zcirc[v] for v in S if v in zcirc]
        if dS != 6 and len(zin) != 3: continue
        if dS == 6 and len(zin) not in (2, 3, 4): continue
        # per component signed count (black - white) under its reference colouring
        sig = []
        for colouring, ends in paths:
            s = sum((1 if colouring[v] else -1) for v in S if v in colouring); sig.append((s, ends))
        se = sum(abs(sum((1 if colouring[v] else -1) for v in S if v in colouring)) for colouring, _ in evens)
        for flips in itertools.product((0, 1), repeat=len(paths)):
            tot = se + sum(abs(s) for s, _ in sig)     # even circuits: best orientation; paths: |.| independent of flip
            # flips matter only through which z's are black: compute the class and the actual signed sum
            tot = se; black = set()
            for (s, ends), f in zip(sig, flips):
                tot += (s if f == 0 else -s)
                for v in ends:
                    colouring = [c for c, e in paths if e == ends][0]
                    if (colouring[v] ^ f): black.add(zcirc[v])
            if 5 * abs(tot) > 3 * dS:
                cl = frozenset(black)
                if len(cl) > 3 or (len(cl) == 3 and min(cl) != min(allc - cl) and (allc - cl) < cl): cl = allc - cl
                if len(cl) == 3: cl = min(cl, allc - cl, key=sorted)
                out[S].add(cl)
    return out
for fn in sys.argv[1:]:
    d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
    M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    odd = [C for C in circuits if len(C) % 2]; circ_of = {v: k for k, C in enumerate(odd) for v in C}
    cands = candidates(G, circuits)
    def analyse(zc):
        col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
        pairing = frozenset(frozenset(circ_of[v] for v in e) for _, e in comps if e)
        return col, comps, pairing, len(comps)
    fails = [list(z) for z in d.get("fails", [d["zero"]])]
    for base in fails:
        col0, comps0, P0, n0 = analyse(base)
        K0 = killers(G, cands, col0, comps0, odd, base)
        classes0 = set().union(*K0.values()) if K0 else set()
        if len(classes0) < 4: print(fn, base, "base not fully killed by structured candidates:", classes0); continue
        for i in range(len(odd)):
            for j in range(len(odd[i])):
                if j == base[i]: continue
                zc = list(base); zc[i] = j; col, comps, P, n = analyse(zc)
                if n != n0 or P == P0 or sum(1 for _, e in comps if e) != 3: continue
                K = killers(G, cands, col, comps, odd, zc)
                def transversal(cl, Pg): return all(len(cl & p) == 1 for p in Pg)
                shared = P & P0
                oldcl = set(); newcl = set(); old_sets = []; new_sets = []
                for S, cls in K.items():
                    (oldcl if S in K0 else newcl).update(cls); (old_sets if S in K0 else new_sets).append((S, cls))
                def ty(cl): return "reach" if transversal(cl, P0) and transversal(cl, P) else "split"
                key = (len(oldcl), len(newcl), len(oldcl | newcl))
                tally[key] += 1
                typet[("old", tuple(sorted(collections.Counter(ty(c) for c in oldcl).items())), "new", tuple(sorted(collections.Counter(ty(c) for c in newcl).items())))] += 1
                if len(oldcl | newcl) >= 4: print("*** FOUR CLASSES at", fn, base, "->", zc)
                if DETAIL and len(oldcl | newcl) >= 3:
                    desc = {S: next(dd for SS, _, _, dd in cands if SS == S) for S in K}
                    print(f"{fn.split('/')[-1]} {base} move circuit {i} -> pos {j}; moved circuit in shared pair: {any(i in p for p in shared)}; P0={sorted(map(sorted,P0))} P={sorted(map(sorted,P))}")
                    for S, cls in old_sets: print("   old ", [sorted(c) for c in cls], [ty(c) for c in cls], desc[S], "contains new z:", odd[i][j] in S, "old z:", odd[i][base[i]] in S)
                    for S, cls in new_sets: print("   NEW ", [sorted(c) for c in cls], [ty(c) for c in cls], desc[S], "contains new z:", odd[i][j] in S, "old z:", odd[i][base[i]] in S, "was tight-shaped at c:", None)
    print(fn.split("/")[-1], "done; running tally (old,new,union):", dict(sorted(tally.items())), flush=True)
print("TYPES:", dict(typet))
