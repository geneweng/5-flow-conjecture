"""Two-position 0-edge shifts (Posa rotations) from failing choices.  For circuit C with z=a0 and visits (a1a2),(a3a4),(a5a6),
the +2 shift replaces the visit a1a2 by a0a1 and makes a2 the new z (the -2 shift is symmetric).  Claims tested:
 (i) colourings with eps(a0)=eps(a2) are admissible at both choices (identical vertex colourings), so their status is unchanged;
 (ii) the non-shared colourings at c are eps(a0)=eps(a1); flipping pi1 (the a1-side of the visit path) maps them onto the
     non-shared colourings at c'; so both fail iff the Klein orbit {eps, eps^R, eps^pi1, eps^{R pi1}} is all bad;
 (iii) print the tight sets of the four orbit colourings with their decomposition over P, pi1, pi2, R.
usage: shift2.py DUMP [DUMP ...]"""
import sys, os, json, itertools, collections
import networkx as nx
from oddness import canonical_colouring, H_components
from coupling_lib import candidates
def paths_of(G, col):
    H = nx.Graph([tuple(e) for e, c in col.items() if c in (1, 2)]); H.add_nodes_from(G.nodes())
    out = []
    for cc in nx.connected_components(H):
        sub = H.subgraph(cc); ends = [v for v in cc if sub.degree(v) == 1]
        if ends: out.append(nx.shortest_path(sub, ends[0], ends[1]))
        else: out.append(list(cc))
    return out
def bad_sets(G, cands, eps):
    res = []
    for S, dS, bnd, desc in cands:
        dlt = sum(1 if eps[v] else -1 for v in S)
        if 5 * abs(dlt) > 3 * dS: res.append((S, dS, desc, dlt))
    return res
def proper_on(eps, col):
    return all(eps[u] != eps[v] for e, c in col.items() if c in (1, 2) for u, v in [tuple(e)])
stats = collections.Counter()
for fn in sys.argv[1:]:
    d = json.load(open(fn)); G = nx.Graph(); G.add_edges_from(map(tuple, d["edges"])); circuits = d["circuits"]
    M = {frozenset(e) for e in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    odd = [C for C in circuits if len(C) % 2]; cands = candidates(G, circuits)
    circ_of = {v: k for k, C in enumerate(odd) for v in C}
    for base in [list(z) for z in d.get("fails", [d["zero"]])]:
        col0 = canonical_colouring(G, M, circuits, base); comps0 = paths_of(G, col0)
        if any(len(c) and c[0] == c[-1] for c in comps0): pass
        for i, C in enumerate(odd):
            L = len(C); j = base[i]
            for sgn in (+1, -1):
                a0, a1, a2 = C[j], C[(j + sgn) % L], C[(j + 2 * sgn) % L]
                zc = list(base); zc[i] = (j + 2 * sgn) % L
                col1 = canonical_colouring(G, M, circuits, zc); comps1 = paths_of(G, col1)
                # identify P (contains a0 as end), P'' (contains edge a1a2), R
                P = next(c for c in comps0 if a0 in c)
                Pp = next(c for c in comps0 if a1 in c)
                if a2 not in Pp or abs(Pp.index(a1) - Pp.index(a2)) != 1: print("unexpected: a1a2 not an H-edge"); continue
                if Pp is P or len(comps0) != 3 or len(comps1) != 3: continue
                if not all(c[0] in {C[k] for C in odd for k in range(len(C))} for c in comps0): pass
                R = next(c for c in comps0 if c is not P and c is not Pp)
                # check pairings differ and no cycles
                pair0 = {frozenset(circ_of[c[0]] // 1 for _ in [0]) for c in []}
                ia1, ia2 = Pp.index(a1), Pp.index(a2)
                pi1 = Pp[:ia1 + 1] if ia1 < ia2 else Pp[ia1:]
                pi2 = Pp[ia2:] if ia1 < ia2 else Pp[:ia2 + 1]
                assert set(pi1) | set(pi2) == set(Pp) and not set(pi1) & set(pi2)
                # reference colouring eps0 admissible at c with eps(a0)=eps(a1)=black, R's end R[0] black
                ref = {}
                for comp in comps0:
                    for k, v in enumerate(comp): ref[v] = k % 2 == 0
                eps0 = dict(ref)
                if not eps0[a0]:
                    for v in P: eps0[v] = not eps0[v]
                if not eps0[a1]:
                    for v in Pp: eps0[v] = not eps0[v]
                assert eps0[a0] and eps0[a1] and not eps0[a2] and proper_on(eps0, col0)
                if not eps0[R[0]]: R = R[::-1]
                orbit = {}
                for X in ((), ("R",), ("pi1",), ("R", "pi1")):
                    e = dict(eps0)
                    if "R" in X:
                        for v in R: e[v] = not e[v]
                    if "pi1" in X:
                        for v in pi1: e[v] = not e[v]
                    orbit[X] = e
                for X, e in orbit.items():
                    adm0, adm1 = proper_on(e, col0), proper_on(e, col1)
                    expect = ("pi1" not in X, "pi1" in X)
                    if (adm0, adm1) != expect: print("ADMISSIBILITY MISMATCH", fn, base, i, sgn, X, adm0, adm1)
                # shared colourings: eps with eps(a0)=eps(a2): check identical admissibility
                bads = {X: bad_sets(G, cands, e) for X, e in orbit.items()}
                nb = sum(1 for X in bads if bads[X]); stats[("orbit bad count", nb)] += 1
                if nb == 4: print("*** ORBIT ALL BAD", fn, base, i, sgn)
                hdr = f"{fn.split('/')[-1]} {base} C{i} shift {sgn:+d}: pieces |P|={len(P)} |pi1|={len(pi1)} |pi2|={len(pi2)} |R|={len(R)}; orbit bad: " + str({"".join(X) or "e": len(bads[X]) for X in bads})
                print(hdr)
                for X, e in orbit.items():
                    for S, dS, desc, dlt in bads[X]:
                        dec = tuple(sum(1 if e[v] else -1 for v in S if v in set(part)) for part in (P, pi1, pi2, R))
                        # arcs on C: boundary circuit edges on C
                        Ev = P[-1] if P[0] == a0 else P[0]; E1v = pi1[-1] if pi1[0] == a1 else pi1[0]; E2v = pi2[-1] if pi2[0] == a2 else pi2[0]
                        role = {i: "C", circ_of[Ev]: "E", circ_of[E1v]: "1", circ_of[E2v]: "2", circ_of[R[0]]: "F", circ_of[R[-1]]: "F'"}
                        Sb = S if dlt > 0 else set(G.nodes()) - S          # black-heavy side
                        parts = []
                        for k, Ck in enumerate(odd):
                            inside = [v in Sb for v in Ck]
                            if all(inside): parts.append(role[k] + ":W")
                            elif any(inside):
                                # arc boundary vertices: inside vertices with a circuit neighbour outside; report piece and colour
                                bv = []
                                for q, v in enumerate(Ck):
                                    if inside[q] and (not inside[q - 1] or not inside[(q + 1) % len(Ck)]):
                                        pc = "P" if v in set(P) else "1" if v in set(pi1) else "2" if v in set(pi2) else "R"
                                        bv.append(pc + ("b" if e[v] else "w"))
                                parts.append(role[k] + f":{sum(inside)}[{','.join(bv)}]")
                        print(f"    {''.join(X) or 'e':5s} d={dS} |delta|={abs(dlt)} |black side|={len(Sb)} {desc} black side: {' '.join(parts)}")
    print(fn.split("/")[-1], "stats so far:", dict(stats), flush=True)
