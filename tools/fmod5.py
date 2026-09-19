"""What determines F(G;5) mod 5?  For all connected cubic graphs of order n (geng), s = F(G;5) * 2^(-n/2) mod 5 lies
in {0,1,4} (exp-D, §6).  Tabulate s against candidate invariants mod 5."""
import sys, subprocess, collections, itertools
import networkx as nx, numpy as np
from flowcount import F
from flowdata import g6_edges
n = int(sys.argv[1])
out = subprocess.run(["geng", "-c", "-d3", "-D3", str(n)], capture_output=True, text=True).stdout.split()
def pm_count(G):
    # perfect matchings by simple recursion
    def rec(vs, adj):
        if not vs: return 1
        v = min(vs); tot = 0
        for w in adj[v]:
            if w in vs: tot += rec(vs - {v, w}, adj)
        return tot
    return rec(frozenset(G.nodes()), {v: set(G[v]) for v in G})
inv2 = pow(2, -1, 5); tabs = collections.defaultdict(collections.Counter)
for s6 in out:
    G = nx.Graph(g6_edges(s6)[1]); f5 = F(G, 5); f4 = F(G, 4); f3 = F(G, 3); f6 = F(G, 6); f7 = F(G, 7)
    s = f5 * pow(inv2, n // 2, 5) % 5
    A0 = np.array([[1.0 if G.has_edge(i, j) else 0.0 for j in range(n)] for i in range(n)]); L = 3 * np.eye(n) - A0; tau = round(np.linalg.det(L[1:, 1:]))
    A = A0; detA = round(np.linalg.det(A)); detA1 = round(np.linalg.det(A - np.eye(n))); detA2 = round(np.linalg.det(A + 2 * np.eye(n)))
    pm = pm_count(G)
    for name, val in [("F4", f4), ("F4/6", f4 // 6), ("tau", tau), ("pm", pm), ("detA", detA), ("det(A-I)", detA1), ("det(A+2I)", detA2), ("F6", f6), ("F7", f7), ("F3", f3)]:
        tabs[name][(val % 5, s)] += 1
print(n, len(out), "graphs")
for name, t in tabs.items():
    det = all(len({s for (v, s) in t if v == val}) == 1 for val in range(5))
    sq = all(s == val * val % 5 for (val, s) in t)
    print(f"{name:10s} determines s: {det}   s == val^2: {sq}   ", dict(sorted(t.items())))
