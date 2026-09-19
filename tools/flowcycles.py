"""Regression of log F(G;5) on short-cycle counts for random 3-connected cubic graphs (experiment D)."""
import sys, os, random, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, networkx as nx
from flowcount import flow_count
n = int(sys.argv[1]) if len(sys.argv) > 1 else 36; N = int(sys.argv[2]) if len(sys.argv) > 2 else 150
rows = []; seed = 0
while len(rows) < N:
    seed += 1; G = nx.random_regular_graph(3, n, seed=77000 + seed)
    if nx.edge_connectivity(G) < 3: continue
    c = {k: 0 for k in range(3, 8)}
    for cy in nx.simple_cycles(G, length_bound=7): c[len(cy)] += 1
    f = flow_count(list(G.edges()), 5, restarts=8, rng=random.Random(1))
    rows.append([1.0] + [c[k] for k in range(3, 8)] + [math.log(f) - n * math.log(1.5)])
A = np.array(rows); X, y = A[:, :-1], A[:, -1]
beta, res, *_ = np.linalg.lstsq(X, y, rcond=None)
pred = X @ beta; r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
print(f'n={n}, {N} random 3-connected cubic graphs: log(F/1.5^n) = {beta[0]:.4f} + sum_k b_k * #(k-cycles),  R^2 = {r2:.4f}')
print(f'   exp(intercept) = {math.exp(beta[0]):.4f}   (prediction q = 5)')
for k, b in zip(range(3, 8), beta[1:]):
    print(f'   k={k}: factor per cycle exp(b_k) = {math.exp(b):.5f}   prediction 1+3(-1/3)^k = {1 + 3 * (-1 / 3) ** k:.5f}   mean count {X[:, k - 2].mean():.2f}')
