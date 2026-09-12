"""Consistency check: tripled flower snarks J_k (k odd) must have modulo 5-orientations
(all snarks up to 36 vertices are known to have nowhere-zero 5-flows).  Also times the SAT approach."""
import time
from mod5 import has_mod_orientation, triple

def flower(k):
    # vertices: a_i (hub), b_i, c_i, d_i for i in Z_k ; 4k vertices
    a = lambda i: 4 * (i % k); b = lambda i: 4 * (i % k) + 1; c = lambda i: 4 * (i % k) + 2; d = lambda i: 4 * (i % k) + 3
    E = []
    for i in range(k):
        E += [(a(i), b(i)), (a(i), c(i)), (a(i), d(i))]
        E += [(b(i), b(i + 1))]
    # the c/d cycle of length 2k: c_0 c_1 ... c_{k-1} d_0 d_1 ... d_{k-1} c_0
    for i in range(k - 1):
        E += [(c(i), c(i + 1)), (d(i), d(i + 1))]
    E += [(c(k - 1), d(0)), (d(k - 1), c(0))]
    return 4 * k, E

for k in (5, 7, 9, 11, 15, 21):
    n, E = flower(k)
    t = time.time()
    ok5 = has_mod_orientation(n, triple(E), 5)
    ok3 = has_mod_orientation(n, E, 3)
    print(f"J_{k} ({n} vertices): nowhere-zero Z5-flow {ok5}, Z3-flow {ok3}  [{time.time()-t:.2f}s]")
