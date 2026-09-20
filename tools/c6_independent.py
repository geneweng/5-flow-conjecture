#!/usr/bin/env python3
"""Independent re-implementation of the count relaxation for (C6).

Written ONLY from notes/c6-model-spec.md (no code shared with count_c6.py,
count_general.py, enum_c6.py, enum_covers.py, validate_c6.py).

Solver: HiGHS MIP (highspy), big-M formulation, integer slack variables for
parities.  Lazy constraints: path/class connectivity (groups 1, 6) and, for
J = 4, the union facts (group 7) for unions of more than 2 atoms.
Every FEASIBLE answer is re-checked by an exact integer checker (`check`)
that is written separately from the model builder and tests ALL constraints
of the spec, including all 2^(2^J) unions.

usage:
  c6_independent.py stats                 items / covers / cache key sanity
  c6_independent.py one NAME NAME ...     solve one configuration, print solution
  c6_independent.py run                   full job (pairs, certificate, sample)
  c6_independent.py report                summary of tools/c6_independent_cache.json
"""
import sys, os, json, time, itertools, random
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, 'enum_c6_cache.json')
MYCACHE = os.path.join(HERE, 'c6_independent_cache.json')
LOG = os.path.join(HERE, 'c6_independent.log')

XUB, UUB, NUB = 11, 5, 500
TLIM = 120.0

# ----------------------------------------------------------------------------
# items (spec section 4)
# ----------------------------------------------------------------------------
CLASSES = [(0, 0), (0, 1), (1, 0), (1, 1)]


def make_items():
    items = {}

    def add(name, memb, kills, c1, c2, nh, cross=None, point_class=None):
        # path status derived from the membership of the ends
        path = []
        for i in range(3):
            s, t = memb[2 * i], memb[2 * i + 1]
            if s != t:
                path.append(('sep', 1 if s else 0))
            elif cross is not None and cross[0] == i:
                path.append(('cross', cross[1]))
            else:
                path.append(('avoid', None))
        items[name] = dict(name=name, memb=tuple(int(v) for v in memb), path=path,
                           c1=c1, c2=c2, nh=nh, kills=frozenset(kills),
                           point_class=point_class)

    for (b2, b3) in CLASSES:
        base = [1, 0, int(b2 == 0), int(b2 == 1), int(b3 == 0), int(b3 == 1)]
        add('P7_%d%d' % (b2, b3), base, [(b2, b3)], 5, 2, False, point_class=(b2, b3))
        add('P641_%d%d' % (b2, b3), base, [(b2, b3)], 4, 1, True, point_class=(b2, b3))
        betas = [1, 1 - b2, 1 - b3]
        for i in range(3):
            for p in range(2):
                m = list(base)
                m[2 * i] = m[2 * i + 1] = p
                add('Pc%dp%d_%d%d' % (i, p, b2, b3), m, [(b2, b3)], 4, 2, False,
                    cross=(i, betas[i]), point_class=(b2, b3))
    for b in range(2):
        for p in range(2):
            add('L01b%dp%d' % (b, p), [1, 0, int(b == 0), int(b == 1), p, p],
                [c for c in CLASSES if c[0] == b], 4, 2, False)
            add('L02b%dp%d' % (b, p), [1, 0, p, p, int(b == 0), int(b == 1)],
                [c for c in CLASSES if c[1] == b], 4, 2, False)
            add('L12b%dp%d' % (b, p), [p, p, 1, 0, int(b == 0), int(b == 1)],
                [c for c in CLASSES if (c[0] == c[1]) == (b == 0)], 4, 2, False)
    return items


ITEMS = make_items()


def compatible(names):
    pc = [ITEMS[n]['point_class'] for n in names if ITEMS[n]['point_class'] is not None]
    return len(pc) == len(set(pc))


def minimal_covers():
    names = sorted(ITEMS)
    out = []
    allc = set(CLASSES)
    for r in range(1, 5):
        for comb in itertools.combinations(names, r):
            if not compatible(comb):
                continue
            ks = [ITEMS[n]['kills'] for n in comb]
            if set().union(*ks) != allc:
                continue
            minimal = True
            for t in range(r):
                rest = set().union(*[ks[s] for s in range(r) if s != t]) if r > 1 else set()
                if rest == allc:
                    minimal = False
                    break
            if minimal:
                out.append(comb)
    return out


def key(names):
    return ' '.join(sorted(names))


# ----------------------------------------------------------------------------
# combinatorial skeleton shared by builder and checker: which variables exist
# ----------------------------------------------------------------------------
class Skeleton:
    def __init__(self, names):
        self.names = list(names)
        self.its = [ITEMS[n] for n in names]
        J = self.J = len(names)
        self.A = 1 << J
        self.NH = [j for j in range(J) if self.its[j]['nh']]
        self.zatom = [sum(self.its[j]['memb'][k] << j for j in range(J)) for k in range(6)]
        self.start = [self.zatom[2 * i] for i in range(3)]
        self.end = [self.zatom[2 * i + 1] for i in range(3)]
        # x arcs
        self.xarcs = {}     # (i,a,b) -> colour
        for i in range(3):
            for a in range(self.A):
                for b in range(self.A):
                    if a == b:
                        continue
                    cols = set()
                    ok = True
                    for j in range(J):
                        if not ((a ^ b) >> j) & 1:
                            continue
                        kind, beta = self.its[j]['path'][i]
                        if kind == 'avoid':
                            ok = False
                            break
                        leaving = (a >> j) & 1
                        cols.add(1 if leaving == beta else 2)
                    if ok and len(cols) == 1:
                        self.xarcs[(i, a, b)] = cols.pop()
        # u arcs: black end a (inside all differing items), white end b
        self.uarcs = []
        for a in range(self.A):
            for b in range(self.A):
                if a != b and (a ^ b) & a == (a ^ b):
                    for c in (1, 2):
                        self.uarcs.append((c, a, b))
        # non-H pairs
        nhmask = sum(1 << j for j in self.NH)
        self.epairs = [(a, b) for a in range(self.A) for b in range(a + 1, self.A)
                       if (a ^ b) & ~nhmask == 0]
        # pairs that can carry colour-2 edges
        c2p = set()
        for (i, a, b), c in self.xarcs.items():
            if c == 2:
                c2p.add((min(a, b), max(a, b)))
        for (c, a, b) in self.uarcs:
            if c == 2:
                c2p.add((min(a, b), max(a, b)))
        self.wpairs = sorted(c2p)


# ----------------------------------------------------------------------------
# MIP model
# ----------------------------------------------------------------------------
class Model:
    def __init__(self):
        self.lb, self.ub, self.integ, self.vname = [], [], [], []
        self.rows = []   # (dict, lo, hi, group)

    def var(self, lb, ub, name, integer=True):
        self.lb.append(lb); self.ub.append(ub); self.integ.append(integer); self.vname.append(name)
        return len(self.lb) - 1

    def row(self, coefs, lo, hi, group):
        d = {}
        for v, c in coefs:
            d[v] = d.get(v, 0) + c
        d = {v: c for v, c in d.items() if c != 0}
        if not d:
            if lo > 0 or hi < 0:
                # trivially infeasible row: keep it via a fixed dummy
                z = self.var(0, 0, 'dummy')
                self.rows.append(({z: 1}, lo, hi, group))
            return
        self.rows.append((d, lo, hi, group))


INF = float('inf')


def lin_add(*terms):
    out = []
    for t in terms:
        out.extend(t)
    return out


def scale(terms, s):
    return [(v, c * s) for v, c in terms]


class C6Model:
    """Builds the system of spec section 3.  `groups` = set of constraint groups to include."""

    def __init__(self, names, groups=frozenset(range(1, 8)), all_unions=None, log=None):
        self.sk = sk = Skeleton(names)
        self.groups = set(groups)
        self.m = m = Model()
        self.log = log
        J, A = sk.J, sk.A
        self.n = [m.var(0, NUB, 'n[%s]' % self.an(a)) for a in range(A)]
        self.x = {k: m.var(0, XUB, 'x%d[%s>%s]c%d' % (k[0], self.an(k[1]), self.an(k[2]), c))
                  for k, c in sk.xarcs.items()}
        self.u = {k: m.var(0, UUB, 'u%d[%s>%s]' % (k[0], self.an(k[1]), self.an(k[2])))
                  for k in sk.uarcs}
        self.e = {(kap, a, b): m.var(0, 1, 'e%d[%s-%s]' % (kap, self.an(a), self.an(b)))
                  for (a, b) in sk.epairs for kap in range(1, 8)}
        self.w = {(kap, a, b): m.var(0, XUB, 'w%d[%s-%s]' % (kap, self.an(a), self.an(b)))
                  for (a, b) in sk.wpairs for kap in range(1, 8)}
        self.union_done = set()
        self.build()
        if 7 in self.groups:
            if all_unions is None:
                all_unions = (J <= 2)
            if all_unions:
                for U in range(1, 1 << A):
                    self.add_union(U)
            else:
                atoms = list(range(A))
                for r in (1, 2):
                    for comb in itertools.combinations(atoms, r):
                        self.add_union(sum(1 << a for a in comb))

    def an(self, a):
        return format(a, '0%db' % self.sk.J)[::-1]   # character j = a_j

    # ---- linear expressions -------------------------------------------------
    def E(self, a, colour):
        t = []
        for (i, p, q), c in self.sk.xarcs.items():
            if c == colour and (p == a or q == a):
                t.append((self.x[(i, p, q)], 1))
        for (c, p, q) in self.sk.uarcs:
            if c == colour and (p == a or q == a):
                t.append((self.u[(c, p, q)], 1))
        return t

    def Enh(self, a):
        return [(v, 1) for (kap, p, q), v in self.e.items() if p == a or q == a]

    def deg(self, kap, a):
        t = [(v, 1) for (k2, p, q), v in self.w.items() if k2 == kap and (p == a or q == a)]
        t += [(v, 1) for (k2, p, q), v in self.e.items() if k2 == kap and (p == a or q == a)]
        return t

    def dU(self, U):
        t = []
        for (i, p, q), v in self.x.items():
            if ((U >> p) & 1) != ((U >> q) & 1):
                t.append((v, 1))
        for (c, p, q), v in self.u.items():
            if ((U >> p) & 1) != ((U >> q) & 1):
                t.append((v, 1))
        for (kap, p, q), v in self.e.items():
            if ((U >> p) & 1) != ((U >> q) & 1):
                t.append((v, 1))
        return t

    # ---- static constraints -------------------------------------------------
    def build(self):
        sk, m, G = self.sk, self.m, self.groups
        J, A = sk.J, sk.A
        if 1 in G:
            for i in range(3):
                for a in range(A):
                    t = []
                    for (i2, p, q), v in self.x.items():
                        if i2 != i:
                            continue
                        if p == a:
                            t.append((v, 1))
                        if q == a:
                            t.append((v, -1))
                    rhs = int(a == sk.start[i]) - int(a == sk.end[i])
                    m.row(t, rhs, rhs, 1)
        if 2 in G:
            for j in range(J):
                for i in range(3):
                    kind, beta = sk.its[j]['path'][i]
                    if kind == 'avoid':
                        continue
                    t = [(v, 1) for (i2, p, q), v in self.x.items() if i2 == i and ((p ^ q) >> j) & 1]
                    s = m.var(0, 10 * XUB * A, 'par2[%d,%d]' % (j, i))
                    off = 1 if kind == 'sep' else 2
                    m.row(t + [(s, -2)], off, off, 2)
        if 3 in G:
            for j in range(J):
                u1 = [(v, 1) for (c, p, q), v in self.u.items() if c == 1 and ((p ^ q) >> j) & 1]
                u2 = [(v, 1) for (c, p, q), v in self.u.items() if c == 2 and ((p ^ q) >> j) & 1]
                m.row(u1 + scale(u2, -1), 0, 0, 3)
                x1 = [(self.x[k], 1) for k, c in sk.xarcs.items() if c == 1 and ((k[1] ^ k[2]) >> j) & 1]
                x2 = [(self.x[k], 1) for k, c in sk.xarcs.items() if c == 2 and ((k[1] ^ k[2]) >> j) & 1]
                m.row(x1 + u1, sk.its[j]['c1'], sk.its[j]['c1'], 3)
                m.row(x2 + u2, sk.its[j]['c2'], sk.its[j]['c2'], 3)
                if sk.its[j]['nh']:
                    t = [(v, 1) for (kap, p, q), v in self.e.items() if ((p ^ q) >> j) & 1]
                    m.row(t, 1, 1, 3)
        if 4 in G:
            for a in range(A):
                t = []
                for (c, p, q), v in self.u.items():
                    sgn = 1 if c == 1 else -1
                    if p == a:
                        t.append((v, sgn))
                    if q == a:
                        t.append((v, -sgn))
                m.row(t, 0, 0, 4)
        if 5 in G:
            for a in range(A):
                na = self.n[a]
                E1, E2, En = self.E(a, 1), self.E(a, 2), self.Enh(a)
                za = sum(1 for k in range(6) if sk.zatom[k] == a)
                s = m.var(0, NUB, 'par5a[%s]' % self.an(a))
                m.row([(na, 1)] + scale(E1, -1) + [(s, -2)], 0, 0, 5)       # n = E1 + 2s (uses E1<=n)
                s = m.var(0, NUB, 'par5b[%s]' % self.an(a))
                m.row([(na, 1)] + scale(E2, -1) + [(s, -2)], za, za, 5)     # n = E2 + z + 2s
                s = m.var(0, NUB, 'par5c[%s]' % self.an(a))
                m.row(E2 + En + [(s, -2)], 0, 0, 5)
                # E2+z<=n and E1<=n are implied by the slack forms above (s>=0); add explicitly anyway
                m.row([(na, 1)] + scale(E2, -1), za, INF, 5)
                m.row([(na, 1)] + scale(E1, -1), 0, INF, 5)
                visits = [(v, 2) for (i, p, q), v in self.x.items() if q == a]
                nstart = sum(1 for i in range(3) if sk.start[i] == a)
                hc = [(v, 1) for (c, p, q), v in self.u.items() if p == a or q == a]
                m.row([(na, 2)] + scale(visits, -1) + scale(hc, -1), 2 * nstart, INF, 5)
                if za > 0:
                    m.row([(na, 1)] + En, 3, INF, 5)
                if nstart == 0:
                    m.row([(na, -1)] + scale(E1 + E2 + En, NUB), 0, INF, 5)
            m.row([(v, 1) for v in self.n], 42, INF, 5)
        if 6 in G:
            for (a, b) in sk.wpairs:
                t = [(self.w[(kap, a, b)], 1) for kap in range(1, 8)]
                for (i, p, q), c in sk.xarcs.items():
                    if c == 2 and {p, q} == {a, b}:
                        t.append((self.x[(i, p, q)], -1))
                for (c, p, q) in sk.uarcs:
                    if c == 2 and {p, q} == {a, b}:
                        t.append((self.u[(c, p, q)], -1))
                m.row(t, 0, 0, 6)
            DM = 8 * A * XUB
            for kap in range(1, 8):
                for a in range(A):
                    d = self.deg(kap, a)
                    if d:
                        s = m.var(0, DM, 'par6[%d,%s]' % (kap, self.an(a)))
                        m.row(d + [(s, -2)], 0, 0, 6)
            for kap in range(1, 7):
                a = sk.zatom[kap - 1]
                na = self.n[a]
                d = self.deg(kap, a)
                t7 = m.var(0, 1, 't7[%d]' % kap)
                m.row([(na, 1), (t7, -7)], 0, INF, 6)          # t7=1 -> n>=7
                m.row(d + [(t7, 2)], 2, INF, 6)                 # t7=0 -> deg>=2
                q = m.var(0, 1, 'q[%d]' % kap)
                g = m.var(0, 1, 'g[%d]' % kap)
                hh = m.var(0, 1, 'h[%d]' % kap)
                m.row([(q, 1), (g, 1), (hh, 1)], 1, INF, 6)
                m.row([(na, 1), (g, -4)], 0, INF, 6)            # g=1 -> n>=4
                m.row(self.Enh(a) + [(hh, -1)], 0, INF, 6)      # h=1 -> Enh>=1
                # q=1 -> deg_kap = 2, other degs 0   (q=0 is possible only if n!=3 or Enh>0,
                # given n+Enh>=3 from group 5; if group 5 is off this is still a relaxation
                # only when n>=3, so add n+Enh>=3 guard through g/h semantics: n<=2,Enh=0 forces q=1,
                # which is stronger than the spec; handle by extra binary)
                if 5 not in G:
                    lo = m.var(0, 1, 'lo[%d]' % kap)            # lo=1 allowed only if n<=2
                    m.rows.pop()  # remove h row, re-add after fixing first row
                    m.rows.pop()
                    m.rows.pop()
                    m.row([(q, 1), (g, 1), (hh, 1), (lo, 1)], 1, INF, 6)
                    m.row([(na, 1), (g, -4)], 0, INF, 6)
                    m.row(self.Enh(a) + [(hh, -1)], 0, INF, 6)
                    m.row([(na, 1), (lo, NUB)], -INF, NUB + 2, 6)
                m.row(d + [(q, DM)], -INF, 2 + DM, 6)
                m.row(d + [(q, -2)], 0, INF, 6)
                for k2 in range(1, 8):
                    if k2 != kap:
                        d2 = self.deg(k2, a)
                        if d2:
                            m.row(d2 + [(q, DM)], -INF, DM, 6)
        if 7 in G:
            # capped sizes mcap[a] = min(n_a, 6)
            self.mcap, self.s6 = [], []
            for a in range(A):
                mc = m.var(0, 6, 'mcap[%s]' % self.an(a))
                s6 = m.var(0, 1, 's6[%s]' % self.an(a))
                na = self.n[a]
                m.row([(mc, 1), (na, -1)], -INF, 0, 7)
                m.row([(mc, 1), (s6, -6)], 0, INF, 7)
                m.row([(na, 1), (s6, -(NUB - 5))], -INF, 5, 7)
                m.row([(mc, 1), (na, -1), (s6, NUB)], 0, INF, 7)
                self.mcap.append(mc); self.s6.append(s6)

    # ---- union facts ----------------------------------------------------------
    def add_union(self, U):
        A = self.sk.A
        full = (1 << A) - 1
        U &= full
        if U == 0 or U == full:
            return False
        canon = min(U, full ^ U)
        if canon in self.union_done:
            return False
        self.union_done.add(canon)
        m = self.m
        d = self.dU(canon)
        ts = {}
        for side, V in (('U', canon), ('C', full ^ canon)):
            atoms = [a for a in range(A) if (V >> a) & 1]
            nV = [(self.n[a], 1) for a in atoms]
            mV = [(self.mcap[a], 1) for a in atoms]
            m.row(scale(nV, 3) + scale(d, -1), 0, INF, 7)       # d <= 3 n_V
            cap = 6 * len(atoms)
            t = []
            for k in range(1, 7):
                tk = m.var(0, 1, 't%d[%s%d]' % (k, side, canon))
                # m_V >= k  ->  t_k = 1
                m.row(mV + [(tk, -(cap - k + 1))], -INF, k - 1, 7)
                t.append(tk)
            # t6 = 1 -> m_V >= 6
            m.row(mV + [(t[5], -6)], 0, INF, 7)
            # 1 <= n_V <= 5 -> d >= n_V + 2
            m.row(d + [(t[0], -3), (t[1], -1), (t[2], -1), (t[3], -1), (t[4], -1), (t[5], 7)], 0, INF, 7)
            ts[side] = t
        m.row(d + [(ts['U'][0], -3), (ts['C'][0], -3)], -3, INF, 7)
        m.row(d + [(ts['U'][3], -6), (ts['C'][3], -6)], -6, INF, 7)
        return True

    # ---- lazy cuts --------------------------------------------------------------
    def add_path_cut(self, i, K):
        """K: set of atoms not containing the start of P_i; arcs inside K need an entry into K."""
        inflow = [(v, XUB) for (i2, p, q), v in self.x.items() if i2 == i and p not in K and q in K]
        for (i2, p, q), v in self.x.items():
            if i2 == i and p in K and q in K:
                self.m.row([(v, -1)] + inflow, 0, INF, 1)

    def add_class_cut(self, kap, K):
        bnd = [(v, XUB) for (k2, p, q), v in self.w.items() if k2 == kap and ((p in K) != (q in K))]
        bnd += [(v, XUB) for (k2, p, q), v in self.e.items() if k2 == kap and ((p in K) != (q in K))]
        for (k2, p, q), v in list(self.w.items()) + list(self.e.items()):
            if k2 == kap and p in K and q in K:
                self.m.row([(v, -1)] + bnd, 0, INF, 6)

    # ---- solve ------------------------------------------------------------------
    def solve_once(self, tlim, seed=0):
        import highspy
        m = self.m
        nv, nr = len(m.lb), len(m.rows)
        h = highspy.Highs()
        h.setOptionValue('output_flag', False)
        h.setOptionValue('threads', 1)
        h.setOptionValue('time_limit', float(max(1.0, tlim)))
        h.setOptionValue('mip_feasibility_tolerance', 1e-9)
        h.setOptionValue('primal_feasibility_tolerance', 1e-9)
        h.setOptionValue('random_seed', seed)
        h.setOptionValue('mip_rel_gap', 1e9)
        h.setOptionValue('mip_abs_gap', 1e9)
        lp = highspy.HighsLp()
        lp.num_col_, lp.num_row_ = nv, nr
        cost = np.zeros(nv)
        if os.environ.get('C6_OBJ', '1') == '1':
            for v in list(self.x.values()) + list(self.w.values()):
                cost[v] = 1.0
        if os.environ.get('C6_PRESOLVE'):
            h.setOptionValue('presolve', os.environ['C6_PRESOLVE'])
        lp.col_cost_ = cost
        lp.col_lower_ = np.array(m.lb, dtype=float)
        lp.col_upper_ = np.array(m.ub, dtype=float)
        lp.row_lower_ = np.array([-highspy.kHighsInf if r[1] == -INF else r[1] for r in m.rows], dtype=float)
        lp.row_upper_ = np.array([highspy.kHighsInf if r[2] == INF else r[2] for r in m.rows], dtype=float)
        starts, idx, val = [0], [], []
        for r in m.rows:
            for v, c in r[0].items():
                idx.append(v); val.append(float(c))
            starts.append(len(idx))
        lp.a_matrix_.format_ = highspy.MatrixFormat.kRowwise
        lp.a_matrix_.start_ = np.array(starts, dtype=np.int32)
        lp.a_matrix_.index_ = np.array(idx, dtype=np.int32)
        lp.a_matrix_.value_ = np.array(val, dtype=float)
        lp.integrality_ = [highspy.HighsVarType.kInteger if t else highspy.HighsVarType.kContinuous
                           for t in m.integ]
        h.passModel(lp)
        h.run()
        st = h.getModelStatus()
        S = highspy.HighsModelStatus
        if st == S.kInfeasible:
            return 'INFEASIBLE', None
        sol = h.getSolution()
        info = h.getInfo()
        if st in (S.kOptimal, S.kObjectiveBound, S.kObjectiveTarget, S.kSolutionLimit) or \
                (st == S.kTimeLimit and info.primal_solution_status == 2):
            vals = [int(round(v)) for v in sol.col_value]
            if max(abs(v - round(v)) for v in sol.col_value) > 1e-5:
                return 'NUMERIC', None
            return 'SOL', vals
        if st == S.kTimeLimit:
            return 'UNKNOWN', None
        return 'STATUS_%s' % st, None

    def extract(self, vals):
        return dict(n=[vals[v] for v in self.n],
                    x={k: vals[v] for k, v in self.x.items() if vals[v]},
                    u={k: vals[v] for k, v in self.u.items() if vals[v]},
                    e={k: vals[v] for k, v in self.e.items() if vals[v]},
                    w={k: vals[v] for k, v in self.w.items() if vals[v]})

    def solve_once_cpsat(self, tlim, seed=0):
        """LAST-RESORT engine: OR-tools CP-SAT on MY row system (the rows are built above,
        independently of the reference implementation; only the search engine is shared)."""
        from ortools.sat.python import cp_model
        m = self.m
        cp = cp_model.CpModel()
        V = [cp.NewIntVar(int(m.lb[i]), int(m.ub[i]), 'v%d' % i) for i in range(len(m.lb))]
        for d, lo, hi, g in m.rows:
            expr = sum(int(c) * V[v] for v, c in d.items())
            if lo != -INF:
                cp.Add(expr >= int(lo))
            if hi != INF:
                cp.Add(expr <= int(hi))
        sv = cp_model.CpSolver()
        sv.parameters.max_time_in_seconds = float(max(1.0, tlim))
        sv.parameters.num_workers = 1
        sv.parameters.random_seed = seed
        r = sv.Solve(cp)
        if r == cp_model.INFEASIBLE:
            return 'INFEASIBLE', None
        if r in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            return 'SOL', [sv.Value(v) for v in V]
        return 'UNKNOWN', None

    def solve(self, tlim=TLIM, engine='highs'):
        if engine == 'cpsat':
            self.solve_once = self.solve_once_cpsat
        t0 = time.time()
        rounds = 0
        while True:
            rounds += 1
            left = tlim - (time.time() - t0)
            if left <= 0:
                return 'UNKNOWN', None, rounds
            st, vals = self.solve_once(left, seed=rounds)
            if st != 'SOL':
                return st, None, rounds
            sol = self.extract(vals)
            viol = check(self.sk, sol, self.groups, want_cuts=True)
            hard = [v for v in viol if v[0] == 'hard']
            if hard:
                return 'NUMERIC:' + hard[0][1], sol, rounds
            if not viol:
                return 'FEASIBLE', sol, rounds
            nun = 0
            for v in viol:
                if v[0] == 'pathconn':
                    self.add_path_cut(v[1], v[2])
                elif v[0] == 'classconn':
                    self.add_class_cut(v[1], v[2])
                elif v[0] == 'union':
                    if nun < 40 and self.add_union(v[1]):
                        nun += 1


# ----------------------------------------------------------------------------
# exact checker, written directly from the spec (independent of the row builder)
# ----------------------------------------------------------------------------
def components(nodes, edges):
    adj = {v: set() for v in nodes}
    for a, b in edges:
        adj[a].add(b); adj[b].add(a)
    seen, comps = set(), []
    for v in nodes:
        if v in seen:
            continue
        st, comp = [v], {v}
        while st:
            y = st.pop()
            for z in adj[y]:
                if z not in comp:
                    comp.add(z); st.append(z)
        seen |= comp
        comps.append(comp)
    return comps


def check(sk, sol, groups=frozenset(range(1, 8)), want_cuts=False):
    """Returns list of violations.  ('hard', text) for a violated static constraint,
    ('pathconn', i, K), ('classconn', kappa, K), ('union', U) for lazy ones."""
    J, A = sk.J, sk.A
    n, x, u, e, w = sol['n'], sol['x'], sol['u'], sol['e'], sol['w']
    out = []

    def hard(msg):
        out.append(('hard', msg))

    # variable domains / existence
    for k, v in x.items():
        if k not in sk.xarcs or not (0 <= v <= XUB):
            hard('x domain %s' % (k,))
    for k, v in u.items():
        if k not in sk.uarcs or not (0 <= v <= UUB):
            hard('u domain %s' % (k,))
    for k, v in e.items():
        if (k[1], k[2]) not in sk.epairs or v not in (0, 1):
            hard('e domain %s' % (k,))
    for k, v in w.items():
        if not (0 <= v <= XUB):
            hard('w domain %s' % (k,))
    if any(v < 0 or v > NUB for v in n):
        hard('n domain')

    def col_edges(a, c):
        tot = 0
        for (i, p, q), v in x.items():
            if sk.xarcs[(i, p, q)] == c and (p == a or q == a):
                tot += v
        for (c2, p, q), v in u.items():
            if c2 == c and (p == a or q == a):
                tot += v
        return tot

    def enh(a):
        return sum(v for (kap, p, q), v in e.items() if p == a or q == a)

    def deg(kap, a):
        return sum(v for (k2, p, q), v in w.items() if k2 == kap and a in (p, q)) + \
            sum(v for (k2, p, q), v in e.items() if k2 == kap and a in (p, q))

    if 1 in groups:
        for i in range(3):
            for a in range(A):
                o = sum(v for (i2, p, q), v in x.items() if i2 == i and p == a)
                ii = sum(v for (i2, p, q), v in x.items() if i2 == i and q == a)
                if o - ii != int(a == sk.start[i]) - int(a == sk.end[i]):
                    hard('flow %d %d' % (i, a))
            arcs = [(p, q) for (i2, p, q), v in x.items() if i2 == i and v > 0]
            nodes = set([sk.start[i]]) | {p for p, q in arcs} | {q for p, q in arcs}
            comps = components(sorted(nodes), arcs)
            reach = [c for c in comps if sk.start[i] in c][0]
            for c in comps:
                if sk.start[i] not in c:
                    out.append(('pathconn', i, frozenset(c)))
            if len(comps) > 1:
                out.append(('pathconn', i, frozenset(set(range(A)) - reach)))
    if 2 in groups:
        for j in range(J):
            for i in range(3):
                kind, beta = sk.its[j]['path'][i]
                cr = sum(v for (i2, p, q), v in x.items() if i2 == i and ((p ^ q) >> j) & 1)
                if kind == 'sep' and cr % 2 != 1:
                    hard('cross sep')
                if kind == 'cross' and (cr % 2 != 0 or cr < 2):
                    hard('cross cross')
                if kind == 'avoid' and cr != 0:
                    hard('cross avoid')
    if 3 in groups:
        for j in range(J):
            u1 = sum(v for (c, p, q), v in u.items() if c == 1 and ((p ^ q) >> j) & 1)
            u2 = sum(v for (c, p, q), v in u.items() if c == 2 and ((p ^ q) >> j) & 1)
            x1 = sum(v for k, v in x.items() if sk.xarcs[k] == 1 and ((k[1] ^ k[2]) >> j) & 1)
            x2 = sum(v for k, v in x.items() if sk.xarcs[k] == 2 and ((k[1] ^ k[2]) >> j) & 1)
            ne = sum(v for (kap, p, q), v in e.items() if ((p ^ q) >> j) & 1)
            if u1 != u2:
                hard('T3 item %d' % j)
            if x1 + u1 != sk.its[j]['c1'] or x2 + u2 != sk.its[j]['c2']:
                hard('c1c2 item %d' % j)
            if ne != (1 if sk.its[j]['nh'] else 0):
                hard('nh item %d' % j)
    if 4 in groups:
        for a in range(A):
            B = {c: sum(v for (c2, p, q), v in u.items() if c2 == c and p == a) for c in (1, 2)}
            W = {c: sum(v for (c2, p, q), v in u.items() if c2 == c and q == a) for c in (1, 2)}
            if B[1] - B[2] != W[1] - W[2]:
                hard('segments atom %d' % a)
    if 5 in groups:
        for a in range(A):
            E1, E2, En = col_edges(a, 1), col_edges(a, 2), enh(a)
            za = sum(1 for k in range(6) if sk.zatom[k] == a)
            if (n[a] - E1) % 2:
                hard('n=E1 mod 2 atom %d' % a)
            if (n[a] - E2 - za) % 2:
                hard('n=E2+z mod 2 atom %d' % a)
            if (E2 + En) % 2:
                hard('E2+Enh even atom %d' % a)
            if E2 + za > n[a] or E1 > n[a]:
                hard('E<=n atom %d' % a)
            visits = sum(v for (i, p, q), v in x.items() if q == a) + sum(1 for i in range(3) if sk.start[i] == a)
            hc = sum(v for (c, p, q), v in u.items() if a in (p, q))
            if 2 * n[a] < 2 * visits + hc:
                hard('visits atom %d' % a)
            if za > 0 and n[a] + En < 3:
                hard('z atom small %d' % a)
            if all(sk.start[i] != a for i in range(3)) and E1 + E2 + En == 0 and n[a] != 0:
                hard('isolated atom %d' % a)
        if sum(n) < 42:
            hard('n>=42')
    if 6 in groups:
        pairs = set((min(p, q), max(p, q)) for (k2, p, q) in w) | set(sk.wpairs)
        for (a, b) in pairs:
            tot = sum(v for k, v in x.items() if sk.xarcs[k] == 2 and {k[1], k[2]} == {a, b}) + \
                sum(v for (c, p, q), v in u.items() if c == 2 and {p, q} == {a, b})
            if tot != sum(w.get((kap, a, b), 0) for kap in range(1, 8)):
                hard('w sum %d %d' % (a, b))
        for kap in range(1, 8):
            for a in range(A):
                if deg(kap, a) % 2:
                    hard('deg parity %d %d' % (kap, a))
        for kap in range(1, 7):
            a = sk.zatom[kap - 1]
            if n[a] <= 6 and deg(kap, a) < 2:
                hard('deg>=2 class %d' % kap)
            if n[a] == 3 and enh(a) == 0:
                if deg(kap, a) != 2 or any(deg(k2, a) for k2 in range(1, 8) if k2 != kap):
                    hard('n=3 rule class %d' % kap)
            edges = [(p, q) for (k2, p, q), v in list(w.items()) + list(e.items()) if k2 == kap and v > 0]
            if edges:
                nodes = set([a]) | {p for p, q in edges} | {q for p, q in edges}
                comps = components(sorted(nodes), edges)
                reach = [c for c in comps if a in c][0]
                for c in comps:
                    if a not in c:
                        out.append(('classconn', kap, frozenset(c)))
                if len(comps) > 1:
                    out.append(('classconn', kap, frozenset(set(range(A)) - reach)))
    if 7 in groups:
        # all unions, vectorised over bitmasks
        NU = 1 << A
        masks = np.arange(NU, dtype=np.int64)
        nU = np.zeros(NU, dtype=np.int64)
        for a in range(A):
            nU += ((masks >> a) & 1) * n[a]
        dU = np.zeros(NU, dtype=np.int64)
        pairw = {}
        for coll in (x, u, e):
            for k, v in coll.items():
                pq = (k[1], k[2])
                pairw[pq] = pairw.get(pq, 0) + v
        for (p, q), v in pairw.items():
            dU += (((masks >> p) & 1) != ((masks >> q) & 1)) * v
        tot = sum(n)
        nC = tot - nU
        req = np.zeros(NU, dtype=np.int64)
        both = (nU >= 1) & (nC >= 1)
        req = np.where(both, 3, req)
        req = np.where(both & (nU <= 5), np.maximum(req, nU + 2), req)
        req = np.where(both & (nC <= 5), np.maximum(req, nC + 2), req)
        req = np.where((nU >= 4) & (nC >= 4), np.maximum(req, 6), req)
        bad = (dU > 3 * nU) | (dU > 3 * nC) | (dU < req)
        bad[0] = bad[NU - 1] = False
        idx = np.nonzero(bad)[0]
        idx = [int(U) for U in idx if U < (NU - 1 - U)]
        idx.sort(key=lambda U: min(bin(U).count('1'), A - bin(U).count('1')))
        for U in idx:
            out.append(('union', U))
    return out


# ----------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------
def solve_config(names, groups=frozenset(range(1, 8)), tlim=TLIM, all_unions=None, engine='highs'):
    t0 = time.time()
    mod = C6Model(sorted(names), groups=groups, all_unions=all_unions)
    st, sol, rounds = mod.solve(tlim, engine=engine)
    if st == 'FEASIBLE':
        # final exact certification against the complete spec (all groups requested)
        assert not check(mod.sk, sol, groups)
    return st, sol, rounds, time.time() - t0, mod


def worker(job):
    names = job
    try:
        st, sol, rounds, dt, mod = solve_config(names)
        sizes = sol['n'] if sol else None
        return key(names), st, sizes, rounds, dt
    except Exception as ex:  # noqa
        import traceback
        return key(names), 'ERROR:' + repr(ex) + traceback.format_exc()[-300:], None, 0, 0.0


def fmt_solution(mod, sol):
    sk = mod.sk
    lines = ['items: ' + ' '.join(sk.names) + '   (atom string: character j = membership in item j)']
    lines.append('path ends: ' + ', '.join('z%d in %s' % (k + 1, mod.an(sk.zatom[k])) for k in range(6)))
    lines.append('n: ' + ', '.join('%s=%d' % (mod.an(a), v) for a, v in enumerate(sol['n']) if v))
    for (i, p, q), v in sorted(sol['x'].items()):
        lines.append('x P%d %s->%s colour %d : %d' % (i, mod.an(p), mod.an(q), sk.xarcs[(i, p, q)], v))
    for (c, p, q), v in sorted(sol['u'].items()):
        lines.append('u colour %d black %s -> white %s : %d' % (c, mod.an(p), mod.an(q), v))
    for (kap, p, q), v in sorted(sol['e'].items()):
        lines.append('e class %d %s-%s : %d' % (kap, mod.an(p), mod.an(q), v))
    for (kap, p, q), v in sorted(sol['w'].items()):
        lines.append('w class %d %s-%s : %d' % (kap, mod.an(p), mod.an(q), v))
    return '\n'.join(lines)


def load_ref():
    c = json.load(open(REF))
    return {k: v for k, v in c.items() if not k.endswith(' T') and not k.startswith('PURE')}


def load_mine():
    if os.path.exists(MYCACHE):
        return json.load(open(MYCACHE))
    return {}


def logmsg(s):
    with open(LOG, 'a') as f:
        f.write(time.strftime('%H:%M:%S ') + s + '\n')


HIGHS_T, CPSAT_T = 45.0, 120.0


def job_proc(names, engine, q):
    try:
        tl = HIGHS_T if engine == 'highs' else CPSAT_T
        st, sol, rounds, dt, mod = solve_config(names, tlim=tl, engine=engine)
        q.put((key(names), engine, st, sol['n'] if sol else None, rounds, dt))
    except Exception as ex:  # noqa
        q.put((key(names), engine, 'ERROR:' + repr(ex), None, 0, 0.0))


def save(mine):
    json.dump(mine, open(MYCACHE + '.tmp', 'w'), indent=0)
    os.replace(MYCACHE + '.tmp', MYCACHE)


def run_jobs(jobs, mine, tag, pool=None, ncpu=4):
    """Own scheduler: one process per solve, hard kill; HiGHS first, CP-SAT (last resort) on UNKNOWN."""
    import multiprocessing as mp
    import queue as qmod
    todo = []
    for j in jobs:
        v = mine.get(key(j))
        if v is None or not v['status'].endswith('FEASIBLE'):
            eng = 'cpsat' if (v is not None and v.get('highs_failed')) else 'highs'
            todo.append((tuple(j), eng))
    logmsg('%s: %d jobs to run' % (tag, len(todo)))
    q = mp.Queue()
    active = {}
    done = 0
    total = len(todo)
    todo.reverse()
    while todo or active:
        while todo and len(active) < ncpu:
            names, eng = todo.pop()
            p = mp.Process(target=job_proc, args=(names, eng, q))
            p.start()
            active[key(names)] = (p, time.time(), names, eng)
        results = []
        try:
            results.append(q.get(timeout=1.0))
            while True:
                results.append(q.get_nowait())
        except qmod.Empty:
            pass
        now = time.time()
        for k, (p, t0, names, eng) in list(active.items()):
            hard = (HIGHS_T if eng == 'highs' else CPSAT_T) + 25
            if now - t0 > hard and not any(r[0] == k for r in results):
                p.kill()
                results.append((k, eng, 'UNKNOWN', None, 0, now - t0))
                logmsg('%s hard-killed %s (%s)' % (tag, k, eng))
        for (k, eng, st, sizes, rounds, dt) in results:
            if k not in active:
                continue
            p, t0, names, _ = active.pop(k)
            p.join(timeout=5)
            if eng == 'highs' and not st.endswith('FEASIBLE'):
                mine[k] = dict(status=st, n=None, rounds=rounds, time=round(dt, 2), stage=tag,
                               engine='highs', highs_failed=True)
                todo.append((names, 'cpsat'))
                logmsg('%s %s -> %s with HiGHS after %.0fs; retry with CP-SAT' % (tag, k, st, dt))
                continue
            rec = dict(status=st, n=sizes, rounds=rounds, time=round(dt, 2), stage=tag, engine=eng)
            if eng == 'cpsat':
                rec['highs_failed'] = True
            mine[k] = rec
            done += 1
            logmsg('%s %d/%d  %s -> %s [%s] rounds=%d  %.1fs  n=%s' % (tag, done, total, k, st, eng, rounds, dt, sizes))
            if done % 10 == 0:
                save(mine)
    save(mine)


def all_pairs():
    names = sorted(ITEMS)
    return [c for c in itertools.combinations(names, 2) if compatible(c)]


def certificate_jobs(mine, ref):
    """For each minimal cover with no INFEASIBLE pair (my results), choose a cached-INFEASIBLE
    triple, else the whole 4-cover."""
    covers = minimal_covers()
    need = set()
    nopair = 0
    unresolved = []
    for cov in covers:
        if any(mine.get(key(p), {}).get('status') == 'INFEASIBLE' for p in itertools.combinations(cov, 2)):
            continue
        if len(cov) == 2:
            unresolved.append(cov)
            continue
        nopair += 1
        subs = [s for s in itertools.combinations(cov, 3) if ref.get(key(s)) == 'INFEASIBLE'
                and mine.get(key(s), {}).get('status', 'INFEASIBLE') == 'INFEASIBLE']
        # prefer one that is already chosen / already solved
        chosen = None
        for s in subs:
            if key(s) in need or mine.get(key(s), {}).get('status') == 'INFEASIBLE':
                chosen = s
                break
        if chosen is None and subs:
            chosen = subs[0]
        if chosen is None:
            if len(cov) == 4 and mine.get(key(cov), {}).get('status', 'INFEASIBLE') == 'INFEASIBLE':
                chosen = cov
            else:
                unresolved.append(cov)
                continue
        need.add(key(chosen))
    return [tuple(k.split()) for k in sorted(need)], nopair, unresolved


def certified(cov, mine):
    for r in range(2, len(cov) + 1):
        for s in itertools.combinations(cov, r):
            if mine.get(key(s), {}).get('status') == 'INFEASIBLE':
                return True
    return False


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'stats'
    if cmd == 'stats':
        print(len(ITEMS), 'items')
        cov = minimal_covers()
        from collections import Counter
        print(len(cov), 'minimal covers', Counter(len(c) for c in cov))
        ref = load_ref()
        pairs = all_pairs()
        print(len(pairs), 'pairs; in ref cache:', sum(key(p) in ref for p in pairs))

        def typ(n):
            return 'L' if n[0] == 'L' else ('Pc' if n.startswith('Pc') else 'P')
        cnt = Counter()
        for p in pairs:
            cnt[(tuple(sorted(typ(n) for n in p)), ref.get(key(p)))] += 1
        for k in sorted(cnt, key=str):
            print(k, cnt[k])
    elif cmd == 'one':
        names = sys.argv[2:]
        groups = frozenset(range(1, 8))
        st, sol, rounds, dt, mod = solve_config(names)
        print(st, 'rounds', rounds, '%.1fs' % dt, 'vars', len(mod.m.lb), 'rows', len(mod.m.rows))
        if sol:
            print(fmt_solution(mod, sol))
        ref = load_ref()
        print('cache says', ref.get(key(names)))
    elif cmd == 'run':
        import multiprocessing as mp
        ref = load_ref()
        mine = load_mine()
        for v in mine.values():
            v.setdefault('engine', 'highs')
        if True:
            pool = None
            run_jobs(all_pairs(), mine, 'pair', pool)
            rnd = random.Random(20260919)
            feas = sorted(k for k, v in ref.items() if v == 'FEASIBLE')
            f3 = [k for k in feas if len(k.split()) == 3]
            sample = rnd.sample(f3, min(150, len(f3)))
            run_jobs([tuple(k.split()) for k in sample], mine, 'sample', pool)
            z3_stage(sizes=(2,))
            mine = load_mine()
            for rnd_pass in range(3):
                jobs, nopair, unresolved = certificate_jobs(mine, ref)
                jobs = [j for j in jobs if key(j) not in mine or not mine[key(j)]['status'].endswith('FEASIBLE')
                        and not mine[key(j)].get('highs_failed')]
                logmsg('certificate pass %d: %d covers without infeasible pair, %d sub-lists to prove, %d unresolved'
                       % (rnd_pass, nopair, len(jobs), len(unresolved)))
                if not jobs:
                    break
                run_jobs(jobs, mine, 'cert', pool)
        logmsg('ALL DONE')
    elif cmd == 'report':
        report()
    elif cmd == 'z3':
        z3_stage(sizes=(2,))
    elif cmd == 'z3one':
        print(z3_solve(sys.argv[2:]))


def report():
    from collections import Counter
    ref = load_ref()
    mine = load_mine()
    print('my cache: %d entries' % len(mine))
    for stage in ('pair', 'cert', 'sample', 'extra'):
        cnt = Counter()
        dis = []
        for k, v in mine.items():
            if v['stage'] != stage:
                continue
            r = ref.get(k, 'ABSENT')
            cnt[(v['status'], r)] += 1
            if v['status'] != r:
                dis.append((k, v['status'], r))
        print(stage, dict(cnt))
        for d in dis:
            print('   DISAGREE', d)
    covers = minimal_covers()
    ok = sum(certified(c, mine) for c in covers)
    print('minimal covers certified from my results alone: %d / %d' % (ok, len(covers)))
    bad = [c for c in covers if not certified(c, mine)]
    for c in bad[:20]:
        print('   not certified:', c)


# ----------------------------------------------------------------------------
# exact-arithmetic second opinion (z3, linear integer arithmetic) on the row system
# ----------------------------------------------------------------------------
def z3_solve(names, tlim=TLIM, groups=frozenset(range(1, 8))):
    """Same lazy loop as C6Model.solve but each solve is done by z3 in exact arithmetic."""
    import z3
    mod = C6Model(sorted(names), groups=groups)
    t0 = time.time()
    rounds = 0
    while True:
        rounds += 1
        m = mod.m
        s = z3.SolverFor('QF_LIA')
        s.set('timeout', int(max(1, tlim - (time.time() - t0)) * 1000))
        V = [z3.Int('v%d' % i) for i in range(len(m.lb))]
        for i, v in enumerate(V):
            s.add(v >= int(m.lb[i]), v <= int(m.ub[i]))
        for d, lo, hi, g in m.rows:
            expr = z3.Sum([int(c) * V[v] for v, c in d.items()])
            if lo == hi:
                s.add(expr == int(lo))
            else:
                if lo != -INF:
                    s.add(expr >= int(lo))
                if hi != INF:
                    s.add(expr <= int(hi))
        r = s.check()
        if r == z3.unsat:
            return 'INFEASIBLE', rounds, time.time() - t0
        if r != z3.sat:
            return 'UNKNOWN', rounds, time.time() - t0
        mdl = s.model()
        vals = [mdl.eval(v, model_completion=True).as_long() for v in V]
        sol = mod.extract(vals)
        viol = check(mod.sk, sol, groups)
        if any(v[0] == 'hard' for v in viol):
            return 'BUG:' + [v for v in viol if v[0] == 'hard'][0][1], rounds, time.time() - t0
        if not viol:
            return 'FEASIBLE', rounds, time.time() - t0
        nun = 0
        for v in viol:
            if v[0] == 'pathconn':
                mod.add_path_cut(v[1], v[2])
            elif v[0] == 'classconn':
                mod.add_class_cut(v[1], v[2])
            elif v[0] == 'union' and nun < 40 and mod.add_union(v[1]):
                nun += 1


def z3_worker(names):
    try:
        st, rounds, dt = z3_solve(names)
    except Exception as ex:  # noqa
        st, rounds, dt = 'ERROR:' + repr(ex), 0, 0.0
    return key(names), st, rounds, dt


def z3_stage(sizes=(2,), limit=None):
    import multiprocessing as mp
    mine = load_mine()
    jobs = [tuple(k.split()) for k, v in sorted(mine.items())
            if len(k.split()) in sizes and 'z3' not in v]
    if limit:
        random.Random(7).shuffle(jobs)
        jobs = jobs[:limit]
    logmsg('z3 stage: %d jobs' % len(jobs))
    with mp.Pool(4) as pool:
        for i, (k, st, rounds, dt) in enumerate(pool.imap_unordered(z3_worker, jobs)):
            mine[k]['z3'] = st
            logmsg('z3 %d/%d %s -> %s (HiGHS %s) rounds=%d %.1fs' % (i + 1, len(jobs), k, st, mine[k]['status'], rounds, dt))
            if i % 10 == 0:
                json.dump(mine, open(MYCACHE + '.tmp', 'w'), indent=0)
                os.replace(MYCACHE + '.tmp', MYCACHE)
    json.dump(mine, open(MYCACHE + '.tmp', 'w'), indent=0)
    os.replace(MYCACHE + '.tmp', MYCACHE)
    logmsg('z3 stage done')


if __name__ == '__main__':
    main()
