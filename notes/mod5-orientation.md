# Route: modulo 5-orientations

*Working notes, started 2026-09-12. Companion to `survey/5-flow-survey.md` §6.6 and §9 (route 3). Code in `tools/`.*

## 0. Why this route

A **modulo $k$-orientation** of a graph is an orientation with $d^+(v)\equiv d^-(v) \pmod k$ at every vertex. For odd $k=2p+1$ this is the same as a circular $(2+\tfrac1p)$-flow, and as a nowhere-zero $\mathbb Z_k$-flow using only the values $\pm1$ (Zhang, Thm 9.2.3). Jaeger's circular flow conjecture (1984) asserted that $4p$-edge-connected graphs have modulo $(2p+1)$-orientations; Han, Li, Wu and Zhang (2018) refuted it for $p\ge 3$, and the cases $p=1$ (the 3-flow conjecture) and $p=2$ are open. The best general theorem is Lovász–Thomassen–Wu–Zhang (2013): $(3k-3)$-edge-connected graphs have modulo $k$-orientations, so **12-edge-connected graphs have modulo 5-orientations.**

The connection to the 5-flow conjecture is tighter than "implies":

> **Lemma A (tripling).** For any graph $G$ let $3G$ be $G$ with every edge replaced by three parallel edges. Then $G$ has a nowhere-zero $\mathbb Z_5$-flow **if and only if** $3G$ has a modulo 5-orientation.
>
> *Proof.* Orient the three copies of $e=uv$; the net contribution at $u$ is $(\#\text{out}-\#\text{in})\in\{3,1,-1,-3\}$, which mod 5 runs over $\{3,1,4,2\}=\mathbb Z_5\setminus\{0\}$ exactly once each. So modulo 5-orientations of $3G$ correspond bijectively (per edge) to nowhere-zero $\mathbb Z_5$-valued functions on $E(G)$ with the same balance condition, i.e. to nowhere-zero $\mathbb Z_5$-flows of $G$. $\square$

Consequences.

- **5FC $\Leftrightarrow$ every $3G$ with $G$ bridgeless has a modulo 5-orientation.** By the survey's ledger it suffices to take $G$ a cyclically 6-edge-connected cubic graph of girth $\ge 11$ and oddness $\ge 6$; then $3G$ is 9-regular and 9-edge-connected.
- Hence (Jaeger) *every 9-edge-connected graph has a modulo 5-orientation* $\Rightarrow$ 5FC, and a fortiori the open $p=2$ case (8-edge-connected) $\Rightarrow$ 5FC.
- The same argument with boundaries: $G$ is $\mathbb Z_5$-connected (JLPT) iff $3G$ is *strongly* $\mathbb Z_5$-connected (every zero-sum $\beta\colon V\to\mathbb Z_5$ is realised as $d^+-d^-$ by some orientation). So the JLPT conjecture "3-edge-connected $\Rightarrow$ $\mathbb Z_5$-connected" is the statement that every $3G$, $G$ 3-edge-connected, is strongly $\mathbb Z_5$-connected.
- A nowhere-zero $\mathbb Z_5$-flow of a cubic graph is the same as an assignment of integer weights $\pm1,\pm3$ to oriented edges such that every vertex has signed sum exactly $\pm5$ (the three weights at a vertex are $\{3,1,1\}$ or $\{3,3,-1\}$ up to a global sign).

So the route is: understand exactly why the LTWZ induction needs 12, and whether the special structure of $3G$ (or of 9-regular graphs whose only small odd cuts are trivial) lets one push through at 9.

## 1. The LTWZ machinery in $\alpha$-form

I could not download the LTWZ paper (the WVU server refuses connections); the following is the refined presentation in Hasanvand, *Modulo orientations with bounded out-degrees* (arXiv:1702.07039, Thm 3.3), which reproduces the LTWZ argument with a set function $\alpha$ ($2\alpha=\tau_{\text{LTWZ}}$ for odd $k$).

**Setup.** $k\ge3$, $p\colon V(G)\to\mathbb Z_k$ a prescription of out-degrees with $\sum_v p(v)\equiv|E|$. For a vertex set $A$ put $p(A)=\sum_{v\in A}p(v)-e(A)$ and let $\alpha(A)\in\{0,\pm\tfrac12,\dots,\pm\tfrac k2\}$ be the rational with $\alpha(A)\equiv p(A)-d(A)/2\pmod k$. So $|\alpha(A)|$ measures how far the prescribed net out-flow of $A$ is from the balanced value $d(A)/2$, on a circle of circumference $k$. Properties: $\alpha(A\cup B)\equiv\alpha(A)+\alpha(B)$ for disjoint $A,B$; $|\alpha(A)|=|\alpha(A^c)|$; $d(A)-2|\alpha(A)|$ is an even integer; lifting two edges at a vertex does not change $|\alpha|$ of any set.

> **Theorem (LTWZ; Hasanvand Thm 3.3).** Let $z_0\in V(G)$ carry a pre-orientation of its incident edges with $d^+(z_0)\equiv p(z_0)$, and let $v_0$ be a minimum-degree vertex with $\alpha(v_0)=0$ if one exists. Suppose
> (i) $d(z_0)\le 2k-2+2|\alpha(z_0)|$, and
> (ii) $d(A)\ge 2k-2+2|\alpha(A)|$ for every $\emptyset\ne A\subsetneq V\setminus z_0$, $A\ne\{v_0\}$.
> Then the pre-orientation extends to a $p$-orientation of $G$ with $|d^+(v)-d(v)/2|\le k-1+|\alpha(v)|$ for all $v$.

Since $|\alpha|\le k/2$ and $|\alpha(A)|=k/2$ forces $d(A)\equiv k\pmod 2$, hypothesis (ii) holds whenever $G$ is $(3k-3)$-edge-connected; that is the source of 12 for $k=5$.

**Skeleton of the proof** (all of it is needed to see where 9-regular graphs break it). Take a counterexample $(G,p,z_0)$ minimising $|V|+|E(G-z_0)|$.

1. *Claim 1.* Every $A\subsetneq V\setminus z_0$ with $|A|\ge2$ has $d(A)\ge 2k+2|\alpha(A)|$ (otherwise contract $A$, extend, then contract $A^c$ into a new $z_0$ and extend inside $A$).
2. *Claim 2.* No vertex other than $z_0$ has $\alpha=0$ (lift a pair of edges at it).
3. *Claim 3.* $G-z_0$ connected, $d(z_0)\ge k$.
4. *Claim 3.A.* **Every vertex $v\ne z_0$ has $d(v)=2k-2+2|\alpha(v)|$ exactly** (otherwise lift a pair of edges at $v$; Claim 1 keeps (ii) true).
5. *Claim 4/5.* All $\alpha(v)$, $v\ne z_0$, have the same sign and $0<|\alpha(v)|<k/2$ (delete an edge $xy$ with $\alpha(x)>0>\alpha(y)$; the deletion lowers $|\alpha|$ at both ends by $\tfrac12$ and Claim 1 absorbs the loss on larger sets).
6. *Claim 6 and the final step.* All edges at $z_0$ point away, $d(z_0)=k+p(z_0)$; replace one arc $z_0x$ by $k-1$ parallel arcs $xz_0$, which flips the sign of $\alpha(x)$ and contradicts Claim 5.

## 2. Where 9-regular graphs break it, exactly

For a **modulo orientation** ($\beta\equiv0$) the prescription is $p(v)\equiv d(v)\cdot\tfrac{k+1}2$, so
$$\alpha(v)\equiv d(v)\cdot\tfrac k2\pmod k=\begin{cases}0,& d(v)\text{ even},\\ \pm k/2,& d(v)\text{ odd},\end{cases}\qquad \alpha(A)=\begin{cases}0,& d(A)\text{ even},\\ \pm k/2,& d(A)\text{ odd}.\end{cases}$$

So for $k=5$ hypothesis (ii) reads: **even cuts need $d(A)\ge 8$, odd cuts need $d(A)\ge 13$**, and (Hasanvand §3.3, generalising LTWZ Thm 4.12) sets with $\alpha(A)=0$ can be dropped entirely: *odd-edge-connectivity $\ge13$ alone gives a modulo 5-orientation.* (Cranston–Li 2018: odd-11 suffices for planar graphs; Cranston–Li–Su–Wang–Wei 2026: 10-edge-connected planar graphs are even strongly $\mathbb Z_5$-connected.)

Now take $G$ a cyclically 6-edge-connected, 3-connected cubic graph and $H=3G$. Cuts of $H$ are triples of cuts of $G$; $d_G(A)\equiv|A|\pmod 2$; and in $G$ the non-cyclic cuts are those with a tree side, of size $|T|+2$. Hence the odd cuts of $H$ have sizes

| $d_H(A)$ | what $A$ is in $G$ |
|---|---|
| 9 | a single vertex |
| 15 | a 3-vertex path (tree side of a 5-cut) |
| $\ge 21$ | cyclic odd cuts (size $\ge7$ by cyclic 6-connectivity) |

and the even cuts have size $\ge 12$ (an edge, $d_G=4$) or $\ge18$ (cyclic). Therefore

> **Observation B.** For $H=3G$ with $G$ a cyclically 6-edge-connected cubic graph, *every* non-singleton vertex set satisfies the LTWZ hypothesis (ii), and even the strengthened Claim-1 level ($d(A)\ge 2k+2|\alpha(A)|$, i.e. $\ge10$ even, $\ge15$ odd). **The only violations are the singletons**: $d(v)=9$ where the induction wants $2k-2+2|\alpha(v)|=13$.

So on the graphs that matter the LTWZ method is short by exactly 4 at every vertex and by nothing anywhere else. Claim 3.A shows this is intrinsic to the induction: minimal counterexamples are driven to have every vertex at degree exactly $13$ (for $|\alpha|=5/2$), and the lifting/deletion steps (Claims 2, 3.A, 4) all consume degree at a vertex. The same phenomenon at $k=3$: 5-regular graphs have $|\alpha(v)|=3/2$ and LTWZ wants $d(v)\ge7$; the 3-flow conjecture lives exactly in that gap of 2. For $k=5$ the gap is 4, but the graphs $3G$ have far more slack on non-trivial cuts than the 5-regular graphs of the 3-flow problem.

**A sharper sub-target.**

> **(T1)** Every 9-regular graph in which every odd cut other than a singleton has $\ge15$ edges and every even cut has $\ge 10$ (or 12) edges admits a modulo 5-orientation.

(T1) implies the 5-flow conjecture (apply it to $3G$). It is strictly stronger, since it also covers non-tripled 9-regular graphs; random 9-regular graphs satisfy its hypothesis a.a.s. and do have modulo 5-orientations a.a.s. (Delcourt–Huq–Prałat 2025), so there is no cheap random counterexample. A counterexample to (T1) that is not of the form $3G$ would be a counterexample to Jaeger's conjecture for $p=2$, which would itself be new.

## 3. A discrepancy reformulation of 9-regular modulo 5-orientations

A 9-regular graph has a modulo 5-orientation iff it has an orientation with every in-degree in $\{2,7\}$. Let $S$ be the set of in-degree-2 vertices; counting edges gives $|S|=n/2$. By Hakimi's orientation theorem such an orientation exists iff $e(X)\le 2|X\cap S|+7|X\setminus S|$ for every $X\subseteq V$. Writing $e(X)=(9|X|-d(X))/2$ this becomes:

> **Lemma C.** A 9-regular graph $H$ has a modulo 5-orientation iff there is $S\subseteq V$ with $|S|=n/2$ such that
> $$\Bigl|\,|X\cap S|-\tfrac{|X|}{2}\,\Bigr|\le \frac{d(X)}{10}\qquad\text{for every } X\subseteq V.$$
> (The lower bound follows from the upper bound applied to $V\setminus X$.) The constraint is vacuous unless $d(X)<5|X|$, i.e. unless $X$ induces average degree $>4$.

For $(4p+1)$-regular graphs and modulo $(2p+1)$ the same computation gives tolerance $d(X)/(2(2p+1))$; for $k=3$ it is $d(X)/6$.

Two readings:

- **Jaeger $p=2$ is tight in an explicit way.** In an 8-edge-connected 9-regular graph, an 8-cut $X$ (necessarily $|X|$ even) must be split *exactly* evenly by $S$; a 9-cut (odd $|X|$) must be split as evenly as an odd set can be. A counterexample to Jaeger $p=2$ is a 9-regular 8-edge-connected multigraph whose family of 8- and 9-cuts admits no simultaneous exact balancing by a half-set. This is a concrete finite search problem (see §5).
- **For $3G$**, $d_H(X)=3d_G(X)$, so the condition is $\bigl||X\cap S|-|X|/2\bigr|\le 0.3\,d_G(X)$: a cyclic 6-cut of $G$ may be unbalanced by at most 1, a 7-, 8- or 9-cut by at most 2, and so on. **The 5-flow conjecture says: every cyclically 6-edge-connected cubic graph has a half-set $S$ that is within $0.3\,d_G(X)$ of balanced on every vertex set $X$.** Sparse-cut large sets are the danger (a random half-set deviates by $\sim\sqrt{|X|}/2$); expanders are easy.

## 4. Experiments (tools/)

`mod5.py` encodes "$(\mathbb Z_k,\beta)$-orientation exists" as CNF (a mod-$k$ counter per vertex) and solves it with CaDiCaL via python-sat. `z5conn.py` decides $\mathbb Z_k$-connectivity of small graphs exactly, as a Minkowski sum over $\mathbb Z_k^V$ (feasible for $k=5$, $|V|\le12$).

Sanity checks (all as expected): Petersen has a nowhere-zero $\mathbb Z_5$-flow via Lemma A and no nowhere-zero $\mathbb Z_3$-flow; $K_{3,3}$ has one; $aK_2$ is strongly $\mathbb Z_5$-connected iff $a\ge4$; $T_{2,2,3}$ and $T_{1,3,3}$ fail with the same bad boundaries as in Cranston–Li; $K_2$ and $C_5$ are not $\mathbb Z_5$-connected.

Results so far (2026-09-12):

| Experiment | Result |
|---|---|
| $\mathbb Z_5$-connectivity of $K_4$, $K_{3,3}$, prism, cube, **Petersen** (exact) | all $\mathbb Z_5$-connected (and $\mathbb Z_6$-connected); $K_2$, $C_5$ correctly rejected |
| all simple 9-regular graphs on 10, 12, 14 vertices (1, 9, 88,193 graphs) | all have modulo 5-orientations |
| all loopless 9-regular multigraphs on 4 and 6 vertices with $\lambda\ge8$ (6 and 453 graphs up to a degree-sequence invariant) | all have modulo 5-orientations |
| tripled flower snarks $3J_k$, $k=5,\dots,21$ (up to 84 vertices) | modulo 5-orientations found in $<0.05$ s each; no $\mathbb Z_3$-flow, as expected |
| **(S1) test**, exhaustive: multigraphs with degrees in $[9,9]$ ($n=4,6$), $[9,11]$ ($n=4$), $[9,10]$ ($n=5$), every set $A$ with $|A|,|A^c|\ge2$ having $d(A)\ge12$ | all strongly $\mathbb Z_5$-connected (1, 22, 31, 31 graphs) |

Here **(S1)** is the statement: *a multigraph with minimum degree $\ge9$ in which every vertex set $A$ with $|A|,|A^c|\ge 2$ has $d(A)\ge12$ is strongly $\mathbb Z_5$-connected.* It is LTWZ's hypothesis (ii) with singletons exempted. (S1) would imply that $3G$ is strongly $\mathbb Z_5$-connected for every cyclically 4-edge-connected cubic $G$, hence the JLPT conjecture for those graphs and the 5-flow conjecture; so it is not going to be easy, but so far no small counterexample exists, and small counterexamples are exactly what a Cranston–Li style proof would need to list.

For the record, one modulo 5-orientation of $3P_{10}$ has in-degree-2 set $S=\{0,2,3,8,9\}$ (outer cycle $0..4$, spokes $i\to i+5$): $S$ and its complement each induce 2 edges, the cut has 11 edges, and Lemma C's constraint $5\le 2.5+0.3\cdot 11$ is met with slack $0.8$.

Running: random 9-regular multigraphs ($n\le 24$, pairing model, 8-edge-connected), Markov-chain samples of 9-regular multigraphs with heavy multiplicities ($n=8,10$), the (S1) test for $n=6$ with degrees in $[9,10]$ and sampled $n=7,8$. Results will be appended below.

## 5. Next steps on this route

1. Finish and record the experiments above; add tripled snarks (all snarks to 36 vertices are known to have 5-flows, so this is only a consistency check of the pipeline).
2. **Targeted search for a $p=2$ counterexample** using Lemma C: build 9-regular 8-edge-connected multigraphs with many crossing 8-cuts (e.g. blow-ups of small base graphs with parallel classes of size 4 and 8-edge "necks") and test with the SAT encoding. Either a counterexample (new: Jaeger $p=2$ false) or a family of hard positive instances that show what a proof must handle.
3. **Attack (T1) theoretically.** The induction has to survive vertices of degree $2k-1$ with $|\alpha|=k/2$. Ideas to try, in order:
   - Find a pre-orientation/lifting step that trades the surplus on non-trivial odd cuts ($\ge15$ vs $13$) for the deficit at singletons; Claim 1 shows contraction of any $A$ with $|A|\ge2$ is "free" on $3G$, so a counterexample to (T1) restricted to $3G$ would have to be minimal in a very strong sense.
   - Work on the Hakimi form of Lemma C directly: choose $S$ greedily/probabilistically on $G$ (cubic), and analyse violated sets $X$ as an uncrossing family; the constraint is only active on sets with $d_G(X) < \tfrac{10}{3}\cdot$ deviation, i.e. sparse cuts, which in a cyclically 6-edge-connected graph are few.
   - Compare with what Cranston–Li do for planar graphs at odd-11: their proof gets from 13 down to 11 by an explicit list of "troublesome partitions" plus discharging; the analogous non-planar list would be the first thing to compute (small 9-regular multigraphs that are not strongly $\mathbb Z_5$-connected, `mod5.strongly_connected`).
4. Obtain the LTWZ and Han–Li–Wu–Zhang papers (WVU server was down); check how the $p\ge3$ counterexamples interact with Lemma C's tolerance $d(X)/(2(2p+1))$ and why the construction does not descend to $p=2$.

## References

See the survey bibliography; additionally: M. Hasanvand, Modulo orientations with bounded out-degrees, arXiv:1702.07039; D. W. Cranston, J. Li, Circular flows in planar graphs, arXiv:1812.09833; D. W. Cranston, J. Li, B. Su, Z. Wang, C. Wei, Orientations of 10-edge-connected planar multigraphs and applications, arXiv:2603.24292; J. Li, B. Su, Z. Wang, C. Wei, Characterization of strongly $\mathbb Z_\ell$-connected graphs of small order, arXiv:2603.21591; M. Delcourt, R. Huq, P. Prałat, Almost all 9-regular graphs have a modulo-5 orientation, arXiv:2210.12103.
