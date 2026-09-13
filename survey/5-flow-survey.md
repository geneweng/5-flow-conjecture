# Tutte's 5-Flow Conjecture: A Survey

*Prepared for the 5-flow-conjecture project, September 2026. Primary sources: C.-Q. Zhang, "Integer Flows and Cycle Covers of Graphs" (1997), the Open Problem Garden entry, and the post-1997 literature cited in Section 10.*

**Conjecture (Tutte, 1954).** Every bridgeless graph admits a nowhere-zero 5-flow.

The conjecture is 72 years old. It is true for planar graphs (where it is the Five Colour Theorem in disguise), it is tight (the Petersen graph has no nowhere-zero 4-flow), and the best general result is Seymour's 6-flow theorem from 1981. Nobody has moved the general bound since. What has moved is the description of what a counterexample would have to look like, the list of graph classes where the conjecture is verified, and the web of stronger and equivalent statements around it. This survey collects all of that, explains the techniques behind each result, identifies why the techniques stop at 6, and ends with a candid list of attack routes.

---

## Contents

1. The conjecture and where it comes from
2. Foundations: what a nowhere-zero flow is, and equivalent forms
3. What a minimal counterexample must look like
4. What has been proved: upper bounds
5. What has been proved: graph classes
6. Stronger, weaker and equivalent statements
7. The toolbox, and why every tool stops at 6
8. The flow polynomial and computation
9. Attack routes for this project
10. Bibliography

---

## 1. The conjecture and where it comes from

### 1.1 Flows as the dual of face colourings

Tutte introduced integer flows in 1949-1954 as the dual of map colouring. Take a bridgeless plane graph $G$, colour its faces properly with colours $\{1,\dots,k\}$, orient every edge so the smaller colour is on its right, and weight each edge by the difference of the two face colours. Around every vertex the signed weights sum to zero, and every weight lies in $\{1,\dots,k-1\}$. That is a nowhere-zero $k$-flow. The construction reverses for plane graphs, so:

> **Theorem (Tutte 1954; Zhang Thm 1.4.5).** A bridgeless plane graph is face-$k$-colourable if and only if it admits a nowhere-zero $k$-flow.

Two consequences frame everything that follows.

- The Four Colour Theorem is exactly the statement "every bridgeless planar graph has a nowhere-zero 4-flow".
- The Five Colour Theorem (Heawood 1890) is "every bridgeless planar graph has a nowhere-zero 5-flow".

The flow statement, unlike the colouring statement, makes sense for non-planar graphs. Tutte asked which colouring theorems survive that generalisation, and proposed three conjectures in the 1950s-60s that are still the spine of the subject (Zhang §1.1):

| Conjecture | Statement | Status (Sept 2026) |
|---|---|---|
| 5-flow (1954) | Every bridgeless graph has a nowhere-zero 5-flow | **Open** |
| 4-flow (1966) | Every bridgeless graph with no Petersen minor has a nowhere-zero 4-flow | Open in general; **cubic case (the edge-3-colouring conjecture) announced proved in 2026** |
| 3-flow (1972) | Every 4-edge-connected graph has a nowhere-zero 3-flow | Open; true for 6-edge-connected graphs (Lovász-Thomassen-Wu-Zhang 2013) |

### 1.2 Why 5

The Petersen graph $P_{10}$ is bridgeless and has no nowhere-zero 4-flow (a cubic graph has a 4-flow iff it is 3-edge-colourable, and $P_{10}$ is not). So 5 is the smallest constant that could possibly work for all bridgeless graphs. Tutte's first conjecture was only that *some* finite $k$ works; Jaeger and Kilpatrick proved $k=8$ (1975-76) and Seymour proved $k=6$ (1981). The gap between 5 and 6 is the whole problem.

### 1.3 Why it matters

The 5-flow conjecture sits at the centre of a cluster of famous problems: the cycle double cover conjecture, the 5-cycle-double-cover and orientable 5-CDC conjectures, the Berge-Fulkerson conjecture, shortest cycle cover bounds, group connectivity, circular flows, and the theory of snarks. Several of these would imply it (Section 6). Its restriction to cubic graphs is a statement about snarks, the same objects that obstruct the Four Colour Theorem's generalisations. A proof would be the first non-trivial *general* upper bound on flow numbers since 1981; a counterexample would have to be a snark of a kind that has never been constructed (Section 3).

---

## 2. Foundations

Throughout, graphs may have loops and parallel edges. A *circuit* is a connected 2-regular graph; a *cycle* is a subgraph in which every vertex has even degree (an even subgraph); a *bridge* is an edge in no circuit.

### 2.1 Definitions

Let $D$ be an orientation of $E(G)$ and $f\colon E(G)\to\Gamma$ a weight into an abelian group. $(D,f)$ is a **flow** if at every vertex $v$ the sum over out-arcs equals the sum over in-arcs. It is an **integer $k$-flow** if $\Gamma=\mathbb Z$ and $|f(e)|<k$ for all $e$, and **nowhere-zero** if $f(e)\ne 0$ for all $e$. A **modular $k$-flow** (or $\mathbb Z_k$-flow) only requires the balance condition modulo $k$.

Elementary facts (Zhang §1.2):

- A bridge always carries weight 0, so bridgelessness is necessary.
- Reversing an arc and negating its weight preserves flows, so orientation is irrelevant: "$G$ has a nowhere-zero $k$-flow" is a property of the undirected graph, and one may always take the flow positive.
- A nowhere-zero $k$-flow is a nowhere-zero $h$-flow for all $h\ge k$ (this monotonicity fails for group connectivity, see §6.5).
- $G$ has a nowhere-zero 2-flow iff $G$ is a cycle (every degree even). The edges of odd weight in any integer flow form a cycle.

### 2.2 Modular flows and group flows: the value set does not matter

> **Theorem (Tutte 1949; Zhang Thm 1.3.3-1.3.4).** $G$ has a nowhere-zero $k$-flow iff $G$ has a nowhere-zero $\mathbb Z_k$-flow. In fact any modular $k$-flow $(D,f)$ lifts to an integer $k$-flow $(D,f')$ with $f'\equiv f \pmod k$ under the same orientation.

The lift is a potential argument: choose representatives in $\{0,\dots,k-1\}$, then repeatedly push $k$ units along a directed path from a source to a sink; the total imbalance strictly decreases (Younger's proof, Zhang Lemma 1.3.5).

> **Theorem (Tutte 1954; Zhang Thm 2.2.3).** For any abelian group $\Gamma$ of order $k$: $G$ has a nowhere-zero $\Gamma$-flow iff $G$ has a nowhere-zero $k$-flow. Moreover the *number* of nowhere-zero $\Gamma$-flows depends only on $|\Gamma|$.

So the 5-flow conjecture is exactly: **every bridgeless graph has a nowhere-zero $\mathbb Z_5$-flow.** This is the form in which almost all modern work is done. Because $\mathbb Z_5$ is the only group of order 5, there is no product decomposition available for 5, unlike $6=2\cdot 3$ and $8=2\cdot2\cdot2$. Section 7 explains why this single fact is the central obstruction.

### 2.3 Products and sums of flows

> **Theorem (Zhang Thm 2.1.2).** $G$ has a nowhere-zero $k_1k_2$-flow iff $G$ has a $k_1$-flow $f_1$ and a $k_2$-flow $f_2$ (same orientation) with $\operatorname{supp}(f_1)\cup\operatorname{supp}(f_2)=E(G)$. Then $k_2f_1+f_2$ is nowhere-zero.

> **Theorem (Matthews 1978; Zhang Thm 2.4.2).** $G$ has a nowhere-zero $2^r$-flow iff $E(G)$ is covered by $r$ cycles.

> **Theorem (Little-Tutte-Younger 1988; Zhang Thm 2.6.2).** $G$ has a positive $k$-flow $(D,f)$ iff $D(G)$ has $k-1$ directed cycles covering every arc $e$ exactly $f(e)$ times. Equivalently every nonnegative $k$-flow is a sum of $k-1$ nonnegative 2-flows.

The first two are the engines of the 8-flow and 6-flow theorems. The third says the conjecture is equivalent to: *every bridgeless graph has an orientation whose arcs are covered by four directed cycles.*

### 2.4 Orientations with bounded cut ratios

> **Theorem (Hoffman's circulation theorem; Zhang Cor 2.3.2).** $G$ has a nowhere-zero $k$-flow iff $G$ has an orientation $D$ such that for every edge cut $(A,B)$,
> $$\frac{1}{k-1}\le\frac{|[A,B]_D|}{|[B,A]_D|}\le k-1.$$

So the 5-flow conjecture says: **every bridgeless graph can be oriented so that no cut is more than 4:1 unbalanced.** This is also the definition that generalises to real $r$ (circular flows, §6.3).

### 2.5 Face colourings and orientable cycle double covers

For a bridgeless graph $G$ and integer $k$ consider (Zhang Thm 2.5.3):

- (P1) $G$ has a nowhere-zero $k$-flow;
- (P2) $G$ has an orientable $k$-cycle double cover (a CDC by at most $k$ cycles, each orientable as a directed cycle, so that the two cycles through each edge traverse it oppositely);
- (P3) $G$ is face-$k$-colourable on some orientable surface.

In general (P3) $\Rightarrow$ (P2) $\Rightarrow$ (P1); for $k\in\{2,3,4\}$ (P1) $\Leftrightarrow$ (P2); for planar graphs all three coincide. Whether (P1) $\Rightarrow$ (P2) for $k=5$ is precisely the orientable 5-CDC conjecture (§6.4), which together with the 5-flow conjecture would make (P1) $\Leftrightarrow$ (P2) universally.

### 2.6 The flow polynomial

$F(G;k)$, the number of nowhere-zero $\Gamma$-flows for $|\Gamma|=k$, is a polynomial in $k$ satisfying $F(G;k)=F(G/e;k)-F(G\setminus e;k)$ for non-loop $e$, $F(G;k)=(k-1)F(G\setminus e;k)$ for a loop, and $F=0$ if $G$ has a bridge (Zhang §2.7). It is an evaluation of the Tutte polynomial: $F(G;k)=(-1)^{|E|-|V|+c}\,T_G(0,1-k)$. The 5-flow conjecture is: **$F(G;5)>0$ for every bridgeless $G$.** For the Petersen graph, $F(P_{10};q)=(q-1)(q-2)(q-3)(q-4)(q^2-5q+10)$, so $F(P_{10};4)=0$ and $F(P_{10};5)=240$.

---

## 3. What a minimal counterexample must look like

The standard strategy is to assume a counterexample minimal with respect to $|V|+|E|$ and derive structure. The classical reductions are in Zhang §2.8; the modern ones are due mainly to Kochol, Steffen and Mazzuoccolo. Together they give the following ledger.

| Property of a smallest counterexample $G$ | Source |
|---|---|
| simple, 3-edge-connected, hence 3-connected | Tutte/Jaeger/Seymour (Zhang Lemma 2.8.4) |
| cubic (so: a 3-connected cubic graph) | Jaeger 1979 (Zhang Lemma 2.8.6), via vertex splitting |
| not 3-edge-colourable, i.e. a snark | cubic + 3-edge-colourable $\Rightarrow$ 4-flow |
| no non-trivial edge cut of size $\le 3$ | Sekine-Zhang 1997 (Zhang Lemma 2.8.8) |
| girth $\ge 7$ | Möller-Carstens-Brinkmann 1988, Celmins 1984, Jensen (Zhang Lemma 2.8.11 gives girth $\ge 2k-3$) |
| cyclically 5-edge-connected | Celmins 1984 |
| cyclically 6-edge-connected | Kochol 2004 |
| girth $\ge 9$ | Kochol 2006 |
| girth $\ge 11$ | Kochol 2010 |
| oddness $\ge 6$ | Mazzuoccolo-Steffen 2017 (cyclically 6-edge-connected + oddness $\le 4$ $\Rightarrow$ 5-flow) |
| cyclic connectivity $\le \tfrac52\,\omega(G)-4$ | Steffen 2010 |
| $\mu_2(G)\ge 3$ (any two perfect matchings share $\ge 3$ edges) | Steffen 2015 |
| order $\ge 38$ | all snarks up to 36 vertices have circular flow number $\le 5$ (Brinkmann-Goedgebeur-Hägglund-Markström 2013; Goedgebeur-Mattiolo-Mazzuoccolo 2021) |
| not planar, not projective-planar, orientable genus $\ge 3$, non-orientable genus $\ge 5$ | Heawood; Steinberg 1984; Möller-Carstens-Brinkmann 1988 |
| no Hamilton path, indeed no 2-factor with $\le 2$ odd circuits | Jaeger 1979/1988 |
| contains a Petersen minor | RST: cubic graphs of girth $\ge 6$ have Petersen minors; also Kochol 1999 |

Here the **oddness** $\omega(G)$ of a cubic graph is the minimum number of odd circuits in a 2-factor (even, and $\ge 2$ for snarks), and $\mu_2(G)$ is the minimum size of the intersection of two perfect matchings.

**How the classical reductions work (Zhang §2.8).**

- *Parallel edges, degree-2 vertices, loops, 2-edge-cuts* are removed by deletion/suppression/contraction; the 2-flow on a digon lets one adjust a flow on the smaller graph.
- *Cubic*: a vertex of degree $\ge 4$ can be split into two vertices without creating a bridge (Fleischner's splitting lemma, Zhang Thm A.5.2); a flow on the split graph collapses back. Note this step does not preserve minor-closed classes, which is why the 4-flow conjecture is *not* known to reduce to the cubic case, while the 5-flow conjecture is.
- *Small cuts*: for a cut $T$ with $|T|\le 3$ separating $M_1,M_2$, the flow polynomial factorises, $F(G;k)=F(G/M_1;k)F(G/M_2;k)/F(K_2^{|T|};k)$, so one of the two contractions is a smaller counterexample. The factorisation fails for cuts of size $\ge 4$, which is why cyclic connectivity beyond 4 needed new ideas.
- *Girth $\ge 2k-3$*: contract alternate edges of a short circuit $C$, split the resulting degree-4 vertices, take a $k$-flow on the smaller graph, lift it to $G$ (zero possibly only on $C$), and add a suitable multiple of the 2-flow on $C$; a counting argument shows some multiple avoids all zeros when $|C|\le 2k-4$.

**How Kochol's reductions work.** Kochol's "Polynomials associated with nowhere-zero flows" (2002) treats a graph with a $k$-edge-cut as two $k$-poles glued together, and studies, for each $k$-pole, the set of boundary value vectors in $\mathbb Z_5^k$ (summing to zero) that extend to a nowhere-zero $\mathbb Z_5$-flow. The 5-flow question for the whole graph becomes whether the two sets intersect. For cyclic $k$-cuts with $k\le 5$ he shows the answer is forced by the smaller sides, giving the cyclically-6-edge-connected reduction (2004). The girth bounds (9 in 2006 by counting, 11 in 2010) compare the rank of a matrix indexed by boundary vectors with the rank of a submatrix; the 2011 papers with Krivoňáková, Smejová and Šranková reduce the size of the matrices so the computation is feasible. The method is in principle iterable to larger girth at increasing computational cost, but it cannot finish the job on its own: Kochol (1996) constructed cyclically 5-edge-connected snarks of arbitrarily large girth, so no finite girth bound excludes all snarks.

**How the oddness results work.** A 2-factor with only two odd circuits $C_1,C_2$ gives a 5-flow directly: add an edge $e$ joining them, contract the 2-factor, the result is Eulerian, so $G+e$ has a nowhere-zero $\mathbb Z_2\times\mathbb Z_2$-flow, i.e. a 4-flow; then Jaeger's lemma (Zhang Ex. 5.2: if $G+e$ has a nowhere-zero 4-flow and $G$ is bridgeless, then $G$ has a nowhere-zero 5-flow) finishes. Steffen (2010) and Mazzuoccolo-Steffen (2017) push this to oddness 4 by pairing odd circuits along paths and exploiting the freedom given by cyclic 6-edge-connectivity to repair the zeros; the extension to oddness 6 is open and is a natural target (Section 9).

**What is not known to exist.** No snark is currently known that is simultaneously cyclically 6-edge-connected, of girth $\ge 11$, and of oddness $\ge 6$. Cyclically 6-edge-connected snarks exist (Kochol 1996, order 118; smaller ones since), snarks of large girth exist (Kochol 1996, cyclic connectivity 5), and snarks of large oddness exist, but the combination has not been built. Jaeger and Swart conjectured in 1980 that no snark has cyclic connectivity $>6$ (open) and that no snark has girth $>6$ (refuted by Kochol). If the first Jaeger-Swart conjecture is true, a minimal counterexample has cyclic connectivity exactly 6.

---

## 4. What has been proved: upper bounds

### 4.1 The 8-flow theorem (Jaeger 1975-79, Kilpatrick 1975)

> Every bridgeless graph has a nowhere-zero 8-flow.

*Proof sketch (Zhang §5.2).* Reduce to 3-edge-connected $G$. Doubling every edge gives a 6-edge-connected graph, which by Nash-Williams-Tutte contains three edge-disjoint spanning trees $T_1,T_2,T_3$. Each spanning tree contains a parity subgraph $P_i$ (a spanning subgraph with $d_{P_i}(v)\equiv d_G(v)$ mod 2 for all $v$), and the complements $G\setminus P_i$ are cycles. The three parity subgraphs have empty common intersection in $G$ (each original edge lies in at most two of the trees), so the three cycles cover $E(G)$, and Matthews' theorem gives a $2^3$-flow. A second proof uses a perfect matching meeting every 3-cut exactly once (Edmonds' matching polytope) plus Jaeger's 4-flow theorem for 4-edge-connected graphs.

### 4.2 The 6-flow theorem (Seymour 1981)

> Every bridgeless graph has a nowhere-zero 6-flow. Equivalently, a nowhere-zero $\mathbb Z_2\times\mathbb Z_3$-flow.

*Proof sketch (Zhang §5.3).* Let $\mathcal C_k$ be the graphs buildable from a single vertex by repeatedly adding a circuit that uses at most $k$ new edges. The key lemma:

> **Lemma (Seymour; Jaeger's form, Zhang Lemma 5.3.3).** If $G\in\mathcal C_{k-1}$ then for *every* prescription $c\colon E(G)\to\mathbb Z_k$ there is a $\mathbb Z_k$-flow $f$ with $f(e)\ne c(e)$ on every edge. In particular $G$ has a nowhere-zero $\mathbb Z_k$-flow.

The proof is induction on the circuits: a new circuit with at most $k-1$ new edges has $k-1$ "forbidden" values on its new edges, so some multiple of its 2-flow avoids all of them. Then:

> **Lemma (Seymour; Zhang Lemma 5.3.5).** Every 3-edge-connected graph $G$ has a cycle $S$ with $G/S\in\mathcal C_2$. (Younger: $S$ may be taken to be a *necklace*; Fan 1992: a union of vertex-disjoint circuits, for all bridgeless $G$.)

Take $S$ maximal so that $G/S \in \mathcal C_2$; if some component $H$ of $G-V(S)$ remains, 3-edge-connectivity gives two edges from a bridgeless piece $H'$ of $H$ to $S$, and two edge-disjoint paths in $H'$ between their ends; adding this "handle" enlarges $S$. Now $G/S\in\mathcal C_2$ has a nowhere-zero $\mathbb Z_3$-flow, which lifts to a 3-flow $f_2$ on $G$ that is nonzero off $S$; a 2-flow $f_1$ supported on $S$ completes a product $3f_1+f_2$, a nowhere-zero 6-flow.

Three recent re-proofs (DeVos-Rollová-Šámal 2017; DeVos 2024; DeVos-Nurse 2025) shorten this but keep its shape: find a $\mathbb Z_3$-flow whose zero set is an even subgraph.

### 4.3 Why neither proof yields 5

Both proofs are *products*: $8=2\cdot2\cdot2$ and $6=2\cdot3$ come from covering $E(G)$ by supports of small flows and combining them via Theorem 2.1.2. Since 5 is prime, $\mathbb Z_5$ has no non-trivial direct-product decomposition, and no covering argument of this type can give 5. The only route to 5 through Seymour's lemma is to show $G\in\mathcal C_4$ directly, but a graph of girth $\ge 5$ is never in $\mathcal C_4$ (its first circuit already needs $\ge 5$ new edges), so this fails for exactly the graphs that matter. Section 7 returns to this.

---

## 5. What has been proved: graph classes

| Class | Result | Reference |
|---|---|---|
| Planar graphs | Five Colour Theorem; also a two-line proof from girth $\ge 7$ vs Euler's formula | Heawood 1890; Zhang Thm 5.1.2 |
| Projective-planar graphs | 5-flow | Steinberg 1984 |
| Orientable genus $\le 2$, non-orientable genus $\le 4$ | 5-flow, computer-assisted; a minimal counterexample on a fixed surface has all faces of length $\ge 7$ | Möller-Carstens-Brinkmann 1988 |
| Graphs with a Hamilton path | 5-flow | Jaeger 1979 |
| $G+e$ has a 4-flow, $G$ bridgeless; or $G-e$ has a 4-flow | $G$ has a 5-flow | Jaeger; Celmins (Zhang Ex. 5.2, 5.3) |
| Cubic graphs of oddness $\le 2$ | 5-flow | Jaeger 1988 |
| Apex graphs (one vertex whose removal leaves a planar graph) | 5-flow | Gerards-Seymour, unpublished (Zhang Ex. 5.4) |
| 4-edge-connected graphs | 4-flow (two edge-disjoint spanning trees) | Jaeger 1979 |
| Bridgeless graphs with no Petersen minor, cubic | 4-flow (edge-3-colourable) | Robertson-Seymour-Thomas 1997 + Sanders-Seymour + Edwards-Sanders-Seymour-Thomas 2016 + Inoue-Kawarabayashi-Matsuo-Miyashita-Mohar-Sonobe 2026 (preprint); 5-flow earlier by Kochol 1999 |
| Cyclically 6-edge-connected cubic, oddness $\le 4$ | 5-flow | Mazzuoccolo-Steffen 2017 |
| Cyclically $k$-edge-connected cubic with $k\ge\tfrac52\omega-3$ | 5-flow | Steffen 2010 |
| Cyclically 6-edge-connected cubic with $\mu_2\le 2$; or cyclically $(5\mu_2-3)$-edge-connected cubic | 5-flow | Steffen 2015 |
| All snarks on $\le 36$ vertices | circular flow number $\le 5$, hence 5-flow | Brinkmann-Goedgebeur-Hägglund-Markström 2013; Goedgebeur-Mattiolo-Mazzuoccolo 2021 |
| 3-edge-connected graphs | $\mathbb Z_6$-connected (stronger than 6-flow) | Jaeger-Linial-Payan-Tarsi 1992 |
| 12-edge-connected graphs | modulo 5-orientation, hence circular flow number $\le 5/2$ | Lovász-Thomassen-Wu-Zhang 2013 |

**Remark on the 2026 edge-colouring announcement.** Robertson, Seymour and Thomas (1997) reduced Tutte's edge-3-colouring conjecture (every bridgeless cubic graph with no Petersen minor is 3-edge-colourable) to apex and doublecross graphs. The doublecross case was published in 2016. The apex case was long announced by Sanders and Seymour but unpublished; a preprint of August 2026 by Inoue, Kawarabayashi, Matsuo, Miyashita, Mohar and Sonobe gives a constructive, partly computer-checked proof and describes itself as the final piece. If it holds up, cubic Petersen-minor-free graphs have nowhere-zero 4-flows, and a minimal counterexample to the 5-flow conjecture must contain a Petersen minor for a second, independent reason.

---

## 6. Stronger, weaker and equivalent statements

### 6.1 Equivalent forms

For a bridgeless graph $G$, the following are all equivalent to "$G$ has a nowhere-zero 5-flow":

1. $G$ has a nowhere-zero $\mathbb Z_5$-flow (Tutte).
2. $G$ has an orientation in which every edge cut has out/in ratio in $[\tfrac14, 4]$ (Hoffman).
3. $G$ has an orientation whose arcs are covered by four directed cycles (Little-Tutte-Younger).
4. $F(G;5)>0$.
5. $G$ has a real-valued flow with all $|f(e)|\in[1,4]$ (circular flow number $\le 5$; Goddyn-Tarsi-Zhang 1998, and Steffen 2001 for integrality when the circular flow number is exactly 5).
6. (cubic $G$) $G$ has an *even $(1,2)$-factor*: a spanning subgraph $F$ with all degrees 1 or 2 such that every circuit of $G$ has an even number of $F$-balanced edges, where an edge is $F$-balanced if it lies in $F$ or its ends have the same $F$-degree (Matamala-Zamora 2013).

There is deliberately no item of the form "a 2-flow and a 3-flow covering $E(G)$": 5 is prime, so no product formulation exists (§4.3).

And the conjecture as a whole is equivalent to each of:

- every 3-connected cubic graph has a nowhere-zero 5-flow (Zhang §2.8);
- every cyclically 6-edge-connected snark of girth $\ge 11$ and oddness $\ge 6$ has a nowhere-zero 5-flow (Kochol 2004, 2010; Mazzuoccolo-Steffen 2017);
- every bridgeless graph has circular flow number at most 5 (§6.3).

### 6.2 Consequences of the conjecture

- **Shortest cycle covers.** If $G$ has a nowhere-zero 5-flow then $G$ has a cycle cover by 3 cycles of total length $\le\tfrac85|E(G)|$ (Jamshy-Raspaud-Tarsi 1989; Zhang Thm 8.7.1). The proof splits a positive 5-flow into four nonnegative 2-flows and averages over the five flows $f, 5f_i-f$. The general bound is $\tfrac53|E|$ (Alon-Tarsi, Bermond-Jackson-Jaeger); the conjectured bound is $\tfrac75|E|$.
- **Complexity.** Kochol (1998): if the 5-flow conjecture is false, deciding whether a cubic graph has a nowhere-zero 5-flow is NP-complete. So the conjecture is either true or computationally hard, with no middle ground.

### 6.3 Circular flows

A nowhere-zero circular $r$-flow ($r$ real) is a real flow with $1\le|f(e)|\le r-1$; the **circular flow number** $\Phi_c(G)$ is the infimum of such $r$. Goddyn, Tarsi and Zhang (1998) showed the infimum is attained and rational, and by Hoffman it equals $\max_{(A,B)} (|[A,B]|+|[B,A]|)/\min(|[A,B]|,|[B,A]|)$ minimised over orientations. Facts:

- $\Phi_c(G)\le\Phi(G)$ (the integer flow number), and if $\Phi_c(G)$ is an integer $k$ then $G$ has a nowhere-zero integer $k$-flow (Steffen 2001). Hence the 5-flow conjecture is exactly $\Phi_c(G)\le 5$ for all bridgeless $G$.
- Cubic: $\Phi_c=3$ iff bipartite; $\Phi_c\le 4$ iff 3-edge-colourable; no cubic graph has $\Phi_c\in(3,4)$; every rational in $[4,5]$ is the circular flow number of some snark (Lukot'ka-Škoviera 2011, answering Pan-Zhu).
- $\Phi_c(P_{10})=5$. Mohar conjectured $P_{10}$ is the only snark with $\Phi_c=5$ (the "strong circular 5-flow conjecture"); Máčajová and Raspaud (2006) refuted it with an infinite family, all of cyclic connectivity 4. Esperet-Mazzuoccolo-Tarsi (2016) characterised the structure of graphs with $\Phi_c\ge 5$ via sets of transferable flow values modulo 5, gave several constructions, and showed recognising $\Phi_c\ge 5$ is NP-complete. Goedgebeur-Mattiolo-Mazzuoccolo (2021) unified the constructions and computed all snarks with $\Phi_c=5$ up to 36 vertices (98 of order 36).
- **Open (Steffen's Problem 5.2):** is $P_{10}$ the only *cyclically 5-edge-connected* snark with circular flow number 5? A positive answer for cyclically 6-edge-connected snarks would prove the 5-flow conjecture in the strong form $\Phi_c<5$ for every snark except Petersen.
- Bounding $\Phi_c$ from above is also how the computer verifications work: an exact algorithm for $\Phi_c$ of snarks (Goedgebeur-Mattiolo-Mazzuoccolo 2020) is fast enough for all snarks up to 36 vertices.

### 6.4 Cycle double covers

- **Orientable 5-CDC conjecture (Archdeacon 1984, Jaeger 1988; Zhang Conj 2.5.4):** every bridgeless graph has an orientable cycle double cover by at most 5 cycles. By Theorem 2.5.3 this **implies** the 5-flow conjecture.
- **5-CDC conjecture (Preissmann 1981, Celmins 1984; Zhang Conj 7.4.1):** every bridgeless graph has a (not necessarily orientable) CDC by at most 5 cycles. Neither this nor the CDC conjecture is known to imply or follow from the 5-flow conjecture.
- **Petersen colouring conjecture (Jaeger 1988)** implied the 5-CDC and Berge-Fulkerson conjectures. It was **refuted in August 2026**: Putman found a 112-vertex counterexample by SAT solving, and Jooken gave a human-checkable proof and announced a 52-vertex counterexample and an infinite cyclically 4-edge-connected family. This does not touch the 5-flow conjecture (there was no known implication in either direction), but it removes one of the "master conjectures" from the landscape and is a reminder that snark conjectures do fail.

### 6.5 Group connectivity (Jaeger-Linial-Payan-Tarsi 1992)

$G$ is **$\Gamma$-connected** if for every orientation and every zero-sum boundary function $b\colon V\to\Gamma$ there is a nowhere-zero $f\colon E\to\Gamma$ with $\partial f=b$; equivalently, every prescribed function $c\colon E\to\Gamma$ can be avoided pointwise by some $\Gamma$-flow. Being $\Gamma$-connected implies having a nowhere-zero $\Gamma$-flow (take $b=0$). Results (Zhang §9.5):

- graphs with two edge-disjoint spanning trees are $\Gamma$-connected for all $|\Gamma|\ge 4$; so 4-edge-connected graphs are $\mathbb Z_5$-connected;
- 3-edge-connected graphs are $\Gamma$-connected for all $|\Gamma|\ge 6$ (the group-connectivity form of Seymour's theorem; the closure lemma is Zhang Lemma 9.5.5);
- monotonicity fails: there is a graph that is $\mathbb Z_5$-connected but not $\mathbb Z_6$-connected (Zhang Fig. 9.1), so the 5 vs 6 gap is not an artefact of a monotone parameter.

**Conjecture (JLPT; Zhang Conj 9.7.6):** every 3-edge-connected graph is $\mathbb Z_5$-connected. This implies the 5-flow conjecture. It is open, and no partial result beyond 4-edge-connectivity is known to the author.

### 6.6 Modulo orientations and Jaeger's circular flow conjecture

A **modulo $(2t+1)$-orientation** has $d^+(v)\equiv d^-(v)$ mod $2t+1$ at every vertex; it is the same as a nowhere-zero $\mathbb Z_{2t+1}$-flow with all values $\pm1$, and also as a circular $(2+\tfrac1t)$-flow and a circular orientable $(2t+1)$-CDC (Zhang Thm 9.2.3).

- **Jaeger's conjecture (1984):** every $4t$-edge-connected graph has a modulo $(2t+1)$-orientation. $t=1$ is the 3-flow conjecture. $t=2$ (every 8-edge-connected graph has a mod-5 orientation) **implies the 5-flow conjecture** (Zhang Ex. 9.2): triple every edge of a 3-edge-connected cubic graph to get a 9-edge-connected graph; a mod-5 orientation of it gives each original edge a value in $\{\pm1,\pm3\}$, nonzero mod 5.
- Thomassen (2012) proved $(2k^2+k)$-edge-connectivity suffices for mod $k$; Lovász-Thomassen-Wu-Zhang (2013) improved this to $(3k-3)$, so **12-edge-connected graphs have mod-5 orientations** and 6-edge-connected graphs have 3-flows.
- Han-Li-Wu-Zhang (2018) **disproved** Jaeger's conjecture for all $t\ge3$ (there are $4t$-edge-connected graphs with no mod $(2t+1)$-orientation, and $(4t+1)$-edge-connected ones for $t\ge5$). The cases $t=1$ and $t=2$ survive. Delcourt et al. (2025) show almost all 9-regular graphs have mod-5 orientations.

So "every 8-, 9-, 10- or 11-edge-connected graph has a modulo 5-orientation" is an open strengthening of the 5-flow conjecture that is squarely in range of the LTWZ machinery, and Section 9 treats it as a live route.

### 6.7 Recent variants (2025-2026)

- **Reconfiguration.** Esperet, Hendrey, Lagoutte, Marseloo, Norin and Steiner (Dec 2025) show that the reconfiguration analogue of the 5-flow conjecture (any two nowhere-zero 5-flows of a 2-edge-connected graph connected by single-cycle changes) is **false**, in both the $\mathbb Z_5$ and integer settings, while all nowhere-zero $\mathbb Z_2^8$-flows are connected. Cranston, Li, Su, Wang and Xu (June 2026) characterise 3-flow-connectedness and reduce $A$-flow-connectedness to cubic graphs.
- **Multidimensional flows.** Gáborik, Kurz, Mazzuoccolo, Rajník and Rieg (Oct 2025) define Manhattan and Chebyshev flow numbers $\Phi^1_d,\Phi^\infty_d$ (rational; in dimension 2 they separate class 1 from class 2 cubic graphs) and propose conjectures built on "$t$-flow-pairs" abstracted from Seymour's proof, which they position as possibly stronger than the 5-flow conjecture.
- **Signed graphs.** Bouchet's conjecture (every flow-admissible signed graph has a nowhere-zero 6-flow) is the signed analogue; 8-flows for 3-edge-connected signed graphs and 6-flows for cyclically 5-edge-connected cubic signed graphs were shown in 2025-26. Signed graphs have their own 5-flow questions (signed ladders, signed Eulerian graphs), but they are a different problem.

---

## 7. The toolbox, and why every tool stops at 6

It is worth stating plainly what each known technique can and cannot do.

**1. Product/covering arguments (Matthews, Seymour).** Produce flows in $\mathbb Z_{k_1}\times\mathbb Z_{k_2}$. Cannot produce 5. This is not a limitation of cleverness but of arithmetic.

**2. Closure/extension lemmas (Seymour's $\mathcal C_k$, JLPT $k$-closure, Zhang Lemma 9.5.5).** Adding a circuit with $\le j$ new edges to a $\Gamma$-connected piece keeps it $\Gamma$-connected as long as $j<|\Gamma|$ forbidden values can be dodged, i.e. $j\le|\Gamma|-1$. With $\Gamma=\mathbb Z_5$ one can add circuits with at most 4 new edges. Seymour's trick is to first contract a cycle $S$ so that the rest is in $\mathcal C_2$; with $\mathbb Z_5$ one may hope to contract *less*. A direct approach: find an even subgraph $S$ and a $\mathbb Z_5$-flow on $G/S$ that is nonzero on $E(G)\setminus E(S)$ *and* whose lift can be corrected on $S$; the correction needs a $\mathbb Z_5$-flow supported on $S$ avoiding one forbidden value per edge of $S$, which exists if $S$ is $\mathbb Z_5$-connected as a graph, e.g. if $S$ is a union of disjoint circuits **plus** enough freedom. The failure point is that a circuit's 2-flow gives only one degree of freedom (the multiplier), so on a circuit of $S$ carrying several distinct forbidden values one cannot always avoid all of them. Fan's structure theorem (Zhang Thm 5.4.1) that $S$ can be chosen as disjoint circuits with a "last expanded edge" is exactly the kind of extra control that a $\mathbb Z_5$ argument would need more of.

**3. Minimal-counterexample surgery (splitting, contraction, girth, small cuts, Kochol's $k$-pole sets).** Very effective at narrowing the target, structurally incapable of closing it: girth and connectivity bounds can be pushed but snarks with arbitrarily large girth and with cyclic connectivity 6 exist.

**4. 2-factor / oddness arguments (Jaeger, Steffen, Mazzuoccolo).** Work by pairing odd circuits of a 2-factor and constructing a 4-flow on a modified graph, then repairing. They give the sharpest "class" results and have a clear next step (oddness 6), but oddness is unbounded on snarks, so this route needs a new idea to become general.

**5. Modulo-orientation / LTWZ lifting.** Proves mod $k$ orientations under high edge-connectivity by a strengthened inductive statement about prescribed boundaries, contracting and lifting. Yields mod-5 orientations at 12-edge-connectivity. Because tripling a cubic graph gives 9-edge-connectivity, the conjecture would follow from a "9-edge-connected" version, and the Han-Li-Wu-Zhang counterexamples do not touch $t=2$.

**6. Computation.** Exhaustive over snarks to 36 vertices; circular flow numbers to 36; Kochol's rank computations to girth 11. Growth of the snark census (there are 60,167,732 snarks on 36 vertices, and tens of billions expected at 38-40) makes brute force beyond ~38 expensive but not impossible with restriction to cyclically 6-edge-connected, girth-$\ge 7$ snarks, which are rare.

**7. Flow polynomial / algebraic.** The conjecture is $F(G;5)>0$. Known real flow roots cluster near 5 from both sides (§8), so no zero-free interval argument is available; the conjecture is "almost false" in the sense of Jacobsen-Salas.

The honest summary: every general upper-bound proof factorises the group; the conjecture needs a $\mathbb Z_5$ argument; the only $\mathbb Z_5$ arguments in the literature are local (extension through small cuts, circuits with $\le 4$ new edges, pairs of odd circuits) and each is blocked by a family of snarks that escapes it.

---

## 8. The flow polynomial and computation

**Real flow roots.** Welsh conjectured that $F(G;q)>0$ for all real $q\ge4$ (a broad generalisation of the 5-flow conjecture through the Birkhoff-Lewis conjecture for planar duals). Haggard, Pearce and Royle found the generalised Petersen graph $G(16,6)$ has real flow roots near 4.025 and 4.233; they then conjectured $F(G;q)>0$ for $q\ge5$. Jacobsen and Salas (2013, "Is the five-flow conjecture almost false?") disproved that too: the families $G(6n,6)$ and $G(7n,7)$ have real flow roots accumulating at 5 from above ($G(119,7)$ has a root at $\approx 5.0000198$), and another accumulation point near 5.2353. The current guesses are $F(G;q)>0$ for $q\ge6$ (Jacobsen-Salas) or merely for $q\ge c$ for some constant (Dong 2020); Jackson proved $F(G;q)>0$ for $q\ge2\log_2 n$. None of the graphs with roots near 5 is a counterexample (they have Hamilton paths, hence 5-flows), but they show that the value 5 has no safety margin.

**Deciding 5-flows on a given graph.** For cubic $G$, a nowhere-zero $\mathbb Z_5$-flow is a choice of values in $\{1,2,3,4\}$ on edges with the mod-5 balance at each vertex; with one edge of each vertex fixed the others are determined, so a search over a spanning tree's complement of $4^{|E|-|V|+1}$ assignments (or a SAT/CSP encoding) is straightforward for small graphs, and dynamic programming over a path or tree decomposition handles moderate size. For circular flow numbers, Goedgebeur-Mattiolo-Mazzuoccolo (2020) give an exact algorithm based on the Hoffman cut characterisation and Kochol-style boundary sets. Snark generation is available from the House of Graphs / snarkhunter (Brinkmann-Goedgebeur).

---

## 9. Attack routes for this project

Ordered from most incremental to most ambitious. Each item lists what a success would mean and what it needs.

1. **Oddness 6.** Extend Mazzuoccolo-Steffen to cyclically 6-edge-connected cubic graphs of oddness 6 (or all oddness with cyclic connectivity $\ge$ some function). The pairing-and-repair method is explicit; the case analysis is the cost. Even a partial result (oddness 6 with girth $\ge 11$) sharpens the ledger.

2. **Circular flow number below 5 for cyclically 5- or 6-edge-connected snarks.** Steffen's Problem 5.2. Computationally: extend the census of snarks with $\Phi_c=5$ to 38-40 vertices restricted to cyclic connectivity $\ge5$; theoretically: use the Esperet-Mazzuoccolo-Tarsi transferable-value structure to show a cyclic 5- or 6-cut cannot carry the obstruction. A theorem "$\Phi_c(G)<5$ for cyclically 6-edge-connected snarks other than $P_{10}$" would prove the conjecture in strengthened form.

3. **Mod-5 orientations at lower edge-connectivity.** Push the LTWZ bound $3k-3=12$ down toward 9 for $k=5$. Any bound $\le 9$ proves the 5-flow conjecture via edge tripling. Jaeger's conjecture is dead for $t\ge3$ but the $t=2$ case has no known obstruction, and the Delcourt et al. random-regular result is encouraging. This is the route with the most recent, live machinery.

4. **$\mathbb Z_5$-connectivity of 3-edge-connected graphs.** JLPT's conjecture. A weaker but new target: $\mathbb Z_5$-connectivity of 3-edge-connected graphs with a specified structure (for example cubic graphs plus one contracted cycle), which is exactly what a Seymour-style proof would need. Study the $\mathbb Z_5$-connected-but-not-$\mathbb Z_6$-connected example to understand what fails.

5. **Kochol's boundary-set calculus with computer assistance.** Automate the $k$-pole boundary-set computations to (a) push the girth bound past 11, (b) more usefully, attempt a reduction to cyclically 7-edge-connected snarks, which combined with the Jaeger-Swart conjecture (cyclic connectivity $\le 6$ for snarks) would be decisive, and independently would tighten the target.

6. **Counterexample search.** Use superposition (Kochol) to build cyclically 6-edge-connected snarks with oddness $\ge6$ and girth $\ge 7$ (girth 11 is much harder), and test $\mathbb Z_5$-flows by SAT. Two outcomes are valuable: a counterexample ends the problem; a large verified family of the "dangerous" type is evidence and may reveal the structural reason they always have 5-flows.

7. **A genuinely new $\mathbb Z_5$ extension lemma.** The structural bottleneck of §7: find a class of "$\mathbb Z_5$-extendable" configurations larger than "circuit with $\le4$ new edges", e.g. two circuits sharing a path, or a theta with $\le 6$ new edges, so that every 3-connected cubic graph decomposes into them. This is the only item on the list that could give a short proof, and the only one with no existing partial results to build on.

---

## 10. Bibliography

**Books and surveys**

- C.-Q. Zhang, *Integer Flows and Cycle Covers of Graphs*, Marcel Dekker, 1997. (Chapters 1-2, 5, 8.7, 9; the reference for everything marked "Zhang".)
- C.-Q. Zhang, *Circuit Double Cover of Graphs*, LMS Lecture Notes 399, Cambridge, 2012.
- F. Jaeger, Nowhere-zero flow problems, in *Selected Topics in Graph Theory 3*, Academic Press, 1988, 71-95.
- P. D. Seymour, Nowhere-zero flows, in *Handbook of Combinatorics*, 1995.
- D. H. Younger, Integer flows, *J. Graph Theory* 7 (1983) 349-357.
- M. A. Fiol, G. Mazzuoccolo, E. Steffen, On measures of edge-uncolorability of cubic graphs: a brief survey and some new results, arXiv:1702.07156 (2017).
- F. Dong, A survey on the study of real zeros of flow polynomials, arXiv:2007.05195 (2020).
- Open Problem Garden, "5-flow conjecture", https://www.openproblemgarden.org/op/5_flow_conjecture.

**Origins and general bounds**

- W. T. Tutte, A contribution to the theory of chromatic polynomials, *Canad. J. Math.* 6 (1954) 80-91.
- W. T. Tutte, On the algebraic theory of graph colorings, *J. Combin. Theory* 1 (1966) 15-50.
- P. A. Kilpatrick, Tutte's first colour-cycle conjecture, PhD thesis, Cape Town, 1975.
- F. Jaeger, Flows and generalized coloring theorems in graphs, *J. Combin. Theory Ser. B* 26 (1979) 205-216.
- P. D. Seymour, Nowhere-zero 6-flows, *J. Combin. Theory Ser. B* 30 (1981) 130-135.
- G. Fan, Integer flows and cycle covers, *J. Combin. Theory Ser. B* 54 (1992) 113-122.
- M. DeVos, E. Rollová, R. Šámal, A new proof of Seymour's 6-flow theorem, *J. Combin. Theory Ser. B* 122 (2017) 187-195; M. DeVos, Another proof of Seymour's 6-flow theorem, *J. Graph Theory* (2024); M. DeVos, J. Nurse, A short proof of Seymour's 6-flow theorem, *Electron. J. Combin.* 32(4) (2025).
- K. R. Matthews, On the eulericity of a graph, *J. Graph Theory* 2 (1978) 143-148.
- C. H. C. Little, W. T. Tutte, D. H. Younger, A theorem on integer flows, *Ars Combin.* 26A (1988) 109-112.

**Minimal counterexamples**

- U. A. Celmins, On cubic graphs that do not have an edge 3-coloring, PhD thesis, Waterloo, 1984.
- M. Möller, H. G. Carstens, G. Brinkmann, Nowhere-zero flows in low genus graphs, *J. Graph Theory* 12 (1988) 183-190.
- K. Sekine, C.-Q. Zhang, Decomposition of the flow polynomial, *Graphs Combin.* 13 (1997) 189-196.
- M. Kochol, Snarks without small cycles, *J. Combin. Theory Ser. B* 67 (1996) 34-47; A cyclically 6-edge-connected snark of order 118, *Discrete Math.* 161 (1996) 297-300.
- M. Kochol, Superposition and constructions of graphs without nowhere-zero k-flows, *European J. Combin.* 23 (2002) 281-306.
- M. Kochol, Polynomials associated with nowhere-zero flows, *J. Combin. Theory Ser. B* 84 (2002) 260-269.
- M. Kochol, Reduction of the 5-flow conjecture to cyclically 6-edge-connected snarks, *J. Combin. Theory Ser. B* 90 (2004) 139-145.
- M. Kochol, Restrictions on smallest counterexamples to the 5-flow conjecture, *Combinatorica* 26 (2006) 83-89.
- M. Kochol, Smallest counterexample to the 5-flow conjecture has girth at least eleven, *J. Combin. Theory Ser. B* 100 (2010) 381-389.
- M. Kochol, N. Krivoňáková, S. Smejová, K. Šranková, Matrix reduction in a combinatorial computation, *Inform. Process. Lett.* 111 (2011) 164-168.
- E. Steffen, Tutte's 5-flow conjecture for highly cyclically connected cubic graphs, *Discrete Math.* 310 (2010) 385-389.
- E. Steffen, Intersecting 1-factors and nowhere-zero 5-flows, *Combinatorica* 35 (2015) 633-640.
- G. Mazzuoccolo, E. Steffen, Nowhere-zero 5-flows on cubic graphs with oddness 4, *J. Graph Theory* 85 (2017) 363-371.
- G. Brinkmann, J. Goedgebeur, J. Hägglund, K. Markström, Generation and properties of snarks, *J. Combin. Theory Ser. B* 103 (2013) 468-488.
- F. Jaeger, T. Swart, Conjectures 1 and 2, in *Combinatorics 79*, Ann. Discrete Math. 9 (1980) 305.

**Graph classes**

- R. Steinberg, Tutte's 5-flow conjecture for the projective plane, *J. Graph Theory* 8 (1984) 277-285.
- M. Kochol, Cubic graphs without a Petersen minor have nowhere-zero 5-flows, *Acta Math. Univ. Comenian.* 68 (1999) 249-252.
- N. Robertson, P. Seymour, R. Thomas, Tutte's edge-colouring conjecture, *J. Combin. Theory Ser. B* 70 (1997) 166-183.
- K. Edwards, D. Sanders, P. Seymour, R. Thomas, Three-edge-colouring doublecross cubic graphs, *J. Combin. Theory Ser. B* 119 (2016) 66-95.
- Y. Inoue, K. Kawarabayashi, R. Matsuo, A. Miyashita, B. Mohar, T. Sonobe, Three-edge-coloring apex cubic graphs, arXiv:2608.22870 (Aug 2026).
- M. Matamala, J. Zamora, Nowhere-zero 5-flows and even (1,2)-factors, *Graphs Combin.* 29 (2013) 609-616.
- M. Kochol, Hypothetical complexity of the nowhere-zero 5-flow problem, *J. Graph Theory* 28 (1998) 1-11.

**Circular flows, group connectivity, orientations**

- L. A. Goddyn, M. Tarsi, C.-Q. Zhang, On (k,d)-colorings and fractional nowhere-zero flows, *J. Graph Theory* 28 (1998) 155-161.
- E. Steffen, Circular flow numbers of regular multigraphs, *J. Graph Theory* 36 (2001) 24-34.
- Z. Pan, X. Zhu, Construction of graphs with given circular flow numbers, *J. Graph Theory* 43 (2003) 304-318.
- E. Máčajová, A. Raspaud, On the strong circular 5-flow conjecture, *J. Graph Theory* 52 (2006) 307-316.
- R. Lukot'ka, M. Škoviera, Snarks with given real flow numbers, *J. Graph Theory* 68 (2011) 189-201.
- L. Esperet, G. Mazzuoccolo, M. Tarsi, The structure of graphs with circular flow number 5 or more, and the complexity of their recognition problem, *J. Comb.* 7 (2016) 453-479.
- J. Goedgebeur, D. Mattiolo, G. Mazzuoccolo, Computational results and new bounds for the circular flow number of snarks, *Discrete Math.* 343 (2020) 112026; A unified approach to construct snarks with circular flow number 5, *J. Graph Theory* 96 (2021) 313-334.
- F. Jaeger, N. Linial, C. Payan, M. Tarsi, Group connectivity of graphs: a nonhomogeneous analogue of nowhere-zero flow properties, *J. Combin. Theory Ser. B* 56 (1992) 165-182.
- F. Jaeger, On circular flows in graphs, in *Finite and Infinite Sets*, Colloq. Math. Soc. János Bolyai 37 (1984) 391-402.
- C. Thomassen, The weak 3-flow conjecture and the weak circular flow conjecture, *J. Combin. Theory Ser. B* 102 (2012) 521-529.
- L. M. Lovász, C. Thomassen, Y. Wu, C.-Q. Zhang, Nowhere-zero 3-flows and modulo k-orientations, *J. Combin. Theory Ser. B* 103 (2013) 587-598.
- M. Han, J. Li, Y. Wu, C.-Q. Zhang, Counterexamples to Jaeger's circular flow conjecture, *J. Combin. Theory Ser. B* 131 (2018) 1-11.
- M. Delcourt et al., Almost all 9-regular graphs have a modulo-5 orientation, *Electron. J. Combin.* 32(2) (2025).

**Cycle covers and related conjectures**

- D. Archdeacon, Face colorings of embedded graphs, *J. Graph Theory* 8 (1984) 387-398.
- U. Jamshy, A. Raspaud, M. Tarsi, Short circuit covers for regular matroids with a nowhere zero 5-flow, *J. Combin. Theory Ser. B* 43 (1987) 354-357.
- J. Jooken, A human-checkable proof of the 112-vertex counterexample to the Petersen coloring conjecture, arXiv:2608.10028 (Aug 2026); the counterexample is due to Putman (2026).

**Flow polynomial**

- J. L. Jacobsen, J. Salas, Is the five-flow conjecture almost false?, *J. Combin. Theory Ser. B* 103 (2013) 532-565.
- G. Haggard, D. J. Pearce, G. Royle, Computing Tutte polynomials, *ACM Trans. Math. Software* 37 (2010).
- B. Jackson, Zeros of chromatic and flow polynomials of graphs, *J. Geom.* 76 (2003) 95-109.

**2025-2026 variants**

- L. Esperet, K. Hendrey, A. Lagoutte, M. Marseloo, S. Norin, R. Steiner, Nowhere-zero flow reconfiguration, arXiv:2512.17342 (Dec 2025).
- D. W. Cranston, J. Li, B. Su, Z. Wang, N. Xu, Reconfiguration of nowhere-zero flows, arXiv:2606.24685 (June 2026).
- L. Gáborik, S. Kurz, G. Mazzuoccolo, J. Rajník, F. Rieg, Manhattan and Chebyshev flows, arXiv:2510.22234 (Oct 2025).

## Addendum (September 2026): the cycle double cover proof and why its method does not reach 5-flows

The cycle double cover conjecture was proved in July 2026 (OpenAI, "A proof of the cycle double cover conjecture"; exposition by S.-i. Oum, arXiv:2607.16356; case study at geneweng.github.io/cdc-case-study). The route: reduce to cubic graphs; take a nowhere-zero $\mathbb F_2^3$-flow $\varphi$ (the 8-flow theorem, or three edge-disjoint spanning trees after doubling); ask for vertex vectors $t_v\in\mathbb F_2^3$ with the coset conditions $t_u+t_v\in d_e+\langle\varphi(e)\rangle$ for every edge; show the linear system is solvable because every certificate of infeasibility $(h_e)$ satisfies $\sum_e h_e\cdot d_e=0$, the sum counting each edge twice, once per end, and $2=0$; read off two-element edge labels $P_e$ (cosets of $\langle\varphi(e)\rangle$) and the seven even subgraphs $H_s=\{e: s\in P_e\}$, each edge lying in exactly two of them.

Three features make this work, and all three fail for Tutte's conjecture.

1. **The target is a parity object.** "Every edge in exactly two cycles" is a statement in the $\mathbb F_2$-cycle space. A nowhere-zero $\mathbb Z_5$-flow is an element of the $\mathbb F_5$-cycle space avoiding all coordinate hyperplanes: an avoidance condition, not an affine one. No linear system over any field has exactly the nowhere-zero flows as its solutions, so there is nothing for a solvability certificate to certify.
2. **Two-element sets are cosets only over $\mathbb F_2$.** The whole construction rests on $x+\langle p\rangle=\{x,x+p\}$. The analogous local constraint for 5-flows, in any of its forms (flow value in $\mathbb F_5^*$; orientation value $\pm1$ in the mod-5 orientation formulation of $3G$; a black/white choice with cut inequalities in the balanced-valuation formulation), is a subset of $\mathbb Z_5$ or of $\mathbb Z$ that is not a coset of a subgroup, since $\mathbb Z_5$ has no proper subgroups. The double-counting finish ("each edge has two ends, twice anything is zero") has no analogue in characteristic 5.
3. **Characteristic 2 forgets orientations, 5-flows are made of them.** The page notes that $-x=x$ makes edge directions irrelevant; that is exactly what a $\mathbb Z_5$-flow cannot afford. The bridge from covers to flows is the *orientable* cycle double cover: an orientable $k$-cycle double cover gives a nowhere-zero $\mathbb Z_k$-flow, so the orientable 5-cover conjecture would imply the 5-flow conjecture. The new proof produces unoriented covers with seven or eight layers and takes an 8-flow as input; it returns nothing to flows that was not already known (an $\mathbb F_2^3$-flow already yields a 7-cycle 4-cover, Jaeger). Even the unoriented 5-cycle double cover conjecture does not imply the 5-flow conjecture.

What can be borrowed is the pattern, not the mechanism: trade the existence question for a rigid algebraic object whose infeasibility certificates can be analysed. The two routes in these notes are already of that kind (mod-5 orientations of $3G$: a $\mathbb Z_5$-linear condition with a non-affine domain $\{\pm1\}$; balanced valuations: an integer system of cut inequalities, whose "certificates" are the bad cuts catalogued in the oddness notes). In both, the certificates do not cancel identically; the oddness-6 experiments show them failing by a single edge in specific colourings, which is the opposite of the CDC situation where every certificate vanishes for free. The nearest live analogue of the CDC mechanism in flow theory is group connectivity (Jaeger–Linial–Payan–Tarsi), where one asks for flows with prescribed boundary over $\mathbb Z_k$ and the strongest tools are the polynomial method and additive bases rather than linear solvability; the known partial results toward 5-flows (Seymour's 6-flow, the 3-flow results of Lovász–Thomassen–Wu–Zhang) come from there and from tree packing, which naturally produces $\mathbb F_2^k$ or $\mathbb Z_2\times\mathbb Z_3$ flows but not a single prime field $\mathbb F_5$.
