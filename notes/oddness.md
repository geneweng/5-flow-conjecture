# Route: oddness 6

*Working notes, started 2026-09-12. Companion to the survey (§3, §9 route 1). Code: `tools/oddness.py`.*

## 0. Target

> **(O6)** Every cyclically 6-edge-connected cubic graph with oddness at most 6 has a nowhere-zero 5-flow.

Known: oddness $\le 2$ (Jaeger 1988), oddness $\le 4$ with cyclic 6-edge-connectivity (Mazzuoccolo–Steffen 2017, "MS"), and cyclic connectivity $\ge \tfrac52\omega-3$ (Steffen 2010). Since a smallest counterexample to the 5-flow conjecture is cyclically 6-edge-connected (Kochol), (O6) would push its oddness to $\ge 8$. It is a natural next step but not obviously within reach of the MS method, as explained in §2.

## 1. The method (Steffen 2010, MS 2017), restated

**Balanced valuations (Bondy, Jaeger).** $w\colon V\to\mathbb R$ is balanced if $|\sum_{v\in X}w(v)|\le|\partial X|$ for all $X$. Jaeger: $G$ has a nowhere-zero $k$-flow iff it has a balanced valuation of the form $w(v)=\tfrac{k}{k-2}(2d^+(v)-d(v))$ for some orientation. For cubic $G$ and $k=5$: **$G$ has a nowhere-zero 5-flow iff some 2-colouring $V=A\cup B$ (white/black) makes $w=\mp\tfrac53$ balanced**, i.e. $5\,|b_X-a_X|\le 3\,|\partial X|$ for all $X$, where $a_X=|X\cap A|$, $b_X=|X\cap B|$. Note this forces $|A|=|B|$ and, for $|\partial X|\le 5$, $|b_X-a_X|\le 3$.

**Flow partitions.** Fix a 2-factor $F_2$ with odd circuits $C_1,\dots,C_{2t}$ and even circuits, and the complementary perfect matching $M$. A *canonical colouring* $c$ gives colour 1 to $M$, colours 2,3 alternately along each circuit, except that on each odd circuit one chosen edge (the *0-edge*) gets colour 0. Exactly one vertex per odd circuit misses colour 2; call them $z_1,\dots,z_{2t}$. The subgraph $H=c^{-1}(1)\cup c^{-1}(2)$ has all degrees $\le2$, degree 1 exactly at the $z_i$; so $H$ is a disjoint union of even circuits and $t$ paths $P_1,\dots,P_t$ of odd length pairing up the $z_i$ (the pairing is dictated by $G$, $F_2$, $c$).

MS build a nowhere-zero 4-flow on $G$ plus two parallel edges per path and read off a balanced $\pm2$ valuation. Unwinding their construction (I checked this from the definitions): the resulting black/white partition is exactly a **proper 2-colouring of $H$**, each component of $H$ coloured independently (their Lemma 2.4 is the statement that recolouring one path is allowed; recolouring an even circuit is equally allowed since the circuit orientations in the 4-flow are arbitrary). The vertex classes of the 2-colouring of $H$ are automatically equal in size on each component (even circuits, and odd-length paths have an even number of vertices), so $|A|=|B|$.

So the **method** is: choose $F_2$, choose the 0-edges, choose a 2-colouring of each component of $H$, set $w=\pm\tfrac53$, and hope $w$ is balanced. The whole freedom is
$$\underbrace{\text{2-factors with }2t\text{ odd circuits}}_{\text{huge}}\times\underbrace{\prod_i|C_i|}_{\text{0-edges}}\times\underbrace{2^{\#\text{components of }H}}_{\text{colourings}} .$$
MS use only the $2^t$ path recolourings for a fixed $F_2$ and $c$ (and only $t=2$).

**What can go wrong.** Suppose $w$ is not balanced, witnessed by $S$ with $5k>3m$, $k=|b_S-a_S|$, $m=|\partial S|$. Let $c_i=|\partial S\cap c^{-1}(i)|$. Since colour-1 and colour-2 edges are bichromatic, $k\le c_1$ and $k\le c_2+q$, where $q$ is the colour imbalance among the $z$'s inside $S$; also $k\equiv m\pmod 2$. Adding, $2k\le m+q$, hence $q>m/5$, and since $q\le t$ we get **$m<5t$**: bad cuts have fewer than $5t$ edges. Steffen's Proposition 1 refines this by residues mod 5, and with cyclic 6-edge-connectivity MS show for $t=2$ that a bad cut has exactly 6 edges, 4 of colour 1 and 2 of colour 2, and separates the $z$'s into two same-coloured pairs.

## 2. What changes at oddness 6 ($t=3$)

- **Bad cut sizes.** $m<15$; Steffen's Proposition 1 with $\omega=6$ gives $m\le11$ if $m\equiv1$, $m\le7$ if $m\equiv2$, $m\le8$ if $m\equiv3$ (mod 5) and nothing for $m\equiv0,4$ beyond $m\le4$. Hence, after excluding the non-cyclic cuts of size $\le5$ by cyclic 6-edge-connectivity: **$m\in\{6,7,8,11\}$** and $q\in\{2,3\}$.
- **The colour imbalance $q$** of $S\cap Z$: paths of $H$ with both ends in $S$ contribute 0; a path $P_i$ with exactly one end in $S$ contributes $\varepsilon_i\sigma_i=\pm1$ where $\varepsilon_i$ is the chosen colouring of $P_i$ and $\sigma_i$ records which end is in $S$. So $q\ge2$ needs at least two paths separated by $\partial S$ with equal signs; $q=3$ needs all three separated with equal signs.
- **Constraint system.** Write $J(S)\subseteq\{1,2,3\}$ for the paths separated by $\partial S$. A cut with $J=\{a,b\}$ is bad for the colourings with $\varepsilon_a\varepsilon_b=\sigma_a\sigma_b=:s(S)$; a cut with $J=\{1,2,3\}$ is bad for the colourings with $\varepsilon_1\sigma_1=\varepsilon_2\sigma_2=\varepsilon_3\sigma_3$. Up to global complement there are 4 colourings $\varepsilon$ of the three paths. The method (with $F_2,c$ fixed and only path recolourings) succeeds iff the system of "forbidden" relations has a solution. It has **no** solution in exactly these situations:
  1. two bad cuts with the same pair $J=\{a,b\}$ and opposite signs $s$ (this is the only obstruction at $t=2$, and MS's whole proof is that it cannot happen);
  2. three bad cuts with $J=\{1,2\},\{1,3\},\{2,3\}$ and $s_{12}s_{13}s_{23}=+1$ (each forbids one parity; the three forced parities are inconsistent);
  3. combinations with $J=\{1,2,3\}$ cuts (each forbids two of the eight $\varepsilon$, i.e. one class up to complement) together with pair cuts.
- **Why the MS argument does not transfer.** Their final contradiction takes two bad cuts $E_3,E_4$ for $J=\{1,2\}$ with opposite signs, uncrosses them into $U_1,\dots,U_4$ using cyclic 6-edge-connectivity to force $|\partial(U_i,U_j)|\in\{0,3\}$, and then observes that $E_3'\cup E_4'$ is a 6-cut inside $H$ crossed an odd number of times by $P_1\cup P_2$ while every circuit of $H$ crosses it evenly. With a third path $P_3$ in $H$, $P_3$ may cross $E_3'\cup E_4'$ an odd number of times, and the parity contradiction disappears. So even obstruction 1 needs a new argument at $t=3$, and obstructions 2 and 3 are new.
- **Extra freedom not used by MS.** The 0-edge on each odd circuit (which moves the $z_i$ around the circuit and changes the paths), the colourings of the even circuits of $H$, and the choice of $F_2$ itself. Any of these may be needed.

## 3. Decisive experiment before theory

Before trying to prove (O6) with the flow-partition method it should be checked that the method is *capable* of it: for cyclically 6-edge-connected cubic graphs with a 2-factor of exactly 6 odd circuits, is there always a balanced flow partition, and if so, is it always reachable by path recolourings alone (as in MS), or does one need the 0-edge/even-circuit/2-factor freedom? `tools/oddness.py` implements the method exactly (canonical colouring, $H$, 2-colourings, the balance test as one min-cut) and, as a check, reproduces Jaeger's oddness-2 theorem: on Petersen, $J_7$ and $J_9$ every sampled 2-factor with 2 odd circuits and every 0-edge choice gives a balanced partition (150 + 22,354 + 71,657 cases).

Test graphs: generalized Petersen graphs $GP(n,3)$ and $GP(n,4)$ with $n$ odd (girth 6 or 7, so cyclic connectivity is plausibly 6 or more; they are 3-edge-colourable but have 2-factors with 4 or 6 odd circuits, which is all the fixed-$F_2$ template sees), the flower snarks $J_7,J_9$ (cyclically 6-edge-connected, oddness 2), and the oddness-6 snark $R_2$ on 40 vertices of Lukoťka–Máčajová–Mazák–Škoviera (cyclic connectivity only 2, so a failure there would not contradict (O6), but a success is informative). Results are recorded in §4 as they come in.

## 4. Experiment log

(pending)

## 5. Plan

1. Run the experiment of §3. Three possible outcomes: (a) path recolourings always suffice at $t=3$ on cyclically 6-connected graphs: then (O6) is a case analysis of obstructions 1–3 above, and the job is to find the replacement for the parity argument; (b) they do not suffice but 0-edge/even-circuit freedom does: then the proof must first choose $c$ well, and the experiment tells which choices work; (c) some cyclically 6-connected graph with a 6-odd 2-factor has no balanced flow partition at all: then (O6), if true, needs a different valuation family, and the template is dead at oddness 6.
2. If (a): prove obstruction 2 impossible via uncrossing three cuts of size 6–11; prove obstruction 1 with $P_3$ present by controlling the parity of $|E(P_3)\cap(E_3'\cup E_4')|$, possibly by choosing $\varepsilon_3$ or the 0-edge of $C_5,C_6$.
3. Either way, record the smallest configurations found; they are the "bad 6-cut" analogues that a proof must enumerate.
