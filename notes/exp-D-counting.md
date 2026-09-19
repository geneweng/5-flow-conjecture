# Experiment D — count instead of decide: $F(G;5)$ on cubic graphs

*2026-09-19. Implements idea D of `oddness.md` §18.4. Scripts: `tools/flowcount.py` (exact counting), `tools/flowdata.py` (families, random graphs, exhaustive enumeration with `geng`), `tools/flowclimb.py` (annealing climber minimising $F$), `tools/flowcycles.py`, `tools/flowstruct.py`. Logs: `tools/flowdata_*.log`, `tools/flowclimb_*.log`, `tools/flowcycles.log`, `tools/flowstruct.log`; best graphs of the climber in `tools/flowclimb_best_*.json`; enumeration summaries in `tools/flowdata_enum_*.json`.*

$F(G;q)$ = number of nowhere-zero $\mathbb Z_q$-flows; $n=|V|$, $m=3n/2$. Everything below is exact integer arithmetic.

## 0. Summary

1. **Tool.** $F(G;q)$ is computed as a tensor-network contraction (one $4\times4\times4$ Kirchhoff tensor per vertex, contraction tree found by randomised greedy, width = max cut $d(S)$ over the tree). About 1 ms per graph for $n\le30$, 0.02–1 s for the project's cyclically 6-connected graphs on 42–78 vertices, 20 s for a random cubic graph on 80 vertices. Validated against deletion–contraction, brute force, the polynomials of $K_4$, $K_{3,3}$, and $F(P_{10};5)=240$, $F(P_{10};4)=0$.
2. **The right normalisation is $1.5^n$, not $4^m/5^{n-1}$.** For every graph of girth $\ge5$ that was computed (snarks, cages, generalised Petersen graphs, random cubic graphs, all 26 of the project's cyclically 6-connected graphs), $F(G;5)/1.5^n\in[2.6,\,7.3]$; for the project's graphs it is $5.0\pm0.2$. (Local heuristic: three independent uniform non-zero values sum to zero with probability $12/64$, so $F\approx4^{3n/2}(3/16)^n=1.5^n$; the naive $4^m/5^{n-1}=5\cdot1.6^n$ overcounts exponentially.) Short cycles correct this multiplicatively: **odd cycles depress the count, even cycles raise it** (triangle $\times0.889$ exactly, 4-cycle $\approx\times1.02$–$1.04$, 5-cycle $\approx\times0.97$–$0.99$, 6-cycle $\approx\times1.005$; regression in §3).
3. **3-connected cubic graphs: $\min F(G;5)=6\cdot2^{n/2}$ exactly, for every even $n\le20$** (exhaustive, 556,000 graphs), attained precisely by the graphs obtained from $K_4$ by inflating vertices to triangles. Reason: $F$ is multiplicative over 3-edge-cuts, $F(G_1\oplus_3G_2)=F(G_1)F(G_2)/12$, and a triangle is a 3-sum with $K_4$ ($F=24$, factor 2 for 2 vertices). So the growth base over 3-connected graphs is $\sqrt2=1.414$, *also for girth 5* (3-sums at triangle vertices destroy the triangles: the girth-5 minimum is $4800\cdot2^{(n-18)/2}$ for $n=18,20,22$, from $P\oplus_3P$). The interesting classes are therefore the cyclically 4- and 5-connected ones.
4. **Cyclically 4-connected / cyclically 5-connected with girth $\ge5$.** Exact minima for $n\le22$ (resp. upper bounds from the climber for $24\le n\le30$) grow by a factor $2.18$–$2.27$ per two vertices (geometric mean over $n=20\to30$: $2.227$ for cyc 4, $2.215$ for cyc 5, i.e. base $1.492$ and $1.488$), slightly **below** the generic $1.5$ ($2.25$ per two vertices) and indistinguishable, at this range, from the base $24^{1/8}=1.4877$ of the Blanuša chain $P\cdot P\cdots P$ (exact: $F$ is multiplied by $24.000$ per dot product with $P_{10}$). The minimisers are **not snarks** in general, they are graphs with an excess of 5-cycles (8–11 against 4–5 for random girth-5 graphs) and a bounded number of Petersen fragments.
5. **Petersen is extremal only among the smallest graphs.** $F^{1/n}$ is decreasing in $n$ along every family ($P_{10}$: 1.730, $n=30$ minimiser: 1.570, $J_{25}$: 1.530), so "is $P_{10}$ the minimiser of $F^{1/n}$" is the wrong question; the right one is the base of $\min_nF$, and there $P_{10}$ enters through the Blanuša chain (cyclic connectivity 4). $P_{10}$ is the unique minimiser of $F$ at $n=10$ among triangle-free cubic graphs (240; next 288).
6. **Arithmetic.** In all data (exhaustive $n\le20$, girth $\ge5$ for $n=22$, all families, all random graphs, $n\le100$): $12\mid F(G;5)$, and **$F(G;5)\bmod5\in\{0,\pm1\}$ if $4\mid n$, $\in\{0,\pm2\}$ if $n\equiv2\pmod4$.** Each of the three residues occurs for about a third of the graphs. Not every $F$ is divisible by 24 ($K_{3,3}$: 60). Every snark in the data has $5\mid F(G;5)$ (§6).
7. **Conjectures supported by the data** (§7): (A) $F(G;5)\ge6\cdot2^{n/2}$ for 3-connected cubic $G$, equality iff $G$ is a triangle-inflation of $K_4$ — equivalent to the same inequality for cyclically 4-connected graphs, where it holds with exponential room; strong evidence (exhaustive to $n=20$). (B) $F(G;5)\ge c\cdot\beta^n$ for cyclically 4-connected cubic $G$ with some $\beta\in[1.47,1.495]$; evidence is moderate ($n\le22$ exact), the value of $\beta$ and the extremal family are **not** determined by the data.

## 1. The counting tool (`tools/flowcount.py`)

Fix an orientation. $F(G;q)=\sum_{f:E\to\mathbb Z_q^*}\prod_v[\partial f(v)=0]$ is the contraction of a tensor network with one tensor $K_v[a_1,a_2,a_3]=[\pm a_1\pm a_2\pm a_3\equiv0]$ of shape $(q-1)^3$ per vertex. Clusters of vertices are merged pairwise along a contraction tree; the tensor of a cluster $S$ has one axis per edge of $d(S)$, so the cost is $\sum(q-1)^{d(S\cup S')+\#\text{shared}}$ and the governing parameter is $\max d(S)$ over the tree (a carving-width type parameter; a vertex order with small frontier, as suggested in the task, is the caterpillar special case and was markedly worse: width 14 against 11 on a random cubic graph with 60 vertices). Plans come from a randomised greedy (score $d(S\cup S')-\alpha\log_4(4^{d(S)}+4^{d(S')})$ plus noise, 40–200 restarts). A plan remains valid when edges are switched on the same vertex set, which the climber uses.

Arithmetic: float64 with BLAS, exact while all entries are $<2^{52}$ (entries are non-negative integers, so partial sums are bounded by the result; checked after every contraction); otherwise int64 modulo primes just below $2^{20}$ and CRT, with the rigorous bound $F\le(q-1)^{m-n+1}$. The CRT path is exercised by $J_{25}$ ($F\approx2.9\cdot10^{18}$) and tested against the float path.

Validation (`python3 tools/flowcount.py --test`): $F(P_{10};3,4,5,6)=0,0,240,1920$; $F(K_4;q)=(q-1)(q-2)(q-3)$ for $q\le7$; $F(K_{3,3};q)=(q-1)(q-2)(q^2-6q+10)$; prism $0,6,48,180$; cube; Heawood; $J_3$; 40 random cubic graphs ($n\le10$, $q\in\{3,4,5,6\}$) against deletion–contraction and under random contraction orders; 30 random multigraphs with loops and parallel edges against deletion–contraction *and* brute force over the cotree; CRT path against float path. Also verified numerically: $F(G_1\oplus_3G_2)=F(G_1)F(G_2)/12$ on four 3-sums, and the climber reproduces the exhaustive minima (§5).

Widths met in practice: $\le8$ for $n\le30$; 8–11 for the project's graphs ($n=42$–$78$); 13 for $GP(30,6)$, $GP(35,7)$; $GP(36,6)$, $GP(42,6)$, $GP(42,7)$ need width 15–16 and were skipped (limit 14 = 2 GB).

## 2. Families

$r=F/1.5^n$; "naive" $=F\big/(4^m/5^{n-1})$; cyc = cyclic edge-connectivity (SAT, exact, capped at 7).

| graph | $n$ | $F(G;5)$ | $F(G;4)$ | $F^{1/n}$ | $r$ | naive | girth | cyc |
|---|---|---|---|---|---|---|---|---|
| $K_4$ | 4 | 24 | 6 | 2.213 | 4.74 | 0.73 | 3 | – |
| $K_{3,3}$ | 6 | 60 | 12 | 1.979 | 5.27 | 0.72 | 4 | – |
| prism | 6 | 48 | 6 | 1.906 | 4.21 | 0.57 | 3 | 3 |
| cube | 8 | 156 | 24 | 1.880 | 6.09 | 0.73 | 4 | 4 |
| **Petersen** | 10 | **240** | 0 | 1.730 | 4.16 | 0.44 | 5 | 5 |
| Heawood | 14 | 1692 | 48 | 1.701 | 5.80 | 0.47 | 6 | 6 |
| dodecahedron | 20 | 16080 | 60 | 1.623 | 4.84 | 0.27 | 5 | 5 |
| Desargues $=GP(10,3)$ | 20 | 18396 | 192 | 1.634 | 5.53 | 0.30 | 6 | 6 |
| McGee | 24 | 83616 | 204 | 1.604 | 4.97 | 0.21 | 7 | 7 |
| Coxeter | 28 | 406560 | 336 | 1.586 | 4.77 | 0.16 | 7 | 7 |
| Tutte–Coxeter | 30 | 1038348 | 864 | 1.587 | 5.42 | 0.16 | 8 | $\ge7$ |
| flower $J_5$ | 20 | 16200 | 0 | 1.624 | 4.87 | 0.27 | 5 | 5 |
| $J_7$ | 28 | 455280 | 0 | 1.593 | 5.34 | 0.18 | 6 | 6 |
| $J_9$ | 36 | 12247320 | 0 | 1.574 | 5.61 | 0.11 | 6 | 6 |
| $J_{11}$ | 44 | 325290240 | 0 | 1.561 | 5.81 | 0.068 | 6 | 6 |
| $J_{13}$ | 52 | 8607117480 | 0 | 1.553 | 6.00 | 0.042 | 6 | 6 |
| $J_{15}$ | 60 | 227528524560 | 0 | 1.546 | 6.19 | 0.026 | 6 | 6 |
| $J_{25}$ | 100 | 2947091530838746200 | 0 | 1.530 | 7.25 | 0.0023 | 6 | 6 |
| Blanuša chain $P^{\cdot2}$ | 18 | 6240 | 0 | 1.625 | 4.22 | 0.26 | 5 | 4 |
| $P^{\cdot3}$ | 26 | 158400 | 0 | 1.585 | 4.18 | 0.16 | 5 | 4 |
| $P^{\cdot5}$ | 42 | 91499520 | 0 | 1.547 | 3.68 | 0.049 | 5 | 4 |
| $P^{\cdot7}$ | 58 | 52818984960 | 0 | 1.531 | 3.23 | 0.015 | 5 | 4 |
| $P^{\cdot10}$ | 82 | 730213908480000 | 0 | 1.518 | 2.65 | – | 5 | 4 |
| $R_2$ (oddness 6, $n=40$) | 40 | 51840000 | 0 | 1.559 | 4.69 | 0.071 | 5 | 2 |
| $GP(7,2)$ | 14 | 1512 | 42 | 1.687 | 5.18 | 0.42 | 5 | 5 |
| $GP(9,2)$ | 18 | 6672 | 6 | 1.631 | 4.52 | 0.28 | 5 | 5 |
| $GP(13,5)$ | 26 | 187200 | 156 | 1.595 | 4.94 | 0.19 | 7 | 7 |
| $GP(16,6)$ | 32 | 2127816 | 480 | 1.577 | 4.93 | 0.13 | 7 | 7 |
| $GP(18,6)$ (6 triangles) | 36 | 6084864 | 504 | 1.543 | 2.79 | 0.055 | 3 | 3 |
| $GP(24,6)$ | 48 | 1929941064 | 44208 | 1.561 | 6.81 | 0.062 | 4 | 4 |
| $GP(28,7)$ | 56 | 52966976316 | 303936 | 1.554 | 7.29 | 0.039 | 4 | 4 |
| $GP(30,6)$ | 60 | 173832272400 | 33780 | 1.539 | 4.73 | 0.020 | 5 | 5 |
| $GP(35,7)$ | 70 | 9958206060000 | 173880 | 1.534 | 4.70 | 0.010 | 5 | 5 |
| prism $C_{30}\times K_2$ | 60 | 205895427061956 | $2^{30}+8$ | 1.732 | 5600 | 23 | 4 | 4 |
| Möbius ladder $M_{60}$ | 60 | 205895427061944 | $2^{30}+2$ | 1.732 | 5600 | 23 | 4 | 4 |
| `data/…_n42` | 42 | 124221120 | 1620 | 1.559 | 4.99 | 0.066 | 6 | 6 |
| `data/…_n50` | 50 | 3181640280 | 5058 | 1.549 | 4.99 | 0.040 | 6 | 6 |
| `data/robust_stage2b_n50` | 50 | 3154628208 | 4068 | 1.549 | 4.95 | 0.039 | 6 | 6 |
| `data/…_n54` | 54 | 16946665752 | 14040 | 1.547 | 5.25 | 0.032 | 6 | 6 |
| `data/…counterexample` (60) | 60 | 187291945224 | 22860 | 1.541 | 5.09 | 0.021 | 6 | 6 |
| `data/…_n60b` | 60 | 187385666016 | 22482 | 1.541 | 5.10 | 0.021 | 6 | 6 |
| 20 × `tools/fullkill_*.json` | 42–78 | — | 1092–452964 | 1.533–1.560 | **4.91–5.48** | — | 6 | (6) |

(Full list with all 20 archive graphs: `tools/flowdata_families.log`. The double-star and Szekeres snarks were not constructed; a from-memory edge list for the double star failed the $F(G;4)=0$ check and was dropped. All snarks of girth $\ge5$ and order $\le22$ are covered by the census of §6: the cyclically 4-connected ones have $F/1.5^n$ between 3.98 and 4.9, e.g. order 20: 14400 (four graphs), 14880, 16200 $=J_5$.)

**Asymptotic growth per vertex** (exact ratios of consecutive members, converged to 4 digits):

| family | limit of $F(G_{k+1})/F(G_k)$ | vertices added | base |
|---|---|---|---|
| triangle inflation | 2 (exact) | 2 | 1.4142 |
| 3-sum with $P_{10}$ | 20 (exact) | 8 | 1.4542 |
| 3-sum with $K_{3,3}$ | 5 (exact) | 4 | 1.4953 |
| **Blanuša chain (dot product with $P_{10}$)** | 24.000 | 8 | **1.4877** |
| flower snarks $J_k$ | 26.47 | 8 | 1.5061 |
| $GP(n,2)$ | 26.92 | 8 | 1.5092 |
| $GP(n,3)$ | $\approx26.9$ | 8 | $\approx1.509$ |
| prisms, Möbius ladders | 3 | 2 | 1.7321 |
| generic (random, high girth) | — | — | $1.50$ |

Observations. (i) The Jacobsen–Salas families $GP(6k,6)$, $GP(7k,7)$, whose flow roots accumulate at 5, are in no way small at $q=5$: $r=6.8$ and $7.3$ for $GP(24,6)$, $GP(28,7)$ (they have 4-cycles, which *raise* the count); $GP(30,6)$, $GP(35,7)$ have $r=4.7$. "Almost false" in the sense of real roots is invisible in the count at 5. (ii) Snarks are not special: $J_k$ has a *larger* $r$ than random graphs. (iii) The project's cyclically 6-connected oddness-6 template failures all sit at $r\approx5$: the count does not see the failure of the template. (iv) The climber for the template (§8, §12 of `oddness.md`) and this count are unrelated objectives, as far as this data goes.

## 3. Random cubic graphs and the effect of short cycles

`tools/flowdata_random.log`: 3-connected random cubic graphs (`networkx.random_regular_graph`).

| $n$ | #graphs | $F^{1/n}$ min / mean / max | $F/1.5^n$ min / mean / max |
|---|---|---|---|
| 20 | 30 | 1.582 / 1.616 / 1.639 | 2.89 / 4.51 / 5.89 |
| 24 | 30 | 1.524 / 1.589 / 1.625 | 1.46 / 4.09 / 6.78 |
| 30 | 30 | 1.540 / 1.574 / 1.591 | 2.19 / 4.30 / 5.81 |
| 36 | 30 | 1.542 / 1.561 / 1.573 | 2.70 / 4.23 / 5.52 |
| 42 | 30 | 1.538 / 1.553 / 1.564 | 2.86 / 4.33 / 5.79 |
| 50 | 30 | 1.534 / 1.544 / 1.553 | 3.02 / 4.22 / 5.59 |
| 56 | 12 | 1.532 / 1.538 / 1.544 | 3.24 / 4.03 / 4.99 |
| 60 | 12 | 1.531 / 1.538 / 1.546 | 3.46 / 4.56 / 6.05 |

The mean of $F/1.5^n$ does not move between $n=20$ and $n=60$: the growth rate of a random cubic graph is $1.500$ to three digits, and the minimum in every sample is a graph with 3–5 triangles.

Regression (`tools/flowcycles.log`, 150 random 3-connected cubic graphs, $n=36$): $\log(F/1.5^n)=\log5.64+\sum_{k=3}^7b_k\,\#\{k\text{-cycles}\}$, $R^2=0.81$, with per-cycle factors $e^{b_k}$:

| $k$ | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|
| fitted factor | 0.861 | 1.020 | 0.972 | 1.007 | 0.996 |
| ansatz $1+3(-\tfrac13)^k$ | 0.889 | 1.037 | 0.988 | 1.004 | 0.999 |

The ansatz is a guess that matches the exact triangle factor ($F(G)=2F(G/T)$ against $1.5^2$) and the sign pattern; it is not derived. The sign pattern is robust: **the count is depressed by odd cycles and raised by even ones**, and the effect decays like $3^{-k}$. Consistent with this, $P_{10}$ (twelve 5-cycles) has $r=4.16$ against a prediction $5.64\cdot0.972^{12}\cdot1.007^{10}=4.3$.

## 4. Exhaustive enumeration ($n\le22$, `geng`)

All connected cubic graphs for $n\le20$ (`-c`), triangle-free for $n\le20$ (`-ct`), girth $\ge5$ for $n\le22$ (`-ctf`); minima per class (cyclic connectivity by SAT on the lowest graphs until five of each class were seen).

| $n$ | #cubic graphs | min, 3-connected | $=6\cdot2^{n/2}$? | min, 3-conn. triangle-free | min, cyc. 4-conn. | min, cyc. 5-conn. & girth $\ge5$ |
|---|---|---|---|---|---|---|
| 4 | 1 | 24 | yes | – | – | – |
| 6 | 2 | 48 | yes | 60 | – | – |
| 8 | 5 | 96 | yes | 144 | 144 ($M_8$) | – |
| 10 | 19 | 192 | yes | 240 | **240** ($P_{10}$) | **240** |
| 12 | 85 | 384 | yes | 576 | 576 | 576 |
| 14 | 509 | 768 | yes | 1200 | 1272 | 1368 |
| 16 | 4060 | 1536 | yes | 2400 | 2784 | 2904 |
| 18 | 41301 | 3072 | yes | 4800 ($P\oplus_3P$) | 6240 (Blanuša and two non-snarks) | 6432 |
| 20 | 510489 | 6144 | yes | 9600 | 13680 | 14256 |
| 22 | (girth $\ge5$: 90938) | – | – | 19200 (girth $\ge5$) | $\le$ 29760 (girth $\ge5$; snark) | 31488 |

- For $10\le n\le16$ the stored lists contain all graphs attaining the minimum (3, 7, 24, 93 of them), for $n=18,20$ the first 100; **every one of them reduces to $K_4$ by contracting triangles**, all have $F(G;4)=6$ (uniquely 3-edge-colourable), and the second smallest value is $1.25\times$ the minimum $=240\cdot2^{(n-10)/2}$ (triangle inflations of $P_{10}$).
- Ratios of consecutive minima in the last column: 2.40, 2.375, 2.12, 2.21, 2.216, 2.209; cyclically 4-connected column: 1.67, 2.40, 2.21, 2.19, 2.24, 2.19, 2.18.
- The minimisers in the last two columns come with many ties (five or more graphs with the same $F$ at $n=20,22$), and are mostly **not** snarks ($F(G;4)=24$–$48$). At $n=18$ the two Blanuša snarks tie with two 3-edge-colourable graphs at 6240; at $n=20$ the flower snark $J_5$ (16200) and the dodecahedron (16080) are 14% above the minimum; at $n=22$ the cyclically 4-connected minimum 29760 is a snark (one of the Loupekine/Blanuša-type snarks of order 22, containing two copies of the dot-product block $P-uv$).

## 5. Climber (`tools/flowclimb.py`)

Simulated annealing on $\log F$ over 2-switches, temperature $0.25\to0.01$, 3 minutes per order (4 for the second seed), class constraints checked by SAT after the Metropolis test. The machine was heavily loaded by other jobs (load 12–14 on 12 cores), so runs have between 4,000 and 120,000 steps.

Check against §4: the climber finds the exact minimum for $n=12$–$18$ and $22$ in class (cyc 5, girth 5) and for $n=12$–$20$ in class (cyc 4); it misses by 0.3% at $n=20$ (14304 vs 14256) and by 2.4% at $n=22$ in class cyc 4 (30480 vs 29760). For $n\ge24$ its values are **upper bounds** of unknown quality (probably within a few percent).

| $n$ | cyc $\ge4$: best $F$ | $F/1.5^n$ | ratio | cyc $\ge5$, girth $\ge5$: best $F$ | $F/1.5^n$ | ratio |
|---|---|---|---|---|---|---|
| 20 | 13680* | 4.11 | | 14256* | 4.29 | |
| 22 | 29760* | 3.98 | 2.175 | 31488* | 4.21 | 2.209 |
| 24 | 67680 | 4.02 | 2.274 | 69960 (both seeds) | 4.16 | 2.222 |
| 26 | 151152 | 3.99 | 2.233 | 154272 (seed 2: 154512) | 4.07 | 2.205 |
| 28 | 333936 | 3.92 | 2.209 | 341616 (seed 2: 343896) | 4.01 | 2.214 |
| 30 | 749280 | 3.91 | 2.244 | 759552 (seed 1: 764088) | 3.96 | 2.223 |

(* exact, from §4.) For comparison at $n=26$: Blanuša $P^{\cdot3}$ has 158400, so the Blanuša chain is *not* the minimiser in the cyclically 4-connected class beyond $n=18$; at $n=42$: Blanuša $P^{\cdot5}$ 91499520 ($r=3.68$), random $\approx1.07\cdot10^8$.

Cyclically 6-connected, girth 6, $n=42$ (start: `data/oddness6_template_counterexample_n42.json`, $F=124221120$, 18 minutes, only 9,900 steps because of an over-eager replanning rule, since fixed): best $F=121237632$, $F/1.5^n=4.87$, a graph with no 5-cycles and two 6-cycles. So under girth 6 the climber moves the count by 2.4% only; the twenty archived 42-vertex graphs of the project lie between 122038368 and 128708352. In class (cyc 5, girth 5) at $n=42$ the climber is not informative: with cheap replanning (6 greedy restarts per move) the contraction width after a switch is 10–12 instead of 7–9, a step costs over a second, and a 6-minute run made 209 steps (best $F=111773160$, $F/1.5^n=4.49$, against $3.68$ for the Blanuša chain $P^{\cdot5}$ of the same order, which has cyclic connectivity 4). **The climber as written is adequate for $n\le30$ only**; for $n\approx42$ it needs incremental plan repair (keep the contraction tree, re-optimise only the subtrees touching the four switched vertices).

**Structure of the minimisers** (`tools/flowstruct.log`; induced copies of $P-v$ (3-pole), $P-uv$ (4-pole, the dot-product block), $P-P_3$ (5-pole, 7 vertices)):

| class | $n$ | #5-cycles | $P-v$ | $P-uv$ | $P-P_3$ | snark? |
|---|---|---|---|---|---|---|
| 3-connected, girth 5 (exact) | 18–22 | 12 | 2 | 6 | 18 | yes ($P\oplus_3P$ and descendants) |
| cyc 4 (exact) | 18 | 8–10 | 0 | 1–2 | 7–10 | 2 of 4 |
| cyc 4 (exact) | 20 | 8 | 0 | 0 | 8 | no |
| cyc 4 (exact, girth 5) | 22 | 8 | 0 | 2 | 8 | yes |
| cyc 4 (climber) | 24–30 | 9–11 | 0 | 1,1,1,0 | 7–8 | $n=24$ only |
| cyc 5 (exact) | 18 | 8 | 0 | 0 | 2 | no |
| cyc 5 (exact) | 20, 22 | 6–9 | 0 | 0 | **4** (all ten graphs) | no |
| cyc 5 (climber) | 24–30 | 7–10 | 0 | 0 | **4** (all four graphs) | no |
| random girth $\ge5$ | 22 / 30 | 4.8 / 4.0 | 0 | 0 / 0.05 | 1.1 / 0.7 | |

So: triangles and 3-sums with $P_{10}$ where allowed; the dot-product block $P-uv$ where cyclic connectivity 4 is allowed; and in the cyclically 5-connected class the minimisers always contain exactly four copies of the Petersen 5-pole $P-P_3$ and about twice the random number of 5-cycles, no 4-cycles even where they are allowed (they raise the count; the class-cyc-4 climber, which may use 4-cycles, ends on girth 5 in 7 of 10 runs and on a single 4-cycle in the other three). What the four copies of $P-P_3$ are (one cluster or two) and whether a family with *many* such clusters has a smaller base was not investigated.

## 6. Arithmetic of $F(G;5)$

- $4\mid F$ (units of $\mathbb Z_5$ act freely). In all data $12\mid F$; the gcd of all non-zero values is exactly 12 at every order $10\le n\le20$. $F/4\bmod6\in\{0,3\}$ always, with $3$ (i.e. $F\equiv12\bmod24$) for 2 of 19, 5 of 85, 13 of 509, 38 of 4060, 149 of 41301, 703 of 510489 graphs ($n=10,\dots,20$).
- **$F(G;5)\bmod5$:** $\{0,2,3\}$ for $n\equiv2\pmod4$ and $\{0,1,4\}$ for $n\equiv0\pmod4$, without exception in 556,000 + 91,000 enumerated graphs, 91 family graphs ($n\le100$) and 204 random graphs; the three residues are roughly equidistributed (e.g. $n=20$: 159670 / 167428 / 170720). Since $F(G;5)\equiv F(G;0)=(-1)^{n/2+1}T_G(0,1)\pmod5$, this is a statement about $T_G(0,1)\bmod5$ for cubic graphs: $(-2)^{n/2}\,T_G(0,1)$ is a square mod 5. I have no proof and did not search the literature.
- **Every snark in the data has $5\mid F(G;5)$**: a census of all bridgeless cubic graphs of girth $\ge5$ with $F(G;4)=0$ (`tools/flowdata_snarks.log`; this includes snarks with 3-cuts) gives 1, 3, 14, 107 graphs of order 10, 18, 20, 22, all with $F\equiv0\pmod5$; so do $J_3,\dots,J_{25}$, the Blanuša chains and $R_2$. Among all cubic graphs only a third have $5\mid F$, and there is no congruence between $F(G;4)$ and $F(G;5)$ mod 5 in general (table of residues for $n=12,14,16$ computed, all combinations occur), so this looks like a genuine property of non-3-edge-colourable graphs. Many of the small snarks are 3-sums or dot products with $P_{10}$ (factors 20 and 24... only the 3-sum factor is proved), which explains part of the sample but not the flower snarks. No proof; literature not searched.

## 7. Conjectures and the strength of the evidence

**Conjecture A.** For every 3-connected cubic graph $G$ on $n$ vertices, $F(G;5)\ge6\cdot2^{n/2}$, with equality iff $G$ arises from $K_4$ by repeatedly inflating vertices into triangles.
*Evidence: strong.* Exhaustive for $n\le20$ (equality cases checked for $n\le16$ completely and for the 100 stored minimisers at $n=18,20$). *Reduction:* with $\varphi(G)=F(G;5)/12$ and $\nu(G)=n-2$, 3-sums are multiplicative in $\varphi$ and additive in $\nu$, so A is equivalent to $\varphi(G)\ge2^{\nu(G)/2}$ with equality only for $K_4$ **on cyclically 4-connected graphs**, where by §4–5 it holds with exponential room ($n=22$: 29760 against 12288; $n=30$: $\le749280$ against 196608). Hence A would follow from any bound $F\ge c\beta^n$ with $\beta>\sqrt2$ on cyclically 4-connected graphs plus a finite check, (The triangle factor for general $q$ is $q-3$, but the analogous inequality is false for $q=4$ because of snarks, so A is specific to $q\ge5$.)

**Conjecture B (weak form).** There are $c>0$, $\beta>\sqrt2$ with $F(G;5)\ge c\,\beta^n$ for all cyclically 4-edge-connected cubic graphs.
*Evidence: moderate.* Exact minima for $n\le22$ grow by $2.18$–$2.24$ per two vertices; climber upper bounds to $n=30$ continue at $2.21$–$2.27$. The data are compatible with $\beta=24^{1/8}=1.4877$ (the Blanuša chain, which however is beaten by 5% at $n=26$ by a non-snark, so the chain is not extremal for the constant $c$ and perhaps not for $\beta$) and with anything in $[1.47,1.495]$; since $F/1.5^n$ of the minima still decreases slowly at $n=30$ (4.29, 4.21, 4.16, 4.07, 4.01, 3.96 for cyc 5), $\beta<1.5$ is likely but $\beta=1.5$ with a slowly varying prefactor is not excluded. **The data do not determine $\beta$ or the extremal family**, and I do not state a sharp form. For cyclically 5-connected graphs of girth $\ge5$ the ratios are $2.21$–$2.24$, the same within the noise; whether cyclic 5- or 6-connectivity pushes the base up to the generic $1.5$ is open — the cyclically 6-connected graphs of the project and the girth-6 climber at $n=42$ stay at $F/1.5^n\approx4.9$–$5.2$, i.e. show no exponential deficit at all.

**What this means for the programme of §18.4.** (1) A counting induction hypothesis must be formulated on cyclically 4-connected graphs, or must carry the 3-cut structure explicitly through the normalisation $\varphi=F/12$, $\nu=n-2$; the hypothesis "$F\ge c^n$ for 3-connected cubic graphs" has the unique natural form A, with the trivial extremal family. (2) The local structures that depress $F$ are exactly the known reducible ones (triangles, and to a smaller extent 5-cycles and Petersen fragments); 4-cycles and 6-cycles help. A Dvořák–Mohar–Šámal style induction therefore has to lose at most a factor 2 per triangle-type reduction and $\approx2.2$ per two vertices elsewhere, and the gap between $\sqrt2$ and $1.487$ is the room available. (3) The count is blind to everything the oddness route struggles with: snarks, oddness 6, template failures and the Jacobsen–Salas graphs all have $F\approx5\cdot1.5^n$ like random graphs.

## 8. Not done / caveats

- Double-star, Szekeres and the House of Graphs snark lists were not used (no local copy); snarks of order $\le22$ and girth $\ge5$ are covered by the census of §6 (`tools/flowdata_snarks_*.json`).
- Exhaustive data for the cyclically 4-connected class include girth 4 up to $n=20$ (`tools/flowdata_below_20ct.log`: of the 573 triangle-free graphs with $F\le13680$ exactly two are cyclically 4-connected, both of girth 5 with $F=13680$); for $n=22$ only girth $\ge5$ was enumerated, and the climber (girth 4 allowed) found nothing lower.
- Climber values for $n\ge24$ are upper bounds from 3–4 minute runs on a loaded machine (two seeds agree exactly at $n=24$ and within 0.7% at $n=26$–$30$ in class cyc 5; class cyc 4 has one seed). The $n=42$ runs are too short to mean anything (§5). The class 3-connected was run by the climber only at $n=12$ (finds 384); it is settled by §4 anyway.
- The per-cycle factors of §3 are a regression with $R^2=0.81$, not a theorem; only the triangle factor is exact.
