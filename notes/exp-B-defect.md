# Experiment B — defect distance from the MS template to a 5-flow (2026-09-19)

Implements §18.2 of `oddness.md`. Script: `tools/defect_min.py` (run with no arguments; about 90 s on 2 cores); raw output: `tools/defect_min_results.json`.

## Set-up

For each archived template failure (graph, 2-factor, failing 0-edge choice; $H$ = colours 1, 2 of the canonical colouring) we minimise, over **all** balanced $\pm\tfrac53$ bipartitions, the number of monochromatic $H$-edges ("defects"). The template family is exactly the bipartitions with 0 defects.

Balancedness is modelled exactly by a flow: black/white is balanced iff there is an integer flow $f$ with $|f(e)|\le3$ and net outflow $+5$ at black and $-5$ at white vertices (max-flow/min-cut on the network of `oddness.balanced`). CP-SAT model: $x_v\in\{0,1\}$, $f_e\in[-3,3]$, $d_e\in\{0,1\}$ per $H$-edge with $x_u\oplus x_v\oplus d_e=1$, conservation $\sum_{\rm out}f-\sum_{\rm in}f=10x_v-5$, minimise $\sum d_e$; `num_workers=2`. Every bipartition the solver returned, and every bipartition produced in the enumeration, was re-checked with `oddness.balanced` (assertions in the script; none failed). Every CP-SAT run ended with status OPTIMAL or INFEASIBLE (no time-outs).

Enumeration: optimal **defect sets** are enumerated with no-good cuts on the $d$-variables until infeasible (complete for every instance; cap 200 never reached). For each defect set all $2^{\#\text{comp}-1}$ bipartitions with exactly that defect set (up to global complement) are tested by min cut.

Killers: for every proper colouring of $H$ (8, or 16 when $H$ has an even circuit) the min-cut witness `badcuts.violating_set` (smaller side) is computed; where the instance JSON has a stored field `S` those sets are added. Stored and computed sets coincide for 4 of the 5 data files that have `S`; for `tc_n54` they differ (both are used). A defect edge is classified **bdry** if it lies in the boundary of at least one killer, else **in** if both ends lie inside some killer's smaller side, else **out**.

Instances: 26 files (`data/oddness6_*.json`, `tools/fullkill_*.json`); 5 of the `fullkill_*` files are the same (graph, 2-factor, 0-edges) as a `data/` file, so there are **21 distinct instances**. Abbreviations: `tc` = `oddness6_template_counterexample`, `fk_` = `fullkill_`.

## Results

### Main table

Columns: **min** = minimum number of defects; **c2-only / c1-only** = minimum when defects are allowed only on colour-2 / only on colour-1 (matching) edges; **working / all c2** = number of colour-2 $H$-edges that work as the single defect / number of colour-2 $H$-edges; **bip.** = number of balanced 1-defect bipartitions up to complement; **working bdry/in/out** and **all c2 bdry/in/out** = killer classification of the working defect edges and, as baseline, of all colour-2 edges; **T/N** = balanced 1-defect bipartitions that are / are not a *template* bipartition for another 0-edge choice on the same 2-factor (see below); **sets all-N** = working defect edges none of whose bipartitions is of type T.

| instance | n | H | killers (size/cut) | min | c2-only | c1-only | working / all c2 | bip. | working: bdry/in/out | all c2: bdry/in/out | T/N | sets all-N |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| robust_stage2b_n50 | 50 | 3P+1C | 13/11 15/11 19/11 21/11 21/11 | 1 | 1 | 2 | 17/22 | 71 | 8/7/2 | 12/8/2 | 60/11 | 3 |
| tc | 60 | 3P+0C | 13/11 13/11 15/11 17/11 | 1 | 1 | 2 | 25/27 | 61 | 9/4/12 | 9/4/14 | 41/20 | 8 |
| tc_n42 | 42 | 3P+0C | 13/11 13/11 13/11 15/11 | 1 | 1 | 2 | 18/18 | 41 | 9/3/6 | 9/3/6 | 39/2 | 1 |
| tc_n50 | 50 | 3P+0C | 13/11 15/11 21/11 23/11 | 1 | 1 | 2 | 18/22 | 43 | 8/7/3 | 10/8/4 | 40/3 | 1 |
| tc_n54 | 54 | 3P+0C | 12/6 16/6 21/11 27/11 | 1 | 1 | 2 | 20/24 | 38 | 9/8/3 | 11/10/3 | 30/8 | 4 |
| tc_n60b | 60 | 3P+0C | 13/11 13/11 13/11 13/11 | 1 | 1 | 2 | 23/27 | 57 | 8/3/12 | 8/3/16 | 40/17 | 7 |
| fk_climb_six7_1110_1643 | 42 | 3P+0C | 13/11 13/11 15/11 15/11 | 1 | 1 | 2 | 17/18 | 44 | 8/5/4 | 8/5/5 | 39/5 | 1 |
| fk_climb_six7_43_1905 (= tc_n42) | 42 | 3P+0C | 13/11 13/11 13/11 15/11 | 1 | 1 | 2 | 18/18 | 41 | 9/3/6 | 9/3/6 | 39/2 | 1 |
| fk_climb_six7_514_1901 | 42 | 3P+0C | 13/11 13/11 13/11 15/11 | 1 | 1 | 2 | 18/18 | 39 | 11/4/3 | 11/4/3 | 36/3 | 0 |
| fk_climb_six7_521_1890 | 42 | 3P+0C | 15/11 15/11 19/11 21/11 | 1 | 1 | 2 | 16/18 | 37 | 11/3/2 | 12/4/2 | 35/2 | 0 |
| fk_climb_six7_531_1801 | 42 | 3P+0C | 13/11 13/11 13/11 17/11 | 1 | 1 | 2 | 17/18 | 40 | 9/4/4 | 10/4/4 | 39/1 | 0 |
| fk_climb_six7_538_2174 | 42 | 3P+0C | 15/11 21/11 21/11 21/11 | 1 | 1 | 2 | 17/18 | 33 | 9/8/0 | 9/9/0 | 29/4 | 1 |
| fk_climb_straddle6A_31_22440 | 78 | 3P+0C | 25/11 30/6 33/11 37/11 | 1 | 1 | 2 | 28/36 | 48 | 8/19/1 | 10/24/2 | 28/20 | 11 |
| fk_climb_straddleA_709_2275 | 54 | 3P+1C | 22/6 27/11 27/11 | 1 | 1 | 2 | 19/24 | 80 | 8/10/1 | 8/15/1 | 64/16 | 3 |
| fk_complete_102 | 60 | 3P+0C | 13/11 13/11 15/11 17/11 | 1 | 1 | 2 | 24/27 | 61 | 9/4/11 | 9/4/14 | 41/20 | 7 |
| fk_complete_103 | 60 | 3P+0C | 13/11 13/11 15/11 17/11 | 1 | 1 | 2 | 24/27 | 57 | 9/4/11 | 9/4/14 | 41/16 | 7 |
| fk_complete_106 | 54 | 3P+0C | 12/6 13/11 16/6 27/11 | 1 | 1 | 2 | 22/24 | 35 | 10/9/3 | 10/11/3 | 25/10 | 4 |
| fk_complete_3 (= tc) | 60 | 3P+0C | 13/11 13/11 15/11 17/11 | 1 | 1 | 2 | 25/27 | 61 | 9/4/12 | 9/4/14 | 41/20 | 8 |
| fk_complete_33 | 42 | 3P+0C | 15/11 15/11 21/11 21/11 | 1 | 1 | 2 | 16/18 | 37 | 12/4/0 | 12/6/0 | 36/1 | 0 |
| fk_complete_34 | 42 | 3P+0C | 15/11 15/11 15/11 21/11 | 1 | 1 | 2 | 16/18 | 40 | 12/3/1 | 12/3/3 | 35/5 | 0 |
| fk_complete_501 | 42 | 3P+0C | 13/11 15/11 15/11 19/11 | 1 | 1 | 2 | 16/18 | 40 | 9/4/3 | 9/4/5 | 39/1 | 0 |
| fk_complete_505 | 42 | 3P+0C | 13/11 15/11 17/11 19/11 | 1 | 1 | 2 | 16/18 | 38 | 11/5/0 | 12/5/1 | 38/0 | 0 |
| fk_complete_512 | 42 | 3P+0C | 13/11 13/11 15/11 21/11 | 1 | 1 | 2 | 18/18 | 40 | 10/5/3 | 10/5/3 | 39/1 | 1 |
| fk_complete_7 (= tc_n54) | 54 | 3P+0C | 12/6 16/6 21/11 27/11 | 1 | 1 | 2 | 20/24 | 38 | 8/9/3 | 10/11/3 | 30/8 | 4 |
| fk_complete_8 (= tc_n50) | 50 | 3P+0C | 13/11 15/11 21/11 23/11 | 1 | 1 | 2 | 18/22 | 43 | 8/7/3 | 10/8/4 | 40/3 | 1 |
| fk_complete_9 (= tc_n60b) | 60 | 3P+0C | 13/11 13/11 13/11 13/11 | 1 | 1 | 2 | 23/27 | 57 | 8/3/12 | 8/3/16 | 40/17 | 7 |

(The two rows for `tc_n54`/`fk_complete_7` differ in the bdry/in split only because the stored `S` of the data file is added to the killer list there.)

### (1) Minimum defect count

**The minimum is 1 on all 26 files (21 distinct instances)**, including the 78-vertex instance and the instances whose killers include pair-type 6-cuts.

### (2) Where the defects are

- **Colour.** In every optimal solution of every instance the single defect is a **colour-2 edge**. No single colour-1 defect is ever feasible. This is forced by parity: an $H$-path starts at a $z$ with a colour-1 edge, so a colour-1 defect cuts the path into two pieces of odd order whose majority colours agree (the two ends of the defect edge have the same colour); with all other components properly coloured this gives $|{\rm black}|-|{\rm white}|=\pm2$, and $|A|=|B|$ fails. A colour-2 defect cuts a path into two pieces of even order.
- **How many edges work.** Over the 21 distinct instances, **405 of the 458 colour-2 edges (88%) work as the single defect**; in 3 of the 11 distinct 42-vertex instances all 18 do. (The denominator includes the colour-2 edges on the even $H$-circuit of the two distinct instances that have one; those cannot work alone, see below.) By parity of the 2-factor circuit carrying the edge: 352/378 (93%) of the colour-2 edges on odd circuits work, 53/80 (66%) of those on even circuits.
- **Which path.** Working defects occur on every $H$-path that has a colour-2 edge (the only path without one is the 1-edge path of `fk_climb_six7_521_1890`), never on an even $H$-circuit (a circuit needs an even number of defects), in numbers following the number of colour-2 edges on the component (e.g. `tc_n42`: paths with 17, 15, 7 edges carry 8, 7, 3 working defects = all their colour-2 edges). Positions along the path are spread over the whole path (always odd positions, being colour 2). Per-edge data (path ends, position, circuit) are in the JSON.
- **Relative to the killers.** Working edges: 197 bdry / 123 in / 85 out; all colour-2 edges: 210 / 143 / 105. So the success rate is 94% on killer boundaries, 86% inside, 81% outside: a mild enrichment only; **the defects are not confined to the killers' boundaries**. A sharper test, per balanced 1-defect bipartition (980 in total): take the nearest proper colouring of $H$ (component-wise majority) and its own killer; the defect edge is on that killer's boundary in 240 cases (24%), inside it in 211 (22%), outside in 529 (54%). Note that a defect at position $i$ of a path is the same thing as the proper colouring with the tail beyond $i$ flipped, so a defect changes $b_S-a_S$ for every $S$ meeting that tail, however far the defect edge itself is from $\partial S$.

### (3) All optimal solutions

Complete enumeration for every instance: 16–28 optimal defect edges per instance (column "working"), 33–84 balanced 1-defect bipartitions per instance up to complement (980 over the 21 distinct instances); a working defect edge carries between 1 and 8 balanced bipartitions (of the 4, or 8 with an even $H$-circuit, available).

**Relation to moving a 0-edge (T/N).** A bipartition is a template bipartition for some other 0-edge choice on the same 2-factor iff all matching edges are bichromatic and each odd circuit has a position $p$ with the edges $p+1,p+3,\dots$ bichromatic (convention of `canonical_colouring`; the mirrored convention was tested too and never made a difference), each even circuit keeping its colour-2 class (swapping the two classes of an even circuit was tested too and never helped). Of the 980 balanced 1-defect bipartitions, **814 (83%) are of this kind, always with exactly one 0-edge moved** (type T): these are the known single-0-edge-move repairs of §12 seen from the other side. **166 (17%) are not template bipartitions for any 0-edge choice of this 2-factor** (type N), and 59 of the 405 working defect edges have only type-N bipartitions. Type N is rare on the 42-vertex instances (0–5 bipartitions each) and common on the instances with even circuits in the 2-factor (up to 20 of 48 on the 78-vertex instance), in line with the defects on even circuits, which can never be of type T.

### (4) Control: non-failing 0-edge choices

30 uniformly random 0-edge choices on each of `tc_n42`, `fk_climb_six7_1110_1643`, `fk_complete_505` (42 vertices; no product of positions materialised): all 90 sampled choices were non-failing, and the CP-SAT minimum was **0 for all 90**, as expected; the model agrees with the template test.

Additional: the 7 failing 0-edge choices of the 0-edge-robust instance `data/oddness6_robust_stage2b_n50.json` (field `fails`):

| zero | H | killers (size/cut) | min | c2-only | c1-only | working / all c2 | bip. | working: bdry/in/out | all c2: bdry/in/out | T/N | sets all-N |
|---|---|---|---|---|---|---|---|---|---|---|---|
| [1,0,2,5,2,1] | 3P+1C | 13/11 15/11 19/11 19/11 | 1 | 1 | 2 | 19/22 | 77 | 9/7/3 | 9/7/6 | 63/14 | 4 |
| [1,0,2,5,2,3] | 3P+0C | 13/11 15/11 19/11 23/11 | 1 | 1 | 2 | 19/22 | 41 | 9/8/2 | 10/8/4 | 33/8 | 4 |
| [5,0,2,5,2,1] | 3P+1C | 13/11 15/11 19/11 21/11 21/11 | 1 | 1 | 2 | 17/22 | 71 | 8/7/2 | 12/8/2 | 60/11 | 3 |
| [5,0,2,5,2,2] | 3P+0C | 13/11 13/11 15/11 21/11 | 1 | 1 | 2 | 18/22 | 43 | 10/7/1 | 11/8/3 | 36/7 | 2 |
| [5,0,2,5,2,3] | 3P+0C | 13/11 15/11 19/11 21/11 | 1 | 1 | 2 | 18/22 | 39 | 7/7/4 | 10/8/4 | 31/8 | 4 |
| [5,0,2,5,2,4] | 3P+1C | 13/11 13/11 15/11 15/11 21/11 | 1 | 1 | 2 | 18/22 | 84 | 11/6/1 | 13/6/3 | 71/13 | 2 |
| [5,0,2,5,2,5] | 3P+0C | 13/11 13/11 19/11 21/11 | 1 | 1 | 2 | 17/22 | 35 | 11/5/1 | 13/6/3 | 30/5 | 3 |

### (5) Variants

- Defects on colour-2 edges only: minimum **1** everywhere (same as unrestricted).
- Defects on colour-1 (matching) edges only: feasible everywhere, minimum **2** everywhere (1 is excluded by the parity argument above).

## Reading

1. On every archived failure the template is at Hamming distance exactly 1 from a nowhere-zero 5-flow: one monochromatic colour-2 edge suffices, and the choice of that edge is almost free (88% of all colour-2 edges work; 100% on three 42-vertex instances). The statement "template + one colour-2 defect always works" holds on all 21 + 6 tested (instance, failing choice) pairs; nothing in the data distinguishes a preferred place for the defect, in particular not the killers' boundaries.
2. Most of this is the known single-0-edge-move repair in another guise: 83% of the balanced 1-defect bipartitions are template bipartitions of a 0-edge choice differing in one circuit. The remaining 17% (166 bipartitions, 59 defect edges with no T-type bipartition at all) are 5-flows outside the template family of *every* 0-edge choice of the 2-factor; they are concentrated on 2-factors with even circuits.
3. Matching edges alone can also absorb the failure, at cost 2 instead of 1.
4. The experiment says nothing about instances where the template fails for *all* 0-edge choices (none is known), and all instances here come from the same search pipeline (11-cut and pair-6-cut killers, $H$ = three paths plus at most one even circuit).
