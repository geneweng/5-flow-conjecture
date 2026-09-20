# Independent check of the (C6) count relaxation

*2026-09-19. Independent verifier report. The model was re-implemented from `notes/c6-model-spec.md` alone; none of `tools/count_c6.py`, `count_general.py`, `enum_c6.py`, `enum_covers.py`, `validate_c6.py` was opened. Only `tools/enum_c6_cache.json` was used, as reference results.*

Files: `tools/c6_independent.py` (implementation), `tools/c6_independent_cache.json` (my results), `tools/c6_independent.log` (progress log).

## 1. Implementation

- **Solver: HiGHS MIP (`highspy` 1.x), single thread per solve**, big-M formulation. Parities use integer slack variables (`n = E1 + 2s`, etc.). Conditionals use binaries: `t7 = [n_a >= 7]`, the rule "n_a = 3 and E_nh = 0" via `q + g + h >= 1` with `g -> n_a >= 4`, `h -> E_nh >= 1`, `q ->` degree conditions. Union facts use exact capped sizes `m_a = min(n_a, 6)` and threshold binaries on `sum m_a`, which keeps the big-M at most 96 instead of 8000. Tolerances 1e-9.
- **Lazy constraints:** path connectivity (group 1), class connectivity (group 6), and union facts for unions of more than two atoms (group 7). Unions of one or two atoms are in the model from the start; for J = 2 all unions are. This deviates from the letter of the task ("all unions up front for J <= 3") but is logically equivalent. Any subset of unions is a valid relaxation, so INFEASIBLE stays rigorous. Every FEASIBLE answer is accepted only after an **exact integer checker** confirms the solution. The checker is written separately from the row builder, directly from the spec, and tests every constraint of groups 1–7, including all 2^(2^J) unions and both connectivity conditions. Up-front unions made triples 3–30 times slower with no change in status on the 4 triples tried.
- **Exact second opinion:** all 834 pairs were re-solved with **z3 (QF_LIA, exact arithmetic)** on the same row system. This rules out MIP numerics for the pairs.
- **CP-SAT as last resort:** the scheduler retries HiGHS-UNKNOWN jobs with OR-tools CP-SAT (1 worker, 120 s) on my own rows. **No INFEASIBLE result relies on CP-SAT.** It produced only 8 FEASIBLE sample answers, each exact-checked (see 2).
- **Limits:** 4 worker processes, 45 s HiGHS limit (120 s in the first pass), 120 s CP-SAT limit, hard kill.
- **Timings to ignore:** the machine went into macOS "Maintenance Sleep" three times before I held a `caffeinate -s` assertion. Logged times of 570–1530 s are sleep artefacts. The affected jobs were re-run.

## 2. Results

| stage | jobs | agree with cache | disagree | undecided (mine) |
|---|---|---|---|---|
| all pairs (no two points of one class) | 834 | 834 (534 FEASIBLE, 300 INFEASIBLE) | 0 | 0 |
| the same pairs under z3, exact | 834 | 834 | 0 | 0 |
| certificate sub-lists (cached INFEASIBLE) | 384 (368 triples, 16 four-covers) | 367 INFEASIBLE (364 triples, 3 four-covers), all by HiGHS | 0 | 17 |
| random cached-FEASIBLE triples | 150 | 150 FEASIBLE, atom sizes in the log | 0 | 0 |

Notes on the table:

- **Counts.** "946 pairs" in the task is C(44,2). Excluding two points of the same class leaves 834, which is exactly the number of pair keys in the cache.
- **Items and covers.** I get 44 items and **5,292 minimal covers** (12 + 1,184 + 4,096), as in the spec.
- **Pair statistics** match §5 of the spec: L–L 36/66, L–Pc 144/288, Pc–Pc 120/216, and every pair with a P7 or P641 is feasible.
- **Sample engines.** 142 of the 150 feasible sample triples were solved by HiGHS. The other 8 were solved by CP-SAT after a sleep-induced kill of the HiGHS job; all 150 passed the exact checker.
- **There is no contradiction anywhere.** I never found FEASIBLE where the cache says INFEASIBLE, or the reverse. Task items (4a) and (4b), the solution dump and the constraint-group ablation, therefore did not apply.

**Certificate status from my results alone: 5,275 of the 5,292 minimal covers are certified infeasible.**

- 4,380 covers contain an INFEASIBLE pair.
- 895 covers are certified by a re-proved triple or four-cover.
- 17 covers are undecided: HiGHS at 45–120 s and then CP-SAT at 120 s on one core both hit the time limit. They are neither confirmed nor contradicted. The cache records all 17 as INFEASIBLE.

The 17 undecided covers are:

- **4 three-covers:** `L01b1p1 L12b0p1 P641_01`, `L01b1p1 P641_00 Pc0p1_01`, `L01b1p1 P641_00 Pc1p0_01`, `L01b1p1 P641_01 Pc1p0_00`.
- **13 four-point covers**, all made of P7/P641 points only:
  - `P641_00 P641_01 P641_10 P641_11`
  - `P641_00 P641_01 P641_10 P7_11`
  - `P641_00 P641_01 P641_11 P7_10`
  - `P641_00 P641_01 P7_10 P7_11`
  - `P641_00 P641_10 P641_11 P7_01`
  - `P641_00 P641_10 P7_01 P7_11`
  - `P641_00 P641_11 P7_01 P7_10`
  - `P641_00 P7_01 P7_10 P7_11`
  - `P641_01 P641_10 P641_11 P7_00`
  - `P641_01 P641_10 P7_00 P7_11`
  - `P641_01 P641_11 P7_00 P7_10`
  - `P641_10 P641_11 P7_00 P7_01`
  - `P641_11 P7_00 P7_01 P7_10`

Four-point covers that HiGHS did prove are `P7_00 P7_01 P7_10 P7_11`, `P641_10 P7_00 P7_01 P7_11` and `P641_01 P7_00 P7_10 P7_11`.

These 17 deserve a longer independent run: a multi-hour limit, or several cores per solve.

## 3. Validity of the constraints (§3 of the spec)

I checked each constraint against the setting of §1 and found **no invalid constraint**.

- **(T1)** follows from the killer definition. Matching edges inside S join black to white, so k <= c1, with equality iff all S-ends are black. Colour-2 edges pair all vertices except the z's, so k <= c2 + #(black separated ends), with equality likewise. The numbers 4 = 4 = 2+2, 5 = 5 = 2+3 and 4 = 4 = 1+3 make both tight.
- **(T2) and (T3)** follow from the alternation of colours and vertex colours along a component of H.
- **Group 4 (segment count):** checked, including one-vertex segments.
- **Group 5:**
  - The parity constraints and "E2 + E_nh even" rest on the fact that the only colour-0/3 edges between atoms are the P641 non-H edges. This holds because m = c1 + c2 for pair and P7 killers.
  - "2n >= 2 visits + H-crossings": visits and H-segments are vertex-disjoint, and each has at least one vertex.
  - "n + E_nh >= 3 at a z", the isolated-atom rule, and n >= 42: fine.
- **Group 6:** a circuit of at least 7 vertices must leave an atom of at most 6 vertices. The n = 3 rule is valid because girth >= 6 excludes a triangle. It also correctly makes two z's in a 3-vertex atom infeasible.
- **Group 7:**
  - d >= 3: G is 3-edge-connected.
  - d >= n + 2 for n <= 5: girth >= 6, so the side is a forest.
  - d >= 6 when both sides have at least 4 vertices: a side with d <= 5 and at least 4 vertices cannot be a forest, so both sides contain a cycle, contradicting cyclic 6-edge-connectivity.
- **Bounds.** x, w <= 11 and u <= 5 are harmless; in fact x <= 7. n_a <= 500 is harmless by the shrinking argument.

## 4. One doubt, about §4 (completeness of the cover list), not about §3

The item list is **not closed under the normalisation it claims**. The L12 lines in the list always have z3 in S and z4 outside S. A killer of the z1-black colouring of class (1,0) or (1,1) whose third path P0 avoids the cut must contain z4, since a killer contains the black ends. Replacing it by its complement is not allowed inside one configuration: complementing a killer complements the colouring of the H-circuits too, and the model assumes a common H-circuit colouring (the u-variables).

I ran a test. Take one killer for each of the four z1-black colourings: 14 choices each, 38,416 assignments.

- **6,512 assignments contain no listed minimal cover**, for example {P7_11, the line with z4,z5 in S and z1,z2,z3,z6 outside (kills class (1,0)), P7_01, P7_00}.
- With path permutations only, 3,456 assignments still fail. With path reversals only, 800 fail.
- The list becomes complete **only modulo the symmetry group generated by permutations of the three paths and reversal of individual paths** (z_{2i+1} <-> z_{2i+2}, beta -> 1 - beta). Under this group of order 48 every one of the 38,416 assignments contains an image of a listed cover (0 failures). Global complement adds nothing: the orbit is 34,096 covers with or without it.

I checked that the system of §3 is invariant under these symmetries.

- Under path reversal, the colour rule is preserved: leaving/entering swaps together with beta.
- Visits become exits + [end in a], which equals entries + [start in a].
- The isolated-atom rule is equivalent by flow conservation.
- Connectivity to the start is equivalent to connectivity to the end.

So the theorem is fine, but **the spec (and the proof text, if it says the same) should state explicitly that the 5,292 covers are representatives up to path permutation and path reversal, and that the model is invariant under them.** As written — "normalise z1 black ... a killer contains the black ends" together with the fixed L12 orientation — the claim "every configuration contains a listed cover" is false.

The script for the main test is `tools/c6_cover_completeness.py` (imports the item list from `tools/c6_independent.py`; runs in about 10 s). The symmetry action on items is (kind, membership, (crossing path, beta)). The subgroup breakdown (permutations only / reversals only) was a one-off variant of the same script.

## 5. Reproduce

```
python3 tools/c6_independent.py stats      # 44 items, 5292 covers, pair statistics
python3 tools/c6_independent.py run        # pairs, feasible sample, z3 on pairs, certificate (resumable)
python3 tools/c6_independent.py report     # agreement table + certificate count
python3 tools/c6_independent.py one L01b0p0 L02b1p1    # prints a full solution
```

Use `caffeinate -s` on macOS.
