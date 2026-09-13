# Counterexamples to the fixed-2-factor flow-partition template at oddness 6

Three cubic graphs, each cyclically 6-edge-connected with girth 6, each with a 2-factor of six odd circuits and a choice of one 0-edge per odd circuit, such that in the canonical colouring (matching = colour 1, circuit edges alternately 2 and 3, the 0-edges colour 0) the subgraph H of colours 1 and 2 consists of three paths and no even circuit, and **all eight proper 2-colourings of H give unbalanced valuations** w = ±5/3 (Mazzuoccolo–Steffen's construction fails for this 2-factor and these 0-edges). Verify with `python3 tools/verify_fullkill.py data/<file>.json` (checks cubicity, girth, cyclic 6-edge-connectivity by SAT, the 2-factor, and every witness cut by arithmetic). Each file records the witness set S per colouring; the files with `_edges.txt` give edge lists with colours.

| file | vertices | circuit lengths | killing cuts |
|---|---|---|---|
| `oddness6_template_counterexample.json` | 60 | 7,6,7,6,7,7,6,7,7 | four 11-cuts of 13, 13, 15, 17 vertices |
| `oddness6_template_counterexample_n60b.json` | 60 | 7,6,7,6,7,7,6,7,7 | four 11-cuts of 13 vertices each |
| `oddness6_template_counterexample_n50.json` | 50 | 7,7,7,7,7,7,8 | four 11-cuts of 13, 15, 21, 23 vertices |

For the first instance all 117,649 other choices of 0-edges admit a balanced colouring (`tools/zero_sweep.py`), so the failure is an isolated point of the 0-edge space; the graphs have nowhere-zero 5-flows.

## The 60-vertex instance with four 13-vertex cuts, in full

Vertices 0..59. 2-factor circuits (in cyclic order):

- C0 (7): 0 1 2 3 4 5 6
- C1 (6): 7 8 9 10 11 12
- C2 (7): 13 14 15 16 17 18 19
- C3 (6): 20 21 22 23 24 25
- C4 (7): 26 27 28 29 30 31 32
- C5 (7): 33 34 35 36 37 38 39
- C6 (6): 40 41 42 43 44 45
- C7 (7): 46 47 48 49 50 51 52
- C8 (7): 53 54 55 56 57 58 59

0-edges (one per odd circuit): 4-5, 19-13, 28-29, 37-38, 47-48, 59-53

Perfect matching (colour 1): 0-9, 1-33, 2-15, 3-57, 4-43, 5-31, 6-22, 7-23, 8-35, 10-28, 11-59, 12-46, 13-50, 14-34, 16-20, 17-49, 18-30, 19-27, 21-51, 24-42, 25-36, 26-48, 29-47, 32-38, 37-54, 39-45, 40-53, 41-56, 44-55, 52-58

Paths of H (end to end): 28..59, 47..19, 4..37

Witness sets S (black-heavy side) per path colouring (bits = which paths are flipped from the base colouring):

- (0, 0, 0): |S|=47, S = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 20, 21, 22, 23, 24, 25, 33, 34, 35, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59]
- (0, 0, 1): |S|=47, S = [0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 20, 21, 22, 23, 24, 25, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59]
- (0, 1, 0): |S|=47, S = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 33, 34, 35, 39, 40, 41, 42, 43, 44, 45, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59]
- (0, 1, 1): |S|=47, S = [0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59]
- (1, 0, 0): |S|=13, S = [3, 4, 5, 26, 27, 28, 29, 30, 31, 32, 46, 47, 48]
- (1, 0, 1): |S|=13, S = [26, 27, 28, 29, 30, 31, 32, 36, 37, 38, 46, 47, 48]
- (1, 1, 0): |S|=13, S = [3, 4, 5, 13, 18, 19, 26, 27, 28, 29, 30, 31, 32]
- (1, 1, 1): |S|=13, S = [13, 18, 19, 26, 27, 28, 29, 30, 31, 32, 36, 37, 38]

For each S: 7 matching edges and 4 circuit edges of colour 2 leave S, |S ∩ Z| = 3 (one end of each path, all black), and 5(b_S − a_S) = 35 > 33 = 3|∂S|.
