"""Covers of the four classes by (4,2)-type 6-cuts only (codim 2 or 3).  A cut is described by its
z-membership pattern: separated pair {a,b}, the in-end of path a normalised to z_{2a-1}, the in-end of
b (either), and the third path both-in or both-out.  Codim-2 cuts kill the two classes of their
relation; codim-3 cuts kill one of them (enumerate both).  Any pair of cuts whose four uncrossing parts
all contain a z gives (MS) part-cuts of size 6 inside H, hence each part separates an even number of paths.
We search for a set of <= 4 cuts covering all 4 classes and passing every parity test."""
import itertools
paths = [(0, 1), (2, 3), (4, 5)]          # z indices
def classes_of(rel):                      # rel: dict pair->sign(+1 same,-1 different); class = (e1e2, e1e3)
    out = []
    for x, y in itertools.product((1, -1), repeat=2):
        e = {0: 1, 1: x, 2: y}            # eps of the three paths (up to complement)
        if all(e[a] * e[b] == s for (a, b), s in rel.items()): out.append((x, y))
    return out
cut_types = []
for a, b in [(0, 1), (0, 2), (1, 2)]:
    c = 3 - a - b
    for bend in (0, 1):
        for third in (0, 1):
            mem = [False] * 6
            mem[paths[a][0]] = True; mem[paths[b][bend]] = True
            if third: mem[paths[c][0]] = mem[paths[c][1]] = True
            same = (bend == 0)            # in-ends z_{2a-1}, z_{2b-1}: "same colour" relation
            rel = {(a, b): 1 if same else -1}
            killed2 = classes_of(rel)
            cut_types.append((tuple(mem), (a, b), killed2, 2))
            for k in killed2:
                cut_types.append((tuple(mem), (a, b), [k], 3))
def parity_ok(m1, m2):
    parts = {}
    for i in range(6): parts.setdefault((m1[i], m2[i]), []).append(i)
    if len(parts) < 4: return True
    for zs in parts.values():
        if sum(1 for p, q in paths if (p in zs) != (q in zs)) % 2: return False
    return True
allc = {(x, y) for x in (1, -1) for y in (1, -1)}
found = 0; checked = 0
for r in range(2, 5):
    for combo in itertools.combinations_with_replacement(range(len(cut_types)), r):
        killed = set()
        for i in combo: killed |= set(cut_types[i][2])
        if killed != allc: continue
        checked += 1
        if all(parity_ok(cut_types[i][0], cut_types[j][0]) for i, j in itertools.combinations(combo, 2)):
            found += 1
            if found <= 8:
                print("SURVIVING COVER:", [(cut_types[i][1], ''.join('1' if m else '0' for m in cut_types[i][0]), cut_types[i][2], f"codim{cut_types[i][3]}") for i in combo])
print(f"covers checked: {checked}, surviving the pairwise parity test: {found}")
