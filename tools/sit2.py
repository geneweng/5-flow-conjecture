"""Situation 2: three codim-2 pair-type cuts S12 (z1,z3 in; z2,z4 out; P3 avoids), S13 (z1,z5 in; z2,z6 out;
P2 avoids), S23 (z3,z5 in; z4,z6 out; P1 avoids).  Enumerate placements of z1..z6 in the 8 cells
(membership in S12,S13,S23).  For each pair of cuts whose four uncrossing parts all contain a z, the
MS claims give part-cuts 3/0, hence each part X has |dX| = 6 inside H, so X separates an even number
of paths.  Report surviving placements."""
import itertools
Z = ['z1','z2','z3','z4','z5','z6']
paths = [('z1','z2'),('z3','z4'),('z5','z6')]
cuts = ['S12','S13','S23']
fixed = {
 'z1': {'S12': True,  'S13': True,  'S23': None},
 'z2': {'S12': False, 'S13': False, 'S23': None},
 'z3': {'S12': True,  'S13': None,  'S23': True},
 'z4': {'S12': False, 'S13': None,  'S23': False},
 'z5': {'S12': None,  'S13': True,  'S23': True},
 'z6': {'S12': None,  'S13': False, 'S23': False},
}
avoid = {'S12': ('z5','z6'), 'S13': ('z3','z4'), 'S23': ('z1','z2')}
survivors = []; total = 0
free = [(z, c) for z in Z for c in cuts if fixed[z][c] is None]
for bits in itertools.product((0, 1), repeat=len(free)):
    mem = {z: dict(fixed[z]) for z in Z}
    for (z, c), b in zip(free, bits): mem[z][c] = bool(b)
    if any(mem[a][c] != mem[b][c] for c, (a, b) in avoid.items()): continue
    total += 1
    ok = True
    for X, Y in itertools.combinations(cuts, 2):
        parts = {(True, True): [], (True, False): [], (False, True): [], (False, False): []}
        for z in Z: parts[(mem[z][X], mem[z][Y])].append(z)
        if any(len(p) == 0 for p in parts.values()):
            continue
        for key, zs in parts.items():
            if sum(1 for a, b in paths if (a in zs) != (b in zs)) % 2:
                ok = False
    if ok: survivors.append(mem)
print("placements consistent with the cut definitions:", total, " surviving the parity constraints:", len(survivors))
for mem in survivors[:20]:
    print("  " + "  ".join(f"{z}:" + "".join('1' if mem[z][c] else '0' for c in cuts) for z in Z))
