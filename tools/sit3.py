"""Two codim-2 pair cuts on different pairs: S12 (z1,z3 in; z2,z4 out; P3 avoids) and S13 (z1,z5 in;
z2,z6 out; P2 avoids).  Same parity test on the uncrossing parts."""
import itertools
Z = ['z1','z2','z3','z4','z5','z6']; paths = [('z1','z2'),('z3','z4'),('z5','z6')]; cuts = ['S12','S13']
fixed = {'z1': {'S12': True, 'S13': True}, 'z2': {'S12': False, 'S13': False},
         'z3': {'S12': True, 'S13': None}, 'z4': {'S12': False, 'S13': None},
         'z5': {'S12': None, 'S13': True}, 'z6': {'S12': None, 'S13': False}}
avoid = {'S12': ('z5','z6'), 'S13': ('z3','z4')}
free = [(z, c) for z in Z for c in cuts if fixed[z][c] is None]; total = 0; surv = []
for bits in itertools.product((0, 1), repeat=len(free)):
    mem = {z: dict(fixed[z]) for z in Z}
    for (z, c), b in zip(free, bits): mem[z][c] = bool(b)
    if any(mem[a][c] != mem[b][c] for c, (a, b) in avoid.items()): continue
    total += 1
    parts = {(True, True): [], (True, False): [], (False, True): [], (False, False): []}
    for z in Z: parts[(mem[z]['S12'], mem[z]['S13'])].append(z)
    if any(len(p) == 0 for p in parts.values()):
        surv.append((mem, "nested (some part has no z)")); continue
    bad = [key for key, zs in parts.items() if sum(1 for a, b in paths if (a in zs) != (b in zs)) % 2]
    if not bad: surv.append((mem, "parity ok"))
print("placements:", total, "survivors:", len(surv))
for mem, why in surv:
    print("  " + "  ".join(f"{z}:" + "".join('1' if mem[z][c] else '0' for c in cuts) for z in Z), "|", why)
