import sys, itertools
sys.path.insert(0,'/Users/gweng/cc_projects/5-flow-conjecture/tools')
import c6_independent as C
def canon(it):
    kind = 'L' if it['name'][0]=='L' else it['name'].split('_')[0][:2] if it['name'].startswith('Pc') else it['name'].split('_')[0]
    cross=None
    for i,(k,b) in enumerate(it['path']):
        if k=='cross': cross=(i,b)
    return (kind, it['memb'], cross)
listed = {n:canon(it) for n,it in C.ITEMS.items()}
assert len(set(listed.values()))==44
def act(item, perm, rev, comp):
    kind, memb, cross = item
    m=[0]*6
    for i in range(3):
        s,t = memb[2*i], memb[2*i+1]
        b = cross[1] if cross and cross[0]==i else None
        if rev[i]: s,t=t,s
        if comp: s,t=1-s,1-t
        m[2*perm[i]],m[2*perm[i]+1]=s,t
    if cross:
        i,b=cross
        if rev[i]: b=1-b
        if comp: b=1-b
        cross=(perm[i],b)
    return (kind, tuple(m), cross)
G=[(p,r,c) for p in itertools.permutations(range(3)) for r in itertools.product((0,1),repeat=3) for c in (0,1)]
covers = C.minimal_covers()
orbit=set()
for cov in covers:
    for g in G:
        orbit.add(frozenset(act(listed[n],*g) for n in cov))
print('orbit of listed covers:', len(orbit))
# full family: killers of each z-colouring pi (pi[i]=1 iff start of P_i black)
def killers(pi):
    base=[]
    for i in range(3): base += [pi[i],1-pi[i]]
    out=[('P7',tuple(base),None),('P641',tuple(base),None)]
    for i in range(3):
        for p in (0,1):
            m=list(base); m[2*i]=m[2*i+1]=p
            out.append(('Pc',tuple(m),(i,pi[i])))
            out.append(('L',tuple(m),None))
    return out
cols=list(itertools.product((0,1),repeat=3))
K={pi:killers(pi) for pi in cols}
# index orbit covers by item for fast superset test
orb=list(orbit); orbset=set(orb)
def contains_cover(image):
    im=list(set(image))
    for r in range(2,len(im)+1):
        for sub in itertools.combinations(im,r):
            if frozenset(sub) in orbset: return True
    return False
def test(colourings):
    bad=0; ex=None; tot=0
    for choice in itertools.product(*[K[pi] for pi in colourings]):
        tot+=1
        if not contains_cover(choice):
            bad+=1; ex=ex or choice
    return tot,bad,ex
black=[pi for pi in cols if pi[0]==1]
print('z1-black colourings only, WITHOUT symmetry (listed covers only):')
orb_save=orbset; orbset=set(frozenset(listed[n] for n in cov) for cov in covers)
print(test(black))
orbset=orb_save
print('z1-black colourings only, with symmetry orbit:')
print(test(black))

print('all 8 colourings (h fixed), DFS for a killer assignment whose image contains no orbit cover:')
sys.setrecursionlimit(10000)
cnt=[0]
def dfs(idx, image):
    cnt[0]+=1
    if idx==len(cols):
        return image
    for k in K[cols[idx]]:
        im=image if k in image else image+[k]
        if k not in image and contains_cover(im):   # any new cover must use k; full test is fine
            continue
        r=dfs(idx+1, im)
        if r: return r
    return None
r=dfs(0,[])
print('nodes',cnt[0],'counterexample:',r)
