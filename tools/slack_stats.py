"""Over generated oddness-6 instances: for every path colouring that survives, the slack
min_S [3 d(S) - 5(b_S - a_S)] of the best cut realizing a forbidden pattern (all three paths
separated with black ends inside, or any pair separated).  Slack 0 = a cut with 3m = 5k exactly."""
import sys, random, itertools, collections
from gen_cyc6 import generate_straddle, generate
from oddness import canonical_colouring, H_components, partition
from badcuts import violating_set, classify, profiles
import networkx as nx
INF = 10 ** 6
def best_forced(G, black, inside, outside):
    D = nx.DiGraph()
    for u, v in G.edges(): D.add_edge(u, v, capacity=3); D.add_edge(v, u, capacity=3)
    for v in G.nodes():
        if v in black: D.add_edge('s', v, capacity=5)
        else: D.add_edge(v, 't', capacity=5)
    for v in inside: D.add_edge('s', v, capacity=INF)
    for v in outside: D.add_edge(v, 't', capacity=INF)
    cut, (X, Y) = nx.minimum_cut(D, 's', 't'); S = set(X) - {'s'}
    return cut - 5 * len(black), S
def sample_zero_choices(odd, k, rng):
    total = 1
    for C in odd: total *= len(C)
    if total <= 4 * k:
        space = list(itertools.product(*[range(len(C)) for C in odd])); return rng.sample(space, min(k, total))
    seen = set()
    while len(seen) < k: seen.add(tuple(rng.randrange(len(C)) for C in odd))
    return sorted(seen)
name = sys.argv[1]; rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 5); ng = int(sys.argv[3]) if len(sys.argv) > 3 else 6
pr = profiles[name]
graphs = generate_straddle(pr[0], pr[1], rng, ng) if isinstance(pr, tuple) else generate(pr, rng, ng)
inst_hist = collections.Counter(); by_killed = collections.defaultdict(collections.Counter); slack_types = collections.Counter()
for G, circuits in graphs:
    M = {frozenset((u, v)) for u, v in G.edges()} - {frozenset((C[i], C[(i + 1) % len(C)])) for C in circuits for i in range(len(C))}
    odd = [C for C in circuits if len(C) % 2]
    for zc in sample_zero_choices(odd, 40, rng):
        col = canonical_colouring(G, M, circuits, zc); comps = H_components(G, col)
        path_idx = [i for i, (_, ends) in enumerate(comps) if ends]; t = len(path_idx)
        if t != 3: continue
        paths = [tuple(comps[i][1]) for i in path_idx]; ref = [rng.randint(0, 1) for _ in comps]
        killed = 0; min_slack = None
        for bits in itertools.product((0, 1), repeat=t):
            f = list(ref)
            for i, b in zip(path_idx, bits): f[i] = b
            black = partition(comps, f)
            if violating_set(G, black) is not None: killed += 1; continue
            best = None
            for sub in [tuple(range(t))] + list(itertools.combinations(range(t), 2)):
                ins = [z for i in sub for z in paths[i] if z in black]; outs = [z for i in sub for z in paths[i] if z not in black]
                sl, S = best_forced(G, black, ins, outs)
                if best is None or sl < best[0]: best = (sl, classify(G, S, col, comps, black)[:5])
            if min_slack is None or best[0] < min_slack[0]: min_slack = best
        inst_hist[killed] += 1
        if min_slack is not None:
            by_killed[killed][min_slack[0]] += 1
            if min_slack[0] == 0: slack_types[min_slack[1]] += 1
print(f"{name}: killed-colourings histogram {dict(sorted(inst_hist.items()))}")
print("  min survivor slack, by number killed:", {k: dict(sorted(v.items())) for k, v in sorted(by_killed.items())})
print("  cut types realizing slack 0:", dict(slack_types.most_common(10)))
