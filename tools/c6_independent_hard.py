"""Longer runs of the independent model (c6_independent.py, HiGHS) on the sub-lists it left UNKNOWN at 120 s."""
import sys, os, json, time, multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c6_independent as C
T = float(os.environ.get("T", "1800"))
def job(k):
    t0 = time.time(); st, sol, rounds, dt, mod = C.solve_config(k.split(), tlim=T)
    return k, st, rounds, time.time() - t0
if __name__ == "__main__":
    mine = json.load(open("c6_independent_cache.json"))
    todo = [k for k, v in mine.items() if v["status"] == "UNKNOWN"]; print(len(todo), "unresolved", flush=True)
    with mp.Pool(4) as P:
        for k, st, rounds, dt in P.imap_unordered(job, todo):
            print(f"{k:50s} {st} rounds={rounds} {dt:.0f}s", flush=True)
