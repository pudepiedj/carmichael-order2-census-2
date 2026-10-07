#!/usr/bin/env python3
r"""
Driver for cn2x/cn2x.c, the compiled exhaustive CN2 enumerator.

Runs the engine, then does everything that needs arbitrary precision or factoring:
  * every LIST survivor (m, pmax, k, s) is factored with sympy and kept only if s is a
    squarefree product of primes in (pmax, T] and p^2 - 1 | n - 1 for all p | n = m s;
  * every hit, from either path, is re-verified from scratch (factorint of n);
  * the hit set is checked against the six CN2s below 10^22 (complete by Goutier).

    python cn2_exhaustive_c.py 1e16 --threads 30
    python cn2_exhaustive_c.py --parity            # counters vs cn2x_python_reference.json
    python cn2_exhaustive_c.py --scan 14 19 --threads 30 --json cn2x_scan.jsonl
"""

import argparse
import json
import math
import os
import subprocess
import sys
import time
from multiprocessing import Pool

from sympy import factorint

from cn2_build import compile_cmd

HERE = os.path.dirname(os.path.abspath(__file__))
ENGINE = os.path.join(HERE, 'cn2x', 'cn2x')
SOURCE = os.path.join(HERE, 'cn2x', 'cn2x.c')
KNOWN_BELOW_1E22 = [443372888629441, 39671149333495681, 842526563598720001,
                    2380296518909971201, 3188618003602886401, 4208895375600667752001]
COUNTERS = ('nodes', 'children', 'last_tests', 'list_nodes', 'list_candidates',
            'fermat_passes', 'factorisations', 'pruned_budget')


def build():
    if not os.path.exists(ENGINE) or os.path.getmtime(ENGINE) < os.path.getmtime(SOURCE):
        subprocess.run(compile_cmd(SOURCE, ENGINE, ('-Wall',)), check=True)


def icbrt(n):
    c = int(round(n ** (1 / 3)))
    while c ** 3 > n:
        c -= 1
    while (c + 1) ** 3 <= n:
        c += 1
    return c


def is_rigid_cn2(n):
    f = factorint(n)
    return len(f) >= 3 and all(e == 1 for e in f.values()) and all((n - 1) % (p * p - 1) == 0 for p in f)


def _verify_survivor(args):
    m, pmax, k, s, T = args
    f = factorint(s)
    if any(e > 1 for e in f.values()):
        return None
    ps = sorted(f)
    if ps[0] <= pmax or ps[-1] > T or k + len(ps) < 3:
        return None
    n = m * s
    return n if all((n - 1) % (p * p - 1) == 0 for p in ps) else None


def run(X, threads=1, ratio=1.0, fermat=1, donate_q=10**10, msplit=1000, sieve=1, workdir=None, verify_procs=8):
    build()
    workdir = workdir or os.path.join(HERE, 'cn2x', 'runs')
    os.makedirs(workdir, exist_ok=True)
    out = os.path.join(workdir, f'cn2x_{X}_t{threads}.txt')
    t0 = time.time()
    proc = subprocess.run([ENGINE, str(X), out, str(threads), str(ratio), str(fermat), str(donate_q), str(msplit), str(sieve)],
                          capture_output=True, text=True, check=True)
    engine_wall = time.time() - t0
    stats = json.loads(proc.stdout.strip().splitlines()[-1])

    hits, survivors = set(), []
    T = icbrt(X - 1)
    with open(out) as fh:
        for line in fh:
            parts = line.split()
            if parts[0] == 'H':
                hits.add(int(parts[1]))
            elif parts[0] == 'S':
                survivors.append((int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]), T))
    t1 = time.time()
    if survivors:
        with Pool(verify_procs) as pool:
            for n in pool.imap_unordered(_verify_survivor, survivors, chunksize=256):
                if n is not None:
                    hits.add(n)
    verified = sorted(n for n in hits if n <= X and is_rigid_cn2(n))
    stats.update(engine_wall=engine_wall, verify_seconds=time.time() - t1,
                 survivors=len(survivors), hits=[str(n) for n in verified], n_hits=len(verified),
                 unverified_hits=len(hits) - len(verified))
    if X <= 10**22:
        expected = [n for n in KNOWN_BELOW_1E22 if n <= X]
        stats['validation'] = 'OK' if verified == expected else f'MISMATCH expected {expected}'
    else:
        stats['validation'] = f'beyond ground truth; known below 1e22 all present: ' \
                              f'{all(str(n) in stats["hits"] for n in KNOWN_BELOW_1E22)}'
    return stats


def parity(threads_list=(1, 30)):
    ref = json.load(open(os.path.join(HERE, 'cn2x_python_reference.json')))
    ok = True
    for r in ref:
        X = int(r['X'])
        for th in threads_list:
            s = run(X, threads=th, sieve=0)          # the Python reference has no small-prime sieve
            diffs = {c: (r[c], s[c]) for c in COUNTERS if r[c] != s[c]}
            same = not diffs and s['n_hits'] == r['n_hits']
            ok &= same
            print(f"  10^{math.log10(X):.0f} threads={th:>2}: counters {'IDENTICAL' if same else 'DIFFER ' + str(diffs)}"
                  f"  python {r['seconds']:.2f}s  C {s['search_seconds']:.3f}s  ({s['validation']})", flush=True)
    print("PARITY", "PASSED" if ok else "FAILED")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('X', nargs='?', type=float)
    ap.add_argument('--threads', type=int, default=1)
    ap.add_argument('--ratio', type=float, default=1.0)
    ap.add_argument('--fermat', type=int, default=1, choices=[0, 1, 2])
    ap.add_argument('--donate-q', type=int, default=10**10)
    ap.add_argument('--msplit', type=int, default=1000)
    ap.add_argument('--no-sieve', action='store_true')
    ap.add_argument('--scan', nargs=2, type=float)
    ap.add_argument('--step', type=float, default=1.0)
    ap.add_argument('--parity', action='store_true')
    ap.add_argument('--json', metavar='PATH')
    args = ap.parse_args(argv)
    if args.parity:
        parity()
        return
    if args.scan:
        lo, hi = args.scan
        bounds = [int(round(10 ** (lo + i * args.step))) for i in range(int(round((hi - lo) / args.step)) + 1)]
    else:
        bounds = [int(args.X)]
    print(f"{'log10 X':>8} {'nodes':>15} {'list cand':>15} {'survivors':>10} {'hits':>5} "
          f"{'search s':>9} {'verify s':>9}  validation")
    for X in bounds:
        s = run(X, args.threads, args.ratio, args.fermat, args.donate_q, args.msplit, 0 if args.no_sieve else 1)
        print(f"{math.log10(X):>8.2f} {s['nodes']:>15,} {s['list_candidates']:>15,} {s['survivors']:>10,} "
              f"{s['n_hits']:>5} {s['search_seconds']:>9.2f} {s['verify_seconds']:>9.2f}  {s['validation']}", flush=True)
        if args.json:
            with open(args.json, 'a') as fh:
                fh.write(json.dumps(s) + '\n')


if __name__ == '__main__':
    main()
