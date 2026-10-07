#!/usr/bin/env python3
"""tests/validate_cert.py -- checks behind the certified census below 10^28 (a few minutes on 30 threads).

    python3 tests/validate_cert.py [--threads 30]

1. results: results/b175531.txt and factorisations.{txt,json} are exactly what cn2_make_results.py rebuilds from the
   certificate's records (every term re-verified, the tail survivors re-verified from scratch).
2. sources: the engine sources in cn2x/ have the sha256 recorded in results/cert_1e28/ENGINES.sha256.
3. pair classes: every certified term lies in the residue class that the pair family builds from its two largest primes.
4. exact count: cn2pair --count-only reproduces an independent Python count (X = 1e20, T = 1000, q <= 6000).
5. the capped half against the proved census below 10^25: a capped cn2xh run with LASTCALC at T = 800 must find every
   term whose second-largest prime is <= 800 (9 of the 14).
6. the pair tail against the proved census below 10^25: cn2pair at T = 10^4, q <= 2*10^4 must find exactly the terms
   whose second-largest prime exceeds 10^4 (2 of the 14), and nothing else.
"""
import argparse, hashlib, os, subprocess, sys, tempfile
from multiprocessing import Pool

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from sympy import factorint                                     # noqa: E402
from cn2_build import compile_cmd                               # noqa: E402
from cn2_order2_defs import verify_candidate, load_reference   # noqa: E402
import cn2_pair_tail                                            # noqa: E402

TERMS = load_reference(os.path.join(HERE, 'results', 'b175531.txt'))
FAILS = []


def check(name, ok, detail=''):
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ''), flush=True)
    if not ok:
        FAILS.append(name)


def survivors(path):
    jobs = []
    for l in open(path):
        t = l.split()
        if len(t) == 2 and t[0] in 'HN':
            jobs.append((int(t[1]), 1, 'howe'))
        elif len(t) == 5 and t[0] == 'S':
            jobs.append((int(t[1]), int(t[4]), 'howe'))
    return jobs


def verified(jobs, X):
    with Pool(8) as pl:
        return {n for n in pl.imap_unordered(verify_candidate, jobs, chunksize=512) if n is not None and n < X}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--threads', type=int, default=30)
    a = ap.parse_args()
    th = str(a.threads)

    r = subprocess.run([sys.executable, os.path.join(HERE, 'cn2_make_results.py'), '--check'], capture_output=True, text=True)
    check("1. results files rebuild exactly from the certificate records", r.returncode == 0, r.stdout.strip().splitlines()[-1])

    want = {}
    for l in open(os.path.join(HERE, 'results', 'cert_1e28', 'ENGINES.sha256')):
        if l.strip() and not l.startswith('#'):
            h, path = l.split()
            if path.endswith('.c'):
                want[os.path.basename(path)] = h
    got = {f: hashlib.sha256(open(os.path.join(HERE, 'cn2x', f), 'rb').read()).hexdigest() for f in want}
    check("2. engine sources match the certificate's recorded sha256", got == want, ', '.join(sorted(want)))

    bad = []
    for n in TERMS:
        f = sorted(factorint(n))
        p, q = f[-2], f[-1]
        lo = q * (q * q - 1)
        if not (n > lo and any(n % M == N0 for N0, M in cn2_pair_tail.pair_classes(p, q, 'howe'))):
            bad.append(n)
    check("3. every certified term lies in its own pair class", not bad, f"{len(TERMS)} terms")

    with tempfile.TemporaryDirectory() as tmp:
        exe = {}
        for name in ('cn2xh', 'cn2pair'):
            exe[name] = os.path.join(tmp, name)
            subprocess.run(compile_cmd(os.path.join(HERE, 'cn2x', name + '.c'), exe[name]), check=True)
        env = dict(os.environ, CN2X_MODE='howe')

        X, T, Q = 10 ** 20, 1000, 6000
        ps = list(cn2_pair_tail.primerange(T + 1, Q + 1))
        py = sum(cn2_pair_tail.pair_count(p, q, X) for i, q in enumerate(ps) for p in ps[:i])
        d = os.path.join(tmp, 'count'); os.mkdir(d)
        out = subprocess.run([exe['cn2pair'], str(X), str(T), str(Q), d, th, '0', '64', '0', '1'], env=env,
                             capture_output=True, text=True).stdout
        c = int(out.split('"candidates":')[1].split(',')[0])
        check("4. cn2pair count-only equals the independent Python count", c == py, f"{c:,} vs {py:,}")

        X, T = 10 ** 25, 800
        below = [n for n in TERMS if n < X]
        need = {n for n in below if sorted(factorint(n))[-2] <= T}
        d = os.path.join(tmp, 'capped'); os.mkdir(d)
        subprocess.run([exe['cn2xh'], str(X), d, th, '0', '1'], env=dict(env, CN2X_TCAP=str(T), CN2X_LASTCALC='1'),
                       capture_output=True, text=True)
        got = verified(survivors(os.path.join(d, 'out.txt')), X)
        check("5. capped + LASTCALC (X = 1e25, T = 800) finds every term it guarantees",
              need <= got and got <= set(below), f"guaranteed {len(need)}, found {len(got)} of {len(below)}")

        X, T, Q = 10 ** 25, 10000, 20000
        need = {n for n in below if sorted(factorint(n))[-2] > T}
        d = os.path.join(tmp, 'pair'); os.mkdir(d)
        subprocess.run([exe['cn2pair'], str(X), str(T), str(Q), d, th, '0', '64', '0', '0'], env=env,
                       capture_output=True, text=True)
        got = verified(survivors(os.path.join(d, 'out.txt')), X)
        check("6. pair tail (X = 1e25, T = 1e4, q <= 2e4) finds exactly the terms it must", got == need,
              f"required {len(need)}, found {len(got)}")

    print(f"\n{'ALL PASS' if not FAILS else 'FAILED: ' + ', '.join(FAILS)}")
    sys.exit(1 if FAILS else 0)


if __name__ == '__main__':
    main()
