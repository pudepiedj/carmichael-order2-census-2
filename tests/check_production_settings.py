#!/usr/bin/env python3
"""tests/check_production_settings.py -- run both halves at the certificate's own settings on the proved census below 10^25.

    python3 tests/check_production_settings.py [--threads 30]

The certificate ran at T = 15,000.  The small-T validations in tests/validate_cert.py exercise the same code, but the
LIST/branch boundary sits elsewhere at a small table; this check removes that difference.  Below 10^25 the census is
proved (14 terms, the first repository), so the expected output of each half is known exactly:
  capped, T = 15,000, LASTCALC: every term with second-largest prime <= 15,000 (13 of the 14);
  pair tail, T = 15,000, q <= 20,000: exactly the terms with second-largest prime > 15,000 and largest <= 20,000 (one:
      2.088e24 = ... 15391 * 17291), and nothing else.
Both engines are built from this repository's sources, with the certificate's other parameters.  About 15 minutes.
"""
import argparse, json, os, subprocess, sys, tempfile, time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from sympy import factorint                                     # noqa: E402
from cn2_build import compile_cmd                               # noqa: E402
from cn2_order2_defs import verify_candidate, load_reference   # noqa: E402

X, T, QMAX = 10 ** 25, 15000, 20000


def survivors(path):
    jobs = []
    for l in open(path):
        t = l.split()
        if len(t) == 2 and t[0] in 'HN':
            jobs.append((int(t[1]), 1, 'howe'))
        elif len(t) == 5 and t[0] == 'S':
            jobs.append((int(t[1]), int(t[4]), 'howe'))
    return jobs


def verified(jobs):
    with Pool(8) as pl:
        return {n for n in pl.imap_unordered(verify_candidate, jobs, chunksize=512) if n is not None and n < X}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--threads', type=int, default=30)
    a = ap.parse_args()
    split = json.load(open(os.path.join(HERE, 'results', 'cert_1e28', 'split.json')))
    terms = [n for n in load_reference(os.path.join(HERE, 'results', 'b175531.txt')) if n < X]
    assert len(terms) == 14
    sec = {n: sorted(factorint(n)) for n in terms}
    fails = []
    with tempfile.TemporaryDirectory() as tmp:
        exe = {}
        for name in ('cn2xh', 'cn2pair'):
            exe[name] = os.path.join(tmp, name)
            subprocess.run(compile_cmd(os.path.join(HERE, 'cn2x', name + '.c'), exe[name]), check=True)
        env = dict(os.environ, CN2X_MODE='howe')

        need = {n for n in terms if sec[n][-2] <= T}
        d = os.path.join(tmp, 'capped'); os.mkdir(d)
        t0 = time.time()
        out = subprocess.run([exe['cn2xh'], str(X), d, str(a.threads), '0', str(split['ratio']), str(split['sieve']),
                              str(split['mshard']), str(split['batch']), str(split['donate_q']), '0', '0'],
                             env=dict(env, CN2X_TCAP=str(T), CN2X_LASTCALC='1'), capture_output=True, text=True).stdout
        rec = json.loads([l for l in out.splitlines() if l.startswith('{')][-1])
        got = verified(survivors(os.path.join(d, 'out.txt')))
        ok = rec['complete'] and rec['lastcalc'] == 1 and need <= got and got <= set(terms)
        print(f"[{'PASS' if ok else 'FAIL'}] capped, X = 1e25, T = {T}, LASTCALC: guaranteed {len(need)}, "
              f"found {len(got)} of 14, missing {sorted(need - got)}, outside the census {sorted(got - set(terms))}  "
              f"({time.time() - t0:.0f} s)", flush=True)
        if not ok:
            fails.append('capped')

        need = {n for n in terms if sec[n][-2] > T and sec[n][-1] <= QMAX}
        d = os.path.join(tmp, 'pair'); os.mkdir(d)
        t0 = time.time()
        out = subprocess.run([exe['cn2pair'], str(X), str(T), str(QMAX), d, str(a.threads), '0', '64', '0', '0'],
                             env=env, capture_output=True, text=True).stdout
        rec = json.loads([l for l in out.splitlines() if l.startswith('{')][-1])
        got = verified(survivors(os.path.join(d, 'out.txt')))
        ok = rec['complete'] and got == need
        print(f"[{'PASS' if ok else 'FAIL'}] pair tail, X = 1e25, T = {T}, q <= {QMAX}: required {len(need)} "
              f"{[f'{n:.4e}' for n in sorted(need)]}, found {len(got)}  ({time.time() - t0:.0f} s)", flush=True)
        if not ok:
            fails.append('pair')
    print('ALL PASS' if not fails else f"FAILED: {fails}")
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
