#!/usr/bin/env python3
"""cn2_make_results.py -- build results/b175531.txt and results/factorisations.{txt,json} from the certificate's records.

    python3 cn2_make_results.py [--check]

The terms are the union of
  * the capped half: results/cert_1e28/result.json (written by `cn2_split_campaign.py verify` on the capped run), and
  * the tail: every survivor in results/cert_1e28/pair/out.txt ("N n") and results/cert_1e28/tail/out.txt ("P/Q q a"),
    re-verified here from scratch,
and every term is checked once more against the definition before anything is written.  With --check nothing is written;
the script only confirms that the files on disk match what it would write.
"""
import argparse, json, math, os, sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sympy import factorint                                   # noqa: E402
from cn2_order2_defs import verify_candidate, is_order2, is_rigid  # noqa: E402
from cn2_split_campaign import _verify_tail                   # noqa: E402

X = 10 ** 28
CERT = os.path.join(HERE, 'results', 'cert_1e28')
RES = os.path.join(HERE, 'results')


def tail_terms():
    jobs_n, jobs_qa = [], []
    for l in open(os.path.join(CERT, 'pair', 'out.txt')):
        t = l.split()
        if len(t) == 2 and t[0] == 'N':
            jobs_n.append((int(t[1]), 1, 'howe'))
    for l in open(os.path.join(CERT, 'tail', 'out.txt')):
        t = l.split()
        if len(t) == 3 and t[0] in 'PQ':
            jobs_qa.append((t[0], int(t[1]), int(t[2]), X, 'howe'))
    with Pool(8) as pl:
        got = {n for n in pl.imap_unordered(verify_candidate, jobs_n, chunksize=256) if n is not None and n < X}
        got |= {n for n in pl.imap_unordered(_verify_tail, jobs_qa, chunksize=64) if n is not None and n < X}
    return got, len(jobs_n), len(jobs_qa)


def lcm(a, b):
    return a // math.gcd(a, b) * b


def describe(i, n):
    f = factorint(n)
    assert all(e == 1 for e in f.values()) and is_order2(n, 'howe', f), n
    ps = sorted(f)
    L = 1
    for p in ps:
        L = lcm(L, p * p - 1)
    assert (n - 1) % L == 0, n                                # rigid: L | n - 1
    Lf = factorint(L)
    return dict(n=i, a=str(n), log10=round(math.log10(n), 3), primes=ps, L=str(L),
                L_factored={str(p): e for p, e in sorted(Lf.items())}, log2_L=round(math.log2(L), 2),
                slack_bits=round(math.log2((n - 1) / L), 2), rigid=is_rigid(n, f))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    capped = {int(v) for v in json.load(open(os.path.join(CERT, 'result.json')))['found']}
    tail, nn, nqa = tail_terms()
    terms = sorted(capped | tail)
    print(f"capped half: {len(capped)} terms; tail: {len(tail)} terms from {nn:,} pair and {nqa:,} (q, a) survivors; "
          f"union: {len(terms)}")
    rows = [describe(i, n) for i, n in enumerate(terms, 1)]
    bfile = ''.join(f"{r['n']} {r['a']}\n" for r in rows)
    txt = ["A175531 below 10^28 (complete, certified 7 Oct 2026): Carmichael numbers of order 2 in Howe's sense.",
           "For each term a: its prime factors; L = lcm(p^2 - 1 : p | a), which divides a - 1 (all terms are rigid);",
           "and the slack log2((a - 1)/L).  The same data, machine-readable: factorisations.json", ""]
    for r in rows:
        Ls = ' * '.join(f"{p}^{e}" if e > 1 else p for p, e in r['L_factored'].items())
        txt += [f"a({r['n']}) = {r['a']}   (10^{r['log10']:.2f})",
                f"    primes ({len(r['primes'])}):  " + '  '.join(map(str, r['primes'])),
                f"    L = {r['L']} = {Ls}",
                f"    log2 L = {r['log2_L']:.2f}    slack log2((a-1)/L) = {r['slack_bits']:.2f}", ""]
    js = dict(sequence='A175531', description='Carmichael numbers of order 2 (Howe): all terms below 10^28',
              L_definition='L = lcm(p^2 - 1 : p | a); L divides a - 1 for every term',
              certificate='docs/CN2_Certificate_1e28.md', note='integers are strings to avoid precision loss in JSON readers',
              terms=rows)
    out = {'b175531.txt': bfile, 'factorisations.txt': '\n'.join(txt), 'factorisations.json': json.dumps(js, indent=1) + '\n'}
    if a.check:
        bad = [k for k, v in out.items() if not os.path.exists(os.path.join(RES, k)) or open(os.path.join(RES, k)).read() != v]
        print("results files match" if not bad else f"MISMATCH: {bad}")
        sys.exit(1 if bad else 0)
    for k, v in out.items():
        open(os.path.join(RES, k), 'w').write(v)
        print(f"wrote results/{k}")


if __name__ == '__main__':
    main()
