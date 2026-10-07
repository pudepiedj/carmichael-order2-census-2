#!/usr/bin/env python3
"""cn2_tail_verify.py -- verify the survivors of the two-large-prime tail (cn2pair "N n" lines, cn2tail "P/Q p a" lines).

    python3 cn2_tail_verify.py --X 1e25 --mode howe --pair RUNDIR [--qa RUNDIR] [--T 10000] [--reference results/b175531.txt]

Every survivor is factored and checked against the definition (cn2_order2_defs.verify_candidate, which also flags any
factor >= 2^64 as UNPROVEN-PRIMALITY on stderr).  With --reference and --T, it also checks the certificate logic on a
range where the answer is known: every reference number below X whose second-largest prime exceeds T must be found,
and nothing outside the reference may be.
"""
import argparse, json, os, sys
from decimal import Decimal
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cn2_order2_defs import verify_candidate, load_reference   # noqa: E402
from cn2_split_campaign import _verify_tail            # noqa: E402  the (p, a) verifier the proved census used


def complete(d):
    """the last session line in stdout.txt / sessions.jsonl, if any, says complete"""
    for name in ('sessions.jsonl', 'stdout.txt'):
        p = os.path.join(d, name)
        if os.path.exists(p):
            js = [l for l in open(p) if l.startswith('{')]
            if js:
                return json.loads(js[-1]).get('complete')
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--X', required=True); ap.add_argument('--mode', default='howe')
    ap.add_argument('--pair', required=True); ap.add_argument('--qa')
    ap.add_argument('--T', type=int); ap.add_argument('--reference'); ap.add_argument('--procs', type=int, default=8)
    a = ap.parse_args()
    X = int(Decimal(a.X))
    jobs, tail = [], []
    for l in open(os.path.join(a.pair, 'out.txt')):      # a run in progress may end in a half-written line: skip it
        t = l.split()
        if len(t) == 2 and t[0] == 'N' and t[1].isdigit() and l.endswith('\n'):
            jobs.append((int(t[1]), 1, a.mode))
    if a.qa and os.path.exists(os.path.join(a.qa, 'out.txt')):
        for l in open(os.path.join(a.qa, 'out.txt')):
            t = l.split()
            if len(t) == 3 and t[0] in 'PQ' and t[1].isdigit() and t[2].isdigit() and l.endswith('\n'):
                tail.append((t[0], int(t[1]), int(t[2]), X, a.mode))
    with Pool(a.procs) as pl:
        got = {n for n in pl.imap_unordered(verify_candidate, jobs, chunksize=256) if n is not None and n < X}
        got |= {n for n in pl.imap_unordered(_verify_tail, tail, chunksize=256) if n is not None and n < X}
    print(f"pair survivors {len(jobs):,}, (q,a) survivors {len(tail):,}; verified order-2 numbers below X: {len(got)}")
    for n in sorted(got):
        print(f"  {n}  ({n:.4e})")
    if a.reference:
        from sympy import factorint
        ref = {n for n in load_reference(a.reference) if n < X}
        need = {n for n in ref if a.T is not None and sorted(factorint(n))[-2] > a.T}
        miss, extra = need - got, got - ref
        print(f"required (reference, second-largest prime > T): {len(need)}; missing: {len(miss)}; outside the reference: {len(extra)}")
        for n in sorted(miss):
            print(f"  MISSING {n}")
        for n in sorted(extra):
            print(f"  NOT IN REFERENCE {n}")
        done = [complete(a.pair)] + ([complete(a.qa)] if a.qa else [])
        if not all(done):
            print(f"RUN IN PROGRESS (complete flags {done}): the check above is interim, PASS/FAIL is only decided at the end")
        else:
            print("PASS" if not miss and not extra else "FAIL")


if __name__ == '__main__':
    main()
