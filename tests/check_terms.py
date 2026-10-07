#!/usr/bin/env python3
"""tests/check_terms.py -- a standalone check of the 45 terms, independent of everything else in this repository.

    python3 tests/check_terms.py [results/b175531.txt]

Standard library only.  It reads nothing but the b-file, ignores the repository's factor lists and factorisation code,
and for each term n checks Howe's definition from scratch:
  * the b-file is indexed 1, 2, 3, ... with strictly increasing terms, all below 10^28;
  * n is odd and composite, and its factorisation (recomputed here by trial division) is squarefree;
  * for every prime p | n:  n mod (p^2 - 1) is 1 or p.
The factorisation is exact: it divides by d = 2, 3, 5, 7, 9, ... (every odd d after 2), so each factor removed is prime
(all smaller factors are already gone), and it stops when d^2 exceeds the remaining cofactor, which is then prime.
It also reports which terms are rigid (n mod (p^2 - 1) = 1 for every p).  The PARI/GP equivalent is tests/check_terms.gp.
"""
import sys


def factor(n):
    f, d = [], 2
    while n > 1:
        if d * d > n:
            f.append(n)                     # remaining cofactor is prime (no divisor <= its square root)
            break
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e:
            f.append(d if e == 1 else (d, e))
        d = 3 if d == 2 else d + 2
    return f


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'results/b175531.txt'
    rows = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    idx = [int(r[0]) for r in rows]
    terms = [int(r[1]) for r in rows]
    ok = idx == list(range(1, len(terms) + 1)) and all(a < b for a, b in zip(terms, terms[1:])) and terms[-1] < 10 ** 28
    print(f"{len(terms)} terms; indexed 1..{len(terms)}, strictly increasing, all below 10^28: {ok}")
    rigid = 0
    for i, n in zip(idx, terms):
        f = factor(n)
        sqfree = all(isinstance(p, int) for p in f)
        good = n % 2 == 1 and len(f) >= 2 and sqfree and all(n % (p * p - 1) in (1, p) for p in f)
        r = sqfree and all(n % (p * p - 1) == 1 for p in f)
        rigid += r
        ok &= good
        print(f"a({i}) = {n}  {'OK' if good else 'FAIL'}{'  rigid' if r else ''}  = {' * '.join(map(str, f))}")
    print(f"\n{'ALL 45 TERMS SATISFY THE DEFINITION' if ok and len(terms) == 45 else 'FAILED'}; rigid: {rigid} of {len(terms)}")
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
