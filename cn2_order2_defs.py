"""The three order-2 definitions shared by the split driver and its verifiers.

    rigid  n == 1        (mod p^2-1)       for every p | n    (the programme's CN2)
    howe   n == 1 or p   (mod p^2-1)                          (Howe 2000; OEIS A175531)
    cheb   n == 1 or p   (mod (p^2-1)/2)                      (Carmichael & Chebyshev; OEIS A299799)
All three: odd, squarefree, composite (>= 3 prime factors is forced for all of them; checked anyway).
"""
import sys

from sympy import factorint

DD = dict(rigid=1, howe=1, cheb=2)          # plus 'testD' (validation only): modulus (p^2-1)/D, D | 24


def _dd(mode):
    return int(mode[4:]) if mode.startswith('test') else DD[mode]


def residues_ok(n, p, mode):
    g = (p * p - 1) // _dd(mode)
    r = n % g
    return r == 1 % g if mode == 'rigid' else r in (1 % g, p % g)


def is_order2(n, mode, fac=None):
    if n < 3 or n % 2 == 0 or n % 3 == 0:
        return False
    f = fac or factorint(n)
    if len(f) < 3 or any(e > 1 for e in f.values()):
        return False
    return all(residues_ok(n, p, mode) for p in f)


def is_rigid(n, fac=None):
    return is_order2(n, 'rigid', fac)


def verify_candidate(args):
    """(n_or_0, m, s, mode): n = m*s; the primes of s are factored, those of m too (small)."""
    m, s, mode = args
    n = m * s
    fm, fs = factorint(m), factorint(s)
    if set(fm) & set(fs):
        return None
    f = dict(fm); f.update(fs)
    big = [p for p in f if p >= 1 << 64]                # BPSW (sympy) is proven only below 2^64; every prime of an
    if big:                                             # order-2 n is < 2^43, so a factor this large is either a true
        print(f"UNPROVEN-PRIMALITY m={m} s={s} factors={big}", file=sys.stderr, flush=True)   # prime (n is not
    return n if is_order2(n, mode, f) else None         # order-2) or a pseudoprime: certification needs none of these


def load_reference(path):
    """known terms from a b-file ('index value' lines) or a JSON file with a 'found' list"""
    if path.endswith('.json'):
        import json
        return [int(v) for v in json.load(open(path))['found']]
    out = []
    for l in open(path):
        t = l.split()
        if len(t) == 2 and not l.startswith('#') and t[1].isdigit():
            out.append(int(t[1]))
    return out
