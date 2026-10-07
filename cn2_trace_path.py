#!/usr/bin/env python3
"""cn2_trace_path.py -- replay the cn2xh search path of one order-2 number, node by node.

    python3 cn2_trace_path.py N [--X 1e28] [--R0 30000] [--mode howe] [--ratio 1.0] [--show 12]

Walks the primes of N in increasing order, exactly as cn2xh's DFS would descend, and prints for
each node the state the engine sees:
    m          product of the primes chosen so far
    L, c       n == c (mod L), L = lcm g(p) over the chosen primes (g(p) = (p^2-1)/D)
    rem        (X-1)/m, the room left for the cofactor s = n/m
    prog       rem/L + 1, the length of the progression s == c m^{-1} (mod L), s <= rem
    nch        table primes t in (pmax, tmax] still eligible as children
    bounds     b_last = isqrt(Dm-1), b_second = min(Dm-1, isqrt(rem)), b_early = icbrt(rem)
The engine switches the node to LIST mode when prog <= ratio * nch: from there it never picks
another prime, it walks the progression and Fermat-tests every term.  The trace marks that node,
lists the progression's survivors, and then (for illustration only) continues the descent past it
to show how the progression collapses to one term -- the point where the largest prime is forced.
Mirrors node_head/child_step/list_mode in cn2x/cn2xh.c; the type of each prime (1 or p) is read
off N itself.
"""
import argparse
from math import gcd, isqrt, prod
from bisect import bisect_right
from sympy import factorint, primerange


def icbrt(n):
    r = round(n ** (1 / 3))
    while r ** 3 > n:
        r -= 1
    while (r + 1) ** 3 <= n:
        r += 1
    return r


def lcm(a, b):
    return a // gcd(a, b) * b


def fermat2(s, nm1):
    return pow(2, nm1, s) == 1


def node_state(m, L, k, pmax, X, T, D, ratio, table):
    """what cn2xh's node_head computes at node (m, L, k, pmax): rem, prog, nch, tmax, eb, and the LIST switch"""
    rem = (X - 1) // m
    b_last = isqrt(D * m - 1) if k >= 2 else 0
    b_second = min(D * m - 1, isqrt(rem)) if k >= 1 else 0
    b_early = icbrt(rem)
    tmax = min(max(b_last, b_second, b_early), T)
    nch = max(0, bisect_right(table, tmax) - bisect_right(table, pmax)) if rem > pmax else 0   # node_head: nothing fits
    prog = rem // L + 1
    return dict(rem=rem, prog=prog, nch=nch, tmax=tmax, eb=max(b_second, b_early),
                is_list=L > 1 and prog <= ratio * nch)


def switch_point(ps, X, T, D=1, ratio=1.0, table=None):
    """(k, p) at the node where the search of n = prod(ps) switches to LIST: p is the largest prime the
    DFS had to pick from the table there.  A capped look with cap T finds n if p <= T.  None if the DFS
    reaches n itself without switching (n is then emitted as a hit on its own path)."""
    table = table or list(primerange(3, T + 1))
    m, L, pmax = 1, 1, 1
    for k in range(len(ps)):
        if node_state(m, L, k, pmax, X, T, D, ratio, table)['is_list']:
            return k, pmax
        t = ps[k]
        m, L, pmax = m * t, lcm(L, (t * t - 1) // D), t
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('N', type=int)
    ap.add_argument('--X', type=float, default=1e28)
    ap.add_argument('--R0', type=int, default=30000, help='table cap T of the capped phase')
    ap.add_argument('--mode', default='howe', choices=['rigid', 'howe', 'cheb'])
    ap.add_argument('--ratio', type=float, default=1.0)
    ap.add_argument('--show', type=int, default=12, help='progression survivors to print at the LIST node')
    a = ap.parse_args()
    N, X, T, D = a.N, int(a.X), a.R0, 2 if a.mode == 'cheb' else 1
    ps = sorted(factorint(N))
    assert prod(ps) == N, 'N must be squarefree'
    table = list(primerange(3, T + 1))

    def g(p):
        return (p * p - 1) // D

    for p in ps:
        assert N % g(p) in (1, p % g(p)), f'{p}: N is not order-2 at this prime'
    print(f"N = {N}  (10^{len(str(N)) - 1}+),  k = {len(ps)},  X = {X:.3g},  T = R0 = {T:,},  mode {a.mode}")
    print(f"primes: {', '.join(map(str, ps))}\n")
    hdr = f"{'k':>2} {'+prime':>7} {'m':>30} {'L':>22} {'rem':>12} {'prog':>10} {'nch':>6} {'tmax':>6}  switch"
    print(hdr)
    print('-' * len(hdr))

    m, L, k, pmax = 1, 1, 0, 1
    listed = None
    for i in range(len(ps) + 1):
        c = N % L
        st = node_state(m, L, k, pmax, X, T, D, a.ratio, table)
        rem, prog, nch, tmax, is_list = st['rem'], st['prog'], st['nch'], st['tmax'], st['is_list']
        added = '' if k == 0 else str(ps[k - 1])
        tag = ''
        if listed is None and is_list:
            tag = '<== LIST: engine stops choosing primes here'
            listed = (k, m, L, c, rem, pmax)
        elif listed is not None:
            tag = '(beyond the engine: illustration)'
        elif nch == 0:
            tag = 'no table child fits'
        print(f"{k:>2} {added:>7} {m:>30} {L:>22} {rem:>12.3e} {prog:>10.3e} {nch:>6} {tmax:>6}  {tag}")
        if i == len(ps):
            break
        t = ps[i]
        if L % t == 0 or gcd(g(t), m % g(t)) != 1:
            print(f"   child {t} is rejected by child_step (t | L or gcd(g(t), m) > 1)")
        eb = st['eb']
        if listed is None and k >= 1 and t > eb:
            print(f"   child {t} > eb = {eb}: handled as a last-prime test, not a node")
        m, L, k, pmax = m * t, lcm(L, g(t)), k + 1, t
        if L >= X and listed is None:
            print(f"   L >= X after {t}: determined node, n resolved on the spot")

    if listed:
        k, m, L, c, rem, pmax = listed
        s0 = c * pow(m, -1, L) % L
        terms = range(s0, rem + 1, L)
        surv = []
        n_terms = n_cand = 0
        for s in terms:
            n_terms += 1
            if s <= pmax or s % 2 == 0 or s % 3 == 0:
                continue
            if any(s % q == 0 for q in primerange(5, min(pmax, 509) + 1)):
                continue
            n_cand += 1
            if fermat2(s, m * s - 1):
                surv.append(s)
        print(f"\nLIST node at k = {k}: s == {s0} (mod {L}), s <= {rem}")
        print(f"  {n_terms:,} progression terms -> {n_cand:,} after parity/3/sieve -> {len(surv)} Fermat survivors")
        for s in surv[:a.show]:
            f = factorint(s)
            hit = N == m * s
            ok = all((m * s) % g(q) in (1, q % g(q)) for q in f) and all(e == 1 for e in f.values())
            print(f"    s = \033[33m{s:>16}\033[0m  = {' * '.join(map(str, sorted(f))):<28} {'order-2' if ok else '':<8}{'  <== N' if hit else ''}")
        if len(surv) > a.show:
            print(f"    ... {len(surv) - a.show} more")


if __name__ == '__main__':
    main()
