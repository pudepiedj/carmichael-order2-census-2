#!/usr/bin/env python3
"""cn2_pair_tail.py -- the two-large-prime tail of a certified order-2 census: prototype, unit test and cost model.

A complete capped cn2xh run with CN2X_LASTCALC=1 and table bound T finds every order-2 n < X whose second-largest
prime is <= T (CN2_LastCalc_Coverage_Proof.md).  What is left is
    R(X, T) = { n < X order-2 : p_{k-1} > T },   i.e. n has at least two primes above T.
For n in R write q = p_k (largest) and p = p_{k-1}, so T < p < q.  For each candidate largest prime q choose ONE of:
  (Q) the (q, a) family of cn2tail: n = q * m with m == r_q q (mod g(q)) -- covers EVERY n divisible by q;
  (P) the pair family: for every prime p in (T, q) and every compatible pair of types, n lies in one residue class
      modulo M = p q lcm(g(p), g(q)); enumerate it -- covers every n whose two largest primes are p < q.
Either choice per q covers every n in R whose largest prime is q, so any per-q mixture is complete for R.

    python3 cn2_pair_tail.py unittest                  # each known n with p_{k-1} > T is found by its own pair
    python3 cn2_pair_tail.py count --X 1e20 --T 1000 --Q 10000      # exact candidate count vs the model
    python3 cn2_pair_tail.py model --X 1e30 --T 30000              # cost of the tail, per-q cheaper family
"""
import argparse, bisect, json, math, os, random, sys
from decimal import Decimal

from sympy import factorint, primerange

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cn2_order2_defs import verify_candidate, DD, load_reference  # noqa: E402

REF = os.path.join(HERE, 'results', 'b175531.txt')     # the certified terms


def crt(r1, m1, r2, m2):
    """x == r1 (mod m1), x == r2 (mod m2) -> (x, lcm) or None if incompatible"""
    g = math.gcd(m1, m2)
    if (r2 - r1) % g:
        return None
    l = m1 // g * m2
    t = ((r2 - r1) // g * pow(m1 // g, -1, m2 // g)) % (m2 // g)
    return (r1 + m1 * t) % l, l


def types(p, mode):
    g = (p * p - 1) // DD[mode]
    return [1 % g] if mode == 'rigid' else sorted({1 % g, p % g})


def pair_classes(p, q, mode):
    """the residue classes (N0, M) of n with p q | n and the type conditions at p and q; [] if p | g(q)"""
    D = DD[mode]
    gp, gq = (p * p - 1) // D, (q * q - 1) // D
    if gq % p == 0:                       # p | g(q): impossible for an order-2 n (F5)
        return []
    out = []
    for rp in types(p, mode):
        for rq in types(q, mode):
            c = crt(rp, gp, rq, gq)
            if c is None:
                continue
            c = crt(c[0], c[1], 0, p * q)
            if c is not None:
                out.append(c)
    return out


def pair_search(p, q, X, mode='howe'):
    """every order-2 n < X whose two largest primes are p < q (verified)"""
    D = DD[mode]
    found = []
    for N0, M in pair_classes(p, q, mode):
        lo = q * ((q * q - 1) // D)       # F2 at q: n = q (n/q) > q g(q)
        n = N0 if N0 > lo else N0 + ((lo - N0) // M + 1) * M
        while n < X:
            if n & 1 and n % 3 and pow(2, n - 1, n) == 1:
                v = verify_candidate((n, 1, mode))
                if v is not None:
                    f = sorted(factorint(v))
                    if f[-1] == q and f[-2] == p:
                        found.append(v)
            n += M
    return found


def pair_count(p, q, X, mode='howe'):
    D = DD[mode]
    lo = q * ((q * q - 1) // D)
    tot = 0
    for N0, M in pair_classes(p, q, mode):
        first = N0 if N0 > lo else N0 + ((lo - N0) // M + 1) * M
        if first < X:
            tot += (X - 1 - first) // M + 1
    return tot


def qa_count(q, X, mode='howe'):
    """(q, a) candidates for q, both type families, as cn2tail counts them: a <= (X-1-q^2) / (q g(q))"""
    D = DD[mode]
    g = (q * q - 1) // D
    A = max(0, (X - 1 - q * q) // (q * g))
    return A * (1 if mode == 'rigid' else 2)


def cmd_unittest(a):
    """every known order-2 number lies in one of the classes of its own top pair (p, q), above the F2 bound;
    where the enumeration is cheap (X = n + 1) it is also run in full and must return n"""
    ref = sorted(load_reference(REF))
    member = full = 0
    for n in ref:
        f = sorted(factorint(n))
        p, q = f[-2], f[-1]
        lo = q * ((q * q - 1) // DD['howe'])
        if not (n > lo and any(n % M == N0 for N0, M in pair_classes(p, q, 'howe'))):
            sys.exit(f"MISSED by the class construction: {n}  p={p} q={q}")
        member += 1
        if pair_count(p, q, n + 1) <= 200_000:
            if n not in pair_search(p, q, n + 1):
                sys.exit(f"MISSED by the enumeration: {n}  p={p} q={q}")
            full += 1
    print(f"{member}/{len(ref)} known numbers lie in a class of their own (p, q); "
          f"{full} also found by the full enumeration to X = n + 1")


def cmd_count(a):
    X, T, Q = int(Decimal(a.X)), a.T, a.Q
    ps = list(primerange(T + 1, Q + 1))
    exact = sum(pair_count(p, q, X) for i, q in enumerate(ps) for p in ps[:i])
    setups = len(ps) * (len(ps) - 1) // 2
    print(f"X={a.X} T={T} Q={Q}: {len(ps)} primes, {setups:,} pair setups, exact pair candidates {exact:,}")
    print(f"model: {model_pairs(X, T, Q, ps):.4g}")


def mean_G(T, Q, mode='howe', samples=20000, seed=1):
    """E[ sum over compatible type pairs of p q / M ] for random primes T < p < q <= Q (the per-pair density factor)"""
    ps = list(primerange(T + 1, Q + 1))
    rnd = random.Random(seed)
    s = 0.0
    for _ in range(samples):
        p, q = sorted(rnd.sample(ps, 2))
        D = DD[mode]
        base = (p * p - 1) // D * ((q * q - 1) // D)
        s += sum(p * q * base / M for _, M in pair_classes(p, q, mode))   # = sum over compatible types of gcd(g(p), g(q))
    return s / samples   # candidates per pair ~ X * Gbar / (p q g(p) g(q))


def model_pairs(X, T, Q, ps=None, mode='howe'):
    """sum over T < p < q <= Q of X * Gbar / (p q g(p) g(q)), Gbar from sampling"""
    D = DD[mode]
    ps = ps or list(primerange(T + 1, Q + 1))
    G = mean_G(T, Q, mode)
    inv = [D / (p * (p * p - 1)) for p in ps]              # 1/(p g(p))
    pref, tot = 0.0, 0.0
    for i, w in enumerate(inv):
        tot += w * pref
        pref += w
    return X * G * tot


def cmd_model(a):
    """per-q min of the (q,a) cost and the pair cost (setups + candidates), for q up to (D X)^(1/3)"""
    X, T, mode = int(Decimal(a.X)), a.T, 'howe'
    D = DD[mode]
    qmax = int(round((D * X) ** (1 / 3))) + 1
    G = mean_G(T, min(qmax, 10 ** 8), mode, samples=a.samples)
    print(f"X={a.X} T={T}: mean pair density factor Gbar = {G:.1f} (sampled)")
    # walk q upwards: pair cost(q) = (#primes in (T,q)) setups + X*G/(q g(q)) * S(q),  S(q) = sum_{T<p<q} 1/(p g(p))
    # (q,a) cost(q) = 2 (X-1-q^2)/(q g(q)).  Once the (q,a) cost is the cheaper it stays so (it falls like q^-3).
    S, npr, tot_pair, tot_setup, tot_qa, switch = 0.0, 0, 0.0, 0, 0.0, None
    seg = 10 ** 7
    lo = T + 1
    done = False
    while lo <= qmax and not done:
        for q in primerange(lo, min(lo + seg, qmax + 1)):
            g = (q * q - 1) // D
            c_pair = npr + X * G * S / (q * g)
            c_qa = 2 * max(0, (X - 1 - q * q) // (q * g))
            if switch is None and c_qa <= c_pair:
                switch = q
            if switch is None:
                tot_pair += X * G * S / (q * g); tot_setup += npr
            else:
                tot_qa += c_qa
            S += D / (q * (q * q - 1)); npr += 1
        lo += seg
        if switch is not None:          # the (q, a) tail above the switch has a closed form: X / (2 Q^2 ln Q) per family
            tot_qa = 2 * X / (2 * switch ** 2 * math.log(switch))
            done = True
    rate = 1.1e8                        # tests per second on 30 threads (cn2tail --bench: 265 ns/test single thread)
    tot = tot_pair + tot_setup + tot_qa
    print(f"  switch to (q,a) at q* = {switch:,}" if switch else "  pairs all the way")
    print(f"  pair candidates {tot_pair:.3g}, pair setups {tot_setup:.3g}, (q,a) candidates above q* {tot_qa:.3g}")
    print(f"  total {tot:.3g} tests = {tot / rate / 3600:.1f} h at {rate:.1e}/s   "
          f"(plain one-prime tail at R0 = T: {X / (T * T * math.log(T)):.3g} tests)")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    sub.add_parser('unittest')
    c = sub.add_parser('count'); c.add_argument('--X', default='1e20'); c.add_argument('--T', type=int, default=1000)
    c.add_argument('--Q', type=int, default=10000)
    m = sub.add_parser('model'); m.add_argument('--X', default='1e30'); m.add_argument('--T', type=int, default=30000)
    m.add_argument('--samples', type=int, default=20000)
    a = ap.parse_args()
    dict(unittest=cmd_unittest, count=cmd_count, model=cmd_model)[a.cmd](a)


if __name__ == '__main__':
    main()
