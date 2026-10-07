/*
 * cn2xh.c -- cn2xc.c generalised to three definitions of an order-2 Carmichael number (CN2X_MODE):
 *   rigid (default)  n == 1 (mod p^2-1) for every p | n          -- identical search and counters to cn2xc
 *   howe             n == 1 or p (mod p^2-1)                      -- Howe (2000); OEIS A175531
 *   cheb             n == 1 or p (mod (p^2-1)/2)                  -- Carmichael & Chebyshev; OEIS A299799
 * Each prime chooses its type ("1" or "p") and a node carries n == c (mod L) built by CRT, instead
 * of n == 1 (mod L).  Incompatible choices die at the merge.  When a merge makes L >= X the number
 * is determined (at most one n < X) and is resolved on the spot.  Role bounds: with g(p) = (p^2-1)/D
 * (D = 1, or 2 for cheb), n/p > g(p) for either type (n/p == p or 1 (mod g(p)), n/p not in {1, p}),
 * so the last prime t has t^2 < D m, the second-to-last t < D m, and every prime is <= (D X)^(1/3).
 *
 * Original header of cn2xc.c follows.
 *
 * cn2xc.c -- checkpointed, restartable version of cn2x.c (same search, same counters).
 *
 * The search tree is cut into deterministic SHARDS: nodes with m < MSHARD are expanded by a
 * single generator; the remaining children of each expanded node (m*t >= MSHARD) are cut, in index
 * order, into consecutive runs of 1, 2, 4, ... BATCH children -- the first (heaviest) children get
 * shards of their own.  No per-child work is needed to generate shards.  Shard ids
 * are assigned in generation order, so the same X/parameters always give the same shards.
 * A finished shard is appended to done.log with its counters; a restart replays the generator
 * and skips done shards.  At a deadline (or RUNDIR/DRAIN) the run DRAINS: no new shards are taken,
 * in-flight shards finish and are recorded, then it exits -- no work is lost.  A hard stop
 * (RUNDIR/STOP, SIGINT/SIGTERM) abandons in-flight shards only.  Final totals = generator counters (recomputed each session) + sum over done shards.
 *
 * Original header of cn2x.c follows.
 *
 * cn2x.c -- exhaustive enumeration of rigid CN2s below X (compiled port of cn2_exhaustive.py).
 *
 * Same search, same counters, same switch rule as the Python reference, so the two
 * can be matched exactly.  Differences are purely mechanical:
 *   - 128-bit integers throughout (X up to ~2^100; L < X by the budget prune);
 *   - base-2 Fermat pre-filter by Montgomery arithmetic (one limb below 2^63, two above);
 *   - LIST survivors are NOT factored here: they are written out as (m, pmax, k, s) and
 *     verified by the Python driver (sympy), which is also where hits are re-verified;
 *   - workers share a LIFO task stack; nodes with m < MSPLIT are expanded into tasks, and whenever a
 *     thread is idle, a busy thread at a node with X/m > DONATE_Q donates its remaining children.
 *     Counters are exact under any split or donation pattern.
 *
 * Output (text, one record per line):
 *   H <n>                      hit found on its own prime path (L | n-1 by construction)
 *   S <m> <pmax> <k> <s>       LIST candidate that passed Fermat; needs factoring
 * Final line on stdout: JSON counters.
 *
 * Build:  clang -O3 -mcpu=native -o cn2x cn2x.c -lpthread -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <stdatomic.h>
#include <signal.h>
#include <unistd.h>

typedef unsigned __int128 u128;
typedef __int128 i128;
typedef uint64_t u64;
typedef uint32_t u32;

static u64 T;
static u32 *PR;
static size_t NPR;
static double RATIO = 1.0;
static int FERMAT = 1;
static int LASTCALC = 0;        /* CN2X_LASTCALC=1: in a capped search, also compute last primes above the cap */
static int SIEVE = 1;
static int MODE = 0;            /* 0 rigid, 1 howe, 2 cheb */
static u64 DD = 1;              /* g(p) = (p^2-1)/DD */
static u64 SKIPMASK = 0;        /* EXPERIMENT ONLY: CN2X_SKIP_PRIMES=5,7 never uses those primes (< 64) as children: the
                                * search is then CONDITIONAL on them not dividing n -- for measuring their subtrees' cost */
static u128 NMIN = 0;           /* CN2X_NMIN: report only n >= NMIN; the search is unchanged */           /* 1: reject s with a prime factor q <= min(pmax, 47) before the Fermat test */
static const u64 SMALLP[] = {5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509};
static int NSMALL = 95;        /* set from X: the crossover is near 1e20 (below it the
                                * Fermat test is cheap enough that 95 divisions cost more than they save) */         /* how many of them to use: CN2X_SIEVEP */          /* 0 off, 1 on, 2 = count candidates but skip the test (timing only) */
static u128 DONATE_Q = (u128)10000000000ull;  /* a node with X/m > DONATE_Q donates its remaining children when a thread is idle */
static atomic_int HUNGRY = 0;
static FILE *OUT;
static pthread_mutex_t out_mu = PTHREAD_MUTEX_INITIALIZER;

/* ---------------------------------------------------------------- u128 helpers */
static u128 parse_u128(const char *s) {
    u128 v = 0;
    for (; *s; s++) if (*s >= '0' && *s <= '9') v = v * 10 + (u64)(*s - '0');
    return v;
}
static void u128_str(u128 v, char *buf) {
    char tmp[64]; int i = 0;
    if (v == 0) { strcpy(buf, "0"); return; }
    while (v) { tmp[i++] = (char)('0' + (int)(v % 10)); v /= 10; }
    int j = 0; while (i) buf[j++] = tmp[--i]; buf[j] = 0;
}
static inline u64 mod128_64(u128 a, u64 b) {
    return (a >> 64) == 0 ? (u64)a % b : (u64)(a % b);
}
static u128 isqrt128(u128 n) {
    if (n == 0) return 0;
    u128 r = (u128)sqrt((double)n);
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return r;
}
static u128 icbrt128(u128 n) {
    if (n == 0) return 0;
    u128 r = (u128)cbrt((double)n);
    while (r * r * r > n) r--;
    while ((r + 1) * (r + 1) * (r + 1) <= n) r++;
    return r;
}
static inline u64 gcd64(u64 a, u64 b) {
    if (!a) return b;
    if (!b) return a;
    int sh = __builtin_ctzll(a | b);
    a >>= __builtin_ctzll(a);
    do {
        b >>= __builtin_ctzll(b);
        if (a > b) { u64 t = a; a = b; b = t; }
        b -= a;
    } while (b);
    return a << sh;
}
static u128 modinv128(u128 a, u128 mod) {          /* gcd(a, mod) = 1, mod > 1 */
    a %= mod;
    i128 x0 = 0, x1 = 1;
    if ((mod >> 64) == 0) {
        u64 r0 = (u64)mod, r1 = (u64)a;
        while (r1) {
            u64 q = r0 / r1, r2 = r0 - q * r1; r0 = r1; r1 = r2;
            i128 x2 = x0 - (i128)q * x1; x0 = x1; x1 = x2;
        }
    } else {
        u128 r0 = mod, r1 = a;
        while (r1) {
            u128 q = r0 / r1, r2 = r0 - q * r1; r0 = r1; r1 = r2;
            i128 x2 = x0 - (i128)q * x1; x0 = x1; x1 = x2;
        }
    }
    if (x0 < 0) x0 += (i128)mod;
    return (u128)x0;
}

static u64 modinv64(u64 a, u64 mod) {             /* gcd(a, mod) = 1, mod > 1 */
    i128 x0 = 0, x1 = 1; u64 r0 = mod, r1 = a % mod;
    while (r1) { u64 q = r0 / r1, r2 = r0 - q * r1; r0 = r1; r1 = r2; i128 x2 = x0 - (i128)q * x1; x0 = x1; x1 = x2; }
    if (x0 < 0) x0 += mod;
    return (u64)x0;
}
static inline int bitlen128(u128 v) {
    u64 hi = (u64)(v >> 64);
    return hi ? 128 - __builtin_clzll(hi) : (v ? 64 - __builtin_clzll((u64)v) : 0);
}
/* a*b mod `mod` for any 1 <= mod < 2^127.  Horner over w-bit limbs of b with w = 127 - bitlen(mod), so that
 * r*2^w < 2^127 and a*limb < 2^127 and their sum fits in 128 bits.  (Until 24 Sep 2026 this used fixed
 * 32-bit limbs, valid only for mod < 2^84, and the howe/cheb modes refused X >= 2^84.) */
static u128 mulmod128(u128 a, u128 b, u128 mod) {
    a %= mod; b %= mod;
    int w = 127 - bitlen128(mod);
    u128 r = 0;
    if (w < 1) {                                         /* 2^126 <= mod < 2^127: double-and-add */
        for (int i = bitlen128(b) - 1; i >= 0; i--) {
            r += r; if (r >= mod) r -= mod;
            if ((b >> i) & 1) { r += a; if (r >= mod) r -= mod; }
        }
        return r;
    }
    int nb = bitlen128(b);
    u128 mask = (((u128)1) << w) - 1;
    for (int sh = ((nb + w - 1) / w - 1) * w; sh >= 0; sh -= w)
        r = ((r << w) + a * ((b >> sh) & mask)) % mod;
    return r;
}
static u128 mulmod_ref(u128 a, u128 b, u128 mod) {     /* reference for --selftest: bitwise double-and-add */
    a %= mod; b %= mod; u128 r = 0;
    for (int i = 127; i >= 0; i--) {
        r += r; if (r >= mod) r -= mod;
        if ((b >> i) & 1) { r += a; if (r >= mod) r -= mod; }
    }
    return r;
}

/* merge n == c (mod L) with n == r (mod gm).  Returns 0 incompatible, 1 node (*Lo, *co; L' < X),
 * 2 determined (L' >= X): *n0 = the unique n < X in the class, or 0 if there is none. */
static u128 X;
static int merge(u128 L, u128 c, u64 gm, u64 r, u128 *Lo, u128 *co, u128 *n0) {
    u64 cm = (u64)(c % gm);
    u64 gg = gcd64(gm, (u64)(L % gm));
    u64 diff = r >= cm ? r - cm : r + (gm - cm);
    if (diff % gg) return 0;
    u64 g2 = gm / gg, k = 0;
    if (g2 > 1) {
        u64 inv = modinv64((u64)((L / gg) % g2), g2);
        k = (u64)(((u128)((diff / gg) % g2) * inv) % g2);
    }
    if (g2 > 1 && L > (X - 1) / g2) {               /* L' = L g2 > X - 1 */
        *n0 = 0;
        if ((u128)k <= (X - 1 - c) / L) *n0 = c + L * (u128)k;   /* else the class has no member below X */
        return 2;
    }
    *Lo = L * g2; *co = c + L * (u128)k;
    return 1;
}

/* ---------------------------------------------------------------- Montgomery Fermat, base 2 */
static inline int topbit128(u128 e) {
    u64 hi = (u64)(e >> 64);
    return hi ? 127 - __builtin_clzll(hi) : 63 - __builtin_clzll((u64)e);
}
static inline u64 inv_mod_2_64(u64 n) {           /* n odd */
    u64 x = n;                                     /* correct to 3 bits */
    for (int i = 0; i < 5; i++) x *= 2 - n * x;
    return x;
}
static inline u128 montmul2(u128 a, u128 b, u128 N, u64 np) {
    u64 a0 = (u64)a, a1 = (u64)(a >> 64), b0 = (u64)b, b1 = (u64)(b >> 64);
    u64 n0 = (u64)N, n1 = (u64)(N >> 64);
    u64 t0, t1, t2, t3, m;
    u128 C;
    /* i = 0 */
    C = (u128)a0 * b0;              t0 = (u64)C; C >>= 64;
    C += (u128)a1 * b0;             t1 = (u64)C; C >>= 64;
    t2 = (u64)C; t3 = 0;
    m = t0 * np;
    C = (u128)t0 + (u128)m * n0;    C >>= 64;
    C += (u128)t1 + (u128)m * n1;   t0 = (u64)C; C >>= 64;
    C += (u128)t2;                  t1 = (u64)C; C >>= 64;
    t2 = t3 + (u64)C;
    /* i = 1 */
    C = (u128)t0 + (u128)a0 * b1;   t0 = (u64)C; C >>= 64;
    C += (u128)t1 + (u128)a1 * b1;  t1 = (u64)C; C >>= 64;
    C += (u128)t2;                  t2 = (u64)C; t3 = (u64)(C >> 64);
    m = t0 * np;
    C = (u128)t0 + (u128)m * n0;    C >>= 64;
    C += (u128)t1 + (u128)m * n1;   t0 = (u64)C; C >>= 64;
    C += (u128)t2;                  t1 = (u64)C; C >>= 64;
    t2 = t3 + (u64)C;
    u128 r = ((u128)t1 << 64) | t0;
    if (t2 || r >= N) r -= N;
    return r;
}
/* returns 1 iff 2^e == 1 (mod s); s odd, 1 < s < 2^126, e >= 1 */
static int fermat2(u128 s, u128 e) {
    int top = topbit128(e);
    if ((s >> 63) == 0) {
        u64 N = (u64)s, np = (u64)0 - inv_mod_2_64(N);
        u64 one = (u64)((((u128)1) << 64) % N), x = one;
        for (int i = top; i >= 0; i--) {
            u128 t = (u128)x * x;
            u64 mm = (u64)t * np;
            x = (u64)((t + (u128)mm * N) >> 64);
            if (x >= N) x -= N;
            if ((e >> i) & 1) { x += x; if (x >= N) x -= N; }
        }
        return x == one;
    } else {
        u128 N = s;
        u64 np = (u64)0 - inv_mod_2_64((u64)N);
        u128 one = ((u128)0 - N) % N, x = one;      /* 2^128 mod N */
        for (int i = top; i >= 0; i--) {
            x = montmul2(x, x, N, np);
            if ((e >> i) & 1) { x += x; if (x >= N) x -= N; }
        }
        return x == one;
    }
}

/* ---------------------------------------------------------------- per-thread context */
typedef struct {
    u64 nodes, children, last_tests, list_nodes, list_candidates, fermat_passes, factorisations, pruned_budget, sieved;
    u64 shard;
} Ctx;

static double now(void);
static atomic_ullong G_NODES = 0, G_CAND = 0, G_HITS = 0;   /* live progress only */
static double PROGRESS_SEC = 60;
static double EXPECT_NODES = 0;     /* predicted node total for X: enables a percentage and an ETA */
static atomic_int STOP = 0;     /* hard stop: abandon in-flight shards */
static atomic_int DRAIN = 0;    /* soft stop: no new shards, finish in-flight ones, then exit */
static void ctx_emit(Ctx *c, const char *line) {
    (void)c;
    pthread_mutex_lock(&out_mu);
    fputs(line, OUT);                     /* unbuffered by design: a shard's records precede its done line */
    pthread_mutex_unlock(&out_mu);
}
static void emit_hit(Ctx *c, u128 n) {
    if (n < NMIN) return;
    atomic_fetch_add(&G_HITS, 1);
    char a[64], line[128]; u128_str(n, a);
    snprintf(line, sizeof line, "H %s\n", a); ctx_emit(c, line);
}
static void emit_survivor(Ctx *c, u128 m, u64 pmax, int k, u128 s) {
    char a[64], b[64], line[200]; u128_str(m, a); u128_str(s, b);
    snprintf(line, sizeof line, "S %s %llu %d %s\n", a, (unsigned long long)pmax, k, b); ctx_emit(c, line);
}

static inline size_t upper(u64 key) {              /* first index with PR[i] > key */
    if (key >= 0xFFFFFFFFull) return NPR;
    size_t lo = 0, hi = NPR;
    while (lo < hi) { size_t mid = (lo + hi) / 2; if (PR[mid] <= key) lo = mid + 1; else hi = mid; }
    return lo;
}

/* ---------------------------------------------------------------- the search */
static int DEBUG_LIST = 0;
static void list_mode(Ctx *c, u128 m, u64 pmax, u128 L, u128 cr, int k, u128 rem) {
    c->list_nodes++;
    if (DEBUG_LIST) { char a[64], b[64]; u128_str(m, a); u128_str(L, b); fprintf(stderr, "LIST %s %s %d\n", a, b, k); }
    u128 a = modinv128(m, L);
    if (cr != 1) a = mulmod128(a, cr, L);           /* s == c m^{-1} (mod L) */
    if (NMIN > 0) {                       /* start the progression at the first s with m*s >= NMIN */
        u128 smin = (NMIN + m - 1) / m;
        if (a < smin) a += ((smin - a + L - 1) / L) * L;
    }
    for (u128 s = a; s <= rem; s += L) {
        if (atomic_load_explicit(&STOP, memory_order_relaxed)) return;
        if (s <= pmax) continue;
        c->list_candidates++;
        if (((u64)s & 1) == 0 || mod128_64(s, 3) == 0) continue;
        if (SIEVE) {
            /* every prime factor of a completion exceeds pmax, so any q <= pmax dividing s rejects it */
            int bad = 0;
            if ((s >> 64) == 0) {
                u64 s64 = (u64)s;
                for (int j = 0; j < NSMALL && SMALLP[j] <= pmax; j++) if (s64 % SMALLP[j] == 0) { bad = 1; break; }
            } else {
                for (int j = 0; j < NSMALL && SMALLP[j] <= pmax; j++) if (mod128_64(s, SMALLP[j]) == 0) { bad = 1; break; }
            }
            if (bad) { c->sieved++; continue; }
        }
        if (FERMAT == 1) {
            int pass = fermat2(s, m * s - 1);
            c->fermat_passes += (u64)pass;
            if (!pass) continue;
        } else if (FERMAT == 2) {
            continue;                                   /* timing mode: candidates counted, tests skipped */
        }
        c->factorisations++;
        emit_survivor(c, m, pmax, k, s);
    }
}

typedef struct { size_t lo, hi, nch; u128 rem, eb; int list; } NodeHead;

/* node header shared by dfs and the splitter: returns 0 if the node has nothing below it */
/* The largest prime r of n = m r is fixed by the rest: n == c (mod L) gives r == c m^{-1} (mod L), and r^2 < D m.
 * With CN2X_LASTCALC=1, a capped search enumerates that progression above the table bound T at every branching
 * node, so a capped look covers every n whose primes OTHER THAN THE LARGEST are all <= the cap (not only every n
 * whose primes are all <= the cap).  LIST nodes need nothing extra: their listing already covers any completion. */
static void last_calc(Ctx *c, u128 m, u64 pmax, u128 L, u128 cr, int k, u128 b_last, u128 rem) {
    u128 lo = (u128)(T > pmax ? T : pmax) + 1, hi = b_last < rem ? b_last : rem;
    if (hi < lo || L <= 1) return;
    u128 a = modinv128(m, L);
    if (cr != 1) a = mulmod128(a, cr, L);
    if (a < lo) a += ((lo - a + L - 1) / L) * L;
    for (u128 r = a; r <= hi; r += L) {
        c->list_candidates++;
        if (((u64)r & 1) == 0 || mod128_64(r, 3) == 0) continue;
        if (!FERMAT || fermat2(r, m * r - 1)) { c->fermat_passes++; c->factorisations++; emit_survivor(c, m, pmax, k, r); }
    }
}

static int node_head(Ctx *c, u128 m, u64 pmax, u128 L, u128 cr, int k, NodeHead *h) {
    c->nodes++;
    if (k >= 3 && m % L == cr) emit_hit(c, m);
    u128 rem = (X - 1) / m;
    if (rem <= pmax) return 0;
    u128 b_last = k >= 2 ? isqrt128(DD * m - 1) : 0;
    u128 b_second = 0;
    if (k >= 1) { u128 sq = isqrt128(rem); b_second = (DD * m - 1 < sq) ? DD * m - 1 : sq; }
    u128 b_early = icbrt128(rem);
    u128 mx = b_last; if (b_second > mx) mx = b_second; if (b_early > mx) mx = b_early;
    u64 tmax = mx < T ? (u64)mx : T;
    h->lo = upper(pmax); h->hi = upper(tmax);
    h->nch = h->hi > h->lo ? h->hi - h->lo : 0;
    int lc = LASTCALC && k >= 2 && b_last > T;         /* last primes beyond the (capped) table are possible */
    if (h->nch == 0) {                                   /* no table prime fits any role */
        if (lc) last_calc(c, m, pmax, L, cr, k, b_last, rem);
        return 0;
    }
    h->rem = rem;
    h->eb = b_second > b_early ? b_second : b_early;
    h->list = 0;
    if (L > 1) {
        u128 nl = rem / L + 1;
        if ((double)nl <= RATIO * (double)h->nch) h->list = 1;
    }
    if (lc && !h->list) last_calc(c, m, pmax, L, cr, k, b_last, rem);
    return 1;
}

/* one child step: up to two children (one per admissible type of t); returns how many should
 * become nodes, filling mt, Lt[], ct[]; sets *stop on 'break'.  Types: residue of n mod g(t) is 1
 * (rigid) or t (non-rigid, howe/cheb only).  In rigid mode this is cn2xc's child_step exactly. */
static inline int child_step(Ctx *c, u128 m, u128 L, u128 cr, int k, u128 eb, u64 t,
                             u128 *mt_out, u128 *Lt_out, u128 *ct_out, int *stop) {
    c->children++;
    if (t < 64 && ((SKIPMASK >> t) & 1)) return 0;   /* CN2X_SKIP_PRIMES (experiment only) */
    if (mod128_64(L, t) == 0) return 0;
    u64 t2 = t * t - 1, gm = t2 / DD;
    if (gcd64(gm, mod128_64(m, gm)) != 1) return 0;
    u128 mt = m * t;
    if (mt > X) { *stop = 1; return 0; }
    u64 res[2] = {1 % gm, t % gm}; int ntypes = MODE ? 2 : 1;
    if (MODE && res[1] == res[0]) ntypes = 1;          /* g(t) | t-1: the two types coincide */
    int nout = 0, counted_last = 0;
    for (int ty = 0; ty < ntypes; ty++) {
        u128 Lt, ct, n0;
        int mr = merge(L, cr, gm, res[ty], &Lt, &ct, &n0);
        if (mr == 0) continue;                                   /* incompatible type choice */
        if (mr == 2) {                                           /* determined: at most one n < X */
            if (n0 == 0 || n0 % mt != 0) { c->pruned_budget++; continue; }
            u128 s = n0 / mt;
            if (s == 1) { if (k + 1 >= 3) emit_hit(c, n0); continue; }
            if (s <= t) { c->pruned_budget++; continue; }
            c->list_candidates++;
            if (!FERMAT || fermat2(s, n0 - 1)) { c->fermat_passes++; c->factorisations++; emit_survivor(c, mt, t, k + 1, s); }
            continue;
        }
        if ((u128)t > eb) {
            if (!counted_last) { c->last_tests++; counted_last = 1; }
            if (k + 1 >= 3 && mt % Lt == ct) emit_hit(c, mt);
            continue;
        }
        Lt_out[nout] = Lt; ct_out[nout] = ct; nout++;
    }
    *mt_out = mt;
    return nout;
}

typedef struct { int type; u128 m, L, c, rem, eb; u64 pmax, shard; int k; size_t lo, hi; } Task;
/* type 0 = dfs node, 1 = list_mode node, 2 = batch of children [lo, hi) of parent (m, L, c, k, eb) */
static void push_batch(Task *batch, size_t n, u64 shard);

/* hand the children PR[from..hi) of node (m, L, c, k) to the shared stack; counters advance exactly as in dfs */
static void donate(Ctx *c, u128 m, u128 L, u128 cr, int k, u128 eb, size_t from, size_t hi) {
    Task *batch = malloc(2 * (hi - from) * sizeof(Task));
    size_t n = 0;
    for (size_t i = from; i < hi; i++) {
        u128 mt, Lt[2], ct[2]; int stop = 0;
        int nc = child_step(c, m, L, cr, k, eb, PR[i], &mt, Lt, ct, &stop);
        for (int j = 0; j < nc; j++) {
            Task *t = &batch[n++];
            t->type = 0; t->m = mt; t->L = Lt[j]; t->c = ct[j]; t->pmax = PR[i]; t->k = k + 1; t->rem = 0; t->shard = c->shard;
        }
        if (stop) break;
    }
    push_batch(batch, n, c->shard);
    free(batch);
}

static void dfs(Ctx *c, u128 m, u64 pmax, u128 L, u128 cr, int k) {
    NodeHead h;
    if (!node_head(c, m, pmax, L, cr, k, &h)) return;
    if (h.list) { list_mode(c, m, pmax, L, cr, k, h.rem); return; }
    int may_donate = (X / m) > DONATE_Q;
    for (size_t i = h.lo; i < h.hi; i++) {
        u128 mt, Lt[2], ct[2]; int stop = 0;
        int nc = child_step(c, m, L, cr, k, h.eb, PR[i], &mt, Lt, ct, &stop);
        for (int j = 0; j < nc; j++) dfs(c, mt, PR[i], Lt[j], ct[j], k + 1);
        if (stop) break;
        if (atomic_load_explicit(&STOP, memory_order_relaxed)) return;
        if (may_donate && i + 1 < h.hi && atomic_load_explicit(&HUNGRY, memory_order_relaxed) > 0) {
            donate(c, m, L, cr, k, h.eb, i + 1, h.hi);
            return;
        }
    }
}

/* ---------------------------------------------------------------- shards */
typedef struct { atomic_long outstanding; Ctx acc; } Shard;
#define CHUNK 65536
static Shard *SHARD_CHUNKS[1 << 16];
static u64 NSHARDS = 0;                       /* generated so far this session (under gen_mu) */
static pthread_mutex_t acc_mu = PTHREAD_MUTEX_INITIALIZER;
static uint8_t *DONE_BITS = NULL;             /* shards finished in earlier sessions */
static u64 DONE_CAP = 0;
static Ctx DONE_PRIOR;                        /* summed counters of those shards */
static u64 DONE_PRIOR_N = 0, DONE_NOW_N = 0;
static FILE *DONE;
static u128 MSHARD = 100000;
static u64 BATCH = 4096;

static inline Shard *shard_ptr(u64 id) { return &SHARD_CHUNKS[id >> 16][id & 0xFFFF]; }
static inline int was_done(u64 id) { return id < DONE_CAP && (DONE_BITS[id >> 3] >> (id & 7)) & 1; }
static u64 new_shard(void) {
    u64 id = NSHARDS++;
    if ((id & 0xFFFF) == 0) SHARD_CHUNKS[id >> 16] = calloc(CHUNK, sizeof(Shard));
    return id;
}
static void add_counts(Ctx *dst, const Ctx *src) {
    dst->nodes += src->nodes; dst->children += src->children; dst->last_tests += src->last_tests;
    dst->list_nodes += src->list_nodes; dst->list_candidates += src->list_candidates;
    dst->fermat_passes += src->fermat_passes; dst->factorisations += src->factorisations;
    dst->pruned_budget += src->pruned_budget; dst->sieved += src->sieved;
}
static void record_done(u64 id) {
    Shard *sh = shard_ptr(id);
    pthread_mutex_lock(&out_mu);
    fflush(OUT);
    fprintf(DONE, "D %llu %llu %llu %llu %llu %llu %llu %llu %llu %llu\n", (unsigned long long)id,
            (unsigned long long)sh->acc.nodes, (unsigned long long)sh->acc.children, (unsigned long long)sh->acc.last_tests,
            (unsigned long long)sh->acc.list_nodes, (unsigned long long)sh->acc.list_candidates,
            (unsigned long long)sh->acc.fermat_passes, (unsigned long long)sh->acc.factorisations,
            (unsigned long long)sh->acc.pruned_budget, (unsigned long long)sh->acc.sieved);
    fflush(DONE);
    DONE_NOW_N++;
    pthread_mutex_unlock(&out_mu);
}
static void load_done(const char *path) {
    FILE *f = fopen(path, "r");
    memset(&DONE_PRIOR, 0, sizeof DONE_PRIOR);
    if (!f) return;
    char line[512];
    unsigned long long id, v[9];
    while (fgets(line, sizeof line, f)) {
        if (sscanf(line, "D %llu %llu %llu %llu %llu %llu %llu %llu %llu %llu",
                   &id, &v[0], &v[1], &v[2], &v[3], &v[4], &v[5], &v[6], &v[7], &v[8]) != 10) continue;  /* torn line */
        if (id >= DONE_CAP) {
            u64 cap = DONE_CAP ? DONE_CAP : 1 << 20;
            while (cap <= id) cap *= 2;
            DONE_BITS = realloc(DONE_BITS, cap / 8 + 1);
            memset(DONE_BITS + DONE_CAP / 8, 0, cap / 8 + 1 - DONE_CAP / 8);
            DONE_CAP = cap;
        }
        if (was_done(id)) continue;             /* duplicate line: count once */
        DONE_BITS[id >> 3] |= (uint8_t)(1 << (id & 7));
        Ctx c = {0};
        c.nodes = v[0]; c.children = v[1]; c.last_tests = v[2]; c.list_nodes = v[3]; c.list_candidates = v[4];
        c.fermat_passes = v[5]; c.factorisations = v[6]; c.pruned_budget = v[7]; c.sieved = v[8];
        add_counts(&DONE_PRIOR, &c);
        DONE_PRIOR_N++;
    }
    fclose(f);
}

/* ---------------------------------------------------------------- shared task stack */
static Task *QUEUE;
static size_t QLEN = 0, QCAP = 0;
static int ACTIVE = 0, GEN_DONE = 0;
static pthread_mutex_t q_mu = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t q_cv = PTHREAD_COND_INITIALIZER;

static void push_batch(Task *batch, size_t n, u64 shard) {
    if (!n) return;
    atomic_fetch_add(&shard_ptr(shard)->outstanding, (long)n);
    pthread_mutex_lock(&q_mu);
    if (QLEN + n > QCAP) {
        while (QLEN + n > QCAP) QCAP = QCAP ? 2 * QCAP : 1 << 16;
        QUEUE = realloc(QUEUE, QCAP * sizeof(Task));
    }
    for (size_t i = n; i-- > 0;) QUEUE[QLEN++] = batch[i];
    pthread_cond_broadcast(&q_cv);
    pthread_mutex_unlock(&q_mu);
}

/* ---------------------------------------------------------------- deterministic shard generator */
typedef struct { u128 m, L, c, eb; u64 pmax; int k; size_t i, hi, frontier, xlim; u64 bsize;
                 int has_pend; u128 pm, pL, pc; u64 pt; } Frame;   /* pend: a second typed child still to expand */
static Frame STACK[256];
static int SP = 0;
static Ctx GEN;
static pthread_mutex_t gen_mu = PTHREAD_MUTEX_INITIALIZER;

/* expand a node in the generator: push a frame, or make it a LIST shard (returns 1 and sets task if a shard is due) */
static int gen_node(u128 m, u64 pmax, u128 L, u128 cr, int k, Task *task) {
    NodeHead h;
    if (!node_head(&GEN, m, pmax, L, cr, k, &h)) return 0;
    if (h.list) {
        u64 id = new_shard();
        if (was_done(id)) return 0;
        task->type = 1; task->m = m; task->pmax = pmax; task->L = L; task->c = cr; task->k = k; task->rem = h.rem; task->shard = id;
        return 1;
    }
    Frame *f = &STACK[SP++];
    f->m = m; f->L = L; f->c = cr; f->has_pend = 0; f->eb = h.eb; f->pmax = pmax; f->k = k; f->i = h.lo; f->hi = h.hi;
    u128 q = (MSHARD - 1) / m;                            /* children t > q have m*t >= MSHARD: frontier */
    f->frontier = upper(q < T ? (u64)q : T);
    if (f->frontier < f->i) f->frontier = f->i;
    u128 xq = X / m;                                      /* children t > xq have m*t > X */
    f->xlim = upper(xq < T ? (u64)xq : T);
    f->bsize = 1;
    return 0;
}
static int next_shard(Task *task) {
    pthread_mutex_lock(&gen_mu);
    while (SP > 0 && !atomic_load(&STOP)) {
        Frame *f = &STACK[SP - 1];
        if (f->has_pend) {                                  /* second typed child of an earlier expansion */
            f->has_pend = 0;
            if (gen_node(f->pm, f->pt, f->pL, f->pc, f->k + 1, task)) { pthread_mutex_unlock(&gen_mu); return 1; }
            continue;
        }
        if (f->i >= f->hi) { SP--; continue; }
        if (f->i < f->frontier) {                           /* expansion child: handled here */
            u64 t = PR[f->i++];
            u128 mt2, Lt[2], ct[2]; int stop = 0;
            int nc = child_step(&GEN, f->m, f->L, f->c, f->k, f->eb, t, &mt2, Lt, ct, &stop);
            if (stop) f->i = f->hi;
            if (nc == 2) { f->has_pend = 1; f->pm = mt2; f->pt = t; f->pL = Lt[1]; f->pc = ct[1]; }
            if (nc >= 1) {
                int kk = f->k;
                if (gen_node(mt2, t, Lt[0], ct[0], kk + 1, task)) { pthread_mutex_unlock(&gen_mu); return 1; }
            }
            continue;
        }
        /* frontier: consecutive runs of 1, 2, 4, ... BATCH children; once past xlim, take the rest */
        size_t lo = f->i, hi = lo + f->bsize;
        if (hi > f->hi) hi = f->hi;
        if (hi > f->xlim) hi = f->hi;
        f->i = hi;
        if (f->bsize < BATCH) f->bsize *= 2;
        u64 id = new_shard();
        if (was_done(id)) continue;
        task->type = 2; task->m = f->m; task->L = f->L; task->c = f->c; task->k = f->k; task->eb = f->eb; task->pmax = f->pmax;
        task->lo = lo; task->hi = hi; task->shard = id;
        pthread_mutex_unlock(&gen_mu);
        return 1;
    }
    int exhausted = SP == 0;
    pthread_mutex_unlock(&gen_mu);
    return exhausted ? 0 : -1;                             /* -1: stopped */
}

static void run_task(const Task *tk) {
    Ctx c = {0};
    c.shard = tk->shard;
    if (tk->type == 0) dfs(&c, tk->m, tk->pmax, tk->L, tk->c, tk->k);
    else if (tk->type == 1) list_mode(&c, tk->m, tk->pmax, tk->L, tk->c, tk->k, tk->rem);
    else {
        for (size_t i = tk->lo; i < tk->hi; i++) {
            if (atomic_load_explicit(&STOP, memory_order_relaxed)) break;
            u128 mt, Lt[2], ct[2]; int stop = 0;
            int nc = child_step(&c, tk->m, tk->L, tk->c, tk->k, tk->eb, PR[i], &mt, Lt, ct, &stop);
            for (int j = 0; j < nc; j++) dfs(&c, mt, PR[i], Lt[j], ct[j], tk->k + 1);
            if (stop) break;
        }
    }
    Shard *sh = shard_ptr(tk->shard);
    pthread_mutex_lock(&acc_mu);
    add_counts(&sh->acc, &c);
    pthread_mutex_unlock(&acc_mu);
    atomic_fetch_add(&G_NODES, c.nodes);
    atomic_fetch_add(&G_CAND, c.list_candidates);
    long left = atomic_fetch_sub(&sh->outstanding, 1) - 1;
    if (left == 0 && !atomic_load(&STOP)) record_done(tk->shard);
}

static void *worker(void *arg) {
    (void)arg;
    Task tk;
    pthread_mutex_lock(&q_mu);
    for (;;) {
        if (atomic_load(&STOP)) break;
        if (QLEN > 0) {
            tk = QUEUE[--QLEN]; ACTIVE++;
            pthread_mutex_unlock(&q_mu);
            run_task(&tk);
            pthread_mutex_lock(&q_mu);
            ACTIVE--;
            if (QLEN == 0 && ACTIVE == 0) pthread_cond_broadcast(&q_cv);
            continue;
        }
        if (!GEN_DONE && !atomic_load(&DRAIN)) {
            ACTIVE++;
            pthread_mutex_unlock(&q_mu);
            int got = next_shard(&tk);
            if (got == 1) {
                atomic_store(&shard_ptr(tk.shard)->outstanding, 1);
                run_task(&tk);
            }
            pthread_mutex_lock(&q_mu);
            ACTIVE--;
            if (got == 0) { GEN_DONE = 1; pthread_cond_broadcast(&q_cv); }
            if (got == -1) break;
            continue;
        }
        if (ACTIVE > 0) {
            atomic_fetch_add(&HUNGRY, 1);
            pthread_cond_wait(&q_cv, &q_mu);
            atomic_fetch_sub(&HUNGRY, 1);
            continue;
        }
        break;
    }
    pthread_cond_broadcast(&q_cv);
    pthread_mutex_unlock(&q_mu);
    return NULL;
}

/* ---------------------------------------------------------------- stop sources */
static double DEADLINE = 0;
static char STOP_FILE[1024] = "", DRAIN_FILE[1024] = "";
static void on_signal(int sig) { (void)sig; atomic_store(&STOP, 1); }
static void *monitor(void *arg) {
    double t_start = *(double *)arg, t_last = t_start;
    unsigned long long n_last = 0;
    for (;;) {
        struct timespec ts = {1, 0};
        nanosleep(&ts, NULL);
        if (PROGRESS_SEC > 0 && now() - t_last >= PROGRESS_SEC) {
            double el = now() - t_start, dt = now() - t_last;
            unsigned long long nn = atomic_load(&G_NODES), cc = atomic_load(&G_CAND);
            size_t ri = 0, rh = 0;
            if (pthread_mutex_trylock(&gen_mu) == 0) {           /* never block the generator for a log line */
                if (SP > 0) { ri = STACK[0].i; rh = STACK[0].hi; }
                pthread_mutex_unlock(&gen_mu);
            }
            time_t wall = time(NULL); struct tm tmv; localtime_r(&wall, &tmv);
            double total_nodes = (double)nn + (double)DONE_PRIOR.nodes;     /* earlier sessions included */
            char pct[96] = "";
            if (EXPECT_NODES > 0) {
                double frac = total_nodes / EXPECT_NODES;
                double rate = nn / (el > 0 ? el : 1);
                double eta = rate > 0 ? (EXPECT_NODES - total_nodes) / rate : 0;
                if (eta < 0) eta = 0;
                if (eta >= 3600) snprintf(pct, sizeof pct, "  ~%.2f%% of predicted, ETA %.1f h", 100 * frac, eta / 3600);
                else snprintf(pct, sizeof pct, "  ~%.2f%% of predicted, ETA %.0f min", 100 * frac, eta / 60);
            }
            fprintf(stderr, "[cn2xc %02d:%02d:%02d] elapsed %6.0fs  shards done %9llu  nodes %.3e (%.2e/s)  "
                            "candidates %.3e  hits %llu  root %zu/%zu%s%s\n",
                    tmv.tm_hour, tmv.tm_min, tmv.tm_sec, el, (unsigned long long)DONE_NOW_N,
                    total_nodes, (nn - n_last) / dt, (double)cc + (double)DONE_PRIOR.list_candidates,
                    (unsigned long long)atomic_load(&G_HITS) + 0ull, ri, rh, pct,
                    atomic_load(&DRAIN) ? "  DRAINING" : "");
            fflush(stderr);
            t_last = now(); n_last = nn;
        }
        pthread_mutex_lock(&q_mu);
        int finished = (GEN_DONE || atomic_load(&DRAIN)) && QLEN == 0 && ACTIVE == 0;
        pthread_mutex_unlock(&q_mu);
        if (finished || atomic_load(&STOP)) break;
        if (DEADLINE > 0 && (double)time(NULL) >= DEADLINE) atomic_store(&DRAIN, 1);
        if (STOP_FILE[0]) { FILE *f = fopen(STOP_FILE, "r"); if (f) { fclose(f); atomic_store(&STOP, 1); } }
        if (DRAIN_FILE[0]) { FILE *f = fopen(DRAIN_FILE, "r"); if (f) { fclose(f); atomic_store(&DRAIN, 1); } }
        if (atomic_load(&DRAIN)) { pthread_mutex_lock(&q_mu); pthread_cond_broadcast(&q_cv); pthread_mutex_unlock(&q_mu); }
        pthread_mutex_lock(&out_mu);
        fflush(OUT); fflush(DONE);
        fsync(fileno(OUT)); fsync(fileno(DONE));
        pthread_mutex_unlock(&out_mu);
    }
    pthread_mutex_lock(&q_mu);
    pthread_cond_broadcast(&q_cv);
    pthread_mutex_unlock(&q_mu);
    return NULL;
}

/* ---------------------------------------------------------------- setup */
static void sieve(u64 limit) {
    NPR = 0; PR = NULL;
    if (limit < 5) return;
    uint8_t *comp = calloc(limit + 1, 1);
    size_t cap = (size_t)(1.3 * limit / (log((double)limit) + 1)) + 1000;
    PR = malloc(cap * sizeof(u32));
    for (u64 p = 2; p * p <= limit; p++)
        if (!comp[p]) for (u64 q = p * p; q <= limit; q += p) comp[q] = 1;
    for (u64 p = 5; p <= limit; p++) if (!comp[p]) {
        if (NPR == cap) { cap *= 2; PR = realloc(PR, cap * sizeof(u32)); }
        PR[NPR++] = (u32)p;
    }
    free(comp);
}
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }

static inline u64 xorshift(u64 *x) { *x ^= *x << 13; *x ^= *x >> 7; *x ^= *x << 17; return *x; }
static inline u128 rnd128(u64 *x) { u128 hi = xorshift(x); return (hi << 64) | xorshift(x); }
static int selftest(long n) {                          /* mulmod128 against the reference, random operands */
    u64 x = 0x9E3779B97F4A7C15ull; long bad = 0;
    for (long i = 0; i < n; i++) {
        int mb = 1 + (int)(xorshift(&x) % 126);           /* modulus bit length 1..126 */
        u128 mod = rnd128(&x) >> (128 - mb); if (mod < 1) mod = 1;
        u128 a = rnd128(&x), b = rnd128(&x);
        if (xorshift(&x) & 1) a %= mod;
        if (xorshift(&x) & 1) b %= mod;
        u128 r1 = mulmod128(a, b, mod), r2 = mulmod_ref(a, b, mod);
        if (r1 != r2) bad++;
        if (i < 200) { char sa[64], sb[64], sm[64], sr[64]; u128_str(a, sa); u128_str(b, sb); u128_str(mod, sm); u128_str(r1, sr);
                       printf("%s %s %s %s\n", sa, sb, sm, sr); }
    }
    fprintf(stderr, "selftest: %ld mulmod128 checks, %ld mismatches\n", n, bad);
    return bad ? 1 : 0;
}

int main(int argc, char **argv) {
    if (argc >= 2 && !strcmp(argv[1], "--selftest")) return selftest(argc > 2 ? atol(argv[2]) : 10000000);
    if (argc < 3) {
        fprintf(stderr, "usage: cn2xc X RUNDIR [threads=30] [deadline_epoch=0] [ratio=1.0] [sieve=1] "
                        "[mshard=100000] [batch=4096] [donate_q=1e10] [progress_sec=60] [expected_nodes=0]\n"
                        "  deadline: at that epoch second, stop taking shards, finish in-flight ones, exit (also: touch RUNDIR/DRAIN)\n"
                        "  hard stop (abandon in-flight shards): touch RUNDIR/STOP, or SIGINT/SIGTERM\n");
        return 2;
    }
    X = parse_u128(argv[1]);
    const char *dir = argv[2];
    int threads = argc > 3 ? atoi(argv[3]) : 30;
    DEADLINE = argc > 4 ? atof(argv[4]) : 0;
    RATIO = argc > 5 ? atof(argv[5]) : 1.0;
    SIEVE = argc > 6 ? atoi(argv[6]) : 1;
    if (argc > 7) MSHARD = parse_u128(argv[7]);
    if (argc > 8) BATCH = strtoull(argv[8], NULL, 10);
    if (argc > 9) DONATE_Q = parse_u128(argv[9]);
    if (argc > 10) PROGRESS_SEC = atof(argv[10]);
    if (argc > 11) EXPECT_NODES = atof(argv[11]);
    FERMAT = 1;
    char path[1200];
    snprintf(STOP_FILE, sizeof STOP_FILE, "%s/STOP", dir);
    snprintf(DRAIN_FILE, sizeof DRAIN_FILE, "%s/DRAIN", dir);
    remove(STOP_FILE); remove(DRAIN_FILE);
    snprintf(path, sizeof path, "%s/done.log", dir);
    load_done(path);
    DONE = fopen(path, "a");
    snprintf(path, sizeof path, "%s/out.txt", dir);
    OUT = fopen(path, "a");
    if (!DONE || !OUT) { fprintf(stderr, "cannot open run files in %s\n", dir); return 2; }
    signal(SIGINT, on_signal);
    signal(SIGTERM, on_signal);

    {
        const char *lcs = getenv("CN2X_LASTCALC");
        if (lcs && lcs[0] == '1') LASTCALC = 1;
        const char *md = getenv("CN2X_MODE");
        if (md && !strcmp(md, "howe")) MODE = 1;
        else if (md && !strcmp(md, "cheb")) { MODE = 2; DD = 2; }
        else if (md && !strcmp(md, "test")) {    /* VALIDATION ONLY: modulus (p^2-1)/CN2X_D, Fermat off */
            MODE = 3; const char *dd = getenv("CN2X_D"); DD = dd ? strtoull(dd, NULL, 10) : 24;
            if (24 % DD) { fprintf(stderr, "CN2X_D must divide 24\n"); return 2; }
        }
        else if (md && strcmp(md, "rigid")) { fprintf(stderr, "CN2X_MODE must be rigid, howe, cheb or test\n"); return 2; }
        if (MODE == 3) FERMAT = 0;
        if ((X >> 126) != 0) { fprintf(stderr, "X must be below 2^126 (Fermat test width)\n"); return 2; }
    }
    T = (u64)icbrt128(DD * (X - 1));                       /* n/p > g(p) = (p^2-1)/DD  =>  p^3 < DD X */
    if (T > 0xFFFFFFFFull && !getenv("CN2X_TCAP")) { fprintf(stderr, "X too large: primes up to %llu exceed the 32-bit prime table\n", (unsigned long long)T); return 2; }
    {   /* EXPERIMENT ONLY: CN2X_TCAP makes the search CONDITIONAL on every factor being <= the cap
         * (the lazy sweep's H-factor); CN2X_NMIN only restricts what is reported.  Both change the
         * shard structure, so a resumed run must use the same values -- they are logged in the run dir. */
        const char *nm = getenv("CN2X_NMIN");
        if (nm) NMIN = parse_u128(nm);
        NSMALL = (X < (u128)10000000000ull * 10000000000ull) ? 13 : 95;   /* 1e20 */
        const char *ns = getenv("CN2X_SIEVEP");
        if (ns) { int v = atoi(ns); if (v > 0 && v <= 95) NSMALL = v; }
        const char *sk = getenv("CN2X_SKIP_PRIMES");
        for (const char *q = sk; q && *q; ) { long v = strtol(q, (char **)&q, 10); if (v > 0 && v < 64) SKIPMASK |= 1ull << v; while (*q == ',' || *q == ' ') q++; if (*q && (*q < '0' || *q > '9')) break; }
        const char *cap = getenv("CN2X_TCAP");
        if (cap) { u64 c = strtoull(cap, NULL, 10); if (c < T) T = c; }
    }
    double t0 = now();
    sieve(T);
    double t_sieve = now() - t0;

    t0 = now();
    Task dummy;
    gen_node(1, 3, 1, 0, 0, &dummy);                          /* root always branches */
    pthread_t mon;
    double t_mon = now();
    pthread_create(&mon, NULL, monitor, &t_mon);
    pthread_t *th = malloc(threads * sizeof(pthread_t));
    for (int i = 0; i < threads; i++) pthread_create(&th[i], NULL, worker, NULL);
    for (int i = 0; i < threads; i++) pthread_join(th[i], NULL);
    atomic_store(&STOP, atomic_load(&STOP));               /* monitor exits on its own */
    pthread_join(mon, NULL);
    double t_search = now() - t0;

    int complete = GEN_DONE && QLEN == 0 && !atomic_load(&STOP);
    fflush(OUT); fflush(DONE); fsync(fileno(OUT)); fsync(fileno(DONE));
    fclose(OUT); fclose(DONE);

    /* totals are exact only when complete: generator counters + every shard, earlier sessions included */
    Ctx total = GEN;
    add_counts(&total, &DONE_PRIOR);
    for (u64 id = 0; id < NSHARDS; id++) {
        if (was_done(id)) continue;
        Shard *sh = shard_ptr(id);
        if (atomic_load(&sh->outstanding) == 0 && (sh->acc.nodes || sh->acc.children || sh->acc.list_nodes))
            add_counts(&total, &sh->acc);
    }
    char xs[64]; u128_str(X, xs);
    printf("{\"skipmask\":%llu,\"lastcalc\":%d,\"mode\":%d,\"X\":\"%s\",\"complete\":%s,\"stopped\":%s,\"drained\":%s,\"threads\":%d,\"shards_generated\":%llu,"
           "\"shards_done_prior\":%llu,\"shards_done_now\":%llu,"
           "\"nodes\":%llu,\"children\":%llu,\"last_tests\":%llu,\"list_nodes\":%llu,\"list_candidates\":%llu,"
           "\"fermat_passes\":%llu,\"factorisations\":%llu,\"pruned_budget\":%llu,\"sieved\":%llu,"
           "\"sieve_seconds\":%.3f,\"search_seconds\":%.3f}\n",
           (unsigned long long)SKIPMASK, LASTCALC, MODE, xs, complete ? "true" : "false", atomic_load(&STOP) ? "true" : "false", atomic_load(&DRAIN) ? "true" : "false", threads,
           (unsigned long long)NSHARDS, (unsigned long long)DONE_PRIOR_N, (unsigned long long)DONE_NOW_N,
           (unsigned long long)total.nodes, (unsigned long long)total.children, (unsigned long long)total.last_tests,
           (unsigned long long)total.list_nodes, (unsigned long long)total.list_candidates,
           (unsigned long long)total.fermat_passes, (unsigned long long)total.factorisations,
           (unsigned long long)total.pruned_budget, (unsigned long long)total.sieved, t_sieve, t_search);
    return complete ? 0 : 3;
}
