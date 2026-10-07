/*
 * cn2pair.c -- the pair family of the two-large-prime tail (see cn2_pair_tail.py and CN2_LastCalc_Coverage_Proof.md).
 *
 * Covers every order-2 n < X whose two largest primes p < q satisfy  T < p < q <= QSTAR.
 * For each such pair and each pair of type residues (r_p, r_q), n lies in ONE residue class:
 *     n == 0 (mod q),  n == r_q (mod g(q)),  n == r_p (mod g(p)),  n == 0 (mod p),     g(x) = (x^2-1)/D,
 * so n == c (mod M), M = p q lcm(g(p), g(q)) (p, q are prime to g(q), g(p) for an order-2 n: F5).  The engine enumerates
 * every member n of that class with q g(q) < n < X (F2: n/q > g(q)), keeps n odd with 3 !| n (F1), passing the
 * small-prime sieve (for each s in 5..59 with s | n: s^2 !| n and n == 1 or s (mod g(s)) -- the definition itself),
 * and with 2^(n-1) == 1 (mod n) (an order-2 n is Carmichael: F6), and writes "N n".  Survivors are factored and
 * checked by the Python driver.  The sieve was added on 5 Oct 2026 (before it, every odd n with 3 !| n was tested).
 * Pairs with p | g(q) are skipped (F5).  Types: rigid r = 1; howe/cheb r in {1, x} (env CN2X_MODE as cn2xh; cheb D = 2).
 *
 * Together with a capped cn2xh run (CN2X_LASTCALC=1, table bound T) and cn2tail with R0 = QSTAR (every n with a prime
 * above QSTAR), this covers every order-2 n < X: an n with p_{k-1} <= T is the capped run's; an n with p_{k-1} > T
 * has its largest prime q either > QSTAR (cn2tail) or <= QSTAR (here, with p = p_{k-1}).
 *
 * Arithmetic: X < 2^126.  Every modulus is formed only after checking it stays <= X - 1; once a modulus would exceed
 * X - 1 the class has at most one member below X, which is computed without forming the product (as cn2xh's merge).
 *
 * Checkpointing: the q range is cut into NSHARDS deterministic shards of roughly equal estimated work; a finished shard
 * is appended to RUNDIR/done.log after its survivors are flushed to RUNDIR/out.txt; a restart skips finished shards.
 * RUNDIR/DRAIN or the deadline: take no new shards; RUNDIR/STOP or SIGINT/SIGTERM: abandon in-flight shards.
 *
 * Build:  clang -O3 -mcpu=native -o cn2pair cn2pair.c -lpthread -lm
 * Usage:  cn2pair X T QSTAR RUNDIR [threads=30] [deadline_epoch=0] [nshards=4096] [progress_sec=60] [count_only=0]
 *         count_only=1: count candidates exactly (no Fermat tests, no output) -- the dry run.
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

static u128 X;
static u64 T, QSTAR, DD = 1;
static int MODE = 1, COUNT_ONLY = 0;
static u32 *PR; static size_t NPR;            /* primes in (T, QSTAR] */
static FILE *OUT, *DONE;
static pthread_mutex_t out_mu = PTHREAD_MUTEX_INITIALIZER;
static atomic_int STOP = 0, DRAIN = 0;
static char STOP_FILE[1100], DRAIN_FILE[1100];
static double DEADLINE = 0, PROGRESS_SEC = 60;

/* ---------------------------------------------------------------- arithmetic (as cn2xh.c) */
static u128 parse_u128(const char *s) { u128 v = 0; for (; *s; s++) if (*s >= '0' && *s <= '9') v = v * 10 + (u64)(*s - '0'); return v; }
static void u128_str(u128 v, char *buf) {
    char tmp[64]; int i = 0;
    if (v == 0) { strcpy(buf, "0"); return; }
    while (v) { tmp[i++] = (char)('0' + (int)(v % 10)); v /= 10; }
    int j = 0; while (i) buf[j++] = tmp[--i]; buf[j] = 0;
}
static inline u64 gcd64(u64 a, u64 b) {
    if (!a) return b;
    if (!b) return a;
    int sh = __builtin_ctzll(a | b);
    a >>= __builtin_ctzll(a);
    do { b >>= __builtin_ctzll(b); if (a > b) { u64 t = a; a = b; b = t; } b -= a; } while (b);
    return a << sh;
}
static u64 modinv64(u64 a, u64 mod) {             /* gcd(a, mod) = 1, mod > 1 */
    i128 x0 = 0, x1 = 1; u64 r0 = mod, r1 = a % mod;
    while (r1) { u64 q = r0 / r1, r2 = r0 - q * r1; r0 = r1; r1 = r2; i128 x2 = x0 - (i128)q * x1; x0 = x1; x1 = x2; }
    if (x0 < 0) x0 += mod;
    return (u64)x0;
}
static inline u64 mulmod64(u64 a, u64 b, u64 m) { return (u64)(((u128)a * b) % m); }
static inline int topbit128(u128 e) { u64 hi = (u64)(e >> 64); return hi ? 127 - __builtin_clzll(hi) : 63 - __builtin_clzll((u64)e); }
static inline u64 inv_mod_2_64(u64 n) { u64 x = n; for (int i = 0; i < 5; i++) x *= 2 - n * x; return x; }
static inline u128 montmul2(u128 a, u128 b, u128 N, u64 np) {
    u64 a0 = (u64)a, a1 = (u64)(a >> 64), b0 = (u64)b, b1 = (u64)(b >> 64);
    u64 n0 = (u64)N, n1 = (u64)(N >> 64);
    u64 t0, t1, t2, t3, m;
    u128 C;
    C = (u128)a0 * b0;              t0 = (u64)C; C >>= 64;
    C += (u128)a1 * b0;             t1 = (u64)C; C >>= 64;
    t2 = (u64)C; t3 = 0;
    m = t0 * np;
    C = (u128)t0 + (u128)m * n0;    C >>= 64;
    C += (u128)t1 + (u128)m * n1;   t0 = (u64)C; C >>= 64;
    C += (u128)t2;                  t1 = (u64)C; C >>= 64;
    t2 = t3 + (u64)C;
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
static int fermat2(u128 s, u128 e) {               /* 2^e == 1 (mod s); s odd, 1 < s < 2^126 */
    int top = topbit128(e);
    if ((s >> 63) == 0) {
        u64 N = (u64)s, np = (u64)0 - inv_mod_2_64(N);
        u64 one = (u64)((((u128)1) << 64) % N), x = one;
        for (int i = top; i >= 0; i--) {
            u128 t = (u128)x * x; u64 mm = (u64)t * np;
            x = (u64)((t + (u128)mm * N) >> 64);
            if (x >= N) x -= N;
            if ((e >> i) & 1) { x += x; if (x >= N) x -= N; }
        }
        return x == one;
    } else {
        u128 N = s; u64 np = (u64)0 - inv_mod_2_64((u64)N);
        u128 one = ((u128)0 - N) % N, x = one;
        for (int i = top; i >= 0; i--) { x = montmul2(x, x, N, np); if ((e >> i) & 1) { x += x; if (x >= N) x -= N; } }
        return x == one;
    }
}

/* ---------------------------------------------------------------- per-thread counters */
typedef struct { u64 pairs, skipped, classes, candidates, tests, passes, sieved; } Ctx;

/* small-prime sieve (definition only, as cn2tail's): if a prime s divides an order-2 n then s^2 !| n and
 * n == 1 or s (mod g(s)) (rigid: 1 only).  A candidate failing this for some s | n is not order-2.  s > 3 since
 * 2, 3 never divide an order-2 n (F1, tested separately); all s here are < 10^4 <= T < p < q, so s is never p or q. */
#define NS 15
static const u64 SP[NS] = {5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59};
static u64 SG[NS], SR1[NS], SR2[NS];                 /* g(s) and the allowed residues mod g(s) */
static void sieve_init(void) {
    for (int j = 0; j < NS; j++) { SG[j] = (SP[j] * SP[j] - 1) / DD; SR1[j] = 1 % SG[j]; SR2[j] = MODE ? SP[j] % SG[j] : SR1[j]; }
}
static inline int sieve_ok(u128 n, int j) {          /* called only when SP[j] | n */
    if ((n / SP[j]) % SP[j] == 0) return 0;          /* s^2 | n: not squarefree */
    u64 r = (u64)(n % SG[j]);
    return r == SR1[j] || r == SR2[j];
}

static void emit(u128 n) {
    char a[64], line[80]; u128_str(n, a); snprintf(line, sizeof line, "N %s\n", a);
    pthread_mutex_lock(&out_mu); fputs(line, OUT); pthread_mutex_unlock(&out_mu);
}

/* the members n of the class c (mod M) with lo < n < X; M <= X - 1 */
static void walk(Ctx *c, u128 cls, u128 M, u128 lo) {
    c->classes++;
    u128 n = cls;
    if (n <= lo) n += ((lo - n) / M + 1) * M;
    if (COUNT_ONLY) { if (n < X) c->candidates += (u64)((X - 1 - n) / M + 1); return; }
    if (n >= X) return;
    u64 r[NS], dm[NS];                                    /* n mod s, advanced by M mod s per step: no divisions in the loop */
    for (int j = 0; j < NS; j++) { r[j] = (u64)(n % SP[j]); dm[j] = (u64)(M % SP[j]); }
    u64 r3 = (u64)(n % 3), d3 = (u64)(M % 3), r2 = (u64)(n & 1), d2 = (u64)(M & 1);
    for (; n < X; n += M) {
        c->candidates++;
        int bad = r2 == 0 || r3 == 0;                     /* F1: n odd, 3 !| n */
        r2 ^= d2; r3 += d3; if (r3 >= 3) r3 -= 3;
        for (int j = 0; j < NS; j++) {
            if (!bad && r[j] == 0 && !sieve_ok(n, j)) bad = 1;
            r[j] += dm[j]; if (r[j] >= SP[j]) r[j] -= SP[j];
        }
        if (bad) { c->sieved++; continue; }
        c->tests++;
        if (fermat2(n, n - 1)) { c->passes++; emit(n); }
    }
}
static void single(Ctx *c, u128 n, u64 p, u128 lo) {     /* a determined class: one candidate, still to check n == 0 (mod p) */
    c->classes++;
    if (n <= lo || n >= X || n % p) return;
    c->candidates++;
    if (COUNT_ONLY) return;
    if (!(n & 1) || n % 3 == 0) return;
    for (int j = 0; j < NS; j++) if (n % SP[j] == 0 && !sieve_ok(n, j)) { c->sieved++; return; }
    c->tests++;
    if (fermat2(n, n - 1)) { c->passes++; emit(n); }
}

/* all pairs (p, q) with q = PR[iq] and p = PR[ip], ip < iq */
static void do_q(Ctx *c, size_t iq) {
    u64 q = PR[iq], gq = (q * q - 1) / DD;
    u128 Lq = (u128)q * gq;                               /* n == 0 (mod q), n == r_q (mod g(q)): one class mod q g(q) */
    u128 lo = Lq;                                         /* F2: n > q g(q) */
    if (Lq >= X - 1) return;                              /* then no n < X has largest prime q (and X - 1 - c below is safe) */
    u64 rq[2] = {1 % gq, q % gq}; int nq = (MODE && rq[1] != rq[0]) ? 2 : 1;
    u128 cq[2];
    u64 igq = modinv64(gq % q, q);
    for (int t = 0; t < nq; t++) {                        /* c = r + g k, k == -r g^{-1} (mod q) */
        u64 k = (q - mulmod64(rq[t] % q, igq, q)) % q;
        cq[t] = (u128)rq[t] + (u128)gq * k;
    }
    for (size_t ip = 0; ip < iq; ip++) {
        if ((ip & 255) == 0 && atomic_load_explicit(&STOP, memory_order_relaxed)) return;   /* hard stop: shard abandoned */
        u64 p = PR[ip];
        c->pairs++;
        if (gq % p == 0) { c->skipped++; continue; }      /* p | g(q): no order-2 n (F5) */
        u64 gp = (p * p - 1) / DD;
        u64 rp[2] = {1 % gp, p % gp}; int np_ = (MODE && rp[1] != rp[0]) ? 2 : 1;
        /* step A: merge (Lq, cq) with (gp, rp).  gg, g2, inv depend only on the moduli. */
        u64 Lq_gp = (u64)(Lq % gp);
        u64 gg = gcd64(gp, Lq_gp), g2 = gp / gg;
        u64 invA = g2 > 1 ? modinv64((u64)((Lq / gg) % g2), g2) : 0;
        int detA = g2 > 1 && Lq > (X - 1) / g2;           /* L1 = Lq g2 > X - 1 */
        u128 L1 = detA ? 0 : Lq * g2;
        u64 invB = 0; int detB = 0; u128 L2 = 0;
        if (!detA) {                                      /* step B: merge (L1, c1) with (p, 0); p !| L1 */
            invB = modinv64((u64)(L1 % p), p);
            detB = L1 > (X - 1) / p;
            if (!detB) L2 = L1 * p;
        }
        for (int tq = 0; tq < nq; tq++) {
            u64 cm = (u64)(cq[tq] % gp);
            for (int tp = 0; tp < np_; tp++) {
                u64 diff = rp[tp] >= cm ? rp[tp] - cm : rp[tp] + (gp - cm);
                if (diff % gg) continue;                  /* incompatible types */
                u64 k = g2 > 1 ? mulmod64((diff / gg) % g2, invA, g2) : 0;
                if (detA) {                               /* at most one n < X: cq + Lq k, if it is below X */
                    if ((u128)k <= (X - 1 - cq[tq]) / Lq) single(c, cq[tq] + Lq * (u128)k, p, lo);
                    continue;
                }
                u128 c1 = cq[tq] + Lq * (u128)k;          /* < L1 <= X - 1 */
                u64 k2 = mulmod64((p - (u64)(c1 % p)) % p, invB, p);
                if (detB) {
                    if ((u128)k2 <= (X - 1 - c1) / L1) single(c, c1 + L1 * (u128)k2, p, lo);
                    continue;
                }
                walk(c, c1 + L1 * (u128)k2, L2, lo);
            }
        }
    }
}

/* ---------------------------------------------------------------- shards */
static size_t NSHARDS = 4096, *SH_LO, *SH_HI;
static uint8_t *DONE_BIT;
static atomic_size_t NEXT = 0;
static Ctx TOTAL; static pthread_mutex_t tot_mu = PTHREAD_MUTEX_INITIALIZER;
static atomic_ullong DONE_NOW = 0, G_PAIRS = 0;
static u64 DONE_PRIOR = 0;

static void plan(void) {                     /* weight(q) ~ pairs + estimated candidates; cut into NSHARDS equal parts */
    double *w = malloc(NPR * sizeof(double)), S = 0, tot = 0, Xd = (double)X;
    for (size_t i = 0; i < NPR; i++) {
        double q = PR[i], g = (q * q - 1) / DD;
        w[i] = (double)i + 1 + 3000.0 * Xd * S / (q * g);
        S += DD / (q * (q * q - 1));
        tot += w[i];
    }
    SH_LO = malloc(NSHARDS * sizeof(size_t)); SH_HI = malloc(NSHARDS * sizeof(size_t));
    size_t s = 0, i = 0; double acc = 0;
    for (s = 0; s < NSHARDS; s++) {
        SH_LO[s] = i;
        double target = tot * (double)(s + 1) / (double)NSHARDS;
        while (i < NPR && (acc + w[i] <= target || s == NSHARDS - 1)) acc += w[i++];
        if (i == SH_LO[s] && i < NPR) acc += w[i++];          /* at least one q per shard while any remain */
        SH_HI[s] = i;
    }
    free(w);
}
static double now(void) { struct timespec ts; clock_gettime(CLOCK_REALTIME, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }
static void on_signal(int sig) { (void)sig; atomic_store(&STOP, 1); }

static void *worker(void *arg) {
    (void)arg;
    for (;;) {
        if (atomic_load(&STOP) || atomic_load(&DRAIN)) break;
        size_t s = atomic_fetch_add(&NEXT, 1);
        if (s >= NSHARDS) break;
        if (DONE_BIT[s]) continue;                        /* an empty shard is still recorded, so `complete` can be reached */
        Ctx c = {0};
        for (size_t iq = SH_LO[s]; iq < SH_HI[s]; iq++) {
            if (atomic_load(&STOP)) return NULL;
            do_q(&c, iq);
            atomic_fetch_add(&G_PAIRS, iq);
        }
        if (atomic_load(&STOP)) return NULL;              /* never record a shard that a stop may have cut short */
        pthread_mutex_lock(&out_mu);
        fflush(OUT);
        fprintf(DONE, "D %zu %llu %llu %llu %llu %llu %llu\n", s, (unsigned long long)c.pairs, (unsigned long long)c.skipped,
                (unsigned long long)c.classes, (unsigned long long)c.candidates, (unsigned long long)c.tests,
                (unsigned long long)c.passes);
        fflush(DONE);
        pthread_mutex_unlock(&out_mu);
        pthread_mutex_lock(&tot_mu);
        TOTAL.pairs += c.pairs; TOTAL.skipped += c.skipped; TOTAL.classes += c.classes;
        TOTAL.candidates += c.candidates; TOTAL.tests += c.tests; TOTAL.passes += c.passes; TOTAL.sieved += c.sieved;
        pthread_mutex_unlock(&tot_mu);
        atomic_fetch_add(&DONE_NOW, 1);
    }
    return NULL;
}
static void *monitor(void *arg) {
    double t0 = *(double *)arg, last = t0;
    for (;;) {
        sleep(1);
        if (access(STOP_FILE, F_OK) == 0) atomic_store(&STOP, 1);
        if (access(DRAIN_FILE, F_OK) == 0 || (DEADLINE > 0 && now() >= DEADLINE)) atomic_store(&DRAIN, 1);
        if (PROGRESS_SEC > 0 && now() - last >= PROGRESS_SEC) {
            last = now();
            fprintf(stderr, "[%.0f s] shards %llu+%llu/%zu  pairs %.3e\n", last - t0, (unsigned long long)DONE_PRIOR,
                    (unsigned long long)atomic_load(&DONE_NOW), NSHARDS, (double)atomic_load(&G_PAIRS));
        }
    }
    return NULL;
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "usage: cn2pair X T QSTAR RUNDIR [threads=30] [deadline_epoch=0] [nshards=4096] [progress_sec=60] [count_only=0]\n");
        return 2;
    }
    X = parse_u128(argv[1]); T = strtoull(argv[2], NULL, 10); QSTAR = strtoull(argv[3], NULL, 10);
    const char *dir = argv[4];
    int threads = argc > 5 ? atoi(argv[5]) : 30;
    DEADLINE = argc > 6 ? atof(argv[6]) : 0;
    if (argc > 7) NSHARDS = strtoull(argv[7], NULL, 10);
    if (argc > 8) PROGRESS_SEC = atof(argv[8]);
    if (argc > 9) COUNT_ONLY = atoi(argv[9]);
    const char *md = getenv("CN2X_MODE");
    if (!md || !strcmp(md, "rigid")) MODE = 0; else if (!strcmp(md, "howe")) MODE = 1;
    else if (!strcmp(md, "cheb")) { MODE = 2; DD = 2; } else { fprintf(stderr, "CN2X_MODE must be rigid, howe or cheb\n"); return 2; }
    if ((X >> 126) != 0) { fprintf(stderr, "X must be below 2^126\n"); return 2; }
    if (T < 5 || QSTAR <= T || QSTAR > 0xFFFFFFFFull) { fprintf(stderr, "need 5 <= T < QSTAR < 2^32\n"); return 2; }

    /* primes in (T, QSTAR] */
    u64 lim = QSTAR;
    uint8_t *comp = calloc(lim + 1, 1);
    for (u64 i = 2; i * i <= lim; i++) if (!comp[i]) for (u64 j = i * i; j <= lim; j += i) comp[j] = 1;
    NPR = 0; for (u64 i = T + 1; i <= lim; i++) if (i >= 2 && !comp[i]) NPR++;
    PR = malloc(NPR * sizeof(u32)); NPR = 0;
    for (u64 i = T + 1; i <= lim; i++) if (i >= 2 && !comp[i]) PR[NPR++] = (u32)i;
    free(comp);
    plan();
    sieve_init();

    char path[1200];
    snprintf(STOP_FILE, sizeof STOP_FILE, "%s/STOP", dir); snprintf(DRAIN_FILE, sizeof DRAIN_FILE, "%s/DRAIN", dir);
    remove(STOP_FILE); remove(DRAIN_FILE);
    DONE_BIT = calloc(NSHARDS, 1);
    snprintf(path, sizeof path, "%s/done.log", dir);
    if (!COUNT_ONLY) {
        FILE *f = fopen(path, "r");
        if (f) {
            char line[512]; unsigned long long id, v[6];
            while (fgets(line, sizeof line, f))
                if (sscanf(line, "D %llu %llu %llu %llu %llu %llu %llu", &id, &v[0], &v[1], &v[2], &v[3], &v[4], &v[5]) == 7
                    && id < NSHARDS && !DONE_BIT[id]) { DONE_BIT[id] = 1; DONE_PRIOR++; }
            fclose(f);
        }
        DONE = fopen(path, "a");
        snprintf(path, sizeof path, "%s/out.txt", dir);
        OUT = fopen(path, "a");
        if (!DONE || !OUT) { fprintf(stderr, "cannot open run files in %s\n", dir); return 2; }
    } else {
        DONE = fopen("/dev/null", "w"); OUT = fopen("/dev/null", "w");
    }
    signal(SIGINT, on_signal); signal(SIGTERM, on_signal);

    double t0 = now();
    pthread_t mon, *th = malloc(threads * sizeof(pthread_t));
    pthread_create(&mon, NULL, monitor, &t0);
    for (int i = 0; i < threads; i++) pthread_create(&th[i], NULL, worker, NULL);
    for (int i = 0; i < threads; i++) pthread_join(th[i], NULL);
    u64 done_total = DONE_PRIOR + atomic_load(&DONE_NOW);
    int complete = done_total == NSHARDS && !atomic_load(&STOP);
    char xs[64]; u128_str(X, xs);
    printf("{\"engine\":\"cn2pair\",\"mode\":%d,\"X\":\"%s\",\"T\":%llu,\"QSTAR\":%llu,\"nprimes\":%zu,\"nshards\":%zu,"
           "\"count_only\":%d,\"complete\":%s,\"stopped\":%s,\"drained\":%s,\"shards_done_prior\":%llu,\"shards_done_now\":%llu,"
           "\"pairs\":%llu,\"skipped\":%llu,\"classes\":%llu,\"candidates\":%llu,\"sieved\":%llu,\"tests\":%llu,\"passes\":%llu,"
           "\"seconds\":%.3f}\n",
           MODE, xs, (unsigned long long)T, (unsigned long long)QSTAR, NPR, NSHARDS, COUNT_ONLY,
           complete ? "true" : "false", atomic_load(&STOP) ? "true" : "false", atomic_load(&DRAIN) ? "true" : "false",
           (unsigned long long)DONE_PRIOR, (unsigned long long)atomic_load(&DONE_NOW),
           (unsigned long long)TOTAL.pairs, (unsigned long long)TOTAL.skipped, (unsigned long long)TOTAL.classes,
           (unsigned long long)TOTAL.candidates, (unsigned long long)TOTAL.sieved, (unsigned long long)TOTAL.tests, (unsigned long long)TOTAL.passes, now() - t0);
    return 0;
}
