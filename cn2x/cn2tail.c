/*
 * cn2tail.c -- the large-prime half of the split exhaustive CN2 search.
 *
 * Identity.  If p | n and n is a rigid CN2, then n/p == p (mod p^2-1) with n/p > p (Lem1), so
 *     n = p * m,   m = p + a (p^2 - 1),   n - 1 = (p^2 - 1)(a p + 1),   a >= 1.
 * Every rigid CN2 n < X with SOME prime factor p > R0 is therefore of this form for a pair (p, a)
 * with R0 < p <= X^(1/3) and 1 <= a <= A(p) = floor((X - 1 - p^2) / (p (p^2 - 1))).  There are about
 * X / (2 R0^2 ln R0) pairs.  This engine tests all of them; the capped DFS (cn2xc with
 * CN2X_TCAP=R0) covers the CN2s whose primes are all <= R0.  Together: every rigid CN2 below X.
 *
 * Per pair, only necessary conditions are applied (a CN2 cannot fail any of them):
 *   - p | a       => p | m => p^2 | n            (not squarefree)
 *   - s | m, s small prime, s^2 - 1 does not divide n - 1       (Korselt for s fails)
 *   - s^2 | m                                                    (not squarefree)
 *   - 2^(n-1) != 1 (mod m)       (a CN2 is Carmichael, hence a base-2 Fermat pseudoprime mod n)
 * Survivors are written as "P <p> <a>" and factored by the Python driver.
 *
 * Checkpointing as in cn2xc: deterministic shards (depend only on X, R0, W), done.log with counters,
 * restart skips finished shards; RUNDIR/DRAIN or the deadline = finish in-flight shards and exit;
 * RUNDIR/STOP or SIGINT/SIGTERM = abandon in-flight shards (redone next session).
 *
 * Modes (env CN2X_MODE, as for cn2xh): with g(p) = (p^2-1)/D the type of p is n == 1 or p (mod g(p)),
 * i.e. m = n/p == p or 1 (mod g(p)).  Family R: m = p + a g (all modes; the rigid family when D = 1);
 * family N: m = 1 + a g (howe: D = 1, cheb: D = 2).  Survivors: "P p a" (family R), "Q p a" (family N).
 * In rigid mode (default) this engine is unchanged: same pairs, same sieve, same survivors.
 *
 * Build:  clang -O3 -mcpu=native -o cn2tail cn2tail.c -lpthread -lm
 * Usage:  cn2tail X R0 RUNDIR [threads=30] [deadline_epoch=0] [W=268435456] [progress_sec=60] [nsmall=95]
 *         cn2tail --bench X R0 [seconds=20] [nsmall=95]        (single-thread rate at the heaviest primes)
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
static u64 R0, T;
static u32 *PR;                 /* primes in (R0, T] */
static size_t NPR;

static const u64 SMALLP[] = {5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509};
static int NSMALL = 95;
static int MODE = 0;            /* 0 rigid, 1 howe, 2 cheb */
static u64 DD = 1;

/* ---------------------------------------------------------------- u128 helpers (as in cn2xc.c) */
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
static u128 icbrt128(u128 n) {
    if (n == 0) return 0;
    u128 r = (u128)cbrt((double)n);
    while (r * r * r > n) r--;
    while ((r + 1) * (r + 1) * (r + 1) <= n) r++;
    return r;
}
static u64 modinv64(u64 a, u64 mod) {             /* gcd(a, mod) = 1, mod > 1 */
    i128 x0 = 0, x1 = 1; u64 r0 = mod, r1 = a % mod;
    while (r1) { u64 q = r0 / r1, r2 = r0 - q * r1; r0 = r1; r1 = r2; i128 x2 = x0 - (i128)q * x1; x0 = x1; x1 = x2; }
    if (x0 < 0) x0 += mod;
    return (u64)x0;
}

/* ---------------------------------------------------------------- Montgomery Fermat, base 2 (verbatim from cn2xc.c) */
static inline int topbit128(u128 e) {
    u64 hi = (u64)(e >> 64);
    return hi ? 127 - __builtin_clzll(hi) : 63 - __builtin_clzll((u64)e);
}
static inline u64 inv_mod_2_64(u64 n) {
    u64 x = n;
    for (int i = 0; i < 5; i++) x *= 2 - n * x;
    return x;
}
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
static int fermat2(u128 s, u128 e) {               /* 1 iff 2^e == 1 (mod s); s odd, 1 < s < 2^126 */
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
        u128 one = ((u128)0 - N) % N, x = one;
        for (int i = top; i >= 0; i--) {
            x = montmul2(x, x, N, np);
            if ((e >> i) & 1) { x += x; if (x >= N) x -= N; }
        }
        return x == one;
    }
}

/* ---------------------------------------------------------------- the per-prime scan */
typedef struct { u64 tests, sieved, fermat_calls, fermat_passes; } Cnt;

static atomic_int STOP = 0, DRAIN = 0;
static FILE *OUT;
static pthread_mutex_t out_mu = PTHREAD_MUTEX_INITIALIZER;
static atomic_ullong G_TESTS = 0, G_PASS = 0;

static u64 amax_f(u64 p, u64 r0) {                  /* largest a with p (r0 + a g) <= X - 1, g = (p^2-1)/DD */
    u128 M = (X - 1) / p, g = ((u128)p * p - 1) / DD;
    if (M < r0) return 0;
    u128 a = (M - r0) / g;
    return a > (u128)UINT64_MAX ? UINT64_MAX : (u64)a;
}
static u64 amax(u64 p) { return amax_f(p, MODE ? 1 : p); }          /* the larger family's range */
static u64 pairs_of(u64 p) { return amax_f(p, p) + (MODE ? amax_f(p, 1) : 0); }

static void emit(u64 p, u64 a, u64 r0) {
    char line[64];
    snprintf(line, sizeof line, "%c %llu %llu\n", r0 == 1 ? 'Q' : 'P', (unsigned long long)p, (unsigned long long)a);
    pthread_mutex_lock(&out_mu); fputs(line, OUT); pthread_mutex_unlock(&out_mu);
}

/* Korselt/squarefree check for one small prime s already known to divide m */
static inline int small_ok(u64 s, u64 p, u128 m) {
    u64 gs = (s * s - 1) / DD;
    u64 t = (p % gs) * mod128_64(m, gs) % gs;                     /* n mod g(s); gs < 2^18 */
    if (t != 1 % gs && !(MODE && t == s % gs)) return 0;
    return mod128_64(m, s * s) != 0;
}

#define BLK 4096
/* scan a in [lo, hi] for prime p; returns 0 if interrupted by STOP */
static int scan(u64 p, u64 r0, u64 lo, u64 hi, Cnt *c, int emit_on) {
    u64 q = (p * p - 1) / DD;                        /* g(p); p <= 2^32, so q < 2^64 */
    if (hi > amax_f(p, r0)) hi = amax_f(p, r0);
    if (lo > hi) return 1;
    int use_sieve = (hi - lo + 1) >= 64 && NSMALL > 0;
    u64 as[95]; int ns = 0; u64 sv[95];
    if (use_sieve) {
        for (int j = 0; j < NSMALL; j++) {
            u64 s = SMALLP[j];
            if (q % s == 0) continue;                /* then m == p (mod s): s never divides m */
            u64 inv = modinv64(q % s, s);
            sv[ns] = s; as[ns] = (s - (r0 % s) * inv % s) % s;  /* s | m  <=>  a == -r0 g^{-1} (mod s) */
            ns++;
        }
    }
    static __thread uint8_t rej[BLK];
    for (u64 b = lo; b <= hi; b += BLK) {
        u64 len = hi - b + 1 < BLK ? hi - b + 1 : BLK;
        if (atomic_load_explicit(&STOP, memory_order_relaxed)) return 0;
        memset(rej, 0, len);
        { u64 ap = (p - (r0 % p) * modinv64(q % p, p) % p) % p;   /* p | m  <=>  a == -r0 g^{-1} (mod p) */
          u64 first = b + ((ap + p - b % p) % p);
          for (u64 a = first; a < b + len; a += p) rej[a - b] = 1; }
        if (use_sieve) {
            for (int j = 0; j < ns; j++) {
                u64 s = sv[j], r = b % s, first = b + ((as[j] + s - r) % s);
                for (u64 a = first; a < b + len; a += s) {
                    if (rej[a - b]) continue;
                    if (!small_ok(s, p, (u128)r0 + (u128)a * q)) rej[a - b] = 1;
                }
            }
        }
        for (u64 i = 0; i < len; i++) {
            u64 a = b + i;
            c->tests++;
            if (rej[i]) { c->sieved++; continue; }
            u128 m = (u128)r0 + (u128)a * q;
            if (!use_sieve) {                        /* short ranges: same conditions, checked directly */
                int bad = 0;
                for (int j = 0; j < NSMALL && !bad; j++) {
                    u64 s = SMALLP[j];
                    if (mod128_64(m, s) == 0 && !small_ok(s, p, m)) bad = 1;
                }
                if (bad) { c->sieved++; continue; }
            }
            u128 e = (u128)p * m - 1;                /* n - 1 */
            c->fermat_calls++;
            if (MODE == 3 || fermat2(m, e)) { c->fermat_passes++; if (emit_on) emit(p, a, r0); }
        }
    }
    return 1;
}

/* ---------------------------------------------------------------- shards */
typedef struct { u32 i_lo, i_hi; u64 a_lo, a_hi; u64 work; } Shard;   /* a_hi == 0: full range of each prime */
static Shard *SH; static size_t NSH = 0, SHCAP = 0;
static u64 W = 1ull << 28;
static u128 TOTAL_TESTS = 0;

static void add_shard(u32 lo, u32 hi, u64 alo, u64 ahi, u64 work) {
    if (NSH == SHCAP) { SHCAP = SHCAP ? 2 * SHCAP : 4096; SH = realloc(SH, SHCAP * sizeof(Shard)); }
    SH[NSH++] = (Shard){lo, hi, alo, ahi, work};
}
static void make_shards(void) {
    u64 acc = 0; long cur = -1;
    for (size_t i = 0; i < NPR; i++) {
        u64 A = amax(PR[i]);
        if (A == 0) break;                           /* A(p) is non-increasing in p */
        TOTAL_TESTS += pairs_of(PR[i]);
        if (A >= W) {
            if (cur >= 0) { add_shard((u32)cur, (u32)i, 1, 0, acc); cur = -1; acc = 0; }
            for (u64 lo = 1; lo <= A; lo += W) {
                u64 hi = lo + W - 1 < A ? lo + W - 1 : A;
                add_shard((u32)i, (u32)i + 1, lo, hi, hi - lo + 1);
                if (hi == A) break;
            }
            continue;
        }
        if (cur < 0) cur = (long)i;
        acc += pairs_of(PR[i]);
        if (acc >= W) { add_shard((u32)cur, (u32)i + 1, 1, 0, acc); cur = -1; acc = 0; }
    }
    if (cur >= 0) {                                  /* trailing run: primes cur .. last with A > 0 */
        size_t end = (size_t)cur;
        while (end < NPR && amax(PR[end]) > 0) end++;
        add_shard((u32)cur, (u32)end, 1, 0, acc);
    }
}

/* done bookkeeping */
static uint8_t *DONE_BITS; static Cnt DONE_PRIOR; static u64 DONE_PRIOR_N = 0, DONE_NOW_N = 0, DONE_PRIOR_WORK = 0;
static FILE *DONE;
static atomic_ullong NEXT = 0;
static atomic_ullong WORK_NOW = 0;

static void load_done(const char *path) {
    DONE_BITS = calloc(NSH / 8 + 1, 1);
    FILE *f = fopen(path, "r");
    if (!f) return;
    char line[256]; unsigned long long id, v[4];
    while (fgets(line, sizeof line, f)) {
        if (sscanf(line, "D %llu %llu %llu %llu %llu", &id, &v[0], &v[1], &v[2], &v[3]) != 5) continue;
        if (id >= NSH) continue;
        if ((DONE_BITS[id >> 3] >> (id & 7)) & 1) continue;
        DONE_BITS[id >> 3] |= (uint8_t)(1 << (id & 7));
        DONE_PRIOR.tests += v[0]; DONE_PRIOR.sieved += v[1]; DONE_PRIOR.fermat_calls += v[2]; DONE_PRIOR.fermat_passes += v[3];
        DONE_PRIOR_N++; DONE_PRIOR_WORK += SH[id].work;
    }
    fclose(f);
}

static void *worker(void *arg) {
    (void)arg;
    for (;;) {
        if (atomic_load(&STOP) || atomic_load(&DRAIN)) break;
        u64 id = atomic_fetch_add(&NEXT, 1);
        if (id >= NSH) break;
        if ((DONE_BITS[id >> 3] >> (id & 7)) & 1) continue;
        Shard *s = &SH[id];
        Cnt c = {0}; int ok = 1;
        for (u32 i = s->i_lo; i < s->i_hi && ok; i++) {
            u64 p = PR[i];
            u64 lo = s->a_hi ? s->a_lo : 1, hi = s->a_hi ? s->a_hi : amax(p);
            Cnt before = c;
            ok = scan(p, p, lo, hi, &c, 1);
            if (ok && MODE) ok = scan(p, 1, lo, hi, &c, 1);
            atomic_fetch_add(&G_TESTS, c.tests - before.tests);
            atomic_fetch_add(&G_PASS, c.fermat_passes - before.fermat_passes);
        }
        if (!ok) break;                              /* STOP mid-shard: abandoned, redone next session */
        pthread_mutex_lock(&out_mu);
        fflush(OUT);
        fprintf(DONE, "D %llu %llu %llu %llu %llu\n", (unsigned long long)id, (unsigned long long)c.tests,
                (unsigned long long)c.sieved, (unsigned long long)c.fermat_calls, (unsigned long long)c.fermat_passes);
        fflush(DONE);
        DONE_NOW_N++;
        pthread_mutex_unlock(&out_mu);
        atomic_fetch_add(&WORK_NOW, s->work);
    }
    return NULL;
}

/* ---------------------------------------------------------------- monitor / stop sources */
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + ts.tv_nsec * 1e-9; }
static double DEADLINE = 0, PROGRESS_SEC = 60;
static char STOP_FILE[1100], DRAIN_FILE[1100];
static atomic_int FINISHED = 0;
static void on_signal(int sig) { (void)sig; atomic_store(&STOP, 1); }
static void *monitor(void *arg) {
    double t0 = *(double *)arg, tl = t0;
    while (!atomic_load(&FINISHED)) {
        struct timespec ts = {1, 0}; nanosleep(&ts, NULL);
        if (PROGRESS_SEC > 0 && now() - tl >= PROGRESS_SEC) {
            double el = now() - t0;
            double done = (double)DONE_PRIOR.tests + (double)atomic_load(&G_TESTS);
            double frac = done / (double)TOTAL_TESTS, rate = atomic_load(&G_TESTS) / (el > 0 ? el : 1);
            double eta = rate > 0 ? ((double)TOTAL_TESTS - done) / rate : 0;
            time_t wall = time(NULL); struct tm tmv; localtime_r(&wall, &tmv);
            fprintf(stderr, "[cn2tail %02d:%02d:%02d] elapsed %6.0fs  shards done %llu+%llu/%zu  tests %.3e/%.3e (%.2f%%, %.2e/s)  "
                            "survivors %llu  ETA %.2f h%s\n",
                    tmv.tm_hour, tmv.tm_min, tmv.tm_sec, el, (unsigned long long)DONE_PRIOR_N, (unsigned long long)DONE_NOW_N, NSH,
                    done, (double)TOTAL_TESTS, 100 * frac, rate,
                    (unsigned long long)(DONE_PRIOR.fermat_passes + atomic_load(&G_PASS)), eta / 3600,
                    atomic_load(&DRAIN) ? "  DRAINING" : "");
            fflush(stderr); tl = now();
        }
        if (DEADLINE > 0 && (double)time(NULL) >= DEADLINE) atomic_store(&DRAIN, 1);
        FILE *f;
        if ((f = fopen(STOP_FILE, "r"))) { fclose(f); atomic_store(&STOP, 1); }
        if ((f = fopen(DRAIN_FILE, "r"))) { fclose(f); atomic_store(&DRAIN, 1); }
        pthread_mutex_lock(&out_mu); fflush(OUT); fflush(DONE); fsync(fileno(OUT)); fsync(fileno(DONE)); pthread_mutex_unlock(&out_mu);
    }
    return NULL;
}

static void sieve_range(u64 lo, u64 hi) {           /* primes p with lo < p <= hi */
    uint8_t *comp = calloc(hi + 1, 1);
    for (u64 p = 2; p * p <= hi; p++) if (!comp[p]) for (u64 q = p * p; q <= hi; q += p) comp[q] = 1;
    size_t cap = 1 << 20; PR = malloc(cap * sizeof(u32)); NPR = 0;
    for (u64 p = (lo < 4 ? 5 : lo + 1); p <= hi; p++) if (!comp[p]) {
        if (NPR == cap) { cap *= 2; PR = realloc(PR, cap * sizeof(u32)); }
        PR[NPR++] = (u32)p;
    }
    free(comp);
}

static int set_mode(void) {
    const char *md = getenv("CN2X_MODE");
    if (md && !strcmp(md, "howe")) MODE = 1;
    else if (md && !strcmp(md, "cheb")) { MODE = 2; DD = 2; }
    else if (md && !strcmp(md, "test")) {    /* VALIDATION ONLY: modulus (p^2-1)/CN2X_D, Fermat off */
        MODE = 3; const char *dd = getenv("CN2X_D"); DD = dd ? strtoull(dd, NULL, 10) : 24;
        if (24 % DD) { fprintf(stderr, "CN2X_D must divide 24\n"); return 1; }
    }
    else if (md && strcmp(md, "rigid")) { fprintf(stderr, "CN2X_MODE must be rigid, howe, cheb or test\n"); return 1; }
    return 0;
}

static int bench(int argc, char **argv) {
    X = parse_u128(argv[2]); R0 = strtoull(argv[3], NULL, 10);
    double secs = argc > 4 ? atof(argv[4]) : 20;
    if (argc > 5) NSMALL = atoi(argv[5]);
    if (set_mode()) return 2;
    T = (u64)icbrt128(DD * (X - 1));
    sieve_range(R0, R0 + 200000 < T ? R0 + 200000 : T);
    OUT = fopen("/dev/null", "w");
    Cnt c = {0}; double t0 = now();
    for (size_t i = 0; i < NPR && now() - t0 < secs; i++) {
        u64 p = PR[i], A = amax(p);
        if (!A) break;
        u64 hi = A < (1u << 22) ? A : (1u << 22);
        scan(p, p, 1, hi, &c, 0);
        if (MODE) scan(p, 1, 1, hi, &c, 0);
    }
    double el = now() - t0;
    printf("{\"tests\":%llu,\"sieved\":%llu,\"fermat_calls\":%llu,\"fermat_passes\":%llu,\"seconds\":%.3f,"
           "\"ns_per_test\":%.1f,\"sieved_frac\":%.4f}\n",
           (unsigned long long)c.tests, (unsigned long long)c.sieved, (unsigned long long)c.fermat_calls,
           (unsigned long long)c.fermat_passes, el, 1e9 * el / (double)c.tests, (double)c.sieved / (double)c.tests);
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 4 && !strcmp(argv[1], "--bench")) return bench(argc, argv);
    if (argc < 4) {
        fprintf(stderr, "usage: cn2tail X R0 RUNDIR [threads=30] [deadline_epoch=0] [W=268435456] [progress_sec=60] [nsmall=95]\n"
                        "       cn2tail --bench X R0 [seconds=20] [nsmall=95]\n");
        return 2;
    }
    X = parse_u128(argv[1]); R0 = strtoull(argv[2], NULL, 10);
    const char *dir = argv[3];
    int threads = argc > 4 ? atoi(argv[4]) : 30;
    DEADLINE = argc > 5 ? atof(argv[5]) : 0;
    if (argc > 6) W = strtoull(argv[6], NULL, 10);
    if (argc > 7) PROGRESS_SEC = atof(argv[7]);
    if (argc > 8) NSMALL = atoi(argv[8]);
    if (NSMALL < 0 || NSMALL > 95) { fprintf(stderr, "nsmall must be 0..95\n"); return 2; }

    if (set_mode()) return 2;
    T = (u64)icbrt128(DD * (X - 1));                  /* n/p > g(p) => p^3 < D X */
    if ((X >> 126) != 0 || T > 0xFFFFFFFFull) { fprintf(stderr, "X too large for this engine (X < 2^126 and primes < 2^32 required)\n"); return 2; }
    double t0 = now();
    sieve_range(R0, T);
    make_shards();
    double t_setup = now() - t0;

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
    signal(SIGINT, on_signal); signal(SIGTERM, on_signal);

    fprintf(stderr, "[cn2tail] primes in (%llu, %llu]: %zu   shards %zu   total (p,a) pairs %.4e   setup %.1fs   prior shards done %llu\n",
            (unsigned long long)R0, (unsigned long long)T, NPR, NSH, (double)TOTAL_TESTS, t_setup, (unsigned long long)DONE_PRIOR_N);
    double ts = now();
    pthread_t mon; pthread_create(&mon, NULL, monitor, &ts);
    pthread_t *th = malloc(threads * sizeof(pthread_t));
    for (int i = 0; i < threads; i++) pthread_create(&th[i], NULL, worker, NULL);
    for (int i = 0; i < threads; i++) pthread_join(th[i], NULL);
    atomic_store(&FINISHED, 1);
    pthread_join(mon, NULL);
    double t_search = now() - ts;

    int complete = (DONE_PRIOR_N + DONE_NOW_N == NSH) && !atomic_load(&STOP);
    fflush(OUT); fflush(DONE); fsync(fileno(OUT)); fsync(fileno(DONE)); fclose(OUT); fclose(DONE);

    /* totals from done.log (exact across sessions) */
    free(DONE_BITS); memset(&DONE_PRIOR, 0, sizeof DONE_PRIOR); DONE_PRIOR_N = 0; DONE_PRIOR_WORK = 0;
    snprintf(path, sizeof path, "%s/done.log", dir);
    load_done(path);
    char xs[64]; u128_str(X, xs);
    char tt[64]; u128_str(TOTAL_TESTS, tt);
    printf("{\"mode\":%d,\"X\":\"%s\",\"R0\":%llu,\"T\":%llu,\"complete\":%s,\"stopped\":%s,\"drained\":%s,\"threads\":%d,"
           "\"shards\":%zu,\"shards_done\":%llu,\"shards_done_now\":%llu,\"total_pairs\":\"%s\",\"tests\":%llu,\"sieved\":%llu,"
           "\"fermat_calls\":%llu,\"fermat_passes\":%llu,\"setup_seconds\":%.3f,\"search_seconds\":%.3f}\n",
           MODE, xs, (unsigned long long)R0, (unsigned long long)T, complete ? "true" : "false",
           atomic_load(&STOP) ? "true" : "false", atomic_load(&DRAIN) ? "true" : "false", threads,
           NSH, (unsigned long long)DONE_PRIOR_N, (unsigned long long)DONE_NOW_N, tt,
           (unsigned long long)DONE_PRIOR.tests, (unsigned long long)DONE_PRIOR.sieved,
           (unsigned long long)DONE_PRIOR.fermat_calls, (unsigned long long)DONE_PRIOR.fermat_passes, t_setup, t_search);
    return complete ? 0 : 3;
}
