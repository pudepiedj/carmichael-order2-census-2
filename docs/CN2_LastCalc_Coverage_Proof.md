---
location: post
title: LASTCALC Coverage Theorem
author: Claude Opus 5.5 and JCP
date: 2026-10-05
---

# What a capped search with LASTCALC provably finds

**Date:** 5 October 2026
**Purpose:** first step towards *certifying* order-2 searches without an unbounded run. A capped search (table of primes up to $T$) is fast but has so far been only *conditionally* complete. This note proves exactly which numbers it is guaranteed to find when `CN2X_LASTCALC=1`. The rest must be covered by a separate, provable "tail", which is the subject of later steps.
**Code proved about:** `cn2x/cn2xh.c` (functions `node_head`, `child_step`, `merge`, `last_calc`, `list_mode`, `dfs`, the shard generator `gen_node`/`next_shard`, `run_task`, `donate`), and the verifier `cn2_order2_defs.verify_candidate`.

---

## 1. The statement

**Theorem.** Let $X < 2^{126}$ and $T \ge 5$. Run `cn2xh` on $X$ with `CN2X_TCAP` $= T$, `CN2X_LASTCALC=1`, `CN2X_NMIN` unset, any mode (rigid, howe or cheb), any switch ratio $> 0$, and any sieve setting, until a session reports `complete: true`. Let $n < X$ be an order-2 number of that mode, with primes $p_1 < p_2 < \dots < p_k$. **If $p_{k-1} \le T$**, so that every prime of $n$ except possibly the largest lies in the table, then the run's `out.txt` contains a record that the verifier turns into $n$. The record is either `H n`, or `S m pmax j s` with $m = p_1\cdots p_j$ and $m\thinspace s = n$.

So, with the tail defined accordingly, one capped run plus one tail gives a certified census:
$$\lbrace n < X \rbrace \;=\; \underbrace{\lbrace n : p_{k-1} \le T\rbrace}_{\text{capped run with LASTCALC}} \;\cup\; \underbrace{\lbrace n : p_{k-1} > T\rbrace}_{\text{numbers with at least two primes above } T}.$$

Two things the theorem does **not** depend on are worth noting.
- **The switch ratio.** Any ratio $> 0$ gives the same coverage. The ratio changes only the cost (§7 of the non-rigid note (not part of this repository)).
- **Where LIST happens to switch.** The proof follows $n$'s own path and shows that every possible decision at every node either emits $n$ or passes it to the next node.

**The existing runs do not meet the hypothesis.** The $10^{30}$ run (and the other capped runs) used `lastcalc: 0`. Without LASTCALC the guaranteed set shrinks to $p_k \le T$: every prime in the table. That run's extra finds came through LIST, which is not guaranteed. A certifying run must set `CN2X_LASTCALC=1`, and its cost has not yet been measured.

## 2. Notation

The mode fixes $D$ and the allowed types: $g(p) = (p^2-1)/D$, with $D = 1$ (rigid, howe) or $D = 2$ (cheb). An **order-2 number** is a squarefree composite $n$ with $n \equiv r_p \pmod{g(p)}$ for every prime $p \mid n$, where the type residue $r_p$ is $1$ (rigid mode) or one of $1, p$ (howe, cheb).

For $n$ with primes $p_1 < \dots < p_k$ write, for $0 \le j \le k$:
- $m_j = p_1 \cdots p_j$ (so $m_0 = 1$, $m_k = n$);
- $L_j = \mathrm{lcm}\big(g(p_1), \dots, g(p_j)\big)$ ($L_0 = 1$);
- $c_j$ = the residue of $n$ modulo $L_j$, with $0 \le c_j < L_j$.

The engine's node is a tuple $(m, p_{\max}, L, c, j)$. Its **path node for $n$** at depth $j$ is $N_j = (m_j, p_j, L_j, c_j, j)$, with $p_0 := 3$ as at the root.

## 3. Facts about every order-2 number

Throughout, $n < X$ is order-2 with primes $p_1 < \dots < p_k$.

**F1 (no 2, no 3).** $n$ is odd and $3 \nmid n$.
*Proof.* If $2 \mid n$, take an odd prime $q \mid n$ (one exists, since $n$ is squarefree and composite). Then $g(q)$ is even, because $q^2-1 \equiv 0 \pmod 8$ and $D \le 2$, and $n \equiv 1$ or $q \pmod{g(q)}$ is odd, a contradiction. If $3 \mid n$, take a prime $q \ge 5$ of $n$. Then $3 \mid g(q)$, so $n \equiv 1$ or $q \not\equiv 0 \pmod 3$, a contradiction. $\square$

**F2 (cofactor bound).** For every $p \mid n$: $n/p > g(p)$, and in fact $D\thinspace (n/p) > p^2$.
*Proof.* Since $g(p) \mid p^2-1$, $p^{-1} \equiv p \pmod{g(p)}$, so $n/p \equiv r_p\thinspace  p \pmod{g(p)}$. That is, $n/p \equiv p$ if $r_p = 1$, and $n/p \equiv p^2 \equiv 1$ if $r_p = p$. Now $n/p \ne 1$ ($n$ is composite) and $n/p \ne p$ ($n$ is squarefree). Since $g(p) > p$ for $p \ge 5$, the least positive members of the two classes other than $1$ and $p$ are $1 + g(p)$ and $p + g(p)$. So $n/p \ge 1 + g(p)$, which gives $n/p > g(p)$. For $D = 1$ this reads $n/p \ge p^2$, and equality is impossible because $p \nmid n/p$ ($n$ is squarefree), so $n/p > p^2$. For $D = 2$ it reads $2(n/p) \ge p^2 + 1 > p^2$. Either way $D\thinspace (n/p) > p^2$. $\square$

**F3 (at least three primes).** $k \ge 3$.
*Proof.* If $n = p_1 p_2$, F2 at $p_2$ gives $p_2^2 < D p_1 \le 2 p_1 < p_2^2$. $\square$

**F4 (role bounds).** Let $0 \le j \le k-1$, $t = p_{j+1}$ and $\mathrm{rem}_j = \lfloor (X-1)/m_j \rfloor$. Then:
1. if $t$ is the last prime ($j = k-1$, so $j \ge 2$): $t^2 \le D m_j - 1$, i.e. $t \le b_{\text{last}} := \lfloor\sqrt{D m_j - 1}\rfloor$;
2. if $t$ is second-to-last ($j = k-2$, so $j \ge 1$): $t \le b_{\text{second}} := \min\big(D m_j - 1,\ \lfloor\sqrt{\mathrm{rem}_j}\rfloor\big)$;
3. if $t$ is earlier ($j \le k-3$): $t \le b_{\text{early}} := \lfloor \mathrm{rem}_j^{1/3}\rfloor$.

*Proof.* (1) is F2 at $t$, with $n/t = m_j$. (2): $n/m_j = t\thinspace p_k > t^2$ and $n/m_j \le \mathrm{rem}\_j$ (an integer at most $(X-1)/m_j$), so $t^2 < \mathrm{rem}_j$. And F2 at $p_k$ gives $t^2 < p_k^2 < D\thinspace m_j t$, so $t < D m_j$. (3): $n/m_j \ge t\thinspace p_{j+2}\thinspace p_{j+3} > t^3$, so $t^3 < \mathrm{rem}_j$. $\square$

Consequently every next prime satisfies $t \le \max(b_{\text{last}}, b_{\text{second}}, b_{\text{early}})$, the engine's `mx`. A next prime that is **not** the last satisfies $t \le e_b := \max(b_{\text{second}}, b_{\text{early}})$, the engine's `eb`. So $t > e_b$ implies that $t$ is the last prime.

**F5 (admissibility).** For distinct primes $p, q \mid n$: $q \nmid g(p)$. Hence $\gcd(m_j, L_j) = 1$ for every $j$.
*Proof.* If $q \mid g(p)$, then $n \equiv r_p \pmod q$ with $r_p \in \lbrace 1, p\rbrace$, while $q \mid n$, so $q \mid 1$ or $q \mid p$. $\square$

**F6 (Fermat).** $2^{n-1} \equiv 1 \pmod n$, so $2^{n-1} \equiv 1 \pmod s$ for every $s \mid n$.
*Proof.* For each $p \mid n$, $p - 1 \mid g(p)$ (for $D = 2$, $g(p) = (p-1)\cdot\tfrac{p+1}{2}$). Both types give $n \equiv 1 \pmod{p-1}$: type 1 directly, and type $p$ because $n - 1 = (n - p) + (p - 1)$. So $n$ is squarefree with $p - 1 \mid n - 1$ throughout, which makes it a Carmichael number (Korselt). It is odd by F1. $\square$

**F7 (size of primes).** Every $p \mid n$ satisfies $p^3 < D n < DX$, so $p < (DX)^{1/3} < 2^{43}$ for $X < 2^{126}$. *Proof:* F2. $\square$

## 4. What the engine does at a node

The relevant code paths, with the line references in `cn2xh.c`, for a node $(m, p_{\max}, L, c, j)$:

- **Head** (`node_head`): if $j \ge 3$ and $m \bmod L = c$, emit `H m`. Let $\mathrm{rem} = \lfloor (X-1)/m\rfloor$; if $\mathrm{rem} \le p_{\max}$, stop. Compute the three bounds of F4 (with $b_{\text{last}} = 0$ for $j < 2$ and $b_{\text{second}} = 0$ for $j < 1$), $t_{\max} = \min(\mathrm{mx}, T)$, and $\mathrm{nch}$ = the number of table primes in $(p_{\max}, t_{\max}]$. Set $\mathrm{lc}$ = (LASTCALC and $j \ge 2$ and $b_{\text{last}} > T$).
  - If $\mathrm{nch} = 0$: run `last_calc` if lc, then stop.
  - Otherwise decide LIST ($L > 1$ and $\lfloor \mathrm{rem}/L\rfloor + 1 \le \text{ratio}\cdot\mathrm{nch}$) or branch. If branching and lc, run `last_calc` first.
- **LIST** (`list_mode`): for every $s \equiv c\thinspace m^{-1} \pmod L$ with $p_{\max} < s \le \mathrm{rem}$: discard $s$ if even, if $3 \mid s$, or (sieve) if some prime $q \le p_{\max}$ from the small-prime list divides $s$. Otherwise emit `S m pmax j s` if $2^{ms-1} \equiv 1 \pmod s$.
- **LASTCALC** (`last_calc`): for every $r \equiv c\thinspace m^{-1} \pmod L$ with $\max(T, p_{\max}) < r \le \min(b_{\text{last}}, \mathrm{rem})$, provided $L > 1$: discard even $r$ and $3 \mid r$, and emit `S m pmax j r` if $2^{mr-1} \equiv 1 \pmod r$.
- **Branch** (`dfs`, and identically the generator, batches and donation): for table primes $t$ in $(p_{\max}, t_{\max}]$ in increasing order, run `child_step`:
  - skip $t$ if $t \mid L$, or if $\gcd(g(t), m) \ne 1$;
  - if $m t > X$, stop the loop;
  - for each allowed type residue $\rho \in \lbrace 1, t\rbrace$ (one if they coincide mod $g(t)$; only $1$ in rigid mode), merge $n \equiv \rho \pmod{g(t)}$ into $(L, c)$ by the CRT (`merge`):
    - **incompatible**: skip this type;
    - **determined** ($L' = \mathrm{lcm}(L, g(t)) > X - 1$): let $n_0$ be the unique member of the class below $X$, if any. If $mt \mid n_0$ and $s = n_0/(mt)$: emit `H n0` if $s = 1$ and $j + 1 \ge 3$; skip if $1 < s \le t$; otherwise emit `S (mt) t (j+1) s` if $2^{n_0 - 1} \equiv 1 \pmod s$;
    - **open**: if $t > e_b$, emit `H (mt)` when $j + 1 \ge 3$ and $mt \equiv c' \pmod{L'}$; otherwise create the child node $(mt, t, L', c', j+1)$.

**Run completeness.** A session reports `complete` when the generator has been exhausted in that session, the task queue is empty, and there was no hard stop (`cn2xh.c:859`). Shards finished in earlier sessions are skipped, and their records are already in `out.txt`, because a shard's records are flushed before its `D` line is written (`record_done`). Work done by the generator itself, including `node_head` (and therefore `last_calc`) on generator nodes, is redone in every session. Donated children carry their parent's shard, and a shard is recorded only when all of its tasks have finished. So in a complete run, every node that the plain recursive `dfs` from the root would visit is visited, in some session, with its records in `out.txt`.

## 5. Proof of the theorem

Let $n$ satisfy the hypothesis. We show by induction on $j = 0, 1, \dots$ that either **(E)** a record for $n$ has been emitted, or **(V)** the path node $N_j$ is visited. (V) holds at $j = 0$: the root is $(1, 3, 1, 0, 0) = N_0$, and $p_1 \ge 5 > 3$ by F1.

**Step.** Suppose $N_j$ is visited, with $0 \le j \le k - 1$, and let $t = p_{j+1}$. First, the node survives the head. $\mathrm{rem}_j \ge n/m_j \ge t > p_j$. By F4, $t \le \mathrm{mx}$. And $\gcd(m_j, L_j) = 1$ (F5), so $c_j m_j^{-1}$ is defined and, because $n \equiv c_j \pmod{L_j}$, the cofactor $s_j = n/m_j$ satisfies $s_j \equiv c_j m_j^{-1} \pmod{L_j}$.

**Case A: $t \le T$.** Then $t \le t_{\max}$, so $t$ is one of the node's table children and $\mathrm{nch} \ge 1$.
- *A1: the node lists.* $s_j$ is in the listed progression: $s_j \equiv c_j m_j^{-1}$, and $p_j < t \le s_j \le \mathrm{rem}_j$. It survives every filter. It is odd and $3 \nmid s_j$ (F1). No prime $q \le p_j$ divides it, because its primes are $p_{j+1}, \dots, p_k > p_j$. And $2^{m_j s_j - 1} = 2^{n-1} \equiv 1 \pmod{s_j}$ (F6). So `S` $m_j\ p_j\ j\ s_j$ is emitted: (E).
- *A2: the node branches.* The loop reaches $t$: a stop needs some $t' \le t$ with $m_j t' > X$, but $m_j t \le n < X$. At $t$, `child_step` does not skip: $t \nmid L_j$ and $\gcd(g(t), m_j) = 1$ (F5). The type $\rho = r_t$ of $n$ at $t$ is compatible with $(L_j, c_j)$, since $n$ itself satisfies both congruences. So the merge for $\rho$ is not "incompatible", and the merged class is exactly $c_{j+1} \bmod L_{j+1}$, because the CRT solution class is unique.
  - *Determined.* The class mod $L_{j+1} > X-1$ has at most one member in $[0, X-1]$, and $n$ is one, so $n_0 = n$. Then $m_j t \mid n$, and $s = n/(m_j t)$ is either $1$ (so $j + 1 = k \ge 3$ and `H n` is emitted) or $s \ge p_{j+2} > t$. In the latter case $2^{n-1} \equiv 1 \pmod s$ (F6), and `S` is emitted: (E).
  - *Open, $t > e_b$.* By F4, $t$ is the last prime, so $m_j t = n$, $j + 1 = k \ge 3$ (F3), and $n \equiv c_{j+1} \pmod{L_{j+1}}$. `H n` is emitted: (E).
  - *Open, $t \le e_b$.* The child node $(m_j t, t, L_{j+1}, c_{j+1}, j+1) = N_{j+1}$ is created and, in a complete run, visited: (V) at $j + 1$.

**Case B: $t > T$.** By hypothesis $p_{k-1} \le T$, so $t$ must be the last prime: $j = k - 1 \ge 2$ (F3). By F4, $t \le b_{\text{last}}$, so $b_{\text{last}} > T$, and lc is true.
- If the node lists, then as in A1 the progression contains $s_j = t$ (with $p_j < t \le \mathrm{rem}_j$), and `S` is emitted: (E).
- If $\mathrm{nch} = 0$, or the node branches, `last_calc` runs. Its range is $\max(T, p_j) < r \le \min(b_{\text{last}}, \mathrm{rem}_j)$, and it contains $r = t$: $t > T$, $t > p_j$, $t \le b_{\text{last}}$, $t \le n/m_j \le \mathrm{rem}_j$. Also $t \equiv c_j m_j^{-1} \pmod{L_j}$ and $L_j \ge g(p_1) > 1$. $t$ is odd and $3 \nmid t$, and $2^{m_j t - 1} = 2^{n-1} \equiv 1 \pmod t$. So `S` $m_j\ p_j\ j\ t$ is emitted: (E).

**Termination.** Each (V) step increases $j$, and (V) at $j = k$ cannot be followed by a further step. When $N_k$ is visited, its head has $m = n$, $j = k \ge 3$ and $n \bmod L_k = c_k$ (since $0 \le c_k < L_k$), so `H n` is emitted. Hence (E) holds after at most $k + 1$ steps. $\square$

**Verification.** The record is then verified (`verify_candidate`): it factorises $m$ and $s$, checks that they share no prime, and checks the definition of order-2 numbers on the combined factorisation. For a genuine $n$, this returns $n$ provided the factorisation is correct. `sympy.factorint` uses the BPSW primality test. Its correctness below $2^{64}$ is established by exhaustive computation, not by a theorem: every base-2 Fermat pseudoprime below $2^{64}$ has been listed (Feitsma and Galway), and none of them passes BPSW. Every true prime of $n$ is below $2^{43}$ (F7). The only way verification could fail on a genuine $n$ is if a **composite** cofactor above $2^{64}$ were declared prime (a BPSW pseudoprime; none is known). The verifier therefore **flags every candidate in which some reported prime factor exceeds $2^{64}$**. It doesn't reject such candidates, because a survivor that is not order-2 can genuinely have a prime factor that large. If a run produces no flags, every factorisation it used involves only factors below $2^{64}$, where BPSW gives no false answers (by the computation above), and verification is exact. If there are flags, each flagged factor needs an independent primality proof before the run can certify anything.

## 6. Assumptions and what remains

The theorem is about the algorithm as written. It also relies on:
1. **Arithmetic.** 128-bit integer arithmetic without overflow for $X < 2^{126}$ (checked at start-up), and the Montgomery and `mulmod128` routines (covered by `cn2xh --selftest`).
2. **The prime table.** The table holds exactly the primes up to $T$ (the sieve).
3. **Factorisation.** It is correct. `verify_candidate` now flags (on stderr, prefix `UNPROVEN-PRIMALITY`) any factor above $2^{64}$, and a certifying run must show no flags or prove each flagged factor prime.
4. **The run.** It is complete in the sense of §4, and `out.txt` is intact (no torn `S`/`H` lines; a truncated last line could only occur after a crash and is re-produced when its shard is redone).

## 7. Empirical check against the proved census

Below $10^{25}$ the 14 Howe-mode terms are proved. So a capped run there with a deliberately tiny table must find **every** term whose second-largest prime is $\le T$, and may find others only by luck (through LIST). Howe mode, ratio 1, `cn2_ratio_experiment.py --X 1e25 --R0 300,800 --lastcalc`:

| $T$ | LASTCALC | time | found | guaranteed by the theorem | guaranteed but missed |
|---|---|---|---|---|---|
| 300 | off | 4 s | 4 / 14 | — | — |
| 300 | **on** | 4 s | 6 / 14 | 2 | **0** |
| 800 | off | 73 s | 9 / 14 | — | — |
| 800 | **on** | 74 s | 12 / 14 | 9 | **0** |

At $T = 300$ the term $4.43\times 10^{14}$ (largest prime 331, second-largest 167) is missed without LASTCALC and found with it, exactly as §5 Case B predicts. At $T = 800$ the two terms not found, $1.16\times 10^{24}$ and $4.11\times 10^{24}$, have second-largest primes 4,049 and 3,761, so the theorem says nothing about them. LASTCALC cost nothing measurable here: the same time, and 0.8% more list candidates at $T = 300$.

**A trap found on the way.** The top-level binary `cn2x/cn2xh` was built on 23 September, before LASTCALC was added to the source on the 24th, so it silently ignores `CN2X_LASTCALC=1`. The first attempt at this check used it and gave identical results with and without LASTCALC. The campaign scripts compile a fresh binary into each run directory, so production runs are unaffected. The experiment driver now does the same. **A certifying run should record the source hash and the `"lastcalc":1` field of its session record.**

**Still to do, towards a certified census:**
- **Measure what LASTCALC costs** at the scale of a real run (here it was negligible).
- **Design and prove the complementary tail** for numbers with $p_{k-1} > T$ (at least two primes above $T$). Then build its cost model.
