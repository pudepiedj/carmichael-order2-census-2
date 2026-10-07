---
location: post
title: Tail Coverage Theorem
author: Claude Opus 5.5 and JCP
date: 2026-10-05
---

# The two-large-prime tail: what it provably finds

**Date:** 5 October 2026
**Companion to:** `docs/CN2_LastCalc_Coverage_Proof.md`, whose numbering of facts (F1–F7) this note uses. That note proves that a capped `cn2xh` run with LASTCALC finds every order-2 $n < X$ whose second-largest prime is $\le T$. This note proves that two further runs find every order-2 $n < X$ whose second-largest prime is $> T$. Together they certify a census below $X$ with no hypothesis about caps or switch primes.
**Code proved about:** `cn2x/cn2pair.c` (`do_q`, `walk`, `single`, `sieve_ok`, `plan`, `worker`), `cn2x/cn2tail.c` (`scan`, `small_ok`, `amax_f`, `make_shards`, `worker`), and the verifiers `cn2_order2_defs.verify_candidate` and `cn2_split_campaign._verify_tail`.

---

## 1. The statement

**Theorem.** Fix a mode (rigid, howe or cheb; $D = 1, 1, 2$), $X < 2^{126}$, and integers $T \ge 5$ and $Q$ with $T < Q < 2^{32}$. Suppose two runs are complete, in the sense of §5:
- `cn2pair X T Q`, the **pair family**;
- `cn2tail X Q` with $R_0 = Q$, the **$(q, a)$ family**.

Let $n < X$ be an order-2 number with primes $p_1 < \dots < p_k$ and **$p_{k-1} > T$**. Then one of the two runs' outputs contains a record that its verifier turns into $n$.

**Corollary (the certificate).** Add a complete capped `cn2xh` run with `CN2X_TCAP` $= T$ and `CN2X_LASTCALC=1`, which covers every $n$ with $p_{k-1} \le T$ (companion note). Then the verified outputs of the three runs contain **every** order-2 number below $X$. Every verified output is an order-2 number by construction, so the union is exactly the set of order-2 numbers below $X$.

## 2. The split

Let $n$ satisfy the hypothesis, and write $p = p_{k-1}$ and $q = p_k$, so $T < p < q$. Exactly one of two cases holds:
- **$q > Q$:** $n$ has a prime above $Q$, and the $(q, a)$ family covers it (§3).
- **$q \le Q$:** the two largest primes satisfy $T < p < q \le Q$, and the pair family covers it (§4).

Nothing needs to be true of the other primes of $n$. The split is by the largest prime alone, and each family is complete for its own case.

## 3. The $(q, a)$ family (`cn2tail`, $R_0 = Q$)

**Lemma A.** Let $q$ be a prime of an order-2 $n < X$, with $q > R_0$. Write $g = g(q)$ and $m = n/q$. Then:
1. $q$ is among the engine's primes, which are all primes in $(R_0, \lfloor (D(X-1))^{1/3}\rfloor]$.
2. $m = r_0 + a\thinspace g$ for some $r_0 \in \lbrace q, 1\rbrace$ (only $r_0 = q$ in rigid mode) and some integer $a$ with $1 \le a \le A_{r_0}(q) := \lfloor ((X-1)/q - r_0)/g \rfloor$.

*Proof.*
1. F7: $q^3 < D n \le D(X-1)$.
2. As in the proof of F2, $m \equiv q$ (type 1) or $m \equiv 1$ (type $q$) modulo $g$. Since $g > q$, the least positive members of these classes are $q$ and $1$. We have $m \ne q$ ($n$ is squarefree) and $m \ne 1$ ($n$ is composite), so $a \ge 1$. Finally, $q m = n \le X - 1$ gives $a \le A_{r_0}(q)$. $\square$

**What the engine does with that pair** (`scan`, for $r_0 = q$ and, in howe and cheb modes, also $r_0 = 1$, for every $a$ in $[1, A_{r_0}(q)]$). It discards the candidate if any of these holds:
- **$q \mid m$**, detected by the residue of $a$ modulo $q$. Then $q^2 \mid n$, so $n$ is not squarefree.
- **For some small prime $5 \le s \le 509$ with $s \mid m$** (why this list, and why the list is immaterial to the proof: see the Remark at the end of §4): $s^2 \mid m$, or $n \bmod g(s) \notin \lbrace 1, s\rbrace$. The engine computes $n \bmod g(s)$ as $(q \bmod g(s))(m \bmod g(s))$. By the definition, both tests reject only non-order-2 numbers. Primes $s$ that divide $g = g(q)$ are not tested: for them $m \equiv r_0 \not\equiv 0 \pmod s$, so $s \nmid m$ anyway.
- **$2^{n-1} \not\equiv 1 \pmod m$.** This would contradict F6, since $m \mid n$.

Otherwise it writes `P q a` ($r_0 = q$) or `Q q a` ($r_0 = 1$). The verifier rebuilds $m = r_0 + a\thinspace g$ and $n = q m$, and checks $n$ against the definition. For our $n$, no test discards it, so its record is written. **Lemma A and these checks prove the case $q > Q$.**

The proved census below $10^{25}$ (the first repository (carmichael-order2-census, `results/howe_1e25`), $R_0 = 10^6$) used this engine and verifier. Its build (`038a9710…`) differs from the one frozen here (`aa96fb59…`) by a single added line: an input check that refuses $X \ge 2^{126}$ or primes above $2^{32}$. The search is otherwise identical.

## 4. The pair family (`cn2pair`)

**Lemma B.** Let $n < X$ be order-2 with $p = p_{k-1}$ and $q = p_k$, where $T < p < q \le Q$. Let $r_p, r_q$ be its type residues ($n \equiv r_p \pmod{g(p)}$, $n \equiv r_q \pmod{g(q)}$). Then:
1. $p \nmid g(q)$ and $q \nmid g(p)$, and $q \nmid g(q)$, $p \nmid g(p)$.
2. The four congruences $n \equiv 0 \pmod q$, $n \equiv r_q \pmod{g(q)}$, $n \equiv r_p \pmod{g(p)}$ and $n \equiv 0 \pmod p$ are simultaneously solvable, and their solutions form one residue class modulo $M = pq\cdot\mathrm{lcm}(g(p), g(q))$. That class contains $n$.
3. $q\thinspace g(q) < n < X$.

*Proof.*
1. F5, and $\gcd(x, x^2 - 1) = 1$.
2. $n$ is a solution. By 1, $p$ and $q$ are prime to the other moduli, so the solution set is a single class modulo the lcm of the moduli, which is $M$.
3. F2 at $q$. $\square$

**What the engine does** (`do_q` for $q = $ `PR[iq]`, then every `PR[ip]` $= p < q$, for every allowed pair of type residues):
- **Skip.** If $p \mid g(q)$, it skips the pair. By Lemma B.1 this never skips our pair.
- **Class modulo $L_q = q\thinspace g(q)$.** It forms the class of $\lbrace n \equiv 0 \ (q),\ n \equiv r_q \ (g(q))\rbrace$, using $\gcd(q, g(q)) = 1$. If $L_q \ge X - 1$ it returns. That cannot happen for our $n$, since $L_q < n < X$.
- **Step A.** It merges in $n \equiv r_p \pmod{g(p)}$ by the CRT. Our types are compatible (Lemma B.2), so the merge succeeds, and the new class modulo $L_1 = \mathrm{lcm}(L_q, g(p))$ contains $n$.
  - If $L_1 > X - 1$, the class has at most one member below $X$. The engine computes it without forming $L_1$, by testing $k \le (X - 1 - c)/L_q$ before forming $c + L_q k$. That member is $n$. `single` then checks $n > q\thinspace g(q)$ and $p \mid n$, both true.
- **Step B.** Otherwise it merges in $n \equiv 0 \pmod p$. Since $p \nmid L_1$, $L_1$ is invertible modulo $p$.
  - If $L_1 p > X - 1$, it handles the class as in Step A.
  - Otherwise `walk` runs through every member of the class modulo $M = L_1 p$ in $(q\thinspace g(q), X)$. That includes $n$.
- **Filters** (`walk` and `single` alike). A candidate is discarded only if:
  - it is even, or divisible by 3 (F1);
  - for some $s \in \lbrace 5, \dots, 59\rbrace$ with $s \mid n$: $s^2 \mid n$, or $n \bmod g(s) \notin \lbrace 1, s\rbrace$ (the definition);
  - $2^{n-1} \not\equiv 1 \pmod n$ (F6).

  In `walk`, $n \bmod s$ is carried by adding $M \bmod s$ at each step. That is the same value, computed without division. Our $n$ passes every filter, and `N n` is written.

The verifier factorises $n$ and checks the definition. **Lemma B and these checks prove the case $q \le Q$.** $\square$

**Remark: the choice of sieving primes.** Both engines sieve with a fixed list of small primes, and the lists differ: $5 \le s \le 509$ (95 primes) in `cn2tail` and $5 \le s \le 59$ (15 primes) in `cn2pair`. **The choice affects speed only, never correctness.** Each test rejects a candidate only when a prime $s \mid n$ contradicts the definition ($s^2 \mid n$, or $n \bmod g(s) \notin \lbrace 1, s\rbrace$). No order-2 number can do that, whatever $s$ is. So any list, including the empty one, leaves every order-2 number in place, and the lists could be changed without touching this proof. Both start at 5 because 2 and 3 are handled separately by F1.

The lengths reflect how each engine pays for a sieving prime:

- **`cn2tail` (95 primes).** The candidates for one prime $q$ are consecutive values $a = 1, 2, 3, \dots$, and $s \mid m$ selects a single residue class of $a$ modulo $s$. The engine marks those $a$ in blocks of 4,096 with one stride per prime, so an extra prime costs only about $4096/s$ operations per block. And since nothing constrains $m$ modulo small primes, the sieve pays off: it rejects about 58% of candidates.
- **`cn2pair` (15 primes).** Each candidate lies in a class modulo $\mathrm{lcm}(g(p), g(q))$, and that modulus is divisible by most small primes. For such an $s$, $n$ is congruent modulo $s$ to $r_p$ or $r_q$, whichever of $g(p)$ and $g(q)$ it divides. Neither residue is $0$ modulo $s$, so $s$ can never divide $n$ and testing it is wasted. The residues $n \bmod s$ are carried with each candidate, so every listed prime costs time on every candidate, against an expected gain of at most about $1/s$. The list stops at 59 on that estimate. It wasn't tuned, and the measured rejection rate is 19%.

## 5. Completeness of the runs

**`cn2pair`.**
- **Shard plan.** `plan` cuts the prime list in $(T, Q]$ into `nshards` consecutive index ranges, and the last range takes every remaining prime. The ranges cover every $q$, and so every pair $p < q$, since `do_q` takes every $p$ below $q$. The plan depends only on $X$, $T$, $Q$, the mode and `nshards`.
- **Recording.** A shard is recorded in `done.log` only after `out.txt` has been flushed. A shard interrupted by a hard stop is never recorded: `do_q` returns early and the worker then checks STOP before recording. Since the fix of 5 October, that holds whichever $q$ the stop interrupts.
- **Resumption.** A restart skips recorded shards. The final JSON reports `"complete":true` exactly when every shard is recorded and there was no stop.

**`cn2tail`.**
- **Shard plan.** `make_shards` covers every prime in $(R_0, \lfloor(D(X-1))^{1/3}\rfloor]$ whose range $A(q)$ is non-empty. $A$ is non-increasing in $q$, so stopping at the first empty range loses nothing.
- **Large primes.** A prime with a large range is split into consecutive sub-ranges of $a$ that cover $[1, A(q)]$.
- **Recording and completeness.** Interruption and recording work as in `cn2pair`, as does the `complete` flag.

Both runs must report `"complete":true`, and the verification must finish with no `UNPROVEN-PRIMALITY` flags (companion note, §5).

## 6. Arithmetic

- **`cn2pair`.** $q < 2^{32}$, so $g(q) < 2^{64}$ and $L_q < 2^{96}$. Each modulus is formed only when it is $\le X - 1 < 2^{126}$. In the determined cases the single member is formed only after checking that it is below $X$ (Lemma B, Steps A and B). Products of two values below $2^{64}$ are formed in 128 bits. The Fermat routine is the capped engine's, covered by its self-test.
- **`cn2tail`.** It uses the same Fermat routine, and $m < X$ throughout.

## 7. Checks so far

- **Pair classes.** For every one of the 95 known numbers below $10^{30}$, $n$ lies in the class that Lemma B builds from its own top two primes (`cn2_pair_tail.py unittest`).
- **Exact count.** In count-only mode, `cn2pair` reproduces exactly the candidate count of an independent Python implementation: 137,505,205 at $X = 10^{20}$, $T = 10^3$, $q \le 6{,}000$.
- **Proved census below $10^{25}$.** At $T = 10^4$ the two terms with second-largest prime above $T$ are $2.088\times 10^{24}$ and $2.196\times 10^{24}$. The pair family with $q \le 2\times 10^4$ finds exactly these two and nothing else, both before and after the sieve was added. The sieved engine's survivors are a subset of the unsieved engine's, as soundness requires.
- **The $10^{28}$ run.** Its pair family found $2.088\times 10^{24}$, the one known term below $10^{28}$ with second-largest prime above $T = 15{,}000$, within its first minutes.
