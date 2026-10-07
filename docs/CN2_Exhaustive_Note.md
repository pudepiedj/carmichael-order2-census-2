# Note on exhaustive CN2 enumeration: a complete census to $10^{24}$, extended to $10^{25}$ under Howe's definition, a new CN2 at $10^{24.32}$, and the cost of being unconditional

**Date:** 18 September 2026, revised 20 September 2026 with the completed $10^{24}$ census, and 24 September 2026 with the split search and the complete Howe census below $10^{25}$ (§9)
**Prepared by:** Claude (Opus 5), from the question "are we anywhere near extending Pinch–Goutier to $10^{24}$?"
**Refers to:** `CN2_Desert_Census.md` (the desert question), `CN2_LAZY_HANDOFF.md` (the conditional sweep), `CN2_Index_Note.md` (the SS engine and the count model)
**Engines:** `cn2_exhaustive.py` (Python reference), `cn2x/cn2x.c` (C engine), `cn2x/cn2xc.c` (checkpointed C engine), `cn2x/cn2xh.c` (the same with Howe and A299799 modes), `cn2x/cn2tail.c` (the large-prime tail)
**Drivers:** `cn2_exhaustive_c.py`, `cn2x_campaign.py`, `cn2x/verify_run.py`, `cn2x/test_checkpoint.py`, `cn2_split_campaign.py`, `cn2_order2_defs.py`
**Data:** `cn2x_scan.jsonl`, `cn2x_capped_scan.jsonl`, `cn2x_python_reference.json`, `cn2x/runs/*/`, `cn2x/runs/howe_1e25/` (the §9 census)

---
[TOC]

## 0. Verdict

1. **Exhaustive CN2 enumeration below $X$ is feasible and is now done to $10^{24}$.** 32.3 h of engine time on 30 threads gives **exactly the six known CN2s below $10^{24}$ and no others** — the first complete count of CN2s obtained independently of a Carmichael table, and two decades beyond Pinch–Goutier. (§3, §4)
2. **Decades 22 and 23 are empty, unconditionally.** No factor cap, no (H-base), no (H-rogue): every squarefree $n < 10^{24}$ consistent with Korselt was examined. The conditional runs had already said so; this proves it. (§4)
3. **A new CN2 was found**: $n = 2088144166339753513992001 \approx 10^{24.32}$, the 14th known below $10^{25}$, missed by every previous campaign for a structural reason. It lies *above* the census bound, so the two results are consistent: the desert ends and the first inhabitant appears within a third of a decade. (§5)
4. **Rigour has a measured price.** Capping factors at 20,000 makes the search $5.4\times$ faster at $10^{18}$; the unconditional run to $10^{24}$ took 32.3 h against 20 minutes capped. (§6)
5. **Added 24 September: the census is complete below $10^{25}$, under Howe's definition.** Splitting the search at a prime bound $R_0 = 10^6$ (a capped search, plus a direct scan of every $n$ with a larger prime factor) made the unconditional run 15.1 h instead of about 132 h. It ran under Howe's definition ($n \equiv 1$ **or** $p \pmod{p^2-1}$), which is OEIS A175531's and strictly wider than the rigid CN2 of §1. **There are exactly 14 order-2 Carmichael numbers below $10^{25}$**: the six known below $10^{22}$ and the eight of decade 24 in §4. All are rigid; none has a prime factor above $10^6$ (the largest is 17,291). Decades 22 and 23 are empty under either definition. (§9)

## 1. What makes the search possible

A rigid CN2 is a squarefree composite $n$ with $p^2 - 1 \mid n-1$ for every prime $p \mid n$. The task is: **given a search bound $X$, list every rigid CN2 $n \le X$, with proof that none is missed.**

Notation used throughout (the search bound $X$ is the only free parameter; everything else is derived):

| symbol | meaning |
|---|---|
| $X$ | the **search bound**: the run enumerates all rigid CN2 $n \le X$ ($10^{22}$, $10^{24}$, $10^{25}$ below) |
| $T = \lfloor X^{1/3}\rfloor$ | the **prime bound**: no factor of any $n \le X$ can exceed it (Lem1), so the prime table stops there ($2.15\times10^{8}$ at $X = 10^{25}$) |
| $k = \omega(n)$ | number of prime factors of $n$, all distinct |
| $p_1 < \cdots < p_k$ | the factors; $q < r$ denote the two largest, $P = n/(qr)$ the rest |
| $m$, $p_{\max}$ | a **prefix** $p_1 \cdots p_j$ built by the search, and its largest prime |
| $L = \mathrm{lcm}_{i \le j}(p_i^2-1)$ | the Korselt modulus of the prefix; $L(n)$ when $j = k$ |
| $s = n/m$ | the **completion** of a prefix: the product of the primes still to be chosen |
| $t$ | a candidate next prime extending a prefix, always $t > p_{\max}$ |

Three facts, all proved and all checked on the knowns, cut the space:

| | statement | consequence |
|---|---|---|
| Lem1 | $n/p \equiv p \pmod{p^2-1}$, and $n/p \neq p$, so $n/p \ge p^2 + p - 1$ | **every prime factor is below $X^{1/3}$** |
| Lem2 | $2$ and $3$ divide $p^2-1$ for every other prime $p$ | $2 \nmid n$, $3 \nmid n$ |
| Lem3 | with $q < r$ the two largest prime factors of $n$ and $P = n/(qr)$: $Pq \equiv r \pmod{r^2-1}$ forces $Pq > r^2$ | $q < r < P$, hence $k \ge 4$ |

**Proof of Lem1.** Let $p \mid n$. Rigidity gives $n \equiv 1 \pmod{p^2-1}$, and $p^2 \equiv 1 \pmod{p^2-1}$ trivially. Multiplying the first congruence by $p$,
$$p\,n \equiv p \pmod{p^2-1}, \qquad p\,n = p^2\cdot\frac{n}{p}=(p^2-1)\frac{n}{p}+\frac{n}{p} \equiv \frac{n}{p} \pmod{p^2-1},$$
so $n/p \equiv p\pmod{p^2-1}$. Since $n$ is squarefree, $n/p \neq p$ (else $n = p^2$), and $n/p > 0$, so the least admissible value is one full period higher: $n/p \ge p + (p^2-1) = p^2 + (p-1) > p^2$, i.e. $n > p^3$. With $n \le X$ this bounds every prime factor by $X^{1/3}$ — the fact that makes $10^{24}$ reachable at all (primes to $10^8$, not $10^{12}$). In fact the empirical evidence suggests that the largest prime divisor of at least a small CN2 is orders of magnitude smaller than this cube-root law implies, but the upper bound makes the calculation tractable. Verified on the six knowns before anything was built on it, and independently confirmed by the engine's exact agreement with ground truth at every bound to $10^{22}$: too small a prime table would have lost CN2s.

Lem2 and Lem3 are equally short.

**Proof of Lem2**: for distinct primes $p, q \mid n$ we have $2, 3 \mid p^2-1 \mid n-1$, so $2, 3 \nmid n$.

**Proof of Lem3**: applying Lem1's congruence to the largest prime $r$ with $P q = n/r$ gives $Pq \equiv r \pmod{r^2-1}$ and $Pq \neq r$, so $Pq > r^2 > qr$, hence $P > r > q$; a three-factor CN2 would need $P = p_1 < q$, so $k \ge 4$.

Lem1 is what makes $10^{24}$ reachable at all: the primes must stop at $10^8$, not $10^{12}$. In practice, for empirical reasons, even this upper bound is unnecessarily high (see the capped runs below).

## 2. The search

Depth-first over prefixes $m = p_1 \cdots p_j$ in increasing prime order, carrying $L = \mathrm{lcm}(p_i^2-1)$.

- **Korselt is the invariant, not a test.** Completions satisfy $s \equiv m^{-1} \pmod L$ by construction, so the congruence is generated, never checked. A prefix is pruned unless $\gcd(m, L) = 1$, and a child prime $t$ needs $t \nmid L$ and $\gcd(m, t^2-1) = 1$.
- **Role bounds.** The next prime must be the last ($t^2 < m$, by Lem1), the second-to-last ($t < m$ and $mt^2 < X$, by Lem3), or earlier ($mt^3 < X$). If no prime fits any role, the node is dead.
- **Two regimes.** When $L > X/m$ the completion is forced into at most one candidate, computed by a modular inverse — this is 98% of all nodes. Otherwise the node branches. The switch is tuned by the ratio of candidates to children.
- **Filters before factoring.** A completion's prime factors all exceed $p_{\max}$, so trial division by small primes rejects 56% of candidates. Then a base-2 Fermat test: every CN2 is Carmichael, so $2^{n-1} \equiv 1 \pmod s$ is necessary, and it removes 99.9988% of the rest. Only the survivors are factored. At $10^{22}$: 301.4 billion candidates, 2.34 million factored, 27 s.
- **Budget prune.** $L(n) \mid n-1 < X$, so a prefix with $L \ge X$ is dead. (Measured: this fires rarely — 3 nodes at $10^{14}$ — because deep prefixes have $X/m < L < X$.)

## 3. Validation

| check | result |
|---|---|
| all 8 counters, C vs the Python reference, $10^{10}$–$10^{15}$, 1 and 30 threads | identical |
| 128-bit helpers (inverse, isqrt, icbrt, mod) against Python, incl. edges at $2^{64}$, $2^{84}$ | 72,288 tests, 0 mismatches |
| Montgomery base-2 Fermat, one- and two-limb, against `pow` | 5,886 tests, 0 mismatches |
| CN2s found vs the known list, $10^{15}$–$10^{24}$ | exact at every bound |
| the $10^{24}$ census itself, 1,812,763 survivors factorised | 6 CN2s, all six known, no false positives |
| checkpoint/restart: deadline drains, STOP-file stops, `kill -9`, up to 11 sessions | exact counters and hits |

Two defects were found by the C/Python comparison and fixed in both: a node whose role bounds cross did pointless work, and (in the old two-list engine) `uint32` subset indices silently truncate at $\lvert P\rvert \ge 66$.

## 4. Results

**Unconditional, complete, verified:**

| range | CN2s | nodes | candidates | time (30 threads) |
|---|---|---|---|---|
| $n < 10^{22}$ | 6 (the known six) | $2.7425\times10^{11}$ | $3.0141\times10^{11}$ | 6,985 s (1 h 56 m) |
| $n < 10^{24}$ | **6 (the same six)** | $4.2335\times10^{12}$ | $4.5013\times10^{12}$ | 116,233 s (32.3 h) |

Both are independent of Goutier's table, which is where the six had come from. The $10^{24}$ census ran as two sessions (18–19 September, 10.1 h; 19–20 September, 22.2 h) and finished at 08:46 on 20 September 2026.

Full counters for the $10^{24}$ run: 4,233,486,169,672 nodes, of which 4,233,473,118,867 (99.9997%) switched to LIST mode; 5,379,542,160,505 children tried; 4,501,252,595,528 candidates, of which 2,561,781,895,567 were removed by the small-prime sieve; 1,812,763 survived the base-2 Fermat test and were factorised, taking 19 s on 24 processes; 4 budget prunes; 1,647,987 shards.

The bound was the exact integer $10^{24}$, not a float (§7).

**On the predicted node total.** `init` extrapolates $4.561\times10^{9} \cdot 3.92^{\log_{10}X - 19}$ from the measured decades; it gave $4.2217\times10^{12}$ against the actual $4.2335\times10^{12}$, low by 0.19%. The observed ratio $10^{22} \to 10^{24}$ is $\times 3.929$ per decade, matching the $\times 3.89$–$3.93$ seen from $10^{15}$ to $10^{19}$. Seven decades of a stable geometric law make the progress percentage trustworthy enough to plan nights around, which is what it is for.

**Conditional on every prime factor $\le$ cap, over $[10^{22}, 10^{25})$:**

| cap | time | nodes | candidates | CN2s found |
|---|---|---|---|---|
| 20,000 | 4,586 s | $1.61\times10^{11}$ | $1.25\times10^{11}$ | 8 |
| 30,000 | 6,746 s | $2.18\times10^{11}$ | $1.71\times10^{11}$ | 8 (identical) |
| 100,000 | 14,364 s | $4.88\times10^{11}$ | $3.96\times10^{11}$ | 8 (identical) |

All three recover the seven previously known in range and add one new. Their upper bound was `1e25` parsed as a float, which rounds *up* to 10000000000000000905969664, so they searched a slight superset of $[10^{22}, 10^{25})$ — harmless, since nothing can be missed by searching too far, and all eight lie below $10^{24.98}$.

**The unconditional run subsumes the $[10^{22}, 10^{24})$ part of these three**, and agrees: no cap, and still nothing. The capped runs therefore stand corroborated where they overlap, and remain conditional only above $10^{24}$ — which is where all eight of their finds actually live. **Decades 22 and 23 are empty in every run**, and **no CN2 below $10^{25}$ has a prime factor between 20,000 and 100,000** — against a largest observed factor of 17,291, so the observed maximum is not sitting just below a cliff. (The 100,000 run used the 95-prime sieve of §6; with the 13-prime sieve it would have been slower, so its time is not directly comparable with the other two.)

The eight in $[10^{22}, 10^{25})$, all in decade 24:

| $n$ | $\log_{10}$ | $k$ | largest prime | geometric mean |
|---|---|---|---|---|
| 1159954316194989017102401 | 24.064 | 9 | 4759 | 472 |
| **2088144166339753513992001** | **24.320** | **8** | **17291** | **1096** |
| 2196407820059694924883201 | 24.342 | 9 | 13441 | 507 |
| 3339611018825185787482801 | 24.524 | 12 | 1429 | 111 |
| 4105879060352839139462401 | 24.613 | 8 | 7237 | 1193 |
| 5002862939121639632040001 | 24.699 | 11 | 2549 | 176 |
| 7865064643837556041286401 | 24.896 | 11 | 15809 | 183 |
| 9400084864021826054720641 | 24.973 | 11 | 4523 | 186 |

## 5. The new CN2

$$n = 2088144166339753513992001 = 113 \cdot 199 \cdot 239 \cdot 263 \cdot 701 \cdot 7919 \cdot 15391 \cdot 17291$$

Verified from scratch: squarefree, composite, $p^2-1 \mid n-1$ for all eight primes. It was not in any earlier record of this work, and was first found on 18 September 2026.

**Why every previous campaign missed it.** Its canonical modulus is
$$L(n) = 21731935068043200 = 2^6 \cdot 3^4 \cdot 5^2 \cdot 7 \cdot 11 \cdot 13 \cdot 17 \cdot 19 \cdot 37 \cdot 107 \cdot 131 = 2^{54.27},$$
carrying **three distinct rogue primes (37, 107, 131)**. The lazy sweep's DFS holds two rogue slots and cuts any branch reaching a third, so this number is outside its frame by construction. Per factor the stated (H-rogue) holds (at most two rogues each, all $\le 463$), as do (H-base) ($2^{35.29} \le 2^{51}$) and (H-factor) ($17291 \le 20000$): it is the *cumulative* two-rogue restriction that fails.

**Its statistics are unremarkable, which is the point.** Pool 51, $\log_2\varphi(L) = 51.659$, fecundity $-0.659$, so a predicted count of 0.63 against 3 actually present for that modulus — an ordinary Poisson arrival in the deep tail (cf. `CN2_Index_Note.md` §4). Index 1, no barred primes, 2-slack $2^{26.5}$. $k = 8$ ties the minimum among all 14 known, and its geometric mean 1096 sits inside the known range 68–1193. What *is* a record is the largest prime factor, 17,291, beating the previous 15,809.


## 6. Costs

Growth with $X$, on 30 threads, measured:

| | per decade | $10^{22}$ | $10^{24}$ | $10^{25}$ |
|---|---|---|---|---|
| unconditional | $\times 4.08$ (measured) | 1.9 h (measured) | 32.3 h (measured) | $\approx$ 5.5 days |
| capped at 20,000 | $\times 2.55$ | 297 s (measured) | $\approx$ 20 min | 1.27 h (measured) |
| capped at 30,000 | $\times 2.53$ | 397 s (measured) | $\approx$ 27 min | 1.87 h (measured) |
| capped at 100,000 | — | — | — | 3.99 h (measured) |

So an unconditional result costs roughly $5.4\times$ a capped one at $10^{18}$, and far more at $10^{24}$ where the proved bound is $10^8$ against an observed maximum factor of 17,291. **The unconditional run to $10^{24}$ settles decades 22–23 only**; decade 24, where every known CN2 above $10^{22}$ lives, needs $10^{25}$: at the measured $\times 4.08$ that is about 132 h, or 13 nights of 10 hours.

Time grows at $\times 4.08$ per decade while nodes grow at $\times 3.93$, so the cost *per node* rises about 4% per decade — wider arithmetic (the 2-limb Montgomery path above $2^{63}$) and longer Fermat exponents.

**GPU:** rejected on measurement, not principle. The base-2 Fermat tests are 44.8% of the run (capped, $10^{21}$: 158.7 s with, 87.6 s without), so by Amdahl a perfect offload gives at most $1.81\times$ — about 2.5 hours off a 6-hour run, against a day of building a Metal kernel with 32-bit-limb Montgomery arithmetic. It becomes worthwhile only for the unconditional $10^{25}$, where 40% is roughly seven nights.

**Small-prime sieve:** extending it from 13 primes to 95 (up to 509) raises rejection from 33.9% to 57.5% and saves 9% at $10^{21}$ capped and 13% at $10^{20}$ uncapped (643.2 s $\to$ 560.4 s) — but *costs* 52% at $10^{19}$, where exponents are short and a Fermat test is cheap. The crossover sits near $10^{20}$, so both engines now choose the depth from $X$ (13 primes below $10^{20}$, 95 above) and report it as `sieve_primes`; `CN2X_SIEVEP` overrides it. The gain is bounded by the Fermat share: at best it removes a third of 45%, so $\approx 15\%$. That predicted $\approx 29$ h for the unconditional $10^{24}$ run against a 33 h baseline; the run used the 95-prime depth and took **32.3 h**, so the saving was well under the bound. The sieve did remove 56.9% of candidates ($2.562\times10^{12}$ of $4.501\times10^{12}$), but the ones it removes are the cheap ones — a candidate killed by trial division against a small prime would mostly have failed the Fermat test quickly too.

## 7. Checkpointing and the nightly campaign

The search is cut into deterministic shards: nodes with $m < 10^5$ are expanded, and the remaining children of each are grouped into runs of $1, 2, 4, \ldots, 4096$, so the heaviest children get shards of their own. Shard ids depend only on $X$ and the sharding parameters, so a restart replays the same split and skips finished shards, each logged with its counters.

At a window's end the run **drains**: no new shards are taken, in-flight shards finish and are recorded, then it exits — nothing is lost, and the overrun is bounded by the largest shard (seconds at test scale, minutes at $10^{24}$). A hard stop (`touch <dir>/STOP`, or a signal) abandons in-flight shards, which are simply redone.

```bash
python cn2x_campaign.py init --X 1e24 --dir cn2x/runs/campaign_1e24
python cn2x_campaign.py run  --dir cn2x/runs/campaign_1e24 --window 22:00-08:00 --threads 30
python cn2x_campaign.py run  --dir cn2x/runs/campaign_1e24 --now --window 22:00-08:00   # start now, no deadline
python cn2x_campaign.py status --dir cn2x/runs/campaign_1e24
python cn2x_campaign.py verify --dir cn2x/runs/campaign_1e24
```

`init` freezes the engine binary and parameters; `run` refuses to resume against anything different. `--now` starts a session immediately with no deadline, for when the machine is free during the day; `touch <dir>/DRAIN` then ends it gracefully at any time, after which the runner waits for the *next* window start rather than restarting inside the current one. The $10^{24}$ census used both: one 22:00–08:00 window, then a `--now` session that ran to completion.

**Parse the bound as an integer.** `--X 1e24` through `float()` gives 999999999999999983222784, because $10^{24}$ needs 54 bits of mantissa — a silently *smaller* search than advertised. `parse_bound()` expands the exponent in integer arithmetic and refuses anything it cannot represent exactly. This was caught in `params.json` at launch; the earlier capped runs had rounded the other way (§4), which is harmless but should still be quoted exactly.

**Reading the progress line.** `nodes`, `candidates` and the percentage are cumulative across sessions; `shards done` counts the current session only. `root i/hi` is the generator's position among the admissible smallest prime factors — `hi` = 5,761,453 = $\pi(10^8) - 2$, the primes 2 and 3 being excluded by (L2) — so `root 34` means every $n$ whose least prime factor is at most $\mathrm{PR}[33] = 139$ has been dispatched. It measures dispatch, not completion. `hits` counts only CN2s found *at* a prefix ($L \mid m-1$ with the prefix already complete); almost all CN2s arrive instead as LIST survivors and appear only after `verify`, so `hits 0` throughout a run is normal. In the $10^{24}$ census all six came that way.

**`done.log`** carries one line per finished shard: `D id` then the nine counters (`nodes children last_tests list_nodes list_candidates fermat_passes factorisations pruned_budget sieved`). The id is a serial number in generation order, so lines are neither sorted nor contiguous — shards finish out of order across 30 threads. A restart replays the generator, skips ids already present, and adds their counters back as the prior total, which is what makes the final figures exact rather than approximate. Duplicate ids are counted once and a torn final line is dropped, so `kill -9` is safe.

## 8. Open

1. ~~**Sweep the three-rogue family.**~~ **Withdrawn — it was done in July.** This item originally claimed that three-rogue generators "are exactly what no campaign has enumerated". That is false: `CN2_Desert_Census.md` records a triple-rogue mini-sweep of 11 July 2026 over 48,053 targets (`harvest_targets_3rogue.json`), which found **zero in 108 s**, the modal outcome at P(0) ≈ 99.6%. The census had also checked the rogue axis to five rogues. What is true is narrower: the *lazy sweep's* two-slot frame could not see 2088144166339753513992001, which is why that number survived to be found here.

   A fecundity measurement of 20 September agrees with the July verdict and explains it: fecundity falls by about 7 per additional rogue, because a rogue inflates φ(L) while adding no pool primes (a fresh large prime in L makes p²−1 | L harder, not easier). Sampling random decade-24 prime products gives generators with a median of 103 bits and 8.3 rogues — utterly sterile — against ≤ 3 rogues and ≤ 54.3 bits for all fourteen CN2s below $10^{25}$ (a complete list, §9). Three rogues is the edge of viability, not an unexplored frontier: the 2088 generator scores −0.659 and produced a family of three.
2. ~~**$10^{25}$ unconditional**~~ **Done 24 September, in 15.1 h under Howe's definition (§9).** The original item read: **$10^{25}$ unconditional** — about 132 h, or 13 nights of 10 hours, and roughly 7 with a working GPU Fermat kernel. It would turn the list of eight in decade 24 into a theorem, and settle whether 2088144166339753513992001 has neighbours that every capped run has been blind to. This is the one remaining computation that changes what is *known* rather than what is *believed*.
3. **Re-tune the LIST switch ratio** for the compiled engine — it was tuned in Python, where the cost balance differs. With 99.9997% of nodes going to LIST at $10^{24}$, the switch is firing almost everywhere, and the threshold has never been tested against the C engine's actual costs.
4. **Cosmetic, next engine build:** print the cumulative shard total with the session's count in brackets, and label `hits` as "path hits" (§7). Both misled a reader of the live log during this census.

---

## 9. The split search and the complete census below $10^{25}$ (24 September 2026)

### 9.1 Two gaps in the census as it stood

1. **Cost.** Open item 2 priced the unconditional $10^{25}$ run at about 132 h. The capped runs show why: at $10^{24}$ the capped search took 20 minutes against 32.3 h, so about 97% of the unconditional cost is spent on branches containing a prime factor above 20,000.
2. **Definition.** OEIS A175531 ("Carmichael numbers of order 2") uses Howe's definition: odd composite $n$ with $n\equiv 1$ **or** $p \pmod{p^2-1}$ for every prime $p\mid n$. A prime with $n\equiv 1$ is *type 1*; the rigid CN2 of §1 has only type-1 primes. A prime with $n \equiv p$ is *type $p$*, i.e. $p-1\mid n-1$ but $p+1\mid n+1$. Every search in this note up to §8 was rigid-only, so "decades 22–23 empty" was proved for rigid numbers only. Non-rigid order-2 numbers do exist: Howe (2000, §5) gives one of about $3.92\times10^{59}$, with 23 prime factors, in which only $p = 1153$ is type $p$. Below $10^{22}$ the definitions agree, since Goutier's table shows A175531 has only the six rigid terms there.

### 9.2 The identity behind the split

For any prime $p \mid n$, Lem1's congruence (and its type-$p$ counterpart) fixes the cofactor $n/p$ modulo $p^2-1$:

- type 1: $n/p \equiv p$, $n/p \ne p$, so $n = p\,\bigl(p + a(p^2-1)\bigr)$ and $n - 1 = (p^2-1)(ap+1)$;
- type $p$: $n/p \equiv 1$, $n/p \ne 1$, so $n = p\,\bigl(1 + a(p^2-1)\bigr)$;

with $a \ge 1$ in both cases. So every order-2 number below $X$ that has *some* prime factor $p > R_0$ is one of the pairs $(p, a)$ with $R_0 < p \le X^{1/3}$ in one of two families, about $X/(R_0^2\ln R_0)$ pairs in all. The search splits into two phases whose union is complete, with no hypothesis:

| phase | engine | covers |
|---|---|---|
| **capped** | `cn2x/cn2xh.c`, prime table cut at $R_0$ | every solution whose primes are all $\le R_0$: the cap only removes table primes, so such a solution's whole path stays in the tree |
| **tail** | `cn2x/cn2tail.c` | every solution with a prime $> R_0$: all pairs $(p, a)$, both families |

In Howe mode each prime of the capped search chooses its type, and a node carries $n \equiv c \pmod L$, built by combining congruences, in place of $n\equiv 1$. Incompatible choices are dropped at once. A merge that pushes $L$ past $X$ determines $n$, which is then resolved on the spot. The role bounds of §2 hold unchanged, because $n/p > p^2-1$ for either type. Both phases apply only necessary conditions before factoring: the definition for small primes, squarefreeness, and base-2 Fermat (every order-2 number is a Carmichael number).

### 9.3 Validation

| test | result |
|---|---|
| `cn2xh` rigid mode against `cn2xc`, $10^{16}$ and $10^{18}$ | all nine counters and every output record identical |
| `cn2tail` rigid mode against its first build, $10^{17}$ | identical pairs, counters and survivors |
| split, rigid, $10^{18}$ with $R_0 = 1000$ | 842526563598720001 (largest prime 2,729) found by the tail |
| A299799 mode ($n \equiv 1$ or $p \pmod{(p^2-1)/2}$) to $10^{20}$ | exactly A299799's eight terms, three of them not rigid CN2s |
| relaxed modulus $(p^2-1)/D$, $D = 24, 12, 8, 6$, Fermat off, $X = 10^{10}$ | every brute-force solution found (23, 6, 4, 3), including 16 with type-$p$ primes; 8 found only through the tail's type-$p$ family |
| Howe mode, split, $10^{22}$ | exactly Goutier's six, none non-rigid, in 20 minutes |
| stop, drain and resume, Howe mode, both engines | counters and records identical to uninterrupted runs |

The relaxed-modulus test exists because no known list below $10^{10}$ contains a type-$p$ prime. That includes A299799, whose terms are all type 1 in its own modulus. Loosening the modulus creates small solutions with type-$p$ primes, which a brute-force factorisation of every $n < 10^{10}$ lists independently.

### 9.4 Result

Run `cn2x/runs/howe_1e25`: Howe mode, $X = 10^{25}$ exactly, $R_0 = 10^6$, 30 threads, engines frozen at `init`.

| phase | work | engine time | candidates to factor |
|---|---|---|---|
| tail | 2,568 shards, both families, primes in $(10^6, 2.15\times10^8]$ | 1.42 h | 698 |
| capped | 1,085,720 shards | 13.71 h | 1,812,046 |

Verification factored all 1,812,744 candidates in 18 s and re-checked every find from scratch. **There are exactly 14 order-2 Carmichael numbers below $10^{25}$**, identical to the list of §4 together with the six below $10^{22}$. All 14 were found by the capped phase, and the tail found none, as it must: no order-2 number below $10^{25}$ has a prime factor above $10^6$. None is non-rigid, and the two phases agree on every number.

This settles open item 2 and more:

- **decade 24 is a theorem**: exactly the eight numbers of §4, with no neighbour of 2088144166339753513992001 that the capped runs could have missed;
- **decades 22 and 23 are empty** under Howe's definition, not only the rigid one;
- **no non-rigid order-2 number exists below $10^{25}$**;
- the list extends OEIS A175531 from six terms to fourteen.

**Cost.** 15.1 h against about 132 h projected for the single search, and against 32.3 h for the rigid census to $10^{24}$, one decade shorter.

### 9.5 Reproduction

```bash
python3 cn2_split_campaign.py init   --X 1e25 --R0 1e6 --dir cn2x/runs/howe_1e25     # --mode howe is the default
python3 cn2_split_campaign.py run    --dir cn2x/runs/howe_1e25 --now --threads 30   # or --window 22:00-08:00
python3 cn2_split_campaign.py status --dir cn2x/runs/howe_1e25
python3 cn2_split_campaign.py verify --dir cn2x/runs/howe_1e25
```

---

### Appendix: what each script does

| file | role |
|---|---|
| `cn2_exhaustive.py` | Python reference implementation; the counter oracle for the C engines |
| `cn2x/cn2x.c` | C engine: 128-bit, Montgomery Fermat, shared task stack with donation. Experiment knobs `CN2X_TCAP` (factor cap: makes results conditional), `CN2X_NMIN` (report floor), `CN2X_SIEVEP` |
| `cn2x/cn2xc.c` | the same search, checkpointed: deterministic shards, `done.log`, drain-at-deadline, progress lines |
| `cn2_exhaustive_c.py` | driver: runs `cn2x`, verifies survivors with sympy, validates against the known CN2s; `--parity` compares counters with the Python reference |
| `cn2x_campaign.py` | windowed campaign runner: `init` / `run` / `status` / `verify` |
| `cn2x/verify_run.py` | verifies a run directory's survivors and compares against the known list (must be run as a file: the worker pool re-imports it) |
| `cn2x/test_checkpoint.py` | checkpoint tests: drains, hard stops, `kill -9`, all compared against an uninterrupted run |
| `cn2x/cn2xh.c` | `cn2xc` with `CN2X_MODE` = `rigid` (identical to `cn2xc`), `howe` (A175531), `cheb` (A299799); `test` with `CN2X_D` is validation only (§9.3) |
| `cn2x/cn2tail.c` | the large-prime tail: every $(p, a)$ with $p > R_0$, both families, checkpointed like `cn2xc`; `--bench` measures its rate |
| `cn2_split_campaign.py` | split driver: `init` / `estimate` / `run` / `status` / `drain` / `stop` / `verify`, phases `tail` then `capped`; `verify` reports rigid and non-rigid finds separately |
| `cn2_order2_defs.py` | the three definitions (rigid, Howe, A299799) and the candidate verifier shared by the drivers |

Run directories under `cn2x/runs/` hold `out.txt` (hits and survivors), `done.log` (finished shards with counters), `progress.log`, `session.json`, `result.json`, and — for conditional runs — `ENV.txt` recording the hypothesis and the environment needed to resume.
