---
location: post
title: Certified CN2 Census Below 10^28
author: Claude Opus 5.5 and JCP
date: 2026-10-07
---

# Certificate: the order-2 Carmichael numbers below $10^{28}$

**Date:** 7 October 2026 (runs 5–6 October 2026)
**Definition:** Howe's order-2 Carmichael numbers, OEIS A175531: squarefree composite $n$ with $n \equiv 1$ or $p \pmod{p^2-1}$ for every prime $p \mid n$.
**Status:** the computations are complete and every check below passes. The result rests on two proofs, `docs/CN2_LastCalc_Coverage_Proof.md` and `docs/CN2_Tail_Coverage_Proof.md`. They have been checked by the authors and by an independent review, whose suggested further checks have all been carried out (§5).
**Run directory:** `results/cert_1e28/`

---

## 1. The result

> **There are exactly 45 order-2 Carmichael numbers below $10^{28}$** (§4). The search makes no assumption about table caps or LIST switch primes.

Before this certificate, the census was proved below $10^{25}$ (14 terms). Above that it was conditional on the capped-search hypothesis: that every switch prime is below the table cap. This certificate removes the hypothesis up to $10^{28}$. **Terms 15–45 of A175531 are now proved**, and the proved bound moves from $10^{25}$ to $10^{28}$.

Two by-products, now proved below $10^{28}$:
- **All 45 terms are rigid.** Every prime $p$ has $n \equiv 1 \pmod{p^2-1}$. The smallest non-rigid term therefore exceeds $10^{28}$; previously it was known only to exceed $10^{25}$. The search allowed both residues at every prime (Howe mode), so this is a finding, not an assumption. It is recorded per term in the `rigid` field of the repository's `results/factorisations.json`, and as `non_rigid: []` in the capped run's `result.json`. The standalone checker (§5) confirms it independently.
- **No term below $10^{28}$ is divisible by 5 or 7.** The search covered both primes like any others.

## 2. How the search is split

Every order-2 $n < X = 10^{28}$, with primes $p_1 < \dots < p_k$, falls into exactly one of three parts. The split uses $T = 15{,}000$ and $q^* = 24{,}500{,}000$.

| part | which $n$ | engine | proved by |
|---|---|---|---|
| capped | $p_{k-1} \le T$ | `cn2xh`, table cap $T$, `CN2X_LASTCALC=1` | LASTCALC coverage theorem |
| pair tail | $p_{k-1} > T$ and $p_k \le q^*$ | `cn2pair` (two largest primes $T < p < q \le q^*$) | tail theorem, Lemma B |
| $(q, a)$ tail | $p_{k-1} > T$ and $p_k > q^*$ | `cn2tail`, $R_0 = q^*$ | tail theorem, Lemma A |

The $(q, a)$ engine in fact covers every $n$ with a prime above $q^*$. The overlap is harmless: only gaps would matter.

## 3. The runs

All three ran on 30 threads, in Howe mode. Each engine was compiled into the run directory before the runs began, and its hash was recorded in `ENGINES.sha256`. After the runs the hashes were rechecked and found identical:

| engine | sha256 |
|---|---|
| `capped/cn2xh` | `1a8157d091f1a5d764c15e954a87ba08b2a6f8e337306069510baa654474d074` |
| `pair/cn2pair` | `bc3f49c6c96aab24d0791bb3a6ad4f7e213296be681225b372feb4a98c02e0dc` |
| `pair/cn2pair.c` | `ee337b6e474cbf33da4f7166238e00d151049b981d052a3b6656504f58f113cd` |
| `tail/cn2tail` | `aa96fb5916d4fa15ba7f7757f96662fa9c3c9c7d1cedeba0c87b127ac71d9836` |
| `tail/cn2tail.c` | `c5379c211c663fb795095cd232b24cd8ed0f727f382108f501b219d911e7b83b` |

| | pair tail | $(q, a)$ tail | capped |
|---|---|---|---|
| shards | 4,096 / 4,096 | 3,542 / 3,542 | 213,275 / 213,275 |
| completion record | `"complete":true`, one session | `"complete":true`, one session | `"complete":true`, `"lastcalc":1`, `skipmask 0`, one session |
| work | $1.18\times 10^{12}$ prime pairs; $4.57\times 10^{11}$ candidates | $9.52\times 10^{11}$ $(q, a)$ tests | $1.91\times 10^{12}$ nodes; $1.43\times 10^{12}$ list candidates |
| rejected before Fermat | $8.5\times 10^{10}$ (19%) | $5.6\times 10^{11}$ (59%) | — |
| Fermat survivors | 41,041 | 159 | 4,201,737 |
| order-2 numbers below $X$ | 1 | 0 | 45 |
| time | 8.0 h | 4.7 h | 12.9 h |

The capped run's `done.log` holds exactly one record for each shard id from 0 to 213,274.

**Verification.** Every survivor was factorised and checked against the definition:
- **capped part:** `cn2_split_campaign.py verify --dir results/cert_1e28`, which wrote `result.json`;
- **tail:** `cn2_tail_verify.py --X 1e28 --pair .../pair --qa .../tail --T 15000`, which gave **PASS**.

No `UNPROVEN-PRIMALITY` flag was raised, so every factor met in verification was below $2^{64}$. There the primality test used (BPSW) is known to give no false answers, by exhaustive computation: every base-2 pseudoprime below $2^{64}$ has been listed, and none passes BPSW.

## 4. The 45 terms

"Guaranteed by" names the part whose theorem covers the term. "Found by" lists the parts that actually produced it. The capped run also found term 8 through LIST, although only the tail guarantees it.

| # | $n$ | $\log_{10} n$ | $k$ | prime factors | guaranteed by | found by |
|---|---|---|---|---|---|---|
| 1 | 443372888629441 | 14.65 | 8 | 17 · 31 · 41 · 43 · 89 · 97 · 167 · 331 | capped | capped |
| 2 | 39671149333495681 | 16.60 | 9 | 17 · 37 · 41 · 71 · 79 · 97 · 113 · 131 · 191 | capped | capped |
| 3 | 842526563598720001 | 17.93 | 8 | 17 · 61 · 71 · 89 · 197 · 311 · 769 · 2729 | capped | capped |
| 4 | 2380296518909971201 | 18.38 | 9 | 19 · 41 · 43 · 71 · 89 · 127 · 199 · 449 · 991 | capped | capped |
| 5 | 3188618003602886401 | 18.50 | 8 | 29 · 37 · 79 · 181 · 191 · 449 · 701 · 3457 | capped | capped |
| 6 | 4208895375600667752001 | 21.62 | 10 | 17 · 29 · 31 · 43 · 71 · 79 · 199 · 389 · 2521 · 5851 | capped | capped |
| 7 | 1159954316194989017102401 | 24.06 | 9 | 31 · 47 · 109 · 137 · 307 · 2927 · 3079 · 4049 · 4759 | capped | capped |
| 8 | 2088144166339753513992001 | 24.32 | 8 | 113 · 199 · 239 · 263 · 701 · 7919 · 15391 · 17291 | tail | capped+tail |
| 9 | 2196407820059694924883201 | 24.34 | 9 | 53 · 109 · 127 · 281 · 307 · 379 · 647 · 10529 · 13441 | capped | capped |
| 10 | 3339611018825185787482801 | 24.52 | 12 | 19 · 41 · 43 · 53 · 67 · 89 · 103 · 131 · 137 · 307 · 389 · 1429 | capped | capped |
| 11 | 4105879060352839139462401 | 24.61 | 8 | 89 · 239 · 401 · 2143 · 2311 · 3571 · 3761 · 7237 | capped | capped |
| 12 | 5002862939121639632040001 | 24.70 | 11 | 31 · 53 · 79 · 89 · 101 · 151 · 181 · 251 · 379 · 647 · 2549 | capped | capped |
| 13 | 7865064643837556041286401 | 24.90 | 11 | 23 · 37 · 67 · 89 · 101 · 109 · 181 · 199 · 433 · 571 · 15809 | capped | capped |
| 14 | 9400084864021826054720641 | 24.97 | 11 | 23 · 31 · 53 · 79 · 103 · 197 · 239 · 379 · 521 · 727 · 4523 | capped | capped |
| 15 | 12071465216556111317952001 | 25.08 | 12 | 37 · 41 · 53 · 59 · 89 · 109 · 127 · 151 · 191 · 233 · 379 · 811 | capped | capped |
| 16 | 31280524488319495535546401 | 25.50 | 11 | 19 · 23 · 43 · 59 · 71 · 109 · 271 · 311 · 701 · 1301 · 47431 | capped | capped |
| 17 | 32492501043239735947612801 | 25.51 | 11 | 11 · 31 · 41 · 79 · 103 · 271 · 313 · 389 · 911 · 2521 · 3769 | capped | capped |
| 18 | 36991212507320386599648001 | 25.57 | 11 | 11 · 31 · 37 · 53 · 97 · 191 · 379 · 883 · 1151 · 2549 · 3041 | capped | capped |
| 19 | 44774022865681509785985601 | 25.65 | 11 | 19 · 41 · 67 · 89 · 101 · 131 · 139 · 233 · 271 · 1301 · 63799 | capped | capped |
| 20 | 50473184681492173584023041 | 25.70 | 11 | 37 · 53 · 67 · 131 · 139 · 181 · 229 · 419 · 577 · 911 · 2311 | capped | capped |
| 21 | 54475112449147532388966721 | 25.74 | 10 | 29 · 97 · 103 · 139 · 181 · 307 · 827 · 1061 · 2393 · 11593 | capped | capped |
| 22 | 76587257660615840312697601 | 25.88 | 11 | 13 · 19 · 31 · 67 · 107 · 139 · 197 · 479 · 1471 · 1801 · 40151 | capped | capped |
| 23 | 111356194710608262218643841 | 26.05 | 10 | 41 · 61 · 79 · 113 · 131 · 379 · 683 · 2729 · 4159 · 12959 | capped | capped |
| 24 | 373851853012084355680131361 | 26.57 | 11 | 37 · 71 · 103 · 109 · 131 · 181 · 239 · 419 · 461 · 3079 · 3761 | capped | capped |
| 25 | 419746116449969951645203201 | 26.62 | 12 | 17 · 37 · 43 · 53 · 71 · 109 · 151 · 191 · 281 · 881 · 1409 · 3761 | capped | capped |
| 26 | 480650413824292471675755841 | 26.68 | 11 | 29 · 41 · 89 · 131 · 197 · 239 · 271 · 727 · 881 · 1429 · 2969 | capped | capped |
| 27 | 507710311713618241797676801 | 26.71 | 11 | 41 · 47 · 71 · 89 · 127 · 151 · 241 · 967 · 1103 · 2089 · 4049 | capped | capped |
| 28 | 571197082724681188312817281 | 26.76 | 11 | 29 · 43 · 67 · 71 · 131 · 211 · 307 · 953 · 1021 · 1429 · 8161 | capped | capped |
| 29 | 598764213578029575774297601 | 26.78 | 13 | 23 · 29 · 43 · 53 · 79 · 109 · 113 · 127 · 181 · 239 · 379 · 433 · 449 | capped | capped |
| 30 | 667882655094744395121816001 | 26.82 | 10 | 67 · 83 · 131 · 163 · 199 · 239 · 461 · 2393 · 5167 · 20747 | capped | capped |
| 31 | 768927403921604814222336001 | 26.89 | 10 | 43 · 67 · 137 · 251 · 271 · 349 · 769 · 1427 · 4159 · 17981 | capped | capped |
| 32 | 922904786964947393227420801 | 26.97 | 10 | 17 · 37 · 41 · 311 · 929 · 1217 · 1301 · 1601 · 3191 · 15313 | capped | capped |
| 33 | 946295515632681859275671041 | 26.98 | 11 | 29 · 43 · 89 · 109 · 131 · 181 · 229 · 991 · 1427 · 1429 · 7129 | capped | capped |
| 34 | 1941211286859256920630460801 | 27.29 | 12 | 31 · 37 · 41 · 43 · 89 · 101 · 151 · 199 · 307 · 463 · 2089 · 11969 | capped | capped |
| 35 | 2119047331560739293908113921 | 27.33 | 12 | 17 · 37 · 43 · 47 · 79 · 89 · 113 · 127 · 571 · 911 · 3541 · 8969 | capped | capped |
| 36 | 2142752298785201165308646401 | 27.33 | 11 | 23 · 37 · 83 · 113 · 163 · 191 · 601 · 701 · 911 · 3079 · 7297 | capped | capped |
| 37 | 2325041445996414342573536161 | 27.37 | 11 | 19 · 29 · 47 · 71 · 79 · 137 · 557 · 619 · 2729 · 10711 · 11593 | capped | capped |
| 38 | 3213322593860118229475166721 | 27.51 | 10 | 37 · 61 · 67 · 71 · 419 · 433 · 617 · 1103 · 11593 · 209089 | capped | capped |
| 39 | 4635427415991948819511107841 | 27.67 | 13 | 23 · 29 · 37 · 53 · 71 · 89 · 113 · 127 · 191 · 271 · 379 · 647 · 3079 | capped | capped |
| 40 | 5466953848346666426169619201 | 27.74 | 11 | 37 · 59 · 79 · 109 · 113 · 151 · 229 · 1217 · 2311 · 4523 · 5851 | capped | capped |
| 41 | 6037031106243073961934109201 | 27.78 | 12 | 23 · 37 · 89 · 109 · 131 · 151 · 181 · 211 · 281 · 397 · 1483 · 5851 | capped | capped |
| 42 | 6473975735726053723736568001 | 27.81 | 13 | 23 · 31 · 53 · 71 · 79 · 97 · 131 · 151 · 181 · 239 · 419 · 593 · 1481 | capped | capped |
| 43 | 7006270310872326698183179201 | 27.85 | 10 | 41 · 71 · 89 · 103 · 199 · 211 · 5851 · 6359 · 8161 · 20593 | capped | capped |
| 44 | 8601689148505906566116203201 | 27.93 | 10 | 59 · 71 · 109 · 241 · 349 · 521 · 701 · 2053 · 2843 · 105071 | capped | capped |
| 45 | 9799274291299994931676646401 | 27.99 | 12 | 31 · 41 · 53 · 71 · 101 · 103 · 199 · 239 · 379 · 449 · 4159 · 5851 | capped | capped |

## 5. Checks

| check | result |
|---|---|
| all three runs complete, every shard recorded once | ✓ |
| engine hashes unchanged from before the runs | ✓ |
| capped run built and run with LASTCALC (`"lastcalc":1`) and no experimental skips (`skipmask 0`) | ✓ |
| no `UNPROVEN-PRIMALITY` flags | ✓ |
| the capped part finds all 44 terms its theorem guarantees | ✓ |
| the tail finds the one term only it guarantees ($2.088\times 10^{24}$) | ✓ |
| the union equals the 45 terms of the earlier capped $10^{30}$ run that lie below $10^{28}$: nothing missing, nothing new. That run was conditional (LASTCALC off), so this is a consistency check only: the certificate does not depend on it | ✓ |
| agreement with the proved census below $10^{25}$ (14 terms) | ✓ |

**Before the run**, each engine was validated against the proved census below $10^{25}$:
- **capped with LASTCALC:** at $T = 300$ and $T = 800$, it found every term its theorem guarantees;
- **pair tail:** at $T = 10^4$, it found exactly the two terms with second-largest prime above $T$, before and after the sieve was added.

These validations are described in §7 of each proof.

**After the run**, an independent review of the proofs (GLM-5.3) found no mathematical error. It placed the remaining risk in the gap between the algorithm the proofs describe and the program that ran, and suggested three cheap closures. All three were carried out on 7 October 2026, with the scripts in the repository's `tests/`:

| closure | what it shows | result |
|---|---|---|
| standalone term check (`check_terms.py`, standard library only; a PARI/GP version in `check_terms.gp`, not yet run) | each of the 45 terms satisfies the definition, with factorisations recomputed from scratch by trial division (no repository code used) | ✓ all 45; all rigid |
| production settings on the proved census (`check_production_settings.py`) | at $X = 10^{25}$ with the certificate's own $T = 15{,}000$: the capped half with LASTCALC finds all 13 terms it guarantees (and the 14th through LIST); the pair tail finds exactly $2.088\times 10^{24}$ and nothing else | ✓ |
| shard replay (`replay_shards.py`) | randomly chosen production shards, re-run alone on one thread (no concurrency, no donation), reproduce their recorded counters exactly | ✓ 100 capped, 60 pair and 60 $(q,a)$ shards (uniform sample); and 30 capped shards sampled in proportion to their work, holding 0.48% of all capped nodes, where donation happens |

The replay targets the one part of the argument the proofs state informally: that the parallel, resumable, work-donating runs do exactly the work of a sequential traversal (LastCalc proof §4, tail proof §5). Every counter of every sampled shard matched. That covers, for example, nodes, LIST candidates and Fermat passes in the capped engine, and pairs, candidates and survivors in the tail.

## 6. Reproducing

```bash
cd cn2-census-2
python3 cn2_split_campaign.py init --X 1e28 --R0 15000 --dir DIR --capped-only --lastcalc
# pair tail (copy and compile cn2x/cn2pair.c into DIR/pair first)
cd DIR/pair && CN2X_MODE=howe ./cn2pair 10000000000000000000000000000 15000 24500000 . 30 0 4096 60 0
cd DIR/tail && CN2X_MODE=howe ./cn2tail 10000000000000000000000000000 24500000 . 30 0 268435456 60 95
cd cn2-census-2 && python3 cn2_split_campaign.py run --dir DIR --now --threads 30
python3 cn2_split_campaign.py verify --dir DIR
python3 cn2_tail_verify.py --X 1e28 --pair DIR/pair --qa DIR/tail --T 15000 --reference results/b175531.txt
```

Every run is resumable: rerun the same command and it skips the shards already recorded. A run counts towards a certificate only once its final record shows `"complete":true`.

## 7. Cost, and what it means for $10^{30}$

The certificate took about 25.6 hours of machine time, in the three runs above. Choosing $T$ trades the two halves against each other:
- **the capped part** grows roughly like $T^2$;
- **the pair tail's candidates** fall like $T^{-4}$;
- **the rest of the tail** does not depend on $T$: the pair setups (about $\pi(q^*)^2/2$) and the $(q, a)$ part above $q^*$. Together they cost about 8 h here, and grow like $\sqrt X$.

Scaling these measured rates to $10^{30}$ gives:
- **$T = 2\times 10^4$:** about 190 h capped, 160 h of pair candidates and 80 h for the rest of the tail, roughly **430 h**;
- **$T = 3\times 10^4$:** about 420 h capped, 30 h of pair candidates and 80 h for the rest of the tail, roughly **530 h**.

So a certificate to $10^{30}$ would take about two and a half to three weeks of continuous running: a possible project. The capped figures assume the $T^2$ scaling, which is the least certain part of the estimate.
