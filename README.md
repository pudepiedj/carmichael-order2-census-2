# Carmichael numbers of order 2: a certified census below 10^28

This repository contains the programs, the run records, the proofs and the validation suite behind this result:

> **There are exactly 45 Carmichael numbers of order 2 below 10^28.**

"Order 2" is Howe's definition, the one used by OEIS [A175531](https://oeis.org/A175531): an odd composite $n$ such that, for every prime $p \mid n$,
$$n \equiv 1 \ \text{ or } \ n \equiv p \pmod{p^2-1}.$$

**Status.** The computation is complete and every check passes (`docs/CN2_Certificate_1e28.md`). The result rests on two proofs, `docs/CN2_LastCalc_Coverage_Proof.md` and `docs/CN2_Tail_Coverage_Proof.md`. They have been checked by the authors and by an independent review, whose suggested further checks have all been carried out (certificate §5).

This repository extends the census below $10^{25}$ (14 terms) of [carmichael-order2-census](https://github.com/pudepiedj/carmichael-order2-census). That repository is the record behind the original OEIS submission. Its first 14 terms are identical to the first 14 here.

## The terms

Also in `results/b175531.txt`, with factorisations in `results/factorisations.txt` and `results/factorisations.json`.

| n | a(n) | log10 | prime factors | largest |
|---|---|---|---|---|
| 1 | 443372888629441 | 14.65 | 8 | 331 |
| 2 | 39671149333495681 | 16.60 | 9 | 191 |
| 3 | 842526563598720001 | 17.93 | 8 | 2729 |
| 4 | 2380296518909971201 | 18.38 | 9 | 991 |
| 5 | 3188618003602886401 | 18.50 | 8 | 3457 |
| 6 | 4208895375600667752001 | 21.62 | 10 | 5851 |
| 7 | 1159954316194989017102401 | 24.06 | 9 | 4759 |
| 8 | 2088144166339753513992001 | 24.32 | 8 | 17291 |
| 9 | 2196407820059694924883201 | 24.34 | 9 | 13441 |
| 10 | 3339611018825185787482801 | 24.52 | 12 | 1429 |
| 11 | 4105879060352839139462401 | 24.61 | 8 | 7237 |
| 12 | 5002862939121639632040001 | 24.70 | 11 | 2549 |
| 13 | 7865064643837556041286401 | 24.90 | 11 | 15809 |
| 14 | 9400084864021826054720641 | 24.97 | 11 | 4523 |
| 15 | 12071465216556111317952001 | 25.08 | 12 | 811 |
| 16 | 31280524488319495535546401 | 25.50 | 11 | 47431 |
| 17 | 32492501043239735947612801 | 25.51 | 11 | 3769 |
| 18 | 36991212507320386599648001 | 25.57 | 11 | 3041 |
| 19 | 44774022865681509785985601 | 25.65 | 11 | 63799 |
| 20 | 50473184681492173584023041 | 25.70 | 11 | 2311 |
| 21 | 54475112449147532388966721 | 25.74 | 10 | 11593 |
| 22 | 76587257660615840312697601 | 25.88 | 11 | 40151 |
| 23 | 111356194710608262218643841 | 26.05 | 10 | 12959 |
| 24 | 373851853012084355680131361 | 26.57 | 11 | 3761 |
| 25 | 419746116449969951645203201 | 26.62 | 12 | 3761 |
| 26 | 480650413824292471675755841 | 26.68 | 11 | 2969 |
| 27 | 507710311713618241797676801 | 26.71 | 11 | 4049 |
| 28 | 571197082724681188312817281 | 26.76 | 11 | 8161 |
| 29 | 598764213578029575774297601 | 26.78 | 13 | 449 |
| 30 | 667882655094744395121816001 | 26.82 | 10 | 20747 |
| 31 | 768927403921604814222336001 | 26.89 | 10 | 17981 |
| 32 | 922904786964947393227420801 | 26.96 | 10 | 15313 |
| 33 | 946295515632681859275671041 | 26.98 | 11 | 7129 |
| 34 | 1941211286859256920630460801 | 27.29 | 12 | 11969 |
| 35 | 2119047331560739293908113921 | 27.33 | 12 | 8969 |
| 36 | 2142752298785201165308646401 | 27.33 | 11 | 7297 |
| 37 | 2325041445996414342573536161 | 27.37 | 11 | 11593 |
| 38 | 3213322593860118229475166721 | 27.51 | 10 | 209089 |
| 39 | 4635427415991948819511107841 | 27.67 | 13 | 3079 |
| 40 | 5466953848346666426169619201 | 27.74 | 11 | 5851 |
| 41 | 6037031106243073961934109201 | 27.78 | 12 | 5851 |
| 42 | 6473975735726053723736568001 | 27.81 | 13 | 1481 |
| 43 | 7006270310872326698183179201 | 27.84 | 10 | 20593 |
| 44 | 8601689148505906566116203201 | 27.93 | 10 | 105071 |
| 45 | 9799274291299994931676646401 | 27.99 | 12 | 5851 |

The census also shows:
- **All 45 terms are rigid:** $n \equiv 1 \pmod{p^2-1}$ for every $p \mid n$. So the smallest term using the residue $p$ (Howe 2000, §5, gives one of about $3.92\times 10^{59}$) exceeds $10^{28}$.
- **No term below $10^{28}$ is divisible by 5 or 7.**
- **The largest prime factor of any term below $10^{28}$ is 209,089.**

## How the search works

Every order-2 $n < 10^{28}$, with primes $p_1 < \dots < p_k$, falls into one of three parts. Each part has its own engine and its own completeness theorem, and the split uses $T = 15{,}000$ and $q^* = 24{,}500{,}000$.

| part | which $n$ | engine | why it misses nothing |
|---|---|---|---|
| capped | $p_{k-1} \le T$ | `cn2x/cn2xh.c`: depth-first search over a prime table cut at $T$, with `CN2X_LASTCALC=1` | `docs/CN2_LastCalc_Coverage_Proof.md` |
| pair tail | $p_{k-1} > T$, $p_k \le q^*$ | `cn2x/cn2pair.c`: for each pair $T < p < q \le q^*$, $n$ lies in one residue class modulo $pq\cdot\mathrm{lcm}(p^2-1, q^2-1)$ | `docs/CN2_Tail_Coverage_Proof.md`, Lemma B |
| $(q, a)$ tail | a prime above $q^*$ | `cn2x/cn2tail.c`: $n = q\,(q + a(q^2-1))$ or $q\,(1 + a(q^2-1))$ | `docs/CN2_Tail_Coverage_Proof.md`, Lemma A |

The capped search needs no assumption about where its LIST step switches in. With LASTCALC on, it is proved to find every $n$ whose primes other than the largest lie in its table. Earlier capped searches lacked this, so their results above the proved bound were conditional. Every candidate from every engine is factorised and checked against the definition (`cn2_order2_defs.py`).

## Checking it

```bash
pip install -r requirements.txt
python3 cn2_make_results.py --check      # the results files rebuild exactly from the run records (tail re-verified)
python3 tests/validate_cert.py           # six checks, about 3 minutes on 30 threads
```

Three further checks answer an independent review of the proofs (`docs/CN2_Certificate_1e28.md` §5):

```bash
python3 tests/check_terms.py                 # the 45 terms against the definition, standard library only (< 1 s)
python3 tests/check_production_settings.py   # both halves at T = 15,000 on the proved census below 10^25 (~50 min)
python3 tests/replay_shards.py               # sampled production shards re-run alone on one thread (~15 min)
```

`tests/validate_cert.py` checks:
1. the results files;
2. the engine sources against the recorded hashes;
3. that every term lies in the residue class of its own top two primes;
4. the pair engine's candidate count against an independent Python count;
5. the capped search against the proved census below $10^{25}$;
6. the pair engine against the proved census below $10^{25}$.

## Reproducing the certificate

`docs/CN2_Certificate_1e28.md` §6 gives the exact commands. It took 25.6 hours on 30 threads: 12.9 h capped, 8.0 h pair tail and 4.7 h $(q, a)$ tail. The run records are in `results/cert_1e28/`:
- parameters (`split.json`);
- engine hashes (`ENGINES.sha256`);
- each run's completion record;
- the shard logs;
- the tail's raw survivors;
- the capped half's verified result.

The capped run's raw survivor file (156 MB) is not included.

## Layout

| path | contents |
|---|---|
| `cn2x/` | the three engines (C), exactly as built for the certificate |
| `cn2_split_campaign.py` | driver for the capped search (`init`, `run`, `verify`); also the original split census |
| `cn2_tail_verify.py`, `cn2_capped_peek.py` | verifiers for the tail and for a capped run in progress |
| `cn2_pair_tail.py` | pair-family prototype: unit test, exact counts, cost model |
| `cn2_make_results.py` | builds `results/` from the run records |
| `cn2_order2_defs.py` | the definitions (rigid, Howe, Chebyshev) and the candidate verifier |
| `docs/` | the certificate, the two proofs, and the method note of the first census |
| `results/` | the terms, factorisations and the certificate's run records |
| `tests/validate_cert.py` | the validation suite |

## Licence

MIT; see `LICENSE`.
