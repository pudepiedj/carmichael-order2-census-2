#!/usr/bin/env python3
"""tests/replay_shards.py -- re-run a random sample of the certificate's shards alone, single-threaded, and compare.

    python3 tests/replay_shards.py [--capped 100] [--pair 60] [--tail 60] [--procs 30] [--seed 2026]

The production runs were parallel (30 threads), resumable, and in the capped engine work-donating.  The proofs argue,
but do not prove, that this machinery does exactly the work of a plain sequential traversal (CN2_LastCalc_Coverage_Proof
§4, CN2_Tail_Coverage_Proof §5).  This test samples that claim.  Each engine skips every shard listed in its done.log, so
replaying shard k alone needs no change to the engines: a fresh run directory gets a done.log naming every OTHER shard,
and the engine runs on ONE thread (no concurrency; no donation, which needs an idle second thread).  It then does shard
k sequentially and appends that shard's counters, which must equal the production record exactly: for the capped engine
nodes, children, last-prime tests, LIST nodes, LIST candidates, Fermat passes, factorisations, budget prunes and sieved
candidates; for the pair engine pairs, skipped pairs, classes, candidates, tests and passes; for cn2tail tests, sieved,
Fermat calls and passes.  The engines are built from this repository's sources with the production parameters
(results/cert_1e28/split.json and the commands in docs/CN2_Certificate_1e28.md section 6).
"""
import argparse, json, os, random, subprocess, sys, tempfile, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
from cn2_build import compile_cmd                                # noqa: E402

CERT = os.path.join(HERE, 'results', 'cert_1e28')
X = '10000000000000000000000000000'


def records(path, nfields):
    out = {}
    for l in open(path):
        t = l.split()
        if len(t) == nfields + 2 and t[0] == 'D':
            out.setdefault(int(t[1]), t[2:])
    return out


def replay(engine, exe, sid, nshards, nfields, args_fn, env, tmp):
    d = os.path.join(tmp, f"{engine}_{sid}")
    os.mkdir(d)
    zeros = ' '.join('0' * nfields)
    with open(os.path.join(d, 'done.log'), 'w') as fh:
        for i in range(nshards):
            if i != sid:
                fh.write(f"D {i} {zeros}\n")
    t0 = time.time()
    subprocess.run([exe] + args_fn(d), env=env, capture_output=True, text=True)
    # the prefilled done.log has no line for sid, so any line for sid is the one the engine appended
    new = [l.split()[2:] for l in open(os.path.join(d, 'done.log')) if l.split()[:2] == ['D', str(sid)]]
    return sid, (new[-1] if new else None), time.time() - t0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--capped', type=int, default=100); ap.add_argument('--pair', type=int, default=60)
    ap.add_argument('--tail', type=int, default=60); ap.add_argument('--procs', type=int, default=30)
    ap.add_argument('--seed', type=int, default=2026)
    ap.add_argument('--weighted', action='store_true',
                    help='sample shards with probability proportional to their work (the first counter: nodes, pairs or '
                         'tests) instead of uniformly; uniform sampling picks mostly tiny shards, while donation and '
                         'concurrency matter most in the heavy ones')
    a = ap.parse_args()
    split = json.load(open(os.path.join(CERT, 'split.json')))
    rnd = random.Random(a.seed)
    spec = {   # engine: (source, nfields in its D lines, production done.log, argument builder, extra env)
        'capped': ('cn2xh', 9, os.path.join(CERT, 'capped', 'done.log'),
                   lambda d: [X, d, '1', '0', str(split['ratio']), str(split['sieve']), str(split['mshard']),
                              str(split['batch']), str(split['donate_q']), '0', '0'],
                   dict(CN2X_TCAP=str(split['R0']), CN2X_LASTCALC='1')),
        'pair': ('cn2pair', 6, os.path.join(CERT, 'pair', 'done.log'),
                 lambda d: [X, '15000', '24500000', d, '1', '0', '4096', '0', '0'], {}),
        'tail': ('cn2tail', 4, os.path.join(CERT, 'tail', 'done.log'),
                 lambda d: [X, '24500000', d, '1', '0', '268435456', '0', '95'], {}),
    }
    want = dict(capped=a.capped, pair=a.pair, tail=a.tail)
    fails = 0
    with tempfile.TemporaryDirectory() as tmp:
        for engine, (src, nf, prod_log, args_fn, extra) in spec.items():
            if not want[engine]:
                continue
            exe = os.path.join(tmp, src)
            subprocess.run(compile_cmd(os.path.join(HERE, 'cn2x', src + '.c'), exe), check=True)
            prod = records(prod_log, nf)
            n = len(prod)
            assert sorted(prod) == list(range(n)), f"{engine}: production done.log is not shards 0..{n - 1}"
            k = min(want[engine], n)
            if a.weighted:                         # without replacement, weight = the shard's first counter
                keys = sorted(range(n), key=lambda i: rnd.random() ** (1.0 / max(int(prod[i][0]), 1e-9)), reverse=True)
                sample = sorted(keys[:k])
            else:
                sample = sorted(rnd.sample(range(n), k))
            share = sum(int(prod[i][0]) for i in sample) / max(1, sum(int(v[0]) for v in prod.values()))
            env = dict(os.environ, CN2X_MODE='howe', **extra)
            t0 = time.time()
            with ThreadPoolExecutor(a.procs) as ex:
                res = list(ex.map(lambda s: replay(engine, exe, s, n, nf, args_fn, env, tmp), sample))
            bad = [(s, prod[s], got) for s, got, _ in res if got != prod[s]]
            fails += len(bad)
            print(f"[{'PASS' if not bad else 'FAIL'}] {engine}: {len(sample)} of {n:,} shards replayed alone on one "
                  f"thread; counters identical for {len(sample) - len(bad)}  ({time.time() - t0:.0f} s; "
                  f"longest shard {max(r[2] for r in res):.0f} s; {'weighted' if a.weighted else 'uniform'} sample holding "
                  f"{share:.2%} of the run's work)", flush=True)
            for s, p, g in bad[:10]:
                print(f"    shard {s}: production {p}, replay {g}")
    print('ALL PASS' if not fails else f"FAILED: {fails} shards differ")
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
