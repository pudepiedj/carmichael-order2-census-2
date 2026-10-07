#!/usr/bin/env python3
r"""
Windowed, restartable campaign runner for cn2xc (checkpointed exhaustive CN2 search).

    python cn2x_campaign.py init   --X 1e24 --dir cn2x/runs/campaign_1e24
    python cn2x_campaign.py run    --dir cn2x/runs/campaign_1e24 --window 23:00-07:00 --threads 30
    python cn2x_campaign.py run    --dir cn2x/runs/campaign_1e24 --until 07:00        # one session, starting now
    python cn2x_campaign.py run    --dir cn2x/runs/campaign_1e24 --now --window 22:00-08:00
                                   # start now, no deadline; after a manual DRAIN resume nightly at 22:00
    python cn2x_campaign.py status --dir cn2x/runs/campaign_1e24
    python cn2x_campaign.py verify --dir cn2x/runs/campaign_1e24                       # after completion

Policy: a session runs until the window closes, then DRAINS -- it takes no new shards, lets every
in-flight shard finish and be recorded, and exits.  Nothing is abandoned, so the only cost of the
window is a short overrun past its end.  To reclaim the machine immediately instead, touch
<dir>/STOP (in-flight shards are abandoned and redone next session; results stay exact).
A manual `touch <dir>/DRAIN` ends a session gracefully at any time; the runner then waits for the
next window start rather than restarting inside the current one.

`init` freezes a copy of the engine binary and the parameters into the run directory; `run`
refuses to continue with anything different, because shard ids depend on them.
"""

import argparse
import datetime as dt
import re
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
SOURCE = os.path.join(HERE, 'cn2x', 'cn2xc.c')
sys.path.insert(0, HERE)
from cn2_exhaustive_c import _verify_survivor, icbrt, is_rigid_cn2, KNOWN_BELOW_1E22  # noqa: E402
from cn2_build import compile_cmd, keep_awake  # noqa: E402

COUNTERS = ('nodes', 'children', 'last_tests', 'list_nodes', 'list_candidates',
            'fermat_passes', 'factorisations', 'pruned_budget', 'sieved')


def parse_bound(text):
    """Exact integer from '1e24' or '1000...0'.  float() would give 999999999999999983222784
    for 1e24 (10^24 is not representable in binary64), silently shrinking the search."""
    t = str(text).strip().lower().replace('_', '')
    m = re.fullmatch(r'(\d+)(?:\.(\d+))?e\+?(\d+)', t)
    if m:
        whole, frac, exp = m.group(1), m.group(2) or '', int(m.group(3))
        if len(frac) > exp:
            raise ValueError(f"{text} is not an integer")
        return int(whole + frac) * 10 ** (exp - len(frac))
    if t.isdigit():
        return int(t)
    raise ValueError(f"cannot parse bound {text!r} exactly")


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def load_params(d):
    with open(os.path.join(d, 'params.json')) as fh:
        return json.load(fh)


def cmd_init(a):
    d = os.path.abspath(a.dir)
    if os.path.exists(os.path.join(d, 'params.json')):
        sys.exit(f"{d} already initialised; refusing to overwrite")
    os.makedirs(d, exist_ok=True)
    engine = os.path.join(d, 'cn2xc')
    subprocess.run(compile_cmd(SOURCE, engine), check=True)
    shutil.copy(SOURCE, os.path.join(d, 'cn2xc.c'))
    X = parse_bound(a.X)
    # node totals measured: 4.561e9 at 1e19, 2.743e11 at 1e22 -> x3.92 per decade
    expected_nodes = 4.561e9 * 3.92 ** (math.log10(X) - 19)
    params = dict(X=str(X), mshard=a.mshard, batch=a.batch, ratio=1.0, sieve=1, donate_q=10**10,
                  expected_nodes=expected_nodes,
                  engine_sha256=sha256(engine), created=dt.datetime.now().isoformat(timespec='seconds'))
    with open(os.path.join(d, 'params.json'), 'w') as fh:
        json.dump(params, fh, indent=1)
    print(f"initialised {d}\n  X = {X:.3e}  mshard = {a.mshard}  batch = {a.batch}\n"
          f"  predicted nodes {expected_nodes:.3e} (for the progress percentage and ETA)\n"
          f"  engine frozen: {params['engine_sha256'][:16]}...")


def next_time(hhmm, after):
    h, m = map(int, hhmm.split(':'))
    t = after.replace(hour=h, minute=m, second=0, microsecond=0)
    return t if t > after else t + dt.timedelta(days=1)


def in_window(now, start, end):
    s = now.replace(hour=int(start[:2]), minute=int(start[3:]), second=0, microsecond=0)
    e = now.replace(hour=int(end[:2]), minute=int(end[3:]), second=0, microsecond=0)
    return (s <= now < e) if s < e else (now >= s or now < e)


def run_session(d, p, threads, deadline_epoch, progress=60):
    engine = os.path.join(d, 'cn2xc')
    if sha256(engine) != p['engine_sha256']:
        sys.exit("engine binary in run directory has changed; refusing to resume")
    cmd = keep_awake([engine, p['X'], d, str(threads), f"{deadline_epoch:.0f}",
           str(p['ratio']), str(p['sieve']), str(p['mshard']), str(p['batch']), str(p['donate_q']), str(progress),
           str(p.get('expected_nodes', 0))])
    started = time.time()
    # engine progress goes to stderr and is inherited, so it appears live on screen (or in the log)
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, text=True)
    lines = [l for l in proc.stdout.splitlines() if l.startswith('{')]
    if not lines:
        print(f"engine produced no summary (exit {proc.returncode})", flush=True)
        return None
    s = json.loads(lines[-1])
    s.update(started=dt.datetime.fromtimestamp(started).isoformat(timespec='seconds'),
             ended=dt.datetime.now().isoformat(timespec='seconds'),
             deadline=dt.datetime.fromtimestamp(deadline_epoch).isoformat(timespec='seconds') if deadline_epoch else None,
             overrun_seconds=(time.time() - deadline_epoch) if deadline_epoch else None, threads=threads)
    with open(os.path.join(d, 'sessions.jsonl'), 'a') as fh:
        fh.write(json.dumps(s) + '\n')
    return s


def cmd_run(a):
    d = os.path.abspath(a.dir)
    p = load_params(d)
    now_first = a.now            # first session starts immediately and runs until DRAIN/STOP/complete
    skip_until = None            # after a manual DRAIN, don't restart before this time
    while True:
        now = dt.datetime.now()
        if now_first:
            end = None
        elif a.until:
            end = next_time(a.until, now)
        else:
            start, stop = a.window.split('-')
            if skip_until is not None or not in_window(now, start, stop):
                why = "manual drain" if skip_until else f"outside window {a.window}"
                wake = skip_until or next_time(start, now)
                skip_until = None
                print(f"{now:%Y-%m-%d %H:%M} {why}; sleeping until {wake:%Y-%m-%d %H:%M}", flush=True)
                while dt.datetime.now() < wake:            # short sleeps survive system sleep/wake
                    time.sleep(min(60, max(1, (wake - dt.datetime.now()).total_seconds())))
                continue
            end = next_time(stop, now)
        print(f"{dt.datetime.now():%Y-%m-%d %H:%M} session start; "
              + (f"drain at {end:%H:%M}" if end else "no deadline (touch DRAIN to end gracefully, STOP to end now)"), flush=True)
        s = run_session(d, p, a.threads, end.timestamp() if end else 0, a.progress)
        if s is None:
            sys.exit(1)
        print(f"{dt.datetime.now():%Y-%m-%d %H:%M} session end: complete={s['complete']} drained={s['drained']} "
              f"stopped={s['stopped']} shards done now {s['shards_done_now']:,} (total {s['shards_done_prior'] + s['shards_done_now']:,}) "
              f"overrun {s['overrun_seconds'] or 0:.0f}s", flush=True)
        if s['complete']:
            print("campaign complete -- run `verify`", flush=True)
            return
        if a.until or s['stopped'] or (now_first and not a.window):
            return
        # A drain that came before the deadline (or with none) was a manual DRAIN: the user wants the
        # machine, so wait for the NEXT window start rather than restarting inside the current window.
        if s['drained'] and (end is None or time.time() < end.timestamp() - 60):
            skip_until = next_time(a.window.split('-')[0], dt.datetime.now())
        now_first = False


def cmd_status(a):
    d = os.path.abspath(a.dir)
    p = load_params(d)
    sessions = [json.loads(l) for l in open(os.path.join(d, 'sessions.jsonl'))] if os.path.exists(os.path.join(d, 'sessions.jsonl')) else []
    done = sum(1 for l in open(os.path.join(d, 'done.log')) if l.startswith('D ')) if os.path.exists(os.path.join(d, 'done.log')) else 0
    hits = sum(1 for l in open(os.path.join(d, 'out.txt')) if l.startswith('H ')) if os.path.exists(os.path.join(d, 'out.txt')) else 0
    surv = sum(1 for l in open(os.path.join(d, 'out.txt')) if l.startswith('S ')) if os.path.exists(os.path.join(d, 'out.txt')) else 0
    hours = sum(s['search_seconds'] for s in sessions) / 3600
    print(f"campaign X = {int(p['X']):.3e}   sessions {len(sessions)}   engine time {hours:.1f} h")
    print(f"shards done {done:,}   path hits {hits}   survivors awaiting verification {surv:,}")
    if sessions:
        last = sessions[-1]
        print(f"last session {last['started']} -> {last['ended']}  complete={last['complete']} "
              f"drained={last['drained']} overrun {last.get('overrun_seconds') or 0:.0f}s  shards generated so far {last['shards_generated']:,}")


def cmd_verify(a):
    d = os.path.abspath(a.dir)
    p = load_params(d)
    X = int(p['X'])
    T = icbrt(X - 1)
    sessions = [json.loads(l) for l in open(os.path.join(d, 'sessions.jsonl'))]
    if not sessions or not sessions[-1]['complete']:
        print("warning: campaign not complete; verifying what exists so far")
    hits, survivors = set(), []
    for line in open(os.path.join(d, 'out.txt')):
        parts = line.split()
        if len(parts) == 2 and parts[0] == 'H':
            hits.add(int(parts[1]))
        elif len(parts) == 5 and parts[0] == 'S':
            survivors.append((int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]), T))
    survivors = list(set(survivors))
    t0 = time.time()
    with Pool(a.procs) as pool:
        for n in pool.imap_unordered(_verify_survivor, survivors, chunksize=512):
            if n is not None:
                hits.add(n)
    verified = sorted(n for n in hits if n <= X and is_rigid_cn2(n))
    by_decade = {}
    for n in verified:
        by_decade.setdefault(int(math.log10(n)), []).append(str(n))
    report = dict(X=p['X'], complete=bool(sessions and sessions[-1]['complete']), n_cn2=len(verified),
                  cn2=[str(n) for n in verified], by_decade={str(k): v for k, v in sorted(by_decade.items())},
                  known_below_1e22_all_found=all(n in verified for n in KNOWN_BELOW_1E22 if n <= X),
                  counters={c: sessions[-1][c] for c in COUNTERS} if sessions else None,
                  survivors_verified=len(survivors), verify_seconds=time.time() - t0)
    with open(os.path.join(d, 'result.json'), 'w') as fh:
        json.dump(report, fh, indent=1)
    print(f"{len(verified)} rigid CN2s <= {X:.3e} (complete={report['complete']}); "
          f"known six below 1e22 all found: {report['known_below_1e22_all_found']}")
    for k, v in sorted(by_decade.items()):
        print(f"  decade 10^{k}: {len(v)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    i = sub.add_parser('init'); i.add_argument('--X', required=True); i.add_argument('--dir', required=True)
    i.add_argument('--mshard', type=int, default=100000); i.add_argument('--batch', type=int, default=4096)
    r = sub.add_parser('run'); r.add_argument('--dir', required=True); r.add_argument('--threads', type=int, default=30)
    g = r.add_mutually_exclusive_group()
    g.add_argument('--window', help="HH:MM-HH:MM, repeated nightly until complete")
    g.add_argument('--until', help="HH:MM, a single session starting now")
    r.add_argument('--now', action='store_true',
                   help="start immediately with no deadline; end it with <dir>/DRAIN (graceful) or <dir>/STOP. "
                        "With --window, nightly sessions resume at the next window start afterwards")
    r.add_argument('--progress', type=float, default=60, help="seconds between progress lines (0 = off)")
    s = sub.add_parser('status'); s.add_argument('--dir', required=True)
    v = sub.add_parser('verify'); v.add_argument('--dir', required=True); v.add_argument('--procs', type=int, default=24)
    a = ap.parse_args()
    if a.cmd == 'run' and not (a.window or a.until or a.now):
        ap.error("run needs --window, --until or --now")
    dict(init=cmd_init, run=cmd_run, status=cmd_status, verify=cmd_verify)[a.cmd](a)


if __name__ == '__main__':
    main()
