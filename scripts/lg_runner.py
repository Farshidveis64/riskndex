#!/usr/bin/env python3
"""Resumable runner of the ablation "w/o need gate" (forward layer with the cost gate only, fltheta unset; label RiskIndex-LG; result directories *_lg).

Every job is one invocation of the binary with the arguments of the corresponding script (run_quality_lg.sh,
run_headroom_lg.sh, run_oracle_lg.sh, run_sweep_lg.sh, run_scale_lg.sh), split into one sketch seed per job where
the script passes several (every seed derives its own random streams, so the split changes no result; the
reproduction checks in merge_ad.py and make_tables_ad.py verify this against the runs without the forward layer, results/*_v2).
Before a job starts, the sizes of its output files are journaled; if the job does not finish (machine restart),
the next start truncates the files back to those sizes and reruns the job.  Finished jobs are listed in
results/lg_jobs/<queue>.done.
Usage: lg_runner.py A|B|C  (queue A: Tier B DBLP, headroom, Tier C YouTube and LiveJournal; queue B: Tier B small
graphs, Tier A, size sweep, Tier C part B; queue C: Tier C DBLP)"""
import os, sys, json, subprocess, time
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
B = os.environ.get('RIBIN', 'bin/riskindex'); D = 'datasets/raw'
J = 'results/lg_jobs'; os.makedirs(J, exist_ok=True)
for d in ['quality_lg', 'headroom_lg', 'oracle_lg', 'sweep_lg', 'scale_lg']: os.makedirs(f'results/{d}', exist_ok=True)
LG = 'hyb=1 val=1 hvz=0 rhomin=0.5 sq=1 fl=1 rawidx=0 direct=0 sg=0 sgc=0 inddirect=0 nomat=0'.split()
P = dict(email=('email-Eu-core.txt', 'directed', 0.01622), wiki=('wiki-Vote.txt', 'directed', 0.02215),
         facebook=('facebook_combined.txt', 'undirected', 0.00616), hepph=('ca-HepPh.txt', 'undirected', 0.00408),
         epinions=('soc-Epinions1.txt', 'directed', 0.00939), slashdot=('soc-Slashdot0902.txt', 'directed', 0.00808),
         dblp=('com-dblp.ungraph.txt', 'undirected', 0.00863))


def quality(part, g, mo, seed, R):
    f, dr, pc = P[g]; p = pc if mo == 'IC' else 0
    out = f'results/quality_lg/quality_{part}.csv'; sysf = f'results/quality_lg/quality_{part}_sys.csv'
    cmd = [B, 'quality', f'graph=snap:{D}/{f}:{dr}', f'name={g}', f'model={mo}', f'p={p}', 'M=1024', 'brule=100', 'bmax=65536', 'H=8',
           'Kmax=100', 'ks=1,2,3,5,10,25,50,100', 'alphas=0.05,0.1,0.2,0.5,1.0', f'seeds={seed}', f'R={R}'] + LG + [f'out={out}', f'sysout={sysf}', 'tag=tierB']
    return dict(key=f'quality_{part}:{g}:{mo}:{seed}', cmd=cmd, files=[out, out + '.cert', out + '.hcert', sysf], log=f'results/quality_lg_{part}.log')


def headroom(g, mo, seed=0):
    f, dr, pc = P[g]; p = pc if mo == 'IC' else 0; out = 'results/headroom_lg/headroom.csv'
    cmd = [B, 'headroom', f'graph=snap:{D}/{f}:{dr}', f'name={g}', f'model={mo}', f'p={p}', 'M=1024', 'brule=100', 'bmax=65536', f'seed={seed}',
           'P=60', 'Pscreen=2000', 'Wsel=1000', 'Weval=3000', 'fl=1', f'out={out}']
    return dict(key=f'headroom:{g}:{mo}:{seed}', cmd=cmd, files=[out], log='results/headroom_lg.log')


def oracle(spec, mo, p, k, tag):
    out = 'results/oracle_lg/oracle_lg.csv'
    cmd = [B, 'oracle', f'graph={spec}', f'model={mo}', f'p={p}', f'k={k}', 'alphas=0.05,0.1,0.2,0.5,1.0', 'M=4096', 'b=64', 'H=8', f'Kmax={k}',
           'seeds=0,1,2,3,4', 'RA=5000', 'RB=5000', 'lite=1', 'hyb=1', 'val=1', 'rhomin=0.5', 'fl=1', f'out={out}', f'tag={tag}']
    return dict(key=f'oracle:{spec}:{mo}:{p}:{k}:{tag}', cmd=cmd, files=[out], log='results/oracle_lg.log')


def sweep(fct, mo, seed):
    out = 'results/sweep_lg/sweep.csv'; sysf = 'results/sweep_lg/sweep_sys.csv'
    cmd = [B, 'quality', f'graph=pad:{D}/ca-HepPh.txt:undirected:{fct}', f'name=hepph-pad{fct}', f'model={mo}', 'M=1024', 'brule=100', 'bmax=65536',
           'H=8', 'Kmax=100', 'ks=1,2,3,5,10', 'alphas=0.05,0.1,0.2,0.5,1.0', f'seeds={seed}', 'R=3000'] + LG + [f'out={out}', f'sysout={sysf}', 'tag=sweep']
    return dict(key=f'sweep:{fct}:{mo}:{seed}', cmd=cmd, files=[out, out + '.cert', out + '.hcert', sysf], log='results/sweep_lg.log')


def scale(part, graph, name, mo, seed, tag, ks='10,25,50,100', alphas='0.05,0.1,0.2,0.5,1.0'):
    out = f'results/scale_lg/scale_{part}.csv'; ho = f'results/scale_lg/scale_hyb_{part}.csv'
    cmd = [B, 'scale', f'graph={graph}', f'name={name}', f'model={mo}', 'M=512', 'brule=100', 'bmax=16384', 'H=8', 'Kmax=100',
           'hyb=1', 'val=1', 'rhomin=0.5', 'sq=1', 'fl=1', 'direct=0', f'ks={ks}', f'alphas={alphas}', f'seeds={seed}', 'Q=1000', 'R=0',
           f'out={out}', 'qout=/dev/null', f'hout={ho}', f'tag={tag}']
    return dict(key=f'scale:{name}:{mo}:{seed}', cmd=cmd, files=[out, ho], log='results/scale_lg.log')


def jobs(q):
    L = []
    if q == 'A':
        for mo, R in [('WC', 2000), ('LT', 2000), ('IC', 2000)]:
            for s in range(3): L.append(quality('dblp', 'dblp', mo, s, R))
        for g in ['email', 'wiki', 'facebook', 'hepph', 'epinions', 'slashdot', 'dblp']:
            for mo in ['WC', 'LT', 'IC']: L.append(headroom(g, mo))
        # (DBLP's Tier C jobs moved to queue C when queue B finished early, to use the idle core)
        for g, f, dr in [('youtube', 'com-youtube.ungraph.txt', 'undirected')]:
            for mo in ['WC', 'LT']:
                for s in range(3): L.append(scale('A', f'snap:{D}/{f}:{dr}', g, mo, s, 'real'))
        L.append(scale('A', f'snap:{D}/soc-LiveJournal1.txt:directed', 'livejournal', 'WC', 0, 'real', ks='10,50,100', alphas='0.05,0.2,1.0'))
    elif q == 'C':
        for mo in ['WC', 'LT']:
            for s in range(3): L.append(scale('C', f'snap:{D}/com-dblp.ungraph.txt:undirected', 'dblp', mo, s, 'real'))
    elif q == 'B':
        for g, n in [('email', 5), ('wiki', 5), ('facebook', 5), ('hepph', 5), ('epinions', 3), ('slashdot', 3)]:
            for mo in ['WC', 'LT', 'IC']:
                for s in range(n): L.append(quality('small', g, mo, s, 5000 if n == 5 else 3000))
        for inst in range(1, 6):
            for spec in [f'er:48:4:{inst}', f'er:48:8:{inst}', f'ba:48:2:{inst}', f'ba:48:4:{inst}', f'sbm:48:4:0.5:0.02:{inst}', f'sbm:48:4:0.3:0.08:{inst}',
                         f'ws:48:4:0.2:{inst}', f'ws:48:6:0.5:{inst}', f'er:64:4:{inst}', f'ba:64:2:{inst}', f'sbm:64:4:0.4:0.02:{inst}']:
                for p in ['0.1', '0.2']: L.append(oracle(spec, 'IC', p, 3, f'IC{p}'))
                L.append(oracle(spec, 'WC', 0, 3, 'WC')); L.append(oracle(spec, 'LT', 0, 3, 'LT'))
            for spec in [f'er:48:4:{inst}', f'ba:48:2:{inst}', f'sbm:48:4:0.5:0.02:{inst}', f'ws:48:4:0.2:{inst}']:
                L.append(oracle(spec, 'IC', '0.2', 4, 'IC0.2'))
        for fct in [1, 5, 20, 100]:
            for mo in ['WC', 'LT']:
                for s in range(3): L.append(sweep(fct, mo, s))
        for n in [10000, 100000, 1000000]:
            for s in range(3): L.append(scale('B', f'ba:{n}:5:1', f'ba{n}', 'WC', s, 'ba'))
        for g, f, dr in [('email', 'email-Eu-core.txt', 'directed'), ('wiki', 'wiki-Vote.txt', 'directed'), ('epinions', 'soc-Epinions1.txt', 'directed'),
                         ('slashdot', 'soc-Slashdot0902.txt', 'directed')]:
            for mo in ['WC', 'LT']:
                for s in range(3): L.append(scale('B', f'snap:{D}/{f}:{dr}', g, mo, s, 'real'))
    return L


def size(f):
    return os.path.getsize(f) if os.path.exists(f) else 0


def main(q):
    donef = f'{J}/{q}.done'; jf = f'{J}/{q}.journal'
    done = set(open(donef).read().split()) if os.path.exists(donef) else set()
    if os.path.exists(jf):  # interrupted job: roll its outputs back
        jr = json.load(open(jf))
        if jr['key'] not in done:
            for f, sz in jr['sizes'].items():
                if os.path.exists(f) and os.path.getsize(f) > sz:
                    with open(f, 'r+b') as fh: fh.truncate(sz)
                    print(f'[runner {q}] rolled back {f} to {sz} bytes (interrupted {jr["key"]})', flush=True)
        os.remove(jf)
    L = jobs(q)
    todo = [j for j in L if j['key'] not in done]
    print(f'[runner {q}] {len(L)} jobs, {len(L) - len(todo)} done, {len(todo)} to run', flush=True)
    for i, j in enumerate(todo):
        with open(jf + '.tmp', 'w') as fh:
            json.dump(dict(key=j['key'], sizes={f: size(f) for f in j['files']}), fh); fh.flush(); os.fsync(fh.fileno())
        os.replace(jf + '.tmp', jf)
        t = time.time()
        with open(j['log'], 'a') as lg:
            lg.write(f'=== {j["key"]}\n'); lg.flush()
            r = subprocess.run(j['cmd'], stdout=lg, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            print(f'[runner {q}] FAILED {j["key"]} (exit {r.returncode}); stopping', flush=True); sys.exit(1)
        with open(donef, 'a') as fh:
            fh.write(j['key'] + '\n'); fh.flush(); os.fsync(fh.fileno())
        os.remove(jf)
        print(f'[runner {q}] {i + 1}/{len(todo)} {j["key"]} {time.time() - t:.0f}s', flush=True)
    print(f'[runner {q}] ALL_DONE', flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
