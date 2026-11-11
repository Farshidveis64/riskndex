"""Development analysis for the hybrid mean chain and the selection margin (dev seeds 100-102 only).

Reference: RIS (ind. RR) on the same seed.  Gain = realized CVaR / reference CVaR - 1 (paired, same
evaluation worlds), averaged over cases with alpha <= 0.2, split into k <= 5 and k >= 10.
Usage: python3 experiments/hyb_dev.py results/dev/hyb_rule.csv results/dev/hyb_b512.csv results/dev/hyb_b128.csv
"""
import sys
import numpy as np
import pandas as pd

COLS = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,sec,mean,cvar,uid,stats'.split(',')


def load(path):
    d = pd.read_csv(path, header=None, names=COLS)
    return d


def gains(d, ref='RIS-IndRR', amax=0.2):
    key = ['graph', 'model', 'seed', 'k', 'alpha']
    r = d[d.method == ref][key + ['cvar', 'mean']].rename(columns={'cvar': 'cref', 'mean': 'mref'})
    x = d.merge(r, on=key)
    x = x[x.alpha <= amax + 1e-9].copy()
    x['g'] = 100 * (x.cvar / x.cref - 1)
    x['win'] = np.sign(np.round(x.cvar - x.cref, 6))
    x['kb'] = np.where(x.k <= 5, 'k<=5', 'k>=10')
    return x


def summarize(path):
    d = load(path)
    x = gains(d)
    tab = x.groupby(['method', 'kb']).g.mean().unstack()
    wtl = x.groupby(['method', 'kb']).win.apply(lambda s: '%d/%d/%d' % ((s > 0).sum(), (s == 0).sum(), (s < 0).sum())).unstack()
    out = tab.join(wtl, rsuffix=' w/t/l')
    return out, x


if __name__ == '__main__':
    res = {}
    for p in sys.argv[1:]:
        name = p.split('_')[-1].replace('.csv', '')
        out, x = summarize(p)
        res[name] = (out, x)
        print('==', p, '(gain over RIS ind., %, alpha<=0.2; w/t/l vs RIS)')
        print(out.round(2).to_string())
        # per graph-model for the main candidates
        piv = x[x.method.isin(['RiskIndex', 'RiskIndex-H', 'RiskIndex-H-z1', 'RiskIndex-HV', 'RiskIndex-HV-z1', 'RiskIndex-V'])].groupby(['graph', 'model', 'kb', 'method']).g.mean().unstack()
        print(piv.round(2).to_string())
    # decision rule (see scripts/dev/README.md)
    if {'rule', 'b512'} <= set(res):
        pre = 'RiskIndex-HV' if 'RiskIndex-HV' in res['rule'][0].index else 'RiskIndex-H'
        cand = {pre: 0.0}
        for z in [0.5, 1, 2, 3]:
            if pre + '-z%g' % z in res['rule'][0].index: cand[pre + '-z%g' % z] = float(z)
        feas = []
        for m, z in cand.items():
            ok = all(res[s][0].loc[m, 'k>=10'] >= -0.1 for s in ['rule', 'b512'])
            print('candidate %-18s z=%.1f  k>=10 rule %.2f b512 %.2f  k<=5 rule %.2f  feasible=%s' % (
                m, z, res['rule'][0].loc[m, 'k>=10'], res['b512'][0].loc[m, 'k>=10'], res['rule'][0].loc[m, 'k<=5'], ok))
            if ok:
                feas.append((res['rule'][0].loc[m, 'k<=5'], z, m))
        if feas:
            best = max(f[0] for f in feas)
            pick = min([f for f in feas if f[0] >= best - 0.1], key=lambda f: f[1])
            print('DECISION:', pick[2], '(z=%.1f, k<=5 gain %.2f)' % (pick[1], pick[0]))
        else:
            print('DECISION: no feasible variant')
