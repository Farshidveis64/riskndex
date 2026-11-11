#!/usr/bin/env python3
"""LargeGraph diagnosis: tail headroom over RIS at small k (development seed 100).
Gain (%) in realized CVaR_alpha on fresh worlds over RIS (ind. RR), averaged over k in {1,2,3,5} and
alpha in {0.05,0.1,0.2}, for RiskIndex-HV and two forward-simulated pool references: Oracle-CVaR (myopic greedy
on the exact empirical CVaR) and Oracle-basis (exact truncated-spread chains + selection by exact CVaR).
cert1: share of alpha at k=1 where no node outside the pool can beat the best found single seed (its CVaR exceeds
the largest mean of the 20 best-covered nodes outside the pool).
Usage: lg_headroom.py <headroom csv>..."""
import sys, pandas as pd, numpy as np
cols = ['graph','model','seed','M','b','method','k','alpha','mean','c05','c10','c20','h','rho','S']
d = pd.concat([pd.read_csv(f, header=None, names=cols) for f in sys.argv[1:]])
acol = {0.05: 'c05', 0.1: 'c10', 0.2: 'c20'}
rows = []
for (g, mo, sd), x0 in d.groupby(['graph','model','seed']):
    edge = x0[x0.method == 'PoolEdge']['mean'].max() if (x0.method == 'PoolEdge').any() else np.nan
    for k, x in x0[x0.method != 'PoolEdge'].groupby('k'):
        ris = x[x.method == 'RIS-IndRR'].iloc[0]
        for al, c in acol.items():
            sel = lambda m: x[(x.method == m) & (np.isclose(x.alpha, al))]
            hv, oc = sel('RiskIndex-HV').iloc[0], sel('Oracle-CVaR').iloc[0]
            ob = sel('Oracle-basis').iloc[0] if len(sel('Oracle-basis')) else oc
            base = ris[c]; bestv = max(hv[c], oc[c], ob[c], base)
            rows.append(dict(graph=g, model=mo, k=k, alpha=al, b=ris.b, tail=ris['c05'] / ris['mean'],
                             hv=100 * (hv[c] / base - 1), oc=100 * (oc[c] / base - 1), ob=100 * (ob[c] / base - 1),
                             best=100 * (bestv / base - 1), rho=hv.rho, thr=hv.S != ris.S,
                             cert=(bestv > edge) if (k == 1 and not np.isnan(edge)) else np.nan))
r = pd.DataFrame(rows)
order = ['email','wiki','facebook','hepph','epinions','slashdot','dblp','youtube','livejournal']
r['go'] = r.graph.map({g: i for i, g in enumerate(order)})
agg = r.groupby(['go','graph','model']).agg(b=('b','first'), tail1=('tail', lambda s: s.iloc[0]), HV=('hv','mean'), OrcCVaR=('oc','mean'),
      OrcBasis=('ob','mean'), best=('best','mean'), rho=('rho','mean'), thr=('thr','mean'), cert1=('cert','mean')).reset_index().drop(columns='go')
pd.set_option('display.width', 220)
print(agg.round(2).to_string(index=False))
for nm, gs in (('small', order[:4]), ('large', order[4:])):
    x = r[r.graph.isin(gs)]
    if len(x): print(f'{nm}: HV {x.hv.mean():.2f}  Oracle-CVaR {x.oc.mean():.2f}  Oracle-basis {x.ob.mean():.2f}  best {x.best.mean():.2f}  (n={len(x)})')
