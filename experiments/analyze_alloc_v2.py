#!/usr/bin/env python3
"""Budget allocation (index without the forward layer): M vs b at a fixed budget Mb on the four small graphs, and the deployed
configuration (Tier B, twice the budget).  Gains of the final index over RIS at the same budget, alpha <= 0.2."""
import os, numpy as np, pandas as pd
R = os.environ.get('RESDIR', 'results')
COLS = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,sec,mean,cvar,uid,stats'.split(',')
HC = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,h,cert,rho,mc_cert,se,se_w,raw,guarded'.split(',')
a = pd.read_csv(f'{R}/alloc_v2/alloc.csv', header=None, names=COLS)
tb = pd.read_csv(f'{R}/quality_v2/quality_small.csv', header=None, names=COLS)
tb = tb[tb.graph.isin(['email', 'wiki', 'facebook', 'hepph']) & tb.seed.isin([0, 1, 2])].assign(tag='deployed')
d = pd.concat([a, tb])
key = ['tag', 'graph', 'model', 'seed', 'k', 'alpha']
ref = d[d.method == 'RIS-IndRR'][key + ['cvar']].rename(columns={'cvar': 'ris'})
d = d.merge(ref, on=key); d['g'] = (d.cvar / d.ris - 1) * 100
x = d[(d.alpha <= 0.2) & (d.method == 'RiskIndex-HV')].assign(kg=lambda z: np.where(z.k <= 5, 'k<=5', 'k>=10'))
print(x.groupby(['tag', 'kg']).g.mean().unstack().round(2))
print(x[(x.graph == 'hepph') & (x.model == 'LT')].groupby(['tag', 'kg']).g.mean().unstack().round(1))
h = pd.concat([pd.read_csv(f'{R}/alloc_v2/alloc.csv.hcert', header=None, names=HC),
               pd.read_csv(f'{R}/quality_v2/quality_small.csv.hcert', header=None, names=HC).assign(tag='deployed')])
h = h[(h.method == 'RiskIndex-HV') & h.graph.isin(['email', 'wiki', 'facebook', 'hepph']) & h.seed.isin([0, 1, 2]) & (h.alpha <= 0.2)]
print('mean rho:', h.groupby('tag').rho.mean().round(2).to_dict())
# the GRR mean chain at each split (it still pays the price of grouping), for reference
y = d[(d.alpha <= 0.2) & (d.method == 'Mean-GRR') & (d.k >= 10)]
print('GRR mean chain vs RIS at k>=10:', y.groupby('tag').g.mean().round(2).to_dict())
