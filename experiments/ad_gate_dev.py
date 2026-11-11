#!/usr/bin/env python3
"""Need-gate development analysis: choose the tail-resolution threshold theta from the LargeGraph development
runs (seeds 100-102; results/lg_dev), which contain both RiskIndex-HV and RiskIndex-LG answers for every k <= 5.
R_tail(k) = b * C_hat_0.05(RIS k-prefix on the validation sketch) / n  (mc_cert of the RiskIndex-HV hcert row at
alpha = 0.05), i.e. the number of resolution units n/b in the mean of the worst 5% scenarios of the RIS answer.
Gated answer for threshold theta: LG where R_tail(k) < theta (and the cost cap allowed the layer), HV otherwise.
Note: an approximation of the gated index (the layer built with kf = 5 and served for k < theta); the exact gated
index rebuilds the layer with the smaller kf.  Usage: ad_gate_dev.py"""
import pandas as pd, numpy as np
QC = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,sec,mean,cvar,uid,stats'.split(',')
HC = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,h,cert,rho,mc_cert,se,se_w,raw,guarded'.split(',')
N = dict(email=1005, wiki=7115, facebook=4039, hepph=12008, epinions=75879, slashdot=82168, dblp=317080)
q = pd.concat([pd.read_csv(f'results/lg_dev/lgdev_{p}.csv', header=None, names=QC) for p in ['small', 'large']])
h = pd.concat([pd.read_csv(f'results/lg_dev/lgdev_{p}.csv.hcert', header=None, names=HC) for p in ['small', 'large']])
key = ['graph', 'model', 'seed', 'k']
rt = h[(h.method == 'RiskIndex-HV') & (h.alpha == 0.05) & (h.k <= 5)][key + ['mc_cert', 'b']].copy()
rt['R'] = rt.b * rt.mc_cert / rt.graph.map(N)
q = q[(q.k <= 5) & (q.alpha <= 0.2)]
ris = q[q.method == 'RIS-IndRR'].set_index(key + ['alpha']).cvar
d = q[q.method.isin(['RiskIndex-HV', 'RiskIndex-LG'])].pivot_table(index=key + ['alpha'], columns='method', values='cvar').reset_index()
d = d.merge(rt[key + ['R']], on=key)
d['ris'] = ris.reindex(pd.MultiIndex.from_frame(d[key + ['alpha']])).values
small = d.graph.isin(['email', 'wiki', 'facebook', 'hepph'])
print('R_tail per cell (dev, mean over seeds):'); print(d.groupby(['graph', 'model', 'k']).R.mean().unstack('k').round(1).to_string())
rows = []
for th in [0, 2, 5, 10, 20, 50, 1e9]:
    g = np.where(d.R < th, d['RiskIndex-LG'], d['RiskIndex-HV']) / d.ris - 1
    rows.append(dict(theta=th, small=100 * g[small].mean(), large=100 * g[~small].mean(), frac_forward=(d.R < th).mean()))
r = pd.DataFrame(rows); print(r.round(2).to_string(index=False))
