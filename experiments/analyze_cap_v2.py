#!/usr/bin/env python3
"""DBLP IC-critical at the memory caps b = 2^14, 2^16 (Tier B) and 2^17: realized CVaR_0.05 of the final index
relative to RIS, absolute CVaR, resolution ratio and certificate calibration."""
import os, numpy as np, pandas as pd
R = os.environ.get('RESDIR', 'results')
COLS = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,sec,mean,cvar,uid,stats'.split(',')
HC = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,h,cert,rho,mc_cert,se,se_w,raw,guarded'.split(',')
cap = pd.read_csv(f'{R}/alloc_v2/dblp_cap.csv', header=None, names=COLS)
tb = pd.read_csv(f'{R}/quality_v2/quality_dblp.csv', header=None, names=COLS); tb = tb[tb.model == 'IC']
d = pd.concat([cap, tb.assign(tag='cap65536')])
hc = pd.concat([pd.read_csv(f'{R}/alloc_v2/dblp_cap.csv.hcert', header=None, names=HC),
                pd.read_csv(f'{R}/quality_v2/quality_dblp.csv.hcert', header=None, names=HC).query("model == 'IC'").assign(tag='cap65536')])
key = ['tag', 'seed', 'k', 'alpha']
ref = d[d.method == 'RIS-IndRR'][key + ['cvar']].rename(columns={'cvar': 'ris'})
d = d.merge(ref, on=key); d['g'] = (d.cvar / d.ris - 1) * 100
a = d[d.alpha == 0.05].assign(kg=lambda x: np.where(x.k <= 5, 'k<=5', 'k>=10'))
print(a.groupby(['tag', 'method', 'kg']).g.mean().unstack().round(2).to_string())
print('absolute CVaR_0.05 (mean over k, seeds):'); print(a.groupby(['tag', 'method']).cvar.mean().unstack().round(2).to_string())
h = hc[hc.method == 'RiskIndex-HV'].merge(d[d.method == 'RiskIndex-HV'][key + ['cvar']], on=key)
h = h.assign(kg=np.where(h.k <= 5, 'k<=5', 'k>=10'), calib=h.cert / h.cvar - 1)
print(h[h.alpha <= 0.2].groupby(['tag', 'kg']).agg(rho=('rho', 'mean'), calib=('calib', 'median'), guard=('guarded', 'mean')).round(3).to_string())
