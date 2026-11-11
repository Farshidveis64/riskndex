#!/usr/bin/env python3
"""rho-guard of the final index on under-resolved dev sketches: rho_min = 0 vs 0.5 (dev seeds 100-102)."""
import os, numpy as np, pandas as pd
R = os.environ.get('RESDIR', 'results')
COLS = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,sec,mean,cvar,uid,stats'.split(',')
HC = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,h,cert,rho,mc_cert,se,se_w,raw,guarded'.split(',')
rows = []
for b in [128, 512]:
    for rm in ['0', '0.5']:
        p = f'{R}/dev/guard_b{b}_r{rm}.csv'
        if not os.path.exists(p): continue
        d = pd.read_csv(p, header=None, names=COLS); key = ['graph', 'model', 'seed', 'k', 'alpha']
        ref = d[d.method == 'RIS-IndRR'][key + ['cvar']].rename(columns={'cvar': 'ris'})
        x = d[d.method == 'RiskIndex-HV'].merge(ref, on=key); x = x[x.alpha <= 0.2]
        x['g'] = (x.cvar / x.ris - 1) * 100
        h = pd.read_csv(p + '.hcert', header=None, names=HC); h = h[(h.method == 'RiskIndex-HV') & (h.alpha <= 0.2)]
        rows.append(dict(b=b, rho_min=rm, k_le5=x[x.k <= 5].g.mean(), k_ge10=x[x.k >= 10].g.mean(), rho=h.rho.mean(), guarded=h.guarded.mean()))
print(pd.DataFrame(rows).round(3).to_string(index=False))
