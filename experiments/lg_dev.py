#!/usr/bin/env python3
"""LargeGraph development summary: gain (%) in realized CVaR over RIS (ind. RR) for RiskIndex-HV and RiskIndex-LG,
alpha <= 0.2, k <= 5 and k >= 10, per graph/model and pooled over the small and the large graphs.
Usage: lg_dev.py <quality csv>..."""
import sys, pandas as pd
c = ['tag','graph','model','M','b','H','Kmax','seed','method','k','alpha','sec','mean','cvar','uid','stats']
d = pd.concat([pd.read_csv(f, header=None, names=c) for f in sys.argv[1:]])
d = d[d.alpha <= 0.2]
key = ['graph','model','seed','k','alpha']
ris = d[d.method == 'RIS-IndRR'].set_index(key).cvar
out = []
for m in ['RiskIndex-HV','RiskIndex-LG']:
    x = d[d.method == m].set_index(key).cvar
    g = (x / ris.reindex(x.index) - 1) * 100
    out.append(g.rename(m))
r = pd.concat(out, axis=1).reset_index()
r['bud'] = r.k.map(lambda k: 'k<=5' if k <= 5 else 'k>=10')
order = ['email','wiki','facebook','hepph','epinions','slashdot','dblp']
t = r.groupby(['graph','model','bud'])[['RiskIndex-HV','RiskIndex-LG']].mean().unstack('bud')
t = t.reindex(sorted(t.index, key=lambda i: (order.index(i[0]) if i[0] in order else 99, i[1])))
pd.set_option('display.width', 200); print(t.round(2).to_string())
for nm, gs in (('small', order[:4]), ('large', order[4:]), ('other', [g for g in r.graph.unique() if g not in order])):
    x = r[r.graph.isin(gs)]
    if len(x): print(nm, x.groupby('bud')[['RiskIndex-HV','RiskIndex-LG']].mean().round(2).to_dict(), 'n=', len(x))
