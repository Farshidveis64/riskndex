#!/usr/bin/env python3
"""Need-gate threshold sensitivity (development data only; theta was frozen at 5 before the reported runs).
Rows: theta in {0 (no layer), 5, 10, 20, 50, inf (no need gate)}.  Columns: Tier A development quality (% of the
exhaustive optimum at alpha = 0.05, seeds 100-102, 720 cases; results/ad_dev/oracle_dev_theta*.csv, exact gated
layer), Tier B development gain over RIS (k <= 5, alpha <= 0.2, pooled over the seven graphs, seeds 100-102;
results/lg_dev, gated answer = ungated layer where R_tail(k) < theta, GRR index otherwise), and the fraction of
Tier B development (graph, model, seed, k <= 5) cells whose budget the layer serves.  Writes results/tables_ad/theta.tex."""
import pandas as pd, numpy as np, os
QC = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,sec,mean,cvar,uid,stats'.split(',')
HC = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,h,cert,rho,mc_cert,se,se_w,raw,guarded'.split(',')
N = dict(email=1005, wiki=7115, facebook=4039, hepph=12008, epinions=75879, slashdot=82168, dblp=317080)
q = pd.concat([pd.read_csv(f'results/lg_dev/lgdev_{p}.csv', header=None, names=QC) for p in ['small', 'large']])
h = pd.concat([pd.read_csv(f'results/lg_dev/lgdev_{p}.csv.hcert', header=None, names=HC) for p in ['small', 'large']])
key = ['graph', 'model', 'seed', 'k']
rt = h[(h.method == 'RiskIndex-HV') & (h.alpha == 0.05) & (h.k <= 5)][key + ['mc_cert', 'b']].copy(); rt['R'] = rt.b * rt.mc_cert / rt.graph.map(N)
q = q[(q.k <= 5) & (q.alpha <= 0.2)]
ris = q[q.method == 'RIS-IndRR'].set_index(key + ['alpha']).cvar
d = q[q.method.isin(['RiskIndex-HV', 'RiskIndex-LG'])].pivot_table(index=key + ['alpha'], columns='method', values='cvar').reset_index().merge(rt[key + ['R']], on=key)
d['ris'] = ris.reindex(pd.MultiIndex.from_frame(d[key + ['alpha']])).values
cols = ['tag', 'graph', 'fam', 'n', 'm', 'model', 'k', 'alpha', 'seed', 'method', 'cvB', 'meanB', 'cvA', 'optA', 'S']
orc = pd.read_csv('results/oracle/oracle.csv', header=None, names=cols); orc = orc[orc.method == 'Oracle'].groupby(['tag', 'graph', 'k', 'alpha']).cvB.first().rename('opt')
def tierA(th):
    f = f'results/ad_dev/oracle_dev_theta{th}.csv'
    if not os.path.exists(f): return None
    x = pd.read_csv(f, header=None, names=cols).merge(orc, on=['tag', 'graph', 'k', 'alpha']); x['r'] = x.cvB / x.opt * 100
    t = x[x.alpha == 0.05].groupby('method').r.mean(); return t.get('RiskIndex-AD', np.nan), t.get('RiskIndex-HV', np.nan)
rows = []
for th, lab in [(0, '$0$ (no layer)'), (5, '$5$ (adopted)'), (10, '$10$'), (20, '$20$'), (50, '$50$'), (1e9, '$\\infty$ (no need gate)')]:
    g = np.where(d.R < th, d['RiskIndex-LG'], d['RiskIndex-HV']) / d.ris - 1
    ta = tierA(th if th < 1e9 else '1e300') if th > 0 else None   # theta = inf: results/ad_dev/oracle_dev_theta1e300.csv (fltheta=1e300)
    a05 = ta[0] if ta else tierA(5)[1]
    rows.append(dict(theta=lab, tierA=a05, gain=100 * g.mean(), frac=100 * (d.R < th).mean()))
r = pd.DataFrame(rows); print(r.round(2).to_string(index=False))
os.makedirs('results/tables_ad', exist_ok=True)
with open('results/tables_ad/theta.tex', 'w') as f:
    for _, x in r.iterrows():
        ta = '--' if np.isnan(x.tierA) else f'{x.tierA:.2f}'
        f.write(f'{x.theta} & {ta} & {x.gain:+.1f} & {x.frac:.0f} \\\\\n')
print(open('results/tables_ad/theta.tex').read())
