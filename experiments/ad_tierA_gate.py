#!/usr/bin/env python3
"""Tier A gate behaviour (analysis only).  Joins the gate diagnostic of scripts/ad_dev/tierA_gate_diag.sh
(results/ad_dev/tierA_gate_diag.log: R_tail(k) of the RIS prefixes and the budgets k_f served by the forward layer, per
instance and seed) with the reported Tier A answers (results/oracle_ad/oracle_ad.csv, reproduced exactly by the
diagnostic run) and the exhaustive optimum (results/oracle/oracle.csv)."""
import pandas as pd, numpy as np, re
cols = ['tag', 'graph', 'fam', 'n', 'm', 'model', 'k', 'alpha', 'seed', 'method', 'cvB', 'meanB', 'cvA', 'optA', 'S']
b = pd.read_csv('results/oracle_ad/oracle_ad.csv', header=None, names=cols)
old = pd.read_csv('results/oracle/oracle.csv', header=None, names=cols)
orc = old[old.method == 'Oracle'].groupby(['tag', 'graph', 'k', 'alpha']).cvB.first().rename('opt')
rows = []
for line in open('results/ad_dev/tierA_gate_diag.log'):
    mm = re.match(r'\[flgate\] (\S+) (\S+) (\S+) seed=(\d+) kf=(\d+) rtail=(\S+)', line); tag, g, mo, seed, kf, rt = mm.groups(); rt = [float(v) for v in rt.split(',')]
    rows.append(dict(tag=tag, graph=g, model=mo, seed=int(seed), k=len(rt), kf=int(kf), rk=rt[-1]))
d = pd.DataFrame(rows); d['served'] = d.kf == d.k
print('instances x seeds:', len(d), '; layer built (k_f >= 1):', round(100 * (d.kf > 0).mean(), 1), '%; queried budget k served (k_f = k):', round(100 * d.served.mean(), 1), '%')
print('served by model:', d.groupby('tag').served.mean().round(3).to_dict())
print('R_tail(k) quantiles:', d.rk.quantile([0, .25, .5, .75, 1]).round(2).to_dict(), '; served: max R_tail(k) =', d[d.served].rk.max().round(2))
x = b[b.method.isin(['RiskIndex-AD', 'RiskIndex-HV'])].merge(orc, on=['tag', 'graph', 'k', 'alpha']); x['r'] = x.cvB / x.opt * 100
x = x.merge(d[['tag', 'graph', 'model', 'seed', 'k', 'served']], on=['tag', 'graph', 'model', 'seed', 'k'])
for al in [0.05, 0.1, 0.2, 0.5, 1.0]:
    y = x[x.alpha == al].pivot_table(index=['tag', 'graph', 'seed', 'k', 'served'], columns='method', values='r').reset_index(); y['diff'] = y['RiskIndex-AD'] - y['RiskIndex-HV']
    s = y[y.served]; u = y[~y.served]
    print(f'alpha={al}: served n={len(s)} AD {s["RiskIndex-AD"].mean():.2f} HV {s["RiskIndex-HV"].mean():.2f} diff {s["diff"].mean():+.3f} (W/T/L {(s["diff"]>1e-9).sum()}/{(abs(s["diff"])<=1e-9).sum()}/{(s["diff"]<-1e-9).sum()}); '
          f'unserved n={len(u)} diff {u["diff"].mean():+.3f}; overall diff {y["diff"].mean():+.3f}')
