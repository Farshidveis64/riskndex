#!/usr/bin/env python3
"""Tier A (RiskIndex = RiskIndex-AD: the GRR index plus the gated forward layer, kf = k; RiskIndex-LG = layer without the need gate; RiskIndex-HV = w/o forward layer): realized CVaR_alpha on the evaluation worlds RB as % of the exhaustive optimum.
The exhaustive optimum and the baselines come from results/oracle/oracle.csv; the index variants from
results/oracle_lg/oracle_lg.csv (identical sketches, selection worlds and evaluation worlds; RiskIndex-HV, Mean-GRR
and RIS-IndRR rows are checked for exact equality against results/oracle_v2/oracle_v2.csv)."""
import os
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
R = os.environ.get('RESDIR', 'results')
cols = ['tag', 'graph', 'fam', 'n', 'm', 'model', 'k', 'alpha', 'seed', 'method', 'cvB', 'meanB', 'cvA', 'optA', 'S']
old = pd.read_csv(f'{R}/oracle/oracle.csv', header=None, names=cols)
v2 = pd.read_csv(f'{R}/oracle_v2/oracle_v2.csv', header=None, names=cols)
lgo = pd.read_csv(f'{R}/oracle_lg/oracle_lg.csv', header=None, names=cols)
new = pd.read_csv(f'{R}/oracle_ad/oracle_ad.csv', header=None, names=cols); new['method'] = new.method.replace({'RiskIndex-LG': 'RiskIndex-AD'})
key = ['tag', 'graph', 'k', 'alpha', 'seed']
for m in ['RiskIndex-HV', 'Mean-GRR', 'RIS-IndRR']:
    j = new[new.method == m].merge(v2[v2.method == m], on=key + ['method'], suffixes=('', '_o'))
    bad = ((j.cvB - j.cvB_o).abs() > 1e-9).sum(); print(f'reproduction {m}: {len(j)} rows, {bad} mismatches'); assert bad == 0
d = pd.concat([old, v2[v2.method.isin(['RiskIndex-H', 'RiskIndex-HV-raw'])], lgo[lgo.method == 'RiskIndex-LG'], new[new.method.isin(['RiskIndex-HV', 'RiskIndex-AD'])]], ignore_index=True)
cases = new[new.method == 'RiskIndex-AD'][key].drop_duplicates(); d = d.merge(cases, on=key)
orc = d[d.method == 'Oracle'][key + ['cvB', 'S']].rename(columns={'cvB': 'opt', 'S': 'Sopt'})
d = d.merge(orc, on=key); d['ratio'] = d.cvB / d.opt * 100; d['regret'] = d.opt - d.cvB
print('cases per alpha', len(cases) // 5)
AL = [0.05, 0.1, 0.2, 0.5, 1.0]
tab = d.groupby(['method', 'alpha']).ratio.mean().unstack().reindex(columns=AL)
print(tab.round(2).to_string())
tab.to_csv(f'{R}/oracle_ad/oracle_by_alpha.csv')
F = 'RiskIndex-AD'
x = d[d.method == F]
print('final: mean regret', round(x.regret.mean(), 3), 'max', x.regret.max(), ' exact optimum', round((x.S.apply(lambda s: tuple(sorted(map(int, s.split())))) == x.Sopt.apply(lambda s: tuple(sorted(map(int, s.split()))))).mean(), 3))
def cmp(a, b, sub=d):
    xa = sub[sub.method == a].set_index(key).cvB; xb = sub[sub.method == b].set_index(key).cvB
    j = pd.concat([xa, xb], axis=1, keys=['a', 'b']).dropna(); df = j.a - j.b
    return f'W/T/L={(df>1e-9).sum()}/{(df.abs()<=1e-9).sum()}/{(df<-1e-9).sum()} p={wilcoxon(j.a, j.b).pvalue if (df.abs()>1e-9).any() else 1:.2g} mean diff {df.mean():+.4f}'
for b in ['DirectCVaR-GRR', 'SG', 'RIS-IndRR', 'Mean-GRR', 'Oracle-Mean', 'RiskIndex-HV', 'RiskIndex-LG', 'RiskIndex', 'RiskIndex-H', 'RiskIndex-HV-raw', 'DirectCVaR-IndRR']:
    print(f'  {F} vs {b:18s}', cmp(F, b))
a05 = d[d.alpha == 0.05]
om = a05[a05.method == 'Oracle-Mean'][key + ['S']].rename(columns={'S': 'Smean'})
z = a05.merge(om, on=key)
dif = z[z.Smean.apply(lambda s: tuple(sorted(map(int, s.split())))) != z.Sopt.apply(lambda s: tuple(sorted(map(int, s.split()))))]
print('a=.05 instances where mean-optimal != CVaR-optimal:', round(dif[dif.method == F][key].drop_duplicates().shape[0] / z[z.method == F].shape[0], 3))
print('   on them: final', round(dif[dif.method == F].ratio.mean(), 2), ' mean-optimal', round(dif[dif.method == 'Oracle-Mean'].ratio.mean(), 2),
      ' max regret final', dif[dif.method == F].regret.max(), ' mean-optimal', dif[dif.method == 'Oracle-Mean'].regret.max())

# the index without the forward layer, for the text
for m in ['RiskIndex-HV', 'RiskIndex-LG']:
    y = d[d.method == m]
    print(m, ': mean regret', round(y.regret.mean(), 3), ' exact optimum', round((y.S.apply(lambda s: tuple(sorted(map(int, s.split())))) == y.Sopt.apply(lambda s: tuple(sorted(map(int, s.split()))))).mean(), 3),
          ' on mean!=CVaR instances', round(dif[dif.method == m].ratio.mean(), 2))
print('by alpha, final vs w/o need gate / w/o forward layer:', tab.loc[[F, 'RiskIndex-LG', 'RiskIndex-HV']].round(2).to_dict('index'))
