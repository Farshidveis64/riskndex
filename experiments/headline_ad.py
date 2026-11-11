#!/usr/bin/env python3
"""Every Tier B number quoted in Section V (RiskIndex = RiskIndex-AD; RiskIndex-LG = forward layer without the need gate; RiskIndex-HV = without the forward layer), from results/quality_ad/merged.csv.
Reference: classical RIS (RIS-IndRR, independent RR sketch of the same size Mb)."""
import os
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
R = os.environ.get('RESDIR', 'results')
d = pd.read_csv(f'{R}/quality_ad/merged.csv', low_memory=False); d = d[d.tag == 'tierB']
key = ['graph', 'model', 'seed', 'k', 'alpha']
FINAL = 'RiskIndex-AD'
S4 = d.graph.isin(['email', 'wiki', 'facebook', 'hepph'])


def cmp(sub, a, b):
    x = sub[sub.method == a].set_index(key).cvar; y = sub[sub.method == b].set_index(key).cvar
    j = pd.concat([x, y], axis=1, keys=['x', 'y']).dropna(); r = (j.x - j.y) / j.y; df = j.x - j.y
    p = wilcoxon(j.x, j.y).pvalue if (df != 0).any() else 1
    return f'n={len(j)} mean={r.mean()*100:+.2f}% median={r.median()*100:+.2f}% W/T/L={(df>0).sum()}/{(df==0).sum()}/{(df<0).sum()} p={p:.2g}'


BASES = ['RIS-IndRR', 'RiskIndex-HV', 'RiskIndex-LG', 'Mean-GRR', 'DirectCVaR-GRR', 'DirectCVaR-GRR-VM', 'DirectCVaR-IndRR', 'SG', 'SG-cached', 'RiskIndex',
         'RiskIndex-H', 'RiskIndex-V', 'RiskIndex-HV-raw', 'RiskIndex-HV-z1']
for nm, s in [('k<=5 all graphs a=.05', d[(d.k <= 5) & (d.alpha == 0.05)]),
              ('k<=5 four small a=.05', d[(d.k <= 5) & (d.alpha == 0.05) & S4]),
              ('k<=5 three large a=.05', d[(d.k <= 5) & (d.alpha == 0.05) & ~S4]),
              ('k>=10 a=.05', d[(d.k >= 10) & (d.alpha == 0.05)]),
              ('k<=5 a<=.2', d[(d.k <= 5) & (d.alpha <= 0.2)]), ('k>=10 a<=.2', d[(d.k >= 10) & (d.alpha <= 0.2)]),
              ('all k, alpha=1', d[d.alpha == 1.0])]:
    print('==', nm)
    for b in BASES:
        if (s.method == b).sum() == 0: continue
        print('  RI vs', f'{b:18s}', cmp(s, FINAL, b))
    print('  pooled gain over RIS:', {m: f'{(s[s.method==m].rel.mean()-1)*100:+.2f}%' for m in [FINAL] + BASES if (s.method == m).any()})
a05 = d[d.alpha == 0.05]
p = a05.pivot_table(index=['graph', 'model', 'k'], columns='method', values='cvar', aggfunc='mean')
for b in ['RIS-IndRR', 'Mean-GRR', 'SG', 'DirectCVaR-GRR']:
    r = p[FINAL] / p[b]; print(f'max ratio RI/{b}', r.idxmax(), round(r.max(), 2), ' min', r.idxmin(), round(r.min(), 3))
print('hepph WC k=1 a=.05 realized CVaR:', p.loc[('hepph', 'WC', 1)][[FINAL, 'RIS-IndRR', 'Mean-GRR', 'SG', 'DirectCVaR-GRR', 'RiskIndex']].round(1).to_dict())
for gm in [('dblp', 'IC'), ('hepph', 'IC'), ('wiki', 'IC'), ('facebook', 'IC')]:
    x = a05[(a05.graph == gm[0]) & (a05.model == gm[1])]
    print('cell', gm, {kb: {m: round((x[(x.k >= lo) & (x.k <= hi) & (x.method == m)].rel.mean() - 1) * 100, 2)
                            for m in [FINAL, 'RiskIndex', 'DirectCVaR-GRR', 'SG', 'Mean-GRR']} for kb, (lo, hi) in {'k<=5': (1, 5), 'k>=10': (10, 100)}.items()})
L = d[(d.k >= 10)]
print('median RI-RIS rel diff k>=10 by alpha', L[L.method == FINAL].groupby('alpha').rel.median().sub(1).mul(100).round(3).to_dict())
print('mean RI-RIS rel diff k>=10 by alpha', L[L.method == FINAL].groupby('alpha').rel.mean().sub(1).mul(100).round(3).to_dict())
print('min cell gain (graph,model,kgroup) a=.05:')
cg = a05[a05.method == FINAL].assign(kb=np.where(a05[a05.method == FINAL].k <= 5, 'k<=5', 'k>=10')).groupby(['graph', 'model', 'kb']).rel.mean().sub(1).mul(100)
print('  ', cg.round(2).sort_values().head(6).to_dict())
mg = a05[(a05.k >= 10) & (a05.method == 'Mean-GRR')]
print('Mean-GRR vs RIS k>=10 a=.05', f'{(mg.rel.mean()-1)*100:+.2f}%', 'worst gm', mg.groupby(['graph', 'model']).rel.mean().sort_values().head(3).sub(1).mul(100).round(2).to_dict())
# identical to the RIS answer
x = d[d.alpha <= 0.2]
r = x[x.method == FINAL].set_index(key); m = x[x.method == 'RIS-IndRR'].set_index(key).reindex(r.index)
same = np.isclose(r.cvar, m.cvar) & np.isclose(r['mean'], m['mean'])
print('same realized result as RIS (a<=.2):', round(same.mean(), 3), ' by kgroup:', {kb: round(same[(r.reset_index().k <= 5).values == (kb == 'k<=5')].mean(), 3) for kb in ['k<=5', 'k>=10']})
# tail ratio of the RIS answer
ris = d[d.method == 'RIS-IndRR'].pivot_table(index=['graph', 'model', 'k'], columns='alpha', values='cvar', aggfunc='mean')
tr = (ris[0.05] / ris[1.0]).reset_index(name='tr')
print('tail ratio RIS k=1', tr[tr.k == 1].tr.agg(['min', 'max']).round(2).to_dict(), 'k=100', tr[tr.k == 100].tr.agg(['min', 'max']).round(2).to_dict())
# raw vs VM selection agreement (final index)
k2 = key
rr = d[d.method == 'RiskIndex-HV-raw'].set_index(k2).uid; rv = d[d.method == 'RiskIndex-HV'].set_index(k2).uid.reindex(rr.index)
print('raw==VM same set on the GRR chains (RiskIndex-HV-raw vs RiskIndex-HV)', round((rr == rv).mean(), 3))
print('  for reference, RiskIndex-HV-raw vs the final index (forward-served budgets differ by construction)', round((rr == d[d.method == FINAL].set_index(k2).uid.reindex(rr.index)).mean(), 3))
g = d[d.alpha <= 0.2].pivot_table(index=['graph'], columns='method', values='rel', aggfunc='mean')
print('raw - VM on the GRR chains (points, a<=.2) by graph', ((g['RiskIndex-HV-raw'] - g['RiskIndex-HV']) * 100).round(2).to_dict())
# ablation table (six graphs other than DBLP, alpha = .05)
six = a05[a05.graph != 'dblp']
abl = {}
for m in ['RIS-IndRR', 'Mean-GRR', 'SG-cached', 'RiskIndex-HV-raw', 'DirectCVaR-GRR', 'RiskIndex', 'RiskIndex-H', 'RiskIndex-V', 'RiskIndex-HV-z1', FINAL]:
    abl[m] = tuple(round((six[(six.k >= lo) & (six.k <= hi) & (six.method == m)].rel.mean() - 1) * 100, 2) for lo, hi in [(1, 5), (10, 100)])
print('ablation (six graphs, a=.05; k<=5, k>=10 gain over RIS):', abl)
print('per-query seconds by method (mean over Tier B):', d.groupby('method').sec.mean().round(3).to_dict())
# certificates of the final index
c = pd.concat([pd.read_csv(f'{R}/quality_ad/certificates.csv'), pd.read_csv(f'{R}/quality_v2/certificates.csv')], ignore_index=True)
for meth in [FINAL, 'RiskIndex', 'RiskIndex-H']:
    cc = c[(c.method == meth) & (c.alpha <= 0.2)]
    if cc.empty: continue
    ex = cc[~((cc.graph == 'dblp') & (cc.model == 'IC'))]
    print(f'cert {meth}: median calib by alpha', cc.groupby('alpha').calib.median().mul(100).round(2).to_dict(),
          ' within 10%:', round((ex.calib.abs() <= 0.10).mean(), 3), ' within 12%:', round((ex.calib.abs() <= 0.12).mean(), 3),
          ' raw median', round(cc.calib_raw.median() * 100, 2), ' by kgrp', cc.groupby('kgrp').calib.median().mul(100).round(2).to_dict(),
          ' alpha=1 median', round(c[(c.method == meth) & (c.alpha == 1.0)].calib.median() * 100, 2))
    print(f'   coverage +-2se {cc.in2se.mean():.3f}  +-2se_w {cc.in2sew.mean():.3f}  by kgrp', cc.groupby('kgrp').in2sew.mean().round(3).to_dict(),
          ' median se/cert', round((cc.se / cc.cert).median() * 100, 2), '% se_w/cert', round((cc.se_w / cc.cert).median() * 100, 2), '%',
          ' guarded', round(cc.guarded.mean(), 4))
    dd = c[(c.method == meth) & (c.graph == 'dblp') & (c.model == 'IC')]
    if not dd.empty: print('   DBLP IC: cert median', round(dd.calib.median() * 100, 2), 'rho k<=5/k>=10', dd.groupby('kgrp').rho.mean().round(2).to_dict())
print('rho by graph/model (final):', c[c.method == FINAL].groupby(['graph', 'model']).rho.mean().round(2).to_dict())
