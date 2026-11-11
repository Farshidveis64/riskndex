#!/usr/bin/env python3
"""Matched-grid control for the ablation row 'w/o tail selection' (Table IV).
SG-geogrid = the SG selection rule (largest surrogate value) on chains built on RiskIndex's own geometric threshold
grid, run by scripts/run_quality_sgg.sh on the sketches and evaluation worlds of the Tier B runs (six graphs, no DBLP).
Checks that RIS-IndRR and RiskIndex-HV reproduce results/quality_ad/merged.csv exactly, then reports the gain in realized
CVaR_0.05 over RIS with the formula of the ablation table (mean of per-case ratios minus one).
Writes results/quality_sgg/analysis.txt"""
import os, glob
import numpy as np, pandas as pd
from scipy.stats import wilcoxon
R = os.environ.get('RESDIR', 'results')
m = pd.read_csv(f'{R}/quality_ad/merged.csv', low_memory=False); m = m[(m.tag == 'tierB') & (m.graph != 'dblp')]
cols = list(m.columns[:16])
s = pd.concat([pd.read_csv(f, header=None, names=cols) for f in sorted(glob.glob(f'{R}/quality_sgg/quality_?.csv'))], ignore_index=True)
key = ['graph', 'model', 'seed', 'k', 'alpha']
out = []
for meth in ['RIS-IndRR', 'RiskIndex-HV', 'Mean-GRR']:
    j = s[s.method == meth].merge(m[m.method == meth], on=key, suffixes=('', '_o'))
    bad = int(((j.cvar - j.cvar_o).abs() > 1e-9).sum()); out.append(f'reproduction {meth}: {len(j)} rows, {bad} mismatches'); assert bad == 0 and len(j) == (m.method == meth).sum()
ris = m[m.method == 'RIS-IndRR'].set_index(key).cvar
def col(df, meth): return df[df.method == meth].set_index(key).cvar.reindex(ris.index)
series = {'SG-geogrid (SG rule, index grid)': col(s, 'SG-geogrid'), 'SG-cached (SG rule, uniform grid; Table IV row)': col(m, 'SG-cached'),
          'SG (per query, uniform grid)': col(m, 'SG'), 'RiskIndex-HV (tail selection, index grid)': col(m, 'RiskIndex-HV')}
idx = ris.index.to_frame(index=False)
for name, x in series.items():
    assert x.notna().all(), name
    for lab, sel in [('k<=5', (idx.k <= 5) & (idx.alpha == 0.05)), ('k>=10', (idx.k >= 10) & (idx.alpha == 0.05))]:
        r = (x.values / ris.values - 1)[sel.values]
        out.append(f'{name:50s} {lab:6s} gain over RIS {100 * r.mean():+.2f}%  (n={len(r)})')
hv = series['RiskIndex-HV (tail selection, index grid)']; g = series['SG-geogrid (SG rule, index grid)']
for lab, sel in [('k<=5', (idx.k <= 5) & (idx.alpha == 0.05)), ('k>=10', (idx.k >= 10) & (idx.alpha == 0.05))]:
    a, b = hv.values[sel.values], g.values[sel.values]; d = a / b - 1
    out.append(f'RiskIndex-HV vs SG-geogrid {lab}: {100 * d.mean():+.2f}%  W/T/L={(d > 1e-12).sum()}/{(np.abs(d) <= 1e-12).sum()}/{(d < -1e-12).sum()}  p={wilcoxon(a, b).pvalue:.2g}')
    nh = (sel & (idx.graph != 'hepph')).values; a, b = hv.values[nh], g.values[nh]; d = a / b - 1
    out.append(f'   without ca-HepPh: {100 * d.mean():+.2f}%  W/T/L={(d > 1e-12).sum()}/{(np.abs(d) <= 1e-12).sum()}/{(d < -1e-12).sum()}  p={wilcoxon(a, b).pvalue:.2g}')
open(f'{R}/quality_sgg/analysis.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
