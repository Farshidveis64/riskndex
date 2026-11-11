#!/usr/bin/env python3
"""Tier B quality table (RiskIndex = RiskIndex-AD, i.e. the GRR index plus the gated forward layer): gain (%) in realized CVaR_0.05 over classical RIS (independent RR
sketch, same budget) on the same case, averaged over seeds and budgets per (graph, model); a pooled row; a row
counting the rows in which each baseline is better than / equal to / worse than RiskIndex (rounded values).
Shading marks the best rounded value per row and budget group.  RiskIndex = the final index (RiskIndex-HV).
Usage: make_quality_table_ad.py <results/quality_ad/merged.csv> <out.tex>"""
import sys, pandas as pd
src, out = sys.argv[1], sys.argv[2]
q = pd.read_csv(src, low_memory=False); q = q[q.tag == 'tierB']
pretty = dict(email='email-Eu-core', wiki='Wiki-Vote', facebook='Facebook', hepph='ca-HepPh', epinions='Epinions', slashdot='Slashdot', dblp='DBLP')
base_methods = ['DirectCVaR-GRR', 'SG', 'DirectCVaR-IndRR', 'Mean-GRR']; FINAL = 'RiskIndex-AD'; allm = base_methods + [FINAL]
SH = '\\cellcolor{gray!22}'
fmt = lambda p: '$0.0$' if p == 0 else (f'$+{p:.1f}$' if p > 0 else f'$-{-p:.1f}$')
def gains(da):
    return {m: (None if (da.method == m).sum() == 0 else round((da[da.method == m].rel.mean() - 1) * 100, 1) + 0.0) for m in allm}
def cells(vals):
    best = max(v for v in vals.values() if v is not None)
    return [('--' if vals[m] is None else (SH if vals[m] == best else '') + fmt(vals[m])) for m in allm]
counts = {(g, m): [0, 0, 0] for g in range(2) for m in base_methods}
lines = []; first = True
for g in ['email', 'wiki', 'facebook', 'hepph', 'epinions', 'slashdot', 'dblp']:
    models = [mo for mo in ['WC', 'LT', 'IC'] if not q[(q.graph == g) & (q.model == mo)].empty]
    if not models: continue
    if not first: lines.append('\\addlinespace[1.5pt]')
    first = False
    for i, mo in enumerate(models):
        d = q[(q.graph == g) & (q.model == mo)]; base = d[d.method == 'RIS-IndRR']
        bm = base[(base.k == 10) & (base.alpha == 1.0)].cvar.mean(); b05 = base[(base.k == 10) & (base.alpha == 0.05)].cvar.mean()
        row = [(f'\\multirow{{{len(models)}}}{{*}}{{{pretty[g]}}}' if i == 0 else ''), mo, f'{bm:,.0f}', f'{b05:,.0f}']
        for gi, (lo, hi) in enumerate([(1, 5), (10, 100)]):
            vals = gains(d[(d.alpha == 0.05) & (d.k >= lo) & (d.k <= hi)])
            for m in base_methods:
                if vals[m] is None: continue
                r = vals[FINAL]; counts[(gi, m)][0 if vals[m] > r else (1 if vals[m] == r else 2)] += 1
            row += cells(vals)
        lines.append(' & '.join(row) + ' \\\\')
pool = ['\\multicolumn{4}{@{}l}{Pooled over all cases}']
for lo, hi in [(1, 5), (10, 100)]:
    pool += cells(gains(q[(q.alpha == 0.05) & (q.k >= lo) & (q.k <= hi)]))
foot = ['\\multicolumn{4}{@{}l}{Rows better/equal/worse}']
for gi in range(2):
    foot += ['/'.join(map(str, counts[(gi, m)])) for m in base_methods] + ['']
lines += ['\\midrule', ' & '.join(pool) + ' \\\\', ' & '.join(foot) + ' \\\\']
open(out, 'w').write('\n'.join(lines) + '\n'); print('wrote', out); print('\n'.join(lines))
