#!/usr/bin/env python3
"""LaTeX bodies of Tables III (Tier A), IV (ablation) and V (Tier C) (Adaptive index).
RiskIndex = RiskIndex-AD (GRR index plus the gated forward layer); 'w/o need gate' = RiskIndex-LG; 'w/o forward layer' = RiskIndex-HV.
Writes results/tables_ad/*.tex"""
import os, io
import numpy as np, pandas as pd
R = os.environ.get('RESDIR', 'results'); OUT = f'{R}/tables_ad'; os.makedirs(OUT, exist_ok=True)
SH = '\\cellcolor{gray!22}'
AL = [0.05, 0.1, 0.2, 0.5, 1.0]
# ------------------------------------------------------------------ Tier A
cols = ['tag', 'graph', 'fam', 'n', 'm', 'model', 'k', 'alpha', 'seed', 'method', 'cvB', 'meanB', 'cvA', 'optA', 'S']
old = pd.read_csv(f'{R}/oracle/oracle.csv', header=None, names=cols)
lgo = pd.read_csv(f'{R}/oracle_lg/oracle_lg.csv', header=None, names=cols)
new = pd.read_csv(f'{R}/oracle_ad/oracle_ad.csv', header=None, names=cols); new['method'] = new.method.replace({'RiskIndex-LG': 'RiskIndex-AD'})
v2 = pd.read_csv(f'{R}/oracle_v2/oracle_v2.csv', header=None, names=cols)
key = ['tag', 'graph', 'k', 'alpha', 'seed']
for m in ['RiskIndex-HV', 'RIS-IndRR', 'Mean-GRR']:  # exact reproduction of the Improvement rows
    j = new[new.method == m].merge(v2[v2.method == m], on=key, suffixes=('', '_o'))
    bad = ((j.cvB - j.cvB_o).abs() > 1e-6).sum(); print(f'Tier A reproduction {m}: {len(j)} rows, {bad} mismatches'); assert bad == 0
cases = new[new.method == 'RiskIndex-AD'][key].drop_duplicates()
d = pd.concat([old, lgo[lgo.method == 'RiskIndex-LG'], new[new.method.isin(['RiskIndex-HV', 'RiskIndex-AD'])]], ignore_index=True).merge(cases, on=key)
orc = d[d.method == 'Oracle'][key + ['cvB']].rename(columns={'cvB': 'opt'})
d = d.merge(orc, on=key); d['ratio'] = d.cvB / d.opt * 100
A = d.groupby(['method', 'alpha']).ratio.mean().unstack().reindex(columns=AL)
rowsA = [('Reference', [('Exhaustive mean-optimal', 'Oracle-Mean')]),
         ('Mean-based', [('RIS (ind.\\ RR)', 'RIS-IndRR'), ('Mean chain (GRR)', 'Mean-GRR')]),
         ('Per-query', [('Direct-CVaR (ind.\\ RR)', 'DirectCVaR-IndRR'), ('Direct-CVaR-VM (GRR)', 'DirectCVaR-GRR-VM'), ('Direct-CVaR (GRR)', 'DirectCVaR-GRR'), ('SG~\\cite{zhou2020} (GRR)', 'SG')]),
         ('Index', [('\\RI{} w/o forward layer', 'RiskIndex-HV'), ('\\RI{} w/o need gate', 'RiskIndex-LG'), ('\\RI{} (Ours)', 'RiskIndex-AD')])]
nonex = [m for _, rs in rowsA[1:] for _, m in rs]
best = {a: A.loc[nonex, a].round(2).max() for a in AL}
lines = []
for ci, (cat, rs) in enumerate(rowsA):
    if ci: lines.append('\\midrule')
    for i, (lab, m) in enumerate(rs):
        cells = []
        for a in AL:
            v = round(A.loc[m, a], 2); cells.append((SH if (m in nonex and v == best[a]) else '') + f'{v:.2f}')
        first = (f'\\multirow{{{len(rs)}}}{{*}}{{{cat}}}' if len(rs) > 1 else cat) if i == 0 else ''
        lines.append(f'{first} & {lab} & ' + ' & '.join(cells) + ' \\\\')
open(f'{OUT}/oracle.tex', 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
print('Tier A cases per alpha:', len(cases) // 5)
for m in ['RiskIndex-AD', 'RiskIndex-LG', 'RiskIndex-HV']:
    x = d[d.method == m]; print(m, 'Tier A mean %% of optimum (all alpha): %.2f, alpha<=0.2: %.2f, min %.1f' % (x.ratio.mean(), x[x.alpha <= 0.2].ratio.mean(), x.ratio.min()))
# ------------------------------------------------------------------ ablation
q = pd.read_csv(f'{R}/quality_ad/merged.csv', low_memory=False); q = q[(q.tag == 'tierB') & (q.alpha == 0.05) & (q.graph != 'dblp')]
def tb(m):
    return tuple(round((q[(q.method == m) & (q.k >= lo) & (q.k <= hi)].rel.mean() - 1) * 100, 1) + 0.0 for lo, hi in [(1, 5), (10, 100)])
d2 = pd.concat([old, pd.read_csv(f'{R}/oracle_v2/oracle_v2.csv', header=None, names=cols).query("method in ['RiskIndex-H', 'RiskIndex-HV-raw', 'RiskIndex-V']"), lgo[lgo.method == 'RiskIndex-LG'], new[new.method.isin(['RiskIndex-HV', 'RiskIndex-AD'])]], ignore_index=True).merge(cases, on=key).merge(orc, on=key)
d2['ratio'] = d2.cvB / d2.opt * 100
A2 = d2.groupby(['method', 'alpha']).ratio.mean().unstack()
ta = lambda m: round(A2.loc[m, 0.05], 1)
lat = {}
s = pd.read_csv(f'{R}/scale_v2/scale_agg.csv'); s = s[(s.model == 'WC') & ~s.graph.isin(['ba30000', 'ba300000'])]
lat['ri'] = '%.1f\\,$\\mu$s' % s.p50us.mean(); lat['nomat'] = '%.0f--%.0f\\,ms' % (s.np50ms.min(), s.np50ms.max())
# exact-layer lookup (median over the Tier C runs that build the layer; results/scale_lg)
_fl = [l for l in open(f'{R}/scale_ad/scale_hyb.csv') if l.startswith('fl,')] if os.path.exists(f'{R}/scale_ad/scale_hyb.csv') else []
_fl = pd.read_csv(io.StringIO(''.join(_fl)), header=None, names=['tag', 'graph', 'model', 'seed', 'fl', 'screen', 'greedy', 'val', 'pool', 'bytes', 'p50', 'p95', 'kf']) if _fl else None
lat['lg'] = ('%.1f--%.1f\\,$\\mu$s' % (s.p50us.mean(), _fl[_fl.kf > 0].groupby(['graph', 'model']).p50.mean().max() * 1e6)) if _fl is not None else lat['ri']
tq = pd.read_csv(f'{R}/quality_ad/merged.csv', low_memory=False); tq = tq[tq.tag == 'tierB']
dsec = tq[tq.method == 'DirectCVaR-GRR'].groupby(['graph']).sec.mean(); sgc = tq[tq.method == 'SG-cached'].sec.mean() * 1e3
Y, N = '\\checkmark', '$\\times$'
# Components G, B, S, M, V, R, O, E; rows 1-9 remove one component from the index without the forward layer (E).
rows = [('w/o grouping', [N, N, N, Y, N, Y, N, N], 'RIS-IndRR', 'RIS-IndRR', '--'),
        ('w/o basis', [Y, N, N, Y, Y, N, N, N], 'Mean-GRR', 'Mean-GRR', lat['ri']),
        ('w/o tail selection', [Y, Y, N, N, N, N, N, N], 'SG', 'SG-cached', '%.0f\\,ms' % sgc),
        ('w/o materialization', [Y, Y, Y, N, Y, Y, Y, N], 'RiskIndex-HV', 'RiskIndex-HV', lat['nomat']),
        ('w/o VM', [Y, Y, Y, Y, N, Y, Y, N], 'RiskIndex-HV-raw', 'RiskIndex-HV-raw', lat['ri']),
        ('w/o index', [Y, N, Y, N, N, N, N, N], 'DirectCVaR-GRR', 'DirectCVaR-GRR', '%.1f--%.0f\\,s' % (dsec.min(), dsec.max())),
        ('w/o RIS chain', [Y, Y, Y, Y, Y, N, Y, N], 'RiskIndex-V', 'RiskIndex-V', lat['ri']),
        ('w/o validation', [Y, Y, Y, Y, Y, Y, N, N], 'RiskIndex-H', 'RiskIndex-H', lat['ri']),
        ('w/o forward layer', [Y, Y, Y, Y, Y, Y, Y, N], 'RiskIndex-HV', 'RiskIndex-HV', lat['ri']),
        ('w/o need gate', [Y, Y, Y, Y, Y, Y, Y, Y], 'RiskIndex-LG', 'RiskIndex-LG', lat['lg']),
        ('\\RI{} (Ours)', [Y, Y, Y, Y, Y, Y, Y, Y], 'RiskIndex-AD', 'RiskIndex-AD', lat['lg'])]
vals = [(ta(a), *tb(b)) for _, _, a, b, _ in rows]
bests = [max(v[i] for v in vals) for i in range(3)]
lines = []
for (lab, comp, a, b, l), v in zip(rows, vals):
    if lab == 'w/o forward layer': lines.append('\\midrule')
    cells = [(SH if v[i] == bests[i] else '') + (f'{v[i]:.1f}' if i == 0 else ('$0.0$' if v[i] == 0 else f'${v[i]:+.1f}$')) for i in range(3)]
    lines.append(f'{lab} & ' + ' & '.join(comp) + ' & ' + ' & '.join(cells) + f' & {l} \\\\')
open(f'{OUT}/ablation.tex', 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
# ------------------------------------------------------------------ Tier C systems table (exact-layer build added)
fl = [l for l in open(f'{R}/scale_ad/scale_hyb.csv') if l.startswith('fl,')] if os.path.exists(f'{R}/scale_ad/scale_hyb.csv') else []
if not fl: print('WARNING: no exact-layer build times yet (results/scale_lg); systems table without them')
fl = pd.read_csv(io.StringIO(''.join(fl) or 'x,x,x,0,0,0,0,0,0,0,0,0,0\n'), header=None, names=['tag', 'graph', 'model', 'seed', 'fl', 'screen', 'greedy', 'val', 'pool', 'bytes', 'p50', 'p95', 'kf'])
fl = fl.drop_duplicates(['graph', 'model', 'seed'], keep='last')
# break-even per seed as in analyze_scale_v2.py (index, or sketch + index, over direct - lookup), now with the exact
# layer's build time of the same seed added, then averaged over seeds; unchanged where the layer is not built
vc = ['tag', 'graph', 'model', 'n', 'm', 'M', 'b', 'H', 'Kmax', 'seed', 'load', 'sketch', 'index', 'L', 'skbytes', 'idxbytes', 'p50', 'p95', 'lookup', 'workload', 'np50', 'np95', 'direct', 'dp50', 'dp95', 'dmax', 'qstar', 'rss0', 'rss1', 'avgrr']
sv = pd.read_csv(f'{R}/scale_v2/scale.csv', header=None, names=vc).drop_duplicates(['graph', 'model', 'seed'], keep='last')
sv = sv.merge(fl[['graph', 'model', 'seed', 'fl']], on=['graph', 'model', 'seed'], how='left'); sv['fl'] = sv.fl.fillna(0.0)
sv['qstar_idx'] = (sv['index'] + sv.fl) / (sv.direct - sv.lookup); sv['qstar_all'] = (sv.sketch + sv['index'] + sv.fl) / (sv.direct - sv.lookup)
qs = sv.groupby(['graph', 'model']).agg(qstar_idx=('qstar_idx', 'mean'), qstar_all=('qstar_all', 'mean')).reset_index()
fl = fl.groupby(['graph', 'model']).agg(fl=('fl', 'mean'), kf=('kf', 'max'), flp50=('p50', 'mean')).reset_index()
s = pd.read_csv(f'{R}/scale_v2/scale_agg.csv').drop(columns=['qstar_idx', 'qstar_all']).merge(fl, on=['graph', 'model'], how='left').merge(qs, on=['graph', 'model'], how='left')
s['fl'] = s.fl.fillna(0.0)
s.to_csv(f'{OUT}/scale_agg_ad.csv', index=False)
order = [('Synthetic', ['ba10000', 'ba100000', 'ba1000000']), ('Real', ['email', 'wiki', 'epinions', 'slashdot', 'dblp', 'youtube', 'livejournal'])]
pretty = dict(ba10000='BA-10K', ba100000='BA-100K', ba1000000='BA-1M', email='email-Eu-core', wiki='Wiki-Vote', epinions='Epinions', slashdot='Slashdot', dblp='DBLP', youtube='YouTube', livejournal='LiveJournal')
lines = []
for ti, (typ, gs) in enumerate(order):
    if ti: lines.append('\\midrule')
    for i, gname in enumerate(gs):
        r = s[(s.graph == gname) & (s.model == 'WC')]
        if r.empty: continue
        r = r.iloc[0]
        first = f'\\multirow{{{len(gs)}}}{{*}}{{{typ}}}' if i == 0 else ''
        fls = f'{r.fl:.1f}' if r.kf > 0 else '--'
        lines.append(f'{first} & {pretty[gname]} & {int(r.b):,} & {r.sketch:.1f} & {r["index"]:.1f} & {fls} & {r.skMB:.0f} & {r.p50us:.2f} & {r.p95us:.2f} & {r.thr/1e6:.1f} & {r.np50ms:.1f} & {r.direct:.2f} & {r.dp95:.1f} & {r.qstar_idx:.1f} & {r.qstar_all:.1f} \\\\')
open(f'{OUT}/systems.tex', 'w').write('\n'.join(lines) + '\n'); print('\n'.join(lines))
print(fl.round(3).to_string())
