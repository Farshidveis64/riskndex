#!/usr/bin/env python3
"""Tier C, index without the forward layer: build cost (with the RR sketch, the RIS chain and the validation sketch charged to
the index), index size, latency, break-even, and realized quality on R fresh worlds."""
import os
import numpy as np, pandas as pd
R = os.environ.get('RESDIR', 'results')
cols = ['tag', 'graph', 'model', 'n', 'm', 'M', 'b', 'H', 'Kmax', 'seed', 'load', 'sketch', 'index', 'L', 'skbytes', 'idxbytes', 'p50', 'p95', 'lookup', 'workload', 'np50', 'np95', 'direct', 'dp50', 'dp95', 'dmax', 'qstar', 'rss0', 'rss1', 'avgrr']
sc = pd.read_csv(f'{R}/scale_v2/scale.csv', header=None, names=cols).drop_duplicates(['graph', 'model', 'seed'], keep='last')
import io
blines = [l for l in open(f'{R}/scale_v2/scale_hyb.csv') if l.startswith('build,')]
hb = pd.read_csv(io.StringIO(''.join(blines)), header=None, names=['tag', 'graph', 'model', 'seed', 'ind', 'ris', 'grr', 'index_total', 'bytes', 'val', 'rss']).drop(columns=['tag'])
hb = hb.drop_duplicates(['graph', 'model', 'seed'], keep='last')
sc = sc.merge(hb[['graph', 'model', 'seed', 'ind', 'ris', 'val']], on=['graph', 'model', 'seed'], how='left')
sc['chains'] = sc['index'] - sc['ind'] - sc['ris'] - sc['val']
sc['qstar_idx'] = sc['index'] / (sc.direct - sc.lookup)
sc['qstar_all'] = (sc.sketch + sc['index']) / (sc.direct - sc.lookup)
sc['throughput'] = 1.0 / sc.lookup
g = sc.groupby(['tag', 'graph', 'model', 'n', 'm', 'M', 'b'])
agg = g.agg(sketch=('sketch', 'mean'), index=('index', 'mean'), ind=('ind', 'mean'), ris=('ris', 'mean'), val=('val', 'mean'), chains=('chains', 'mean'),
            L=('L', 'mean'), skMB=('skbytes', lambda x: x.mean() / 2**20), idxMB=('idxbytes', lambda x: x.mean() / 2**20),
            p50us=('p50', lambda x: x.mean() * 1e6), p95us=('p95', lambda x: x.mean() * 1e6), np50ms=('np50', lambda x: x.mean() * 1e3),
            direct=('direct', 'mean'), dp95=('dp95', 'mean'), dmax=('dmax', 'mean'), thr=('throughput', 'mean'),
            qstar_idx=('qstar_idx', 'mean'), qstar_all=('qstar_all', 'mean'), rss=('rss1', 'max'), nseeds=('seed', 'count')).reset_index()
agg.to_csv(f'{R}/scale_v2/scale_agg.csv', index=False)
pd.set_option('display.width', 250)
print(agg.round(3).to_string())
# quality on fresh worlds
qcols = ['tag', 'graph', 'model', 'seed', 'method', 'k', 'alpha_label', 'mean'] + ['c%d' % i for i in range(5)]
rows = []
for line in open(f'{R}/scale_v2/scale_quality.csv'):
    f = line.strip().split(',')
    tag, gname, model, seed, method, k, alab, mean = f[:8]; cv = list(map(float, f[8:]))
    rows.append(dict(tag=tag, graph=gname, model=model, seed=int(seed), method=method, k=int(k), alpha=float(alab), mean=float(mean), cvars=cv))
q = pd.DataFrame(rows)
# the per-row cvars list is over the run's alpha grid; pick CVaR at the record's own alpha
def grid(gname): return [0.05, 0.2, 1.0] if gname == 'livejournal' else [0.05, 0.1, 0.2, 0.5, 1.0]
q['cvar'] = [r.cvars[grid(r.graph).index(min(grid(r.graph), key=lambda a: abs(a - r.alpha)))] for r in q.itertuples()]
q = q.drop(columns=['cvars']).drop_duplicates(['graph', 'model', 'seed', 'method', 'k', 'alpha'])
p = q.pivot_table(index=['graph', 'model', 'seed', 'k', 'alpha'], columns='method', values='cvar').reset_index()
p['vsDirect'] = (p['RiskIndex'] / p['DirectCVaR-GRR'] - 1) * 100; p['vsRIS'] = (p['RiskIndex'] / p['RIS-IndRR'] - 1) * 100
p['DirvsRIS'] = (p['DirectCVaR-GRR'] / p['RIS-IndRR'] - 1) * 100
a = p[p.alpha == 0.05]
print('Tier C quality at alpha=0.05 (mean over seeds), RiskIndex vs Direct-CVaR / vs RIS, by graph and k:')
print(a.groupby(['graph', 'model', 'k'])[['vsDirect', 'vsRIS', 'DirvsRIS']].mean().round(2).to_string())
print('range over graphs/k: vsDirect', a.groupby(['graph', 'model', 'k']).vsDirect.mean().agg(['min', 'max']).round(2).to_dict(),
      ' vsRIS', a.groupby(['graph', 'model', 'k']).vsRIS.mean().agg(['min', 'max']).round(2).to_dict())
print('pooled alpha=0.05: vsDirect %.2f%%  vsRIS %.2f%%' % (a.vsDirect.mean(), a.vsRIS.mean()))
x = p[p.alpha <= 0.2]
print('pooled alpha<=0.2: vsDirect %.2f%%  vsRIS %.2f%%' % (x.vsDirect.mean(), x.vsRIS.mean()))
