#!/usr/bin/env python3
"""Numbers quoted in Section V-D (small budgets on large graphs) (Adaptive index).

Inputs (RESDIR, default results):
  headroom_lg/headroom.csv     tail headroom (seed 0): RIS, RiskIndex-HV, RiskIndex-LG, Oracle-CVaR, Oracle-basis on
                               the same 3000 evaluation worlds (k in {1,2,3,5}; alpha in {0.05, 0.1, 0.2})
  sweep_lg/sweep.csv(+_sys)    size sweep: ca-HepPh padded with isolated nodes, seeds 0-2
  quality_ad/merged.csv        Tier B (gain over RIS, tail ratio of the RIS answer)
  quality_{ad,lg}/quality_*_sys.csv forward-layer build times and served budgets (kf) per graph/model/seed
Gains: realized CVaR_alpha of a method / that of the RIS answer on the same case - 1, mean over k <= 5 and alpha <= 0.2
(the convention of the development log), and over alpha = 0.05 where stated (Table II)."""
import os, io
import numpy as np, pandas as pd
R = os.environ.get('RESDIR', 'results')
SMALL = ['email', 'wiki', 'facebook', 'hepph']; LARGE = ['epinions', 'slashdot', 'dblp']
pd.set_option('display.width', 220)

# ------------------------------------------------------------------ headroom
cols = ['graph', 'model', 'seed', 'M', 'b', 'method', 'k', 'alpha', 'mean', 'c05', 'c10', 'c20', 'h', 'rho', 'S']
h = pd.read_csv(f'{R}/headroom_ad/headroom.csv', header=None, names=cols); h['method'] = h.method.replace({'RiskIndex-LG': 'RiskIndex-AD'})
hl = pd.read_csv(f'{R}/headroom_lg/headroom.csv', header=None, names=cols); h = pd.concat([h, hl[hl.method == 'RiskIndex-LG']])
h = h[h.method != 'PoolEdge']
acol = {0.05: 'c05', 0.1: 'c10', 0.2: 'c20'}
rows = []
for (g, mo, sd, k), x in h.groupby(['graph', 'model', 'seed', 'k']):
    ris = x[x.method == 'RIS-IndRR'].iloc[0]
    for al, c in acol.items():
        r = dict(graph=g, model=mo, k=k, alpha=al, ris=ris[c], ris_mean=ris['mean'])
        for m in ['RiskIndex-HV', 'RiskIndex-LG', 'RiskIndex-AD', 'Oracle-CVaR', 'Oracle-basis']:
            y = x[(x.method == m) & np.isclose(x.alpha, al)]
            r[m] = (y[c].iloc[0] / ris[c] - 1) * 100 if len(y) else np.nan
        rows.append(r)
hr = pd.DataFrame(rows)
hr['ref_case'] = hr[['Oracle-CVaR', 'Oracle-basis']].max(axis=1)   # per-case max (optimistic)
t = hr.groupby(['graph', 'model'])[['RiskIndex-HV', 'RiskIndex-LG', 'RiskIndex-AD', 'Oracle-CVaR', 'Oracle-basis', 'ref_case']].mean()
t['ref'] = t[['Oracle-CVaR', 'Oracle-basis']].max(axis=1)   # the better reference per cell (Fig. 5a)
hr = hr.merge(t['ref'].reset_index(), on=['graph', 'model'])
print('== headroom (seed 0; k<=5, alpha<=0.2; gain over RIS, %)'); print(t.round(2).to_string())
for nm, gs in [('small', SMALL), ('large', LARGE), ('Epinions+Slashdot', ['epinions', 'slashdot'])]:
    x = hr[hr.graph.isin(gs)]
    print(f'  {nm:18s}', x[['RiskIndex-HV', 'RiskIndex-LG', 'RiskIndex-AD', 'Oracle-CVaR', 'Oracle-basis', 'ref']].mean().round(2).to_dict())
es = t.loc[['epinions', 'slashdot']]
print('  Epinions/Slashdot: max reference gain per cell %.2f%%' % es.ref.max(), ' max Oracle-basis %.2f%%' % es['Oracle-basis'].max())
print('  share of the reference gain captured (large graphs, cells with ref > 2%):')
for (g, mo), r in t.iterrows():
    if g in LARGE and r.ref > 2: print(f'    {g} {mo}: HV {r["RiskIndex-HV"]:.2f} LG {r["RiskIndex-LG"]:.2f} AD {r["RiskIndex-AD"]:.2f} ref {r.ref:.2f}')
# tail ratio of the RIS 1-prefix on the headroom worlds (CVaR_0.05 / mean)
tr = hr[(hr.k == 1) & (hr.alpha == 0.05)].assign(tr=lambda z: z.ris / z.ris_mean).groupby(['graph', 'model']).tr.mean()
print('  tail ratio CVaR_0.05/mean of the RIS 1-seed answer:', tr.round(2).to_dict())
# by k on the large graphs
print('  large graphs by k:', hr[hr.graph.isin(LARGE)].groupby('k')[['RiskIndex-HV', 'RiskIndex-LG', 'RiskIndex-AD', 'ref']].mean().round(2).to_dict())

# ------------------------------------------------------------------ size sweep
qc = ['tag', 'graph', 'model', 'M', 'b', 'H', 'Kmax', 'seed', 'method', 'k', 'alpha', 'sec', 'mean', 'cvar', 'uid', 'stats']
s = pd.read_csv(f'{R}/sweep_ad/sweep.csv', header=None, names=qc); s['method'] = s.method.replace({'RiskIndex-LG': 'RiskIndex-AD'})
sl = pd.read_csv(f'{R}/sweep_lg/sweep.csv', header=None, names=qc); s = pd.concat([s, sl[sl.method == 'RiskIndex-LG']])
s = s[(s.alpha <= 0.2) & (s.k <= 5)]
key = ['graph', 'model', 'seed', 'k', 'alpha']
ris = s[s.method == 'RIS-IndRR'].set_index(key).cvar
sw = []
for m in ['RiskIndex-HV', 'RiskIndex-LG', 'RiskIndex-AD']:
    x = s[s.method == m].set_index(key).cvar
    sw.append(((x / ris.reindex(x.index) - 1) * 100).rename(m))
sw = pd.concat(sw, axis=1).reset_index()
sw['factor'] = sw.graph.str.extract(r'pad(\d+)').astype(int)
st = sw.groupby(['model', 'factor'])[['RiskIndex-HV', 'RiskIndex-LG', 'RiskIndex-AD']].mean()
print('== size sweep (k<=5, alpha<=0.2; seeds 0-2)'); print(st.round(2).to_string())
seedsd = sw.groupby(['model', 'factor', 'seed'])[['RiskIndex-HV', 'RiskIndex-LG', 'RiskIndex-AD']].mean().groupby(['model', 'factor']).std()
print('  seed s.d.:'); print(seedsd.round(2).to_string())
for mo in ['WC', 'LT']:
    x = st.loc[mo]
    print(f'  {mo}: LG {x["RiskIndex-LG"].min():.1f}-{x["RiskIndex-LG"].max():.1f}%  HV {x["RiskIndex-HV"].min():.1f}-{x["RiskIndex-HV"].max():.1f}%  '
          f'HV x1 {x.loc[1, "RiskIndex-HV"]:.1f} -> x100 {x.loc[100, "RiskIndex-HV"]:.1f}; LG x1 {x.loc[1, "RiskIndex-LG"]:.1f} -> x100 {x.loc[100, "RiskIndex-LG"]:.1f}')
ss = [l for l in open(f'{R}/sweep_ad/sweep_sys.csv') if l.startswith('fl,')]
if ss:
    f = pd.read_csv(io.StringIO(''.join(ss)), header=None, names=['tag', 'graph', 'model', 'M', 'b', 'H', 'Kmax', 'seed', 'build', 'screen', 'greedy', 'val', 'pool', 'cands', 'bytes', 'kf'])
    print('  forward layer build (s) by size:', f.groupby(['model', 'graph']).build.mean().round(1).to_dict())
sb = [l for l in open(f'{R}/sweep_ad/sweep_sys.csv') if not l.startswith(('fl,', 'hyb,', 'val,', 'rho,'))]
if sb:  # tag, graph, model, n, m, M, b, ...
    b = pd.read_csv(io.StringIO(''.join(sb)), header=None)
    print('  n, b by size/model:', b.groupby([1, 2])[[3, 6]].first().apply(tuple, axis=1).to_dict())
    print('  sketch resolution n/b by size/model:', (b.groupby([1, 2])[3].first() / b.groupby([1, 2])[6].first()).round(1).to_dict())

# ------------------------------------------------------------------ Tier B: exact-layer cost and served budgets
fls = []
for p in ['quality_small_sys.csv', 'quality_dblp_sys.csv']:
    fn = f'{R}/quality_ad/{p}'
    if os.path.exists(fn): fls += [l for l in open(fn) if l.startswith('fl,')]
f = pd.read_csv(io.StringIO(''.join(fls)), header=None, names=['tag', 'graph', 'model', 'M', 'b', 'H', 'Kmax', 'seed', 'build', 'screen', 'greedy', 'val', 'pool', 'cands', 'bytes', 'kf'])
ft = f.groupby(['graph', 'model']).agg(build=('build', 'mean'), kf=('kf', 'max'), kfmin=('kf', 'min'), pool=('pool', 'mean'), bytes=('bytes', 'max'))
print('== Tier B forward layer (mean over seeds)'); print(ft.round(1).to_string())
print('  cells without the layer (kf=0):', ft[ft.kf == 0].index.tolist(), ' partial:', ft[(ft.kf > 0) & (ft.kf < 5)].index.tolist())
q = pd.read_csv(f'{R}/quality_ad/merged.csv', low_memory=False); q = q[q.tag == 'tierB']
x = q[(q.k <= 5) & (q.alpha <= 0.2)]
g = x[x.method.isin(['RiskIndex-HV', 'RiskIndex-LG', 'RiskIndex-AD'])].groupby(['graph', 'model', 'method']).rel.mean().unstack().sub(1).mul(100)
print('== Tier B gain over RIS (k<=5, alpha<=0.2)'); print(g.round(2).to_string())
for nm, gs in [('small', SMALL), ('large', LARGE)]:
    y = x[x.graph.isin(gs)]
    print(f'  {nm}: ', y[y.method.isin(['RiskIndex-HV', 'RiskIndex-LG', 'RiskIndex-AD'])].groupby('method').rel.mean().sub(1).mul(100).round(2).to_dict())
y = q[(q.k <= 5) & (q.alpha == 0.05)]
for nm, gs in [('small', SMALL), ('large', LARGE)]:
    z = y[y.graph.isin(gs)]
    print(f'  alpha=0.05 {nm}: ', z[z.method.isin(['RiskIndex-HV', 'RiskIndex-LG', 'RiskIndex-AD', 'DirectCVaR-GRR', 'SG'])].groupby('method').rel.mean().sub(1).mul(100).round(2).to_dict())
