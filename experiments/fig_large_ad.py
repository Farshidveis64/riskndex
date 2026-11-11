#!/usr/bin/env python3
"""LargeGraph figure: (a) tail headroom at k <= 5 per graph/model: gain over RIS of the forward-simulated reference
(x) vs. RiskIndex (filled) and RiskIndex without the forward layer (hollow) on the same 3000 fresh worlds;
(b) size sweep at fixed influence structure (ca-HepPh padded with isolated nodes): gain over RIS at k <= 5.
Gains: realized CVaR_alpha, mean over k in {1,2,3,5} and alpha in {0.05,0.1,0.2} (and seeds in (b)).
Usage: fig_large_lg.py  (HEADROOM, SWEEP, FIGDIR environment variables)"""
import os, numpy as np, pandas as pd
import figstyle as fs
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
C = fs.C
H = os.environ.get('HEADROOM', '../results/headroom_ad/headroom.csv'); SW = os.environ.get('SWEEP', '../results/sweep_ad/sweep.csv')
HL = os.environ.get('HEADROOM_LG', '../results/headroom_lg/headroom.csv'); SWL = os.environ.get('SWEEP_LG', '../results/sweep_lg/sweep.csv')
cols = ['graph','model','seed','M','b','method','k','alpha','mean','c05','c10','c20','h','rho','S']
d = pd.concat([pd.read_csv(f, header=None, names=cols) for f in H.split(',')]); d['method'] = d.method.replace({'RiskIndex-LG': 'RiskIndex-AD'})
dl = pd.read_csv(HL, header=None, names=cols); d = pd.concat([d, dl[dl.method == 'RiskIndex-LG']])
d = d[d.method != 'PoolEdge']
acol = {0.05: 'c05', 0.1: 'c10', 0.2: 'c20'}
rows = []
for (g, mo, sd, k), x in d.groupby(['graph','model','seed','k']):
    ris = x[x.method == 'RIS-IndRR'].iloc[0]
    for al, c in acol.items():
        def v(m):
            y = x[(x.method == m) & np.isclose(x.alpha, al)][c]
            return y.iloc[0] if len(y) else v('RiskIndex-HV')  # (development files have no LG rows)
        rows.append(dict(graph=g, model=mo, oc=v('Oracle-CVaR') / ris[c] - 1, ob=v('Oracle-basis') / ris[c] - 1, lg=v('RiskIndex-AD') / ris[c] - 1, ng=v('RiskIndex-LG') / ris[c] - 1, hv=v('RiskIndex-HV') / ris[c] - 1))
r = pd.DataFrame(rows).groupby(['graph','model'])[['oc','ob','lg','ng','hv']].mean().reset_index()
for c in ['oc','ob','lg','ng','hv']: r[c] *= 100
r['ref'] = r[['oc', 'ob']].max(axis=1)   # the better of the two forward-simulated references (cell means)
small = ['email','wiki','facebook','hepph']
print(r.round(2).to_string(index=False))
fig, axes = plt.subplots(1, 2, figsize=(fs.COLW, 0.92), gridspec_kw=dict(wspace=0.40))
ax = axes[0]
MK = {'WC': 'o', 'LT': 's', 'IC': '^'}
large = r[~r.graph.isin(small)]
ABB = {'epinions': 'Ep', 'slashdot': 'Sl', 'dblp': 'DBLP'}
for _, p in large.iterrows():
    ax.plot([p.ref, p.ref], [p.hv, p.lg], color=C['gray'], lw=0.5, zorder=2)
    ax.plot(p.ref, p.hv, marker=MK[p.model], ms=3.2, mfc='white', mec=C['salmon'], mew=0.8, ls='none', zorder=3)
    ax.plot(p.ref, p.ng, marker=MK[p.model], ms=2.0, color=C['gray'], ls='none', zorder=3)
    ax.plot(p.ref, p.lg, marker=MK[p.model], ms=3.2, color=C['purple'], ls='none', zorder=4)
    top = p.ref == large.ref.max()
    if p.ref > 2: ax.annotate(f'{ABB[p.graph]} {p.model}', (p.ref, p.lg), xytext=(-5, 0) if top else (-4, 1.5), textcoords='offset points', ha='right', va='center' if top else 'bottom', fontsize=5.4, color=C['ink2'], bbox=dict(boxstyle='round,pad=0.08', fc='white', ec='none', alpha=0.85), zorder=6)
ax.text(0.97, 0.04, 'Epinions, Slashdot:\nall six near 0', transform=ax.transAxes, ha='right', va='bottom', fontsize=5.4, color=C['ink2'], linespacing=1.0)
lim = [-2, 40]
ax.plot(lim, lim, color=C['gray'], lw=0.6, ls=(0, (3, 2)), zorder=1)
for a in (ax.set_xscale, ax.set_yscale): a('symlog', linthresh=1, linscale=0.5)
ax.set_xlim(-0.5, 36); ax.set_ylim(-1, 36)
ax.set_xticks([0, 1, 5, 20]); ax.set_xticklabels(['0', '1', '5', '20']); ax.set_yticks([0, 1, 5, 20]); ax.set_yticklabels(['0', '1', '5', '20'])
ax.set_xlabel('reference gain (%)', labelpad=1); ax.set_ylabel('gain over RIS (%)', labelpad=1)
fs.sub(ax, '(a) Epinions, Slashdot, DBLP', y=-0.46); fs.style(ax)
# (b) size sweep
cq = ['tag','graph','model','M','b','H','Kmax','seed','method','k','alpha','sec','mean','cvar','uid','stats']
q = pd.read_csv(SW, header=None, names=cq); q['method'] = q.method.replace({'RiskIndex-LG': 'RiskIndex-AD'})
ql = pd.read_csv(SWL, header=None, names=cq); q = pd.concat([q, ql[ql.method == 'RiskIndex-LG']]); q = q[(q.alpha <= 0.2) & (q.k <= 5)]
key = ['graph','model','seed','k','alpha']
ris = q[q.method == 'RIS-IndRR'].set_index(key).cvar
ax = axes[1]
for m, ls, mk, col, lab in (('RiskIndex-AD', '-', 'P', C['purple'], 'RiskIndex'), ('RiskIndex-LG', ':', 'o', C['gray'], 'w/o need gate'), ('RiskIndex-HV', '--', 'o', C['salmon'], 'w/o forward layer')):
    x = q[q.method == m].set_index(key).cvar
    gg = ((x / ris.reindex(x.index) - 1) * 100).reset_index()
    gg['f'] = gg.graph.str.extract(r'pad(\d+)').astype(int); gg['n'] = gg.f * 12008
    for mo in ('WC', 'LT'):
        s = gg[gg.model == mo].groupby('n').cvar.mean()
        print(m, mo, s.round(1).to_dict())
        ax.plot(s.index, s.values, color=col, ls=ls, marker=MK[mo], ms=3.0 if m != 'RiskIndex-LG' else 2.0, mfc=col if m != 'RiskIndex-HV' else 'white', mec=col, lw=0.9 if m != 'RiskIndex-LG' else 0.7)
ax.set_xscale('log'); ax.set_xticks([12008, 60040, 240160, 1200800]); ax.set_xticklabels(['12K', '60K', '240K', '1.2M'])
ax.xaxis.set_minor_locator(plt.NullLocator())
ax.set_xlabel('nodes (isolated-node padding)', labelpad=1); ax.set_ylabel('gain over RIS (%)', labelpad=1)
ax.set_ylim(0, None)
fs.sub(ax, '(b) graph size (ca-HepPh)', y=-0.46); fs.style(ax)
hd = [Line2D([], [], color=C['purple'], marker='P', ms=3.0, lw=0.9, label='RiskIndex'),
      Line2D([], [], color=C['gray'], marker='o', ms=2.0, lw=0.7, ls=':', label='w/o need gate'),
      Line2D([], [], color=C['salmon'], marker='o', ms=3.0, mfc='white', lw=0.9, ls='--', label='w/o forward layer'),
      Line2D([], [], color=C['ink2'], marker='o', ms=3.0, ls='none', label='WC'),
      Line2D([], [], color=C['ink2'], marker='s', ms=3.0, ls='none', label='LT'),
      Line2D([], [], color=C['ink2'], marker='^', ms=3.0, ls='none', label='IC')]
fig.legend(handles=hd, loc='lower center', bbox_to_anchor=(0.5, 0.955), ncol=6, fontsize=5.3, frameon=False, handlelength=1.8, columnspacing=0.7, borderaxespad=0.1)
fs.save(fig, 'fig_large', figdir=os.environ.get('FIGDIR', '../figures'))
