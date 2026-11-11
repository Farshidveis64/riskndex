#!/usr/bin/env python3
"""Fig. 7: sensitivity of the GRR index (RIS chain + validation; forward layer not involved) to H, M, b, and the
price of grouping.  Gains are relative to classical RIS (independent RR sketch of the same size) on the same case;
Direct-CVaR (GRR) of the H sweep comes from results/sens (identical sketches; it does not depend on H)."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import figstyle as fs
from figstyle import C
import os
R = os.environ.get('RESDIR', 'results')
COLS = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,sec,mean,cvar,uid,stats'.split(',')

q = pd.read_csv(f'{R}/sens_v2/sens.csv', header=None, names=COLS)
old = pd.read_csv(f'{R}/sens/sens.csv', header=None, names=COLS)
kk = ['tag', 'graph', 'model', 'M', 'b', 'H', 'Kmax', 'seed', 'k', 'alpha']
for m in ['Mean-GRR']:  # identical sketches: shared rows must reproduce
    j = q[q.method == m].merge(old[old.method == m], on=kk, suffixes=('', '_o')); assert ((j.cvar - j.cvar_o).abs() < 1e-6).all(), 'reproduction failed'
q = pd.concat([q, old[(old.method == 'DirectCVaR-GRR') & (old.tag == 'Hsweep')]], ignore_index=True)
ref = q[q.method == 'RIS-IndRR'][kk + ['cvar']].rename(columns={'cvar': 'ris'})
q = q.merge(ref, on=kk, how='left'); q['rel'] = q.cvar / q.ris
hc = pd.read_csv(f'{R}/sens_v2/sens.csv.hcert', header=None, names='tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,h,cert,rho,mc_cert,se,se_w,raw,guarded'.split(','))
rho = hc[hc.method == 'RiskIndex-HV']

def ci95(x):
    x = np.asarray(x, float); x = x[~np.isnan(x)]
    return 1.96 * x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else 0.0

# one line per graph/model setting: colour = graph, IC dashed with open markers, WC solid filled
SET = [('wiki', 'WC'), ('wiki', 'IC'), ('email', 'WC'), ('email', 'IC'), ('epinions', 'WC')]
GCOL = {'wiki': C['blue'], 'email': C['salmon'], 'epinions': C['sage']}
GMK = {'wiki': 'o', 'email': '^', 'epinions': 's'}
GNAME = {'wiki': 'Wiki-Vote', 'email': 'email', 'epinions': 'Epinions'}
def sty(key):
    col = GCOL[key[0]]; ic = key[1] == 'IC'
    return dict(color=col, marker=GMK[key[0]], ls=(0, (3, 1.5)) if ic else '-', mfc='white' if ic else col, mec=col, mew=0.6, ms=2.6, lw=0.8)
fig, axes = plt.subplots(1, 4, figsize=(fs.TEXTW, fs.ROWH), gridspec_kw=dict(wspace=0.38))

def draw(ax, df, xcol, ycol):
    """one thin line per setting (mean over seeds, budgets and alpha <= 0.2) and their mean (black dotted)."""
    df = df.copy(); df['pct'] = (df[ycol] - 1) * 100
    for key in SET:
        d = df[(df.graph == key[0]) & (df.model == key[1])].groupby(xcol).pct.mean()
        ax.plot(d.index, d.values, zorder=4, **sty(key))
    per = df.groupby(['graph', 'model', xcol]).pct.mean().groupby(level=xcol).mean()
    ax.plot(per.index, per.values, color=C['black'], lw=0.9, ls=(0, (1, 1)), zorder=5)
    ax.axhline(0, color=C['ink2'], lw=0.5, zorder=1)
    return per
tail = q[q.alpha <= 0.2]

# (a) H sweep: RiskIndex vs Direct-CVaR (same sketch), tail levels
ax = axes[0]
h = tail[(tail.tag == 'Hsweep')]
base = h[h.method == 'DirectCVaR-GRR'][['graph', 'model', 'seed', 'k', 'alpha', 'H', 'cvar']].rename(columns={'cvar': 'direct'})
hh = h[h.method == 'RiskIndex-HV'].merge(base, on=['graph', 'model', 'seed', 'k', 'alpha', 'H']); hh['r'] = hh.cvar / hh.direct
print('(a)', draw(ax, hh, 'H', 'r').round(2).to_dict())
ax.set_xscale('log', base=2); ax.set_xticks([1, 2, 4, 8, 16, 32]); ax.set_xticklabels(['1', '2', '4', '8', '16', '32'])
ax.set_xlabel('threshold chains $H$', labelpad=1); ax.set_ylabel('CVaR vs. Direct-CVaR (%)', labelpad=1)
fs.sub(ax, '(a) threshold chains $H$'); fs.style(ax)

# (b) M sweep: RiskIndex relative to RIS-Mean(GRR) at same M; and build time
ax = axes[1]
mm = tail[(tail.tag == 'Msweep') & (tail.method == 'RiskIndex-HV')]
print('(b)', draw(ax, mm, 'M', 'rel').round(2).to_dict())
ax.set_xscale('log', base=2); ax.set_xticks([64, 128, 256, 512, 1024]); ax.set_xticklabels(['64', '128', '256', '512', '1K'])
ax.set_xlabel('worlds $M$', labelpad=1); ax.set_ylabel('CVaR vs. RIS (%)', labelpad=1)
fs.sub(ax, '(b) scenario count $M$'); fs.style(ax)

# (c) b sweep: quality and tail-resolution ratio rho
ax = axes[2]
bb = tail[(tail.tag == 'bsweep') & (tail.method == 'RiskIndex-HV')]
print('(c)', draw(ax, bb, 'b', 'rel').round(2).to_dict())
ax.set_xscale('log', base=2); ax.set_xticks([64, 256, 1024, 4096, 16384]); ax.set_xticklabels(['64', '256', '1K', '4K', '16K'])
ax.set_xlabel('targets per world $b$', labelpad=1); ax.set_ylabel('CVaR vs. RIS (%)', labelpad=1)
rr = rho[(rho.tag == 'bsweep') & (rho.alpha <= 0.2)].groupby('b').rho.mean()
print('rho by b:', rr.round(2).to_dict())
fs.sub(ax, '(c) targets per world $b$'); fs.style(ax)

# (d) price of grouping: Mean-GRR vs RIS-IndRR realized CVaR_1 (= mean) as M grows
ax = axes[3]
pm = q[(q.tag == 'Mprice') & (q.alpha == 1.0)]
base = pm[pm.method == 'RIS-IndRR'][['graph', 'model', 'seed', 'k', 'M', 'cvar']].rename(columns={'cvar': 'ind'})
pp = pm[pm.method == 'Mean-GRR'].merge(base, on=['graph', 'model', 'seed', 'k', 'M']); pp['r'] = pp.cvar / pp.ind
print('(d)', draw(ax, pp, 'M', 'r').round(2).to_dict())
ph = pm[pm.method == 'RiskIndex-HV'].merge(base, on=['graph', 'model', 'seed', 'k', 'M']); ph['pct'] = (ph.cvar / ph.ind - 1) * 100
phm = ph.groupby(['graph', 'model', 'M']).pct.mean().groupby(level='M').mean()
ax.plot(phm.index, phm.values, color=C['purple'], lw=1.1, ls='--', marker='P', ms=3.0, zorder=6)
print('(d) RiskIndex at alpha=1 vs RIS', phm.round(2).to_dict())
ax.set_xscale('log', base=2); ax.set_xticks([128, 256, 512, 1024, 2048, 4096]); ax.set_xticklabels(['128', '256', '512', '1K', '2K', '4K'])
ax.set_xlabel('worlds $M$ (equal budget $Mb$)', labelpad=1); ax.set_ylabel('E[$\\sigma$] vs. RIS (%)', labelpad=1)
fs.sub(ax, '(d) price of grouping'); fs.style(ax)
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
hd = [Line2D([], [], label='%s %s' % (GNAME[k[0]], k[1]), **sty(k)) for k in SET]
hd.append(Line2D([], [], color=C['black'], lw=0.9, ls=(0, (1, 1)), label='mean'))
hd.append(Line2D([], [], color=C['purple'], lw=1.1, ls='--', marker='P', ms=3.0, label='RiskIndex, $\\alpha{=}1$ (d)'))
for a in axes: a.xaxis.set_minor_locator(plt.NullLocator())
fig.legend(handles=hd, loc='lower center', bbox_to_anchor=(0.5, 0.98), ncol=7, fontsize=5.8, frameon=False, handlelength=2.6, columnspacing=1.4)
fs.save(fig, 'fig5_sens', figdir=os.environ.get('FIGDIR', 'figures'), bbox=fs.ROW_BBOX)
