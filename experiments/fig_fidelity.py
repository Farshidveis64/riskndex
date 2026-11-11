#!/usr/bin/env python3
"""Figure 3 (RQ1): GRR vs independent RR at equal budget M*b (two-column, 4 panels).
Thin lines: the six graph/model settings; bold: their median.  One shared legend row on top,
no text inside the plot areas (same layout as Figs. 4-6)."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import figstyle as fs
from figstyle import C

a = pd.read_csv('results/fidelity/fidelity_agg.csv')
a = a[(a.group == 'random') & (a.tag == 'bsweep')]
graphs = [('email', 'WC'), ('email', 'IC'), ('wiki', 'WC'), ('wiki', 'IC'), ('wiki', 'LT'), ('epinions', 'WC')]
fig, axes = plt.subplots(1, 4, figsize=(fs.TEXTW, fs.ROWH), gridspec_kw=dict(wspace=0.38))

def panel(ax, typ, col, ycol, ls, mk='P', band=True):
    """median over the six settings, shaded interquartile range (no per-setting lines)"""
    allv = []
    for (g, mo) in graphs:
        dd = a[(a.graph == g) & (a.model == mo) & (a.type == typ)].sort_values('b')
        if dd.empty: continue
        allv.append(dd.set_index('b')[ycol])
    D = pd.concat(allv, axis=1)
    med, lo, hi = D.median(axis=1), D.quantile(0.25, axis=1), D.quantile(0.75, axis=1)
    if band: ax.fill_between(med.index, lo.values, hi.values, color=col, alpha=0.15, lw=0, zorder=2)
    ax.plot(med.index, med.values, color=col, ls=ls, lw=0.9, marker=mk, ms=2.6, mfc=col if mk != 'D' else 'white', mec=col, mew=0.5, zorder=4)

ax = axes[0]
panel(ax, 'GRR', C['purple'], 'nmae_rawcvar0.05', '-'); panel(ax, 'IndRR', C['salmon'], 'nmae_rawcvar0.05', (0, (3, 1.2)), 's')
ax.set_yscale('log'); ax.set_ylim(0.02, 5); ax.set_ylabel(r'NMAE of $\mathrm{CVaR}_{0.05}$', labelpad=1); fs.sub(ax, '(a) tail error')
ax = axes[1]
ax.axhline(0, color=C['black'], lw=0.5)
panel(ax, 'GRR', C['purple'], 'bias_rawcvar0.05', '-'); panel(ax, 'IndRR', C['salmon'], 'bias_rawcvar0.05', (0, (3, 1.2)), 's')
ax.set_ylim(-1.1, 2.2); ax.set_ylabel(r'signed bias of $\mathrm{CVaR}_{0.05}$', labelpad=1); fs.sub(ax, '(b) tail bias')
ax = axes[2]
panel(ax, 'GRR', C['purple'], 'nmae_mean', '-'); panel(ax, 'IndRR', C['salmon'], 'nmae_mean', (0, (3, 1.2)), 's')
ax.set_yscale('log'); ax.set_ylabel(r'NMAE of $\mathbb{E}[\sigma]$', labelpad=1); fs.sub(ax, '(c) mean error')
ax = axes[3]
panel(ax, 'GRR', C['gray'], 'nmae_rawcvar0.05', (0, (1, 1.2)), 'D', band=False); panel(ax, 'GRR', C['purple'], 'nmae_cvar0.05', '-')
ax.set_yscale('log'); ax.set_ylim(0.02, 5); ax.set_ylabel(r'NMAE of $\mathrm{CVaR}_{0.05}$', labelpad=1); fs.sub(ax, '(d) variance matching (GRR)')
for ax in axes:
    ax.set_xscale('log', base=2); ax.set_xlabel('targets per world $b$', labelpad=1)
    ax.set_xticks([1, 4, 16, 64, 256, 1024, 4096]); ax.set_xticklabels(['1', '4', '16', '64', '256', '1K', '4K'])
    ax.xaxis.set_minor_locator(plt.NullLocator()); fs.style(ax)
h = [Line2D([], [], color=C['purple'], lw=0.9, marker='P', ms=2.6, label='GRR (variance-matched in d)'),
     Line2D([], [], color=C['salmon'], lw=0.9, ls=(0, (3, 1.2)), marker='s', ms=2.6, label='independent RR'),
     Line2D([], [], color=C['gray'], lw=0.9, ls=(0, (1, 1.2)), marker='D', ms=2.6, mfc='white', label='GRR raw profile (d)'),
     Patch(facecolor=C['ink2'], alpha=0.18, label='band: interquartile range')]
fig.legend(handles=h, loc='lower center', bbox_to_anchor=(0.5, 0.98), ncol=4, fontsize=5.8, frameon=False, handlelength=2.6, columnspacing=1.4)
fs.save(fig, 'fig3_fidelity', bbox=fs.ROW_BBOX)
