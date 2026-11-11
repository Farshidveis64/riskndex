#!/usr/bin/env python3
"""Fig. 4 (RiskIndex = RiskIndex-AD) (reference: classical RIS on an independent RR sketch of the same size):
one line per method, identical colour/marker/dash per method in every panel, one legend row.
(a) Tier A: realized CVaR as % of the exhaustive out-of-sample optimum vs alpha.
(b) Tier B: realized CVaR_0.05 gain over RIS (ind. RR) vs budget k.
(c) Tier B, k <= 5: the same gain vs risk level alpha (realized CVaR_alpha).
(d) Tier B: the gain vs quartile of the tail ratio CVaR_0.05/E[sigma] of the RIS answer (Q1 dispersed).
Means with 95% CIs.  The previous grouped-bar version is fig_quality_bars.py."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import figstyle as fs
from figstyle import C, METHOD, LABEL
import os
R = os.environ.get('RESDIR', 'results')
METHOD = dict(METHOD); LABEL = dict(LABEL)
METHOD['RiskIndex-AD'] = METHOD['RiskIndex']; LABEL['RiskIndex-AD'] = 'RiskIndex'
LABEL['RIS-IndRR'] = 'RIS (ind. RR)'; LABEL['Mean-GRR'] = 'Mean chain (GRR)'

def ci95(x):
    x = np.asarray(x, float); x = x[~np.isnan(x)]
    return 1.96 * x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else 0.0

def sty(m, main=False):
    col, mk, ls = METHOD[m]
    return dict(color=col, marker=mk, ls=ls, lw=1.0 if main else 0.75, ms=3.2 if main else 2.4,
                mfc=col, mec=col, mew=0.5, capsize=1.2, elinewidth=0.5, zorder=4 if main else 3)

def band(ax, x, mu, ci, m, both=False):
    main = m in ('RiskIndex', 'RiskIndex-AD')
    st = sty(m, main); col = st['color']
    if main or (both and m == 'DirectCVaR-GRR'):
        ax.fill_between(x, mu - ci, mu + ci, color=col, alpha=0.18 if main else 0.14, lw=0, zorder=st['zorder'] - 1)
    ax.plot(x, mu, color=col, marker=st['marker'], ls=st['ls'], lw=st['lw'], ms=st['ms'], mfc=st['mfc'], mec=st['mec'], mew=st['mew'], zorder=st['zorder'])

fig, axes = plt.subplots(1, 4, figsize=(fs.TEXTW, fs.ROWH), gridspec_kw=dict(wspace=0.36))

# ---------------------------------------------------------------- (a) Tier A
cols = ['tag', 'graph', 'fam', 'n', 'm', 'model', 'k', 'alpha', 'seed', 'method', 'cvB', 'meanB', 'cvA', 'optA', 'S']
df = pd.concat([pd.read_csv(f'{R}/oracle/oracle.csv', header=None, names=cols),
                pd.read_csv(f'{R}/oracle_ad/oracle_ad.csv', header=None, names=cols).replace({'method': {'RiskIndex-LG': 'RiskIndex-AD'}}).query("method == 'RiskIndex-AD'")])
orc = df[df.method == 'Oracle'][['tag', 'graph', 'k', 'alpha', 'seed', 'cvB']].rename(columns={'cvB': 'opt'})
d = df.merge(orc, on=['tag', 'graph', 'k', 'alpha', 'seed']); d['ratio'] = d.cvB / d.opt * 100
AL = [0.05, 0.1, 0.2, 0.5, 1.0]; xa = np.arange(len(AL))
ax = axes[0]
for m in ['Oracle-Mean', 'Mean-GRR', 'RIS-IndRR', 'SG', 'DirectCVaR-GRR', 'RiskIndex-AD']:
    g = d[d.method == m].groupby('alpha').ratio
    band(ax, xa, g.mean().reindex(AL).values, g.apply(ci95).reindex(AL).values, m, both=True)
ax.set_xticks(xa); ax.set_xticklabels(['.05', '.1', '.2', '.5', '1'])
ax.set_ylim(97.4, 100.25); ax.set_yticks([98, 99, 100])
ax.set_xlabel(r'risk level $\alpha$', labelpad=1); ax.set_ylabel('% of optimum', labelpad=1)
fs.sub(ax, r'(a) Tier A: vs. exhaustive optimum'); fs.style(ax)

# ---------------------------------------------------------------- Tier B
q = pd.read_csv(f'{R}/quality_ad/merged.csv', low_memory=False)
q = q[(q.tag == 'tierB') & np.isfinite(q.rel) & (q.rel > 0)].copy(); q['g'] = (q.rel - 1) * 100
TB = ['Mean-GRR', 'SG', 'DirectCVaR-GRR', 'RiskIndex-AD']
def lines(ax, df, xcol, xs):
    x = np.arange(len(xs)); out = {}
    for m in TB:
        g = df[df.method == m].groupby(xcol).g
        mu = g.mean().reindex(xs); out[m] = mu.round(1).tolist()
        band(ax, x, mu.values, g.apply(ci95).reindex(xs).values, m)
    c, _, ls = METHOD['RIS-IndRR']
    ax.axhline(0, color=c, lw=0.75, ls=ls, zorder=2)
    ax.set_xticks(x); fs.style(ax); return out

ax = axes[1]
K = [1, 2, 3, 5, 10, 25, 50, 100]
ob = lines(ax, q[q.alpha == 0.05], 'k', K)
ax.set_xticklabels([str(k) for k in K]); ax.tick_params(axis='x', labelsize=5.0); ax.axvspan(3.5, 7.5, color=C['gray'], alpha=0.10, lw=0, zorder=0)
ax.set_ylim(-8, 60); ax.set_yticks([0, 20, 40, 60])
ax.set_xlabel('budget $k$', labelpad=1); ax.set_ylabel(r'CVaR$_{0.05}$ gain (%)', labelpad=1)
fs.sub(ax, '(b) Tier B: by budget')

ax = axes[2]
oc = lines(ax, q[q.k <= 5], 'alpha', AL)
ax.set_xticklabels(['.05', '.1', '.2', '.5', '1']); ax.set_ylim(-3, 20); ax.set_yticks([0, 5, 10, 15, 20])
ax.set_xlabel(r'risk level $\alpha$ ($k\leq5$)', labelpad=1); ax.set_ylabel(r'CVaR$_\alpha$ gain (%)', labelpad=1)
fs.sub(ax, '(c) Tier B: by risk level')

ris = q[q.method == 'RIS-IndRR'].pivot_table(index=['graph', 'model', 'k'], columns='alpha', values='cvar', aggfunc='mean')
reg = (ris[0.05] / ris[1.0]).reset_index(name='tail_ratio05')
QS = ['Q1', 'Q2', 'Q3', 'Q4']; reg['Q'] = pd.qcut(reg.tail_ratio05, 4, labels=QS).astype(str)
qq = q[q.alpha == 0.05].merge(reg, on=['graph', 'model', 'k'])
ax = axes[3]
od = lines(ax, qq, 'Q', QS)
ax.set_xticklabels(['Q1', 'Q2', 'Q3', 'Q4']); ax.set_ylim(-6, 42); ax.set_yticks([0, 10, 20, 30, 40])
ax.set_xlabel(r'dispersion quartile', labelpad=1); ax.set_ylabel(r'CVaR$_{0.05}$ gain (%)', labelpad=1)
fs.sub(ax, '(d) Tier B: by dispersion')

for ax in axes: ax.xaxis.set_minor_locator(plt.NullLocator())
order = ['RiskIndex-AD', 'DirectCVaR-GRR', 'SG', 'Mean-GRR', 'RIS-IndRR', 'Oracle-Mean']
h = [Line2D([], [], **{k: v for k, v in sty(m, m == 'RiskIndex-AD').items() if k not in ('capsize', 'elinewidth', 'zorder')}, label=LABEL[m] + (' (a)' if m == 'Oracle-Mean' else '')) for m in order]
fig.legend(handles=h, loc='lower center', bbox_to_anchor=(0.5, 0.98), ncol=6, fontsize=5.8, frameon=False, handlelength=2.6, columnspacing=1.6)
print('(b)', ob); print('(c)', oc); print('(d)', od)
fs.save(fig, 'fig4_quality', figdir=os.environ.get('FIGDIR', 'figures'))
