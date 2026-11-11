#!/usr/bin/env python3
"""RQ1: GRR vs independent RR estimation fidelity (equal budget M*b).

Input : results/fidelity/fidelity.csv (raw per-set rows written by `riskindex fidelity`)
Output: results/fidelity/fidelity_summary.csv (per config x seed metrics) and
        results/fidelity/fidelity_agg.csv (mean/sd over seeds)
"""
import sys, numpy as np, pandas as pd
from scipy.stats import spearmanr, kendalltau

ALPHAS = [0.05, 0.1, 0.2]
cols = ['tag', 'graph', 'model', 'M', 'b', 'type', 'seed', 'set', 'kind', 'true_mean', 'est_mean']
for a in ALPHAS: cols += [f'true_cvar{a}', f'raw_cvar{a}', f'est_cvar{a}']
cols += ['rho', 'build_sec', 'L']
df = pd.read_csv('results/fidelity/fidelity.csv', header=None, names=cols)

def partial_spearman(x, y, z):
    """rank correlation of x and y after removing the rank-linear effect of z"""
    rx, ry, rz = pd.Series(x).rank(), pd.Series(y).rank(), pd.Series(z).rank()
    def resid(a):
        A = np.vstack([rz, np.ones(len(rz))]).T
        beta = np.linalg.lstsq(A, a, rcond=None)[0]; return a - A @ beta
    ex, ey = resid(rx.values), resid(ry.values)
    if ex.std() == 0 or ey.std() == 0: return 0.0
    return float(np.corrcoef(ex, ey)[0, 1])

rows = []
df['group'] = np.where(df.kind == 'random', 'random', 'family')
for key, d in df.groupby(['tag', 'graph', 'model', 'M', 'b', 'type', 'seed', 'group']):
    r = dict(zip(['tag', 'graph', 'model', 'M', 'b', 'type', 'seed', 'group'], key))
    r['nmae_mean'] = np.mean(np.abs(d.est_mean - d.true_mean)) / np.mean(d.true_mean)
    r['rho_mean'] = spearmanr(d.est_mean, d.true_mean).correlation
    r['rho'] = d.rho.mean()
    for a in ALPHAS:
        t, e, raw = d[f'true_cvar{a}'], d[f'est_cvar{a}'], d[f'raw_cvar{a}']
        r[f'nmae_cvar{a}'] = np.mean(np.abs(e - t)) / np.mean(t)
        r[f'bias_cvar{a}'] = np.mean(e - t) / np.mean(t)
        r[f'nmae_rawcvar{a}'] = np.mean(np.abs(raw - t)) / np.mean(t)
        r[f'bias_rawcvar{a}'] = np.mean(raw - t) / np.mean(t)
        r[f'rho_cvar{a}'] = spearmanr(e, t).correlation if e.std() > 0 else 0.0
        r[f'tau_cvar{a}'] = kendalltau(e, t).correlation if e.std() > 0 else 0.0
        # tail-gap fidelity: does the sketch rank sets by their *excess* downside
        # (mean - CVaR) rather than merely by their mean?
        tg, eg = d.true_mean - t, d.est_mean - e
        r[f'rho_gap{a}'] = spearmanr(eg, tg).correlation if eg.std() > 0 else 0.0
        r[f'nmae_gap{a}'] = np.mean(np.abs(eg - tg)) / np.mean(tg)
        # partial rank correlation controlling for the true mean (random sets)
        r[f'prho_cvar{a}'] = partial_spearman(e.values, t.values, d.true_mean.values)
        # within-budget ranking for family sets: average Spearman over k blocks
        if key[-1] == 'family':
            rs = []
            for kk, dd in d.groupby('kind'):
                if len(dd) >= 4 and dd[f'est_cvar{a}'].std() > 0 and dd[f'true_cvar{a}'].std() > 0:
                    rs.append(spearmanr(dd[f'est_cvar{a}'], dd[f'true_cvar{a}']).correlation)
            r[f'wrho_cvar{a}'] = np.mean(rs) if rs else np.nan
            # top-1 regret: true CVaR of the sketch-argmax candidate vs the true best, per k
            regs = []
            for kk, dd in d.groupby('kind'):
                best_true = dd[f'true_cvar{a}'].max(); pick = dd.loc[dd[f'est_cvar{a}'].idxmax(), f'true_cvar{a}']
                regs.append(1 - pick / best_true if best_true > 0 else 0)
            r[f'top1reg{a}'] = np.mean(regs)
    r['build_sec'] = d.build_sec.iloc[0]; r['L'] = d.L.iloc[0]
    rows.append(r)
s = pd.DataFrame(rows)
s.to_csv('results/fidelity/fidelity_summary.csv', index=False)
g = s.groupby(['tag', 'graph', 'model', 'M', 'b', 'type', 'group'])
agg = g.mean(numeric_only=True).drop(columns=['seed']).join(g.std(numeric_only=True).drop(columns=['seed']), rsuffix='_sd').reset_index()
agg.to_csv('results/fidelity/fidelity_agg.csv', index=False)
pd.set_option('display.width', 250)
print(agg[agg.group=='random'][['tag', 'graph', 'model', 'M', 'b', 'type', 'rho', 'nmae_mean', 'nmae_rawcvar0.05', 'bias_rawcvar0.05', 'nmae_cvar0.05', 'bias_cvar0.05', 'prho_cvar0.05', 'nmae_cvar0.1', 'bias_cvar0.1', 'nmae_cvar0.2']].round(3).to_string())
