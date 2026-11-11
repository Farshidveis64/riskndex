#!/usr/bin/env python3
"""Dev analysis: relative realized CVaR of every method vs Mean-GRR, paired by
(graph, model, cfg, seed, k, alpha).  Usage: an.py file.csv [file2.csv ...] [--ref METHOD] [--alphas 0.05,0.1,0.2]"""
import sys, pandas as pd, numpy as np
from scipy.stats import wilcoxon
files = [f for f in sys.argv[1:] if not f.startswith('--')]
opts = dict(a.split('=') for a in sys.argv[1:] if a.startswith('--') and '=' in a)
ref = opts.get('--ref', 'Mean-GRR')
alphas = [float(x) for x in opts.get('--alphas', '0.05,0.1,0.2').split(',')]
cols = ['tag','graph','model','M','b','H','Kmax','seed','method','k','alpha','sec','mean','cvar','setid','stats']
d = pd.concat([pd.read_csv(f, header=None, names=cols, usecols=range(16)) for f in files])
d['alpha'] = d['alpha'].astype(float)
key = ['tag','graph','model','M','b','H','Kmax','seed','k','alpha']
r = d[d.method == ref][key + ['cvar','setid']].rename(columns={'cvar':'ref','setid':'refid'})
m = d.merge(r, on=key)
m['rel'] = m.cvar / m.ref
m['same'] = (m.setid == m.refid)
m = m[m.alpha.isin(alphas)]
def summ(x):
    diff = x.rel - 1
    w = (diff > 1e-9).sum(); l = (diff < -1e-9).sum(); t = len(x) - w - l
    try: p = wilcoxon(x.cvar, x.ref).pvalue if w + l > 0 else 1.0
    except Exception: p = float('nan')
    return pd.Series({'n': len(x), 'mean_rel%': 100 * diff.mean(), 'gmean_rel%': 100 * (np.exp(np.log(x.rel.clip(1e-9)).mean()) - 1),
                      'W/T/L': f'{w}/{t}/{l}', 'p': p, 'same%': 100 * x.same.mean(), 'sec': x.sec.mean()})
m['kgrp'] = np.where(m.k <= 5, 'k<=5', 'k>=10')
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500)
print('== by method x k-group (alphas', alphas, ') rel. to', ref)
print(m.groupby(['kgrp','method']).apply(summ).round(4))
if '--bygraph' in sys.argv:
    print('== by graph/model x method, k<=5')
    print(m[m.k <= 5].groupby(['graph','model','method']).apply(summ)[['n','mean_rel%','W/T/L']].round(2))
    print('== by graph/model x method, k>=10')
    print(m[m.k >= 10].groupby(['graph','model','method']).apply(summ)[['n','mean_rel%','W/T/L']].round(2))
if '--byalpha' in sys.argv:
    print('== by alpha x method')
    print(m.groupby(['alpha','method']).apply(summ)[['n','mean_rel%','W/T/L','p']].round(4))
if '--byk' in sys.argv:
    print('== by k x method')
    print(m.groupby(['k','method']).apply(summ)[['n','mean_rel%','W/T/L']].round(3))
