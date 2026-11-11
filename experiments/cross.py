#!/usr/bin/env python3
"""Cross-variant comparison with a FIXED reference: Mean-GRR of variant A (M=512, in-sample).
Evaluation worlds are identical across variants for a given (graph, model, seed), so absolute
realized CVaR is comparable.  Usage: cross.py A.csv B.csv ... [--alphas ...] [--bygraph]"""
import sys, pandas as pd, numpy as np
from scipy.stats import wilcoxon
files = [f for f in sys.argv[1:] if not f.startswith('--')]
opts = dict(a.split('=') for a in sys.argv[1:] if a.startswith('--') and '=' in a)
alphas = [float(x) for x in opts.get('--alphas', '0.05,0.1,0.2').split(',')]
cols = ['tag','graph','model','M','b','H','Kmax','seed','method','k','alpha','sec','mean','cvar','setid','stats']
d = pd.concat([pd.read_csv(f, header=None, names=cols, usecols=range(16)) for f in files])
d['alpha'] = d['alpha'].astype(float)
key = ['graph','model','seed','k','alpha']
ref = d[(d.tag == d.tag.iloc[0]) & (d.method == 'Mean-GRR')][key + ['cvar']].rename(columns={'cvar': 'ref'})
m = d.merge(ref, on=key); m['rel'] = m.cvar / m.ref; m = m[m.alpha.isin(alphas)]
m['kgrp'] = np.where(m.k <= 5, 'k<=5', 'k>=10')
m['variant'] = m.tag + ':' + m.method
def summ(x):
    diff = x.rel - 1; w = (diff > 1e-9).sum(); l = (diff < -1e-9).sum()
    return pd.Series({'n': len(x), 'mean_rel%': 100 * diff.mean(), 'W/T/L': f'{w}/{len(x)-w-l}/{l}'})
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500)
print('== vs Mean-GRR of', d.tag.iloc[0], 'alphas', alphas)
print(m.groupby(['kgrp', 'variant']).apply(summ).round(3))
if '--bygraph' in sys.argv:
    print(m.groupby(['graph', 'model', 'kgrp', 'variant']).apply(summ)['mean_rel%'].round(2).unstack('variant'))
if '--byalpha' in sys.argv:
    print(m.groupby(['alpha', 'kgrp', 'variant']).apply(summ)['mean_rel%'].round(2).unstack('variant'))
