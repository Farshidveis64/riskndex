#!/usr/bin/env python3
"""Certificate calibration: (cert - realized)/realized, coverage of cert +- z*SE, rho diagnostics.
Usage: cert.py quality.csv  (expects quality.csv.cert alongside)"""
import sys, pandas as pd, numpy as np
f = sys.argv[1]
cols = ['tag','graph','model','M','b','H','Kmax','seed','method','k','alpha','sec','mean','cvar','setid','stats']
q = pd.read_csv(f, header=None, names=cols, usecols=range(16)); q['alpha'] = q.alpha.astype(float)
c = pd.read_csv(f + '.cert', header=None, names=['tag','graph','model','M','b','H','Kmax','seed','k','alpha','h','cert','rho','meancert','se','raw'])
c['alpha'] = c.alpha.astype(float)
key = ['tag','graph','model','M','b','H','Kmax','seed','k','alpha']
m = c.merge(q[q.method == 'RiskIndex'][key + ['cvar']], on=key)
m['err'] = (m.cert - m.cvar) / m.cvar
m['rawerr'] = (m.raw - m.cvar) / m.cvar
m['z'] = (m.cert - m.cvar) / m.se.replace(0, np.nan)
m['in1'] = (m.cert - m.cvar).abs() <= 1 * m.se
m['in2'] = (m.cert - m.cvar).abs() <= 2 * m.se
m['kgrp'] = np.where(m.k <= 5, 'k<=5', 'k>=10')
pd.set_option('display.width', 250)
def s(x): return pd.Series({'n': len(x), 'mean_err%': 100 * x.err.mean(), 'med_err%': 100 * x.err.median(), 'p90_abs%': 100 * x.err.abs().quantile(0.9),
                            'raw_med%': 100 * x.rawerr.median(), 'cov1se': x.in1.mean(), 'cov2se': x.in2.mean(), 'se/cert': (x.se / x.cert).median(), 'rho': x.rho.mean()})
print('== by alpha'); print(m.groupby('alpha').apply(s).round(3))
print('== by alpha x kgrp'); print(m.groupby(['alpha', 'kgrp']).apply(s).round(3))
print('== by graph/model (alpha<=0.2)'); print(m[m.alpha <= 0.2].groupby(['graph', 'model']).apply(s).round(3))
