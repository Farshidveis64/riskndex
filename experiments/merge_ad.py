#!/usr/bin/env python3
"""Merge the final runs (results/quality_ad; the gated layer, recorded by the driver as
"RiskIndex-AD"; older runs used the label "RiskIndex-LG", renamed here) and the ungated-layer ablation (results/quality_lg, "RiskIndex-LG") with the Tier B table of the index without the forward layer (results/quality_v2/merged.csv).

The reruns use the identical sketches (seeds) and evaluation worlds, so every baseline of results/quality_v2 stays
paired with the new rows.  Rows that both runs contain (RiskIndex-HV, Mean-GRR, RIS-IndRR) are checked for exact
equality.  Output: results/quality_ad/merged.csv, which is results/quality_v2/merged.csv plus the RiskIndex-AD and RiskIndex-LG
rows (rel = realized CVaR / that of RIS on the same case), and results/quality_ad/certificates.csv for RiskIndex-AD.
"""
import os
import numpy as np, pandas as pd

R = os.environ.get('RESDIR', 'results')
COLS = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,sec,mean,cvar,uid,stats'.split(',')
HC = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,h,cert,rho,mc_cert,se,se_w,raw,guarded'.split(',')
key = ['graph', 'model', 'seed', 'k', 'alpha']
old = pd.read_csv(f'{R}/quality_v2/merged.csv', low_memory=False)
parts = [p for p in [f'{R}/quality_ad/quality_small.csv', f'{R}/quality_ad/quality_dblp.csv'] if os.path.exists(p) and os.path.getsize(p) > 0]
new = pd.concat([pd.read_csv(p, header=None, names=COLS) for p in parts]); new['method'] = new.method.replace({'RiskIndex-LG': 'RiskIndex-AD'})
lgparts = [f'{R}/quality_lg/quality_small.csv', f'{R}/quality_lg/quality_dblp.csv']
lgrows = pd.concat([pd.read_csv(p, header=None, names=COLS) for p in lgparts]); lgrows = lgrows[lgrows.method == 'RiskIndex-LG']
for m in ['RiskIndex-HV', 'Mean-GRR', 'RIS-IndRR']:
    j = new[new.method == m].merge(old[old.method == m], on=key, suffixes=('', '_o'))
    bad = ((j.cvar - j.cvar_o).abs() > 1e-6).sum() + ((j['mean'] - j['mean_o']).abs() > 1e-6).sum()
    print(f'reproduction {m:12s}: {len(j)} rows, {bad} mismatches'); assert bad == 0 and len(j) > 0
lg = pd.concat([new[new.method == 'RiskIndex-AD'][COLS], lgrows[COLS]])
cases = old[old.method == 'RIS-IndRR'][key + ['cvar', 'mean', 'cvar_mg', 'mean_mg']].rename(columns={'cvar': 'cvar_ris', 'mean': 'mean_ris'})
lg = lg.merge(cases, on=key)
lg['rel'] = lg.cvar / lg.cvar_ris; lg['rel_mg'] = lg.cvar / lg.cvar_mg; lg['relmean'] = lg['mean'] / lg.mean_ris
d = pd.concat([old, lg], ignore_index=True)
d.to_csv(f'{R}/quality_ad/merged.csv', index=False)
print('merged rows', len(d), 'AD+LG rows', len(lg), 'methods', sorted(d.method.unique()))
hc = pd.concat([pd.read_csv(p + '.hcert', header=None, names=HC) for p in parts])
hc['method'] = hc.method.replace({'RiskIndex-LG': 'RiskIndex-AD'})
hc = hc[hc.method == 'RiskIndex-AD'].merge(d[['method'] + key + ['cvar', 'mean']], on=['method'] + key)
# Answers of the forward layer are recorded with rho = 1, raw = cert and se_w = 0: they have no within-world term, so
# their full standard error is the between-world one (se).  Answers of k not served by the layer are RiskIndex-HV's.
hc['exact'] = np.isclose(hc.rho, 1.0) & (hc.se_w == 0) & np.isclose(hc.raw, hc.cert)
hc.loc[hc.exact, 'se_w'] = hc.loc[hc.exact, 'se']
hc['calib'] = hc.cert / hc.cvar - 1
hc['calib_raw'] = hc.raw / hc.cvar - 1
hc['in2se'] = (hc.cert - hc.cvar).abs() <= 2 * hc.se
hc['in2sew'] = (hc.cert - hc.cvar).abs() <= 2 * hc.se_w
hc['kgrp'] = np.where(hc.k <= 5, 'k<=5', 'k>=10')
hc.to_csv(f'{R}/quality_ad/certificates.csv', index=False)
print('certificate rows', len(hc), 'served by the forward layer', int(hc.exact.sum()))
