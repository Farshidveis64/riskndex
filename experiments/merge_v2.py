#!/usr/bin/env python3
"""Merge the runs of the index without the forward layer (results/quality_v2) with the Tier B baseline runs (results/quality).

The reruns use the identical sketches (seeds) and evaluation worlds, so every baseline row of
results/quality/quality_merged.csv stays paired with the new index rows; rows that both runs contain
(RiskIndex, Mean-GRR, RIS-IndRR) are checked for exact equality.  Output: results/quality_v2/merged.csv with
  rel      = realized CVaR / that of RIS (ind. RR) on the same case   (reference of the paper)
  rel_mg   = realized CVaR / that of the GRR mean chain (Mean-GRR)     (reference of an earlier draft; unused)
and results/quality_v2/certificates.csv for RiskIndex-HV.
"""
import glob, os, sys
import numpy as np, pandas as pd

R = os.environ.get('RESDIR', 'results')
COLS = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,sec,mean,cvar,uid,stats'.split(',')
HC = 'tag,graph,model,M,b,H,Kmax,seed,method,k,alpha,h,cert,rho,mc_cert,se,se_w,raw,guarded'.split(',')
key = ['graph', 'model', 'seed', 'k', 'alpha']

old = pd.read_csv(f'{R}/quality/quality_merged.csv', low_memory=False)
old = old[old.tag == 'tierB'][COLS]
parts = [p for p in [f'{R}/quality_v2/quality_small.csv', f'{R}/quality_v2/quality_dblp.csv'] if os.path.exists(p) and os.path.getsize(p) > 0]
new = pd.concat([pd.read_csv(p, header=None, names=COLS) for p in parts])
ind = f'{R}/quality_v2/dblp_inddirect.csv'
extra = pd.read_csv(ind, header=None, names=COLS) if os.path.exists(ind) and os.path.getsize(ind) > 0 else pd.DataFrame(columns=COLS)
extra = extra[extra.method == 'DirectCVaR-IndRR']

# exact reproduction of shared rows
for m in ['RiskIndex', 'Mean-GRR', 'RIS-IndRR']:
    j = new[new.method == m].merge(old[old.method == m], on=key, suffixes=('', '_o'))
    bad = ((j.cvar - j.cvar_o).abs() > 1e-6).sum() + ((j['mean'] - j['mean_o']).abs() > 1e-6).sum()
    print(f'reproduction {m:10s}: {len(j)} rows, {bad} mismatches'); assert bad == 0
# a new run may also have produced DirectCVaR-IndRR rows that the submitted runs lack (DBLP)
have = old[old.method == 'DirectCVaR-IndRR'][key]
extra = extra.merge(have, on=key, how='left', indicator=True); extra = extra[extra._merge == 'left_only'][COLS]

newm = new[~new.method.isin(['RiskIndex', 'Mean-GRR', 'RIS-IndRR'])]
cases = new[new.method == 'RIS-IndRR'][key].drop_duplicates()
d = pd.concat([old.merge(cases, on=key), newm, extra.merge(cases, on=key)], ignore_index=True)
ref = d[d.method == 'RIS-IndRR'][key + ['cvar', 'mean']].rename(columns={'cvar': 'cvar_ris', 'mean': 'mean_ris'})
mg = d[d.method == 'Mean-GRR'][key + ['cvar', 'mean']].rename(columns={'cvar': 'cvar_mg', 'mean': 'mean_mg'})
d = d.merge(ref, on=key).merge(mg, on=key)
d['rel'] = d.cvar / d.cvar_ris; d['rel_mg'] = d.cvar / d.cvar_mg; d['relmean'] = d['mean'] / d.mean_ris
d.to_csv(f'{R}/quality_v2/merged.csv', index=False)
print('merged rows', len(d), 'cases', len(cases), 'methods', sorted(d.method.unique()))

# certificates of the final index (and of the in-sample variants, for comparison)
hc = pd.concat([pd.read_csv(p + '.hcert', header=None, names=HC) for p in parts])
hc = hc.merge(d[['method'] + key + ['cvar', 'mean']], on=['method'] + key)
hc['calib'] = hc.cert / hc.cvar - 1
hc['calib_raw'] = hc.raw / hc.cvar - 1
hc['in2se'] = (hc.cert - hc.cvar).abs() <= 2 * hc.se
hc['in2sew'] = (hc.cert - hc.cvar).abs() <= 2 * hc.se_w
hc['kgrp'] = np.where(hc.k <= 5, 'k<=5', 'k>=10')
hc.to_csv(f'{R}/quality_v2/certificates.csv', index=False)
print('certificate rows', len(hc))
