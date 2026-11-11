#!/usr/bin/env python3
"""Need-gate development, Tier A: mean realized CVaR_alpha as % of the exhaustive optimum for RiskIndex-HV and
the gated layer (RiskIndex-AD rows of results/ad_dev/oracle_dev_theta<theta>.csv), dev sketch seeds 100-102.
The optimum is joined from results/oracle/oracle.csv (same selection/evaluation worlds)."""
import sys, pandas as pd, numpy as np
cols = ['tag', 'graph', 'fam', 'n', 'm', 'model', 'k', 'alpha', 'seed', 'method', 'cvB', 'meanB', 'cvA', 'optA', 'S']
old = pd.read_csv('results/oracle/oracle.csv', header=None, names=cols)
orc = old[old.method == 'Oracle'].groupby(['tag', 'graph', 'k', 'alpha']).cvB.first().rename('opt')   # optimum does not depend on the sketch seed
out = []
for th in sys.argv[1:]:
    d = pd.read_csv(f'results/ad_dev/oracle_dev_theta{th}.csv', header=None, names=cols)
    d = d.merge(orc, on=['tag', 'graph', 'k', 'alpha']); d['ratio'] = d.cvB / d.opt * 100
    t = d.groupby(['method', 'alpha']).ratio.mean().unstack()
    hv, lg = t.loc['RiskIndex-HV'], t.loc['RiskIndex-AD']
    n05 = len(d[(d.method == 'RiskIndex-AD') & (d.alpha == 0.05)])
    out.append(dict(theta=th, cases=n05, HV_05=hv[0.05], gated_05=lg[0.05], diff_05=lg[0.05] - hv[0.05], HV_all=hv.mean(), gated_all=lg.mean()))
print(pd.DataFrame(out).round(3).to_string(index=False))
