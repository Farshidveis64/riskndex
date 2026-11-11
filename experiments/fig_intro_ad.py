"""Fig. 1 (introduction): the build in (b) charges sketch, index and the
forward layer (results/tables_ad/scale_agg_ad.csv, written by make_tables_ad.py).


(a) Spread distribution of two seed sets on one graph (ER, n=48, IC p=0.2) over 20,000 sampled
    worlds: results/fig1_dist.csv (row = seed set, then 20,000 realized spreads).
    A = {9,15,43}, B = {4,21,27}.  Step histograms (bin width 2); the worst 10% of each set is
    shaded; triangles: CVaR_0.1 of each set; dotted line: the (common) mean.
(b) Cumulative time of Q queries on DBLP (WC, Tier C, results/scale/scale_agg.csv and the 60 measured
    per-query times of results/scale/scale_queries.csv: 20 (k, alpha) queries x 3 seeds).  Per-query
    CVaR greedy: mean cumulative time Q * mean cost (line), and the 25-75% and 5-95% ranges over
    20,000 random workloads that draw their queries uniformly from the measured ones (bands).
    RiskIndex: build (sketch + index, both charged to RiskIndex) plus one P50 lookup per query.
    The lines cross at the pessimistic break-even Q* = (sketch + index)/(direct - lookup) of Table V.
Usage: python3 experiments/fig_intro.py [results_dir] [figure_dir]
"""
import sys, numpy as np, pandas as pd
import figstyle as fs
from figstyle import plt, C
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

RES = sys.argv[1] if len(sys.argv) > 1 else 'results'
OUT = sys.argv[2] if len(sys.argv) > 2 else 'figures'

# ---------------------------------------------------------------- data (a)
rows = [l.strip().split(',') for l in open(f'{RES}/fig1_dist.csv') if l.strip()]
dist = {('A' if r[0].startswith('9') else 'B'): np.array(r[1:], float) for r in rows}
def cvar(v, a):
    s = np.sort(v); m = len(s) * a; r = int(np.ceil(m))
    return (s[:r - 1].sum() + (m - r + 1) * s[r - 1]) / m
stat = {k: dict(mean=v.mean(), cvar=cvar(v, 0.1), var=np.sort(v)[int(np.ceil(0.1 * len(v))) - 1])
        for k, v in dist.items()}
for k, d in stat.items():
    print(f"{k}: mean {d['mean']:.2f}  CVaR_0.1 {d['cvar']:.2f}  VaR_0.1 {d['var']:.0f}")

# ---------------------------------------------------------------- data (b)
sc = pd.read_csv(f'{RES}/tables_ad/scale_agg_ad.csv')     # scale_v2 aggregates plus the exact-layer build (scale_lg)
row = sc[(sc.graph == 'dblp') & (sc.model == 'WC')].iloc[0]
build = row.sketch + row['index'] + row.fl  # s: sketch + index + forward layer, all charged to RiskIndex
direct = row.direct                         # s per query (mean over the (k, alpha) grid)
lookup = row.p50us * 1e-6                   # s per query (P50)
qd = pd.read_csv(f'{RES}/scale_v2/scale_queries.csv', header=None)
qt = qd[(qd[1] == 'dblp') & (qd[2] == 'WC')][6].values          # measured per-query times (s)
assert len(qt) == 60 and abs(qt.mean() - direct) < 1e-6 * max(1.0, direct) * 1e3
rng = np.random.default_rng(2027)
QMAX, NW = 1000, 20000
Cw = np.cumsum(rng.choice(qt, size=(NW, QMAX), replace=True), axis=1)
Qi = np.arange(1, QMAX + 1)
lo, hi = np.quantile(Cw, 0.05, axis=0), np.quantile(Cw, 0.95, axis=0)
be = np.argmax(Cw >= build + Qi * lookup, axis=1) + 1                # break-even of each random workload
print('random workloads: break-even median %d, 5-95%%: %d-%d; band at Q=1000: %.0f-%.0f s'
      % (np.median(be), np.quantile(be, 0.05), np.quantile(be, 0.95), lo[-1], hi[-1]))
qstar = build / (direct - lookup)
ratio = 1e3 * direct / (build + 1e3 * lookup)
print(f'DBLP: build {build:.1f}s (sketch {row.sketch:.1f} + index {row["index"]:.1f} + forward layer {row.fl:.1f}), direct {direct:.2f}s/query, '
      f'lookup {row.p50us:.2f}us, Q*={qstar:.2f}, ratio at 1000 queries {ratio:.0f}x')

# ---------------------------------------------------------------- figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fs.COLW, 1.20), gridspec_kw=dict(wspace=0.44))

# (a) step histograms with the worst 10% shaded
edges = np.arange(0, 50, 2)
col = {'A': (C['sky'], C['blue']), 'B': (C['salmon'], C['red'])}
for k in ['A', 'B']:
    h, _ = np.histogram(dist[k], bins=edges, density=True)
    fc, ec = col[k]
    # worst 10%: the part of the histogram at or below the 10% quantile, filled
    cut = stat[k]['var'] + 1
    ht = np.where(edges[:-1] < cut, h, 0.0)
    et = np.minimum(edges, cut)
    ax1.stairs(ht, et, fill=True, color=fc, alpha=0.55 if k == 'A' else 0.50, lw=0, zorder=2 if k == 'A' else 3)
    ax1.stairs(h, edges, color=ec, lw=1.0, zorder=4 if k == 'A' else 5, ls='-' if k == 'A' else (0, (3, 1.3)))
for k in ['A', 'B']:
    ax1.plot([stat[k]['cvar']], [0], marker='^', ms=4.6, color=col[k][1], mec='white', mew=0.4,
             ls='none', zorder=6, clip_on=False)
ax1.axvline(np.mean([stat['A']['mean'], stat['B']['mean']]), color=C['ink2'], lw=0.7, ls=(0, (1, 1.2)), zorder=4)
ax1.set_xlim(0, 48); ax1.set_ylim(0, 0.083)
ax1.set_xticks([0, 10, 20, 30, 40])
ax1.set_yticks([0, 0.02, 0.04, 0.06, 0.08]); ax1.set_yticklabels(['0', '.02', '.04', '.06', '.08'])
ax1.set_xlabel('realized spread (nodes)', labelpad=1.5)
ax1.set_ylabel('density', labelpad=1.5)
fs.style(ax1)
h = [Line2D([], [], color=C['blue'], lw=1.0, label='set A'),
     Line2D([], [], color=C['red'], lw=1.0, ls=(0, (3, 1.3)), label='set B'),
     Patch(facecolor='#d9d9d9', edgecolor='none', label='worst 10%'),
     Line2D([], [], color=C['ink2'], lw=0.7, ls=(0, (1, 1.2)), label='mean'),
     Line2D([], [], color=C['ink2'], marker='^', ms=4, ls='none', label=r'CVaR$_{0.1}$')]
lg = ax1.legend(handles=h, loc='upper left', frameon=True, framealpha=0.95, edgecolor='#cccccc', fancybox=False,
                fontsize=5.3, handlelength=1.4, handleheight=0.8, borderpad=0.35, labelspacing=0.22, borderaxespad=0.3)
lg.get_frame().set_linewidth(0.4)
fs.sub(ax1, '(a) Same mean, different tail', y=-0.37)

# (b) cumulative time on DBLP (log-log): mean line with markers, fan of measured-workload quantiles,
#     break-even point with a drop line and a triangle on the x-axis (as the CVaR triangles of (a))
q25, q75 = np.quantile(Cw, 0.25, axis=0), np.quantile(Cw, 0.75, axis=0)
qm = np.array([1, 10, 100, 1000])
cD, mD, lD = fs.METHOD['DirectCVaR-GRR']
cR, mR, lR = fs.METHOD['RiskIndex']
YLO = 0.5
ax2.fill_between(Qi, np.maximum(lo, YLO * 0.5), hi, color=C['salmon'], alpha=0.16, lw=0, zorder=1)
ax2.fill_between(Qi, np.maximum(q25, YLO * 0.5), q75, color=C['salmon'], alpha=0.30, lw=0, zorder=1)
ax2.plot(Qi, Qi * direct, color=cD, ls=lD, lw=1.0, zorder=3)
ax2.plot(qm, qm * direct, color=cD, marker=mD, ms=3.0, ls='none', zorder=3)
ax2.plot(Qi, build + Qi * lookup, color=cR, ls=lR, lw=1.1, zorder=4)
ax2.plot(qm, build + qm * lookup, color=cR, marker=mR, ms=3.4, ls='none', zorder=4)
# break-even: point, dotted drop line, triangle on the axis, label
ax2.plot([qstar, qstar], [YLO, build], color=C['ink2'], lw=0.7, ls=(0, (1, 1.2)), zorder=5)
ax2.plot([qstar], [build + qstar * lookup], marker='o', ms=3.4, color=C['black'], ls='none', zorder=6)
ax2.plot([qstar], [YLO], marker='^', ms=4.6, color=cR, mec='white', mew=0.4, ls='none', zorder=7, clip_on=False)
ax2.text(qstar * 1.22, YLO * 1.55, r'$Q^{\star}\!=\!%.1f$' % qstar, fontsize=5.6, color=C['black'], ha='left', va='bottom',
         zorder=7, bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.85))
# gap after 1000 queries
ax2.annotate('', xy=(1e3, build * 1.12), xytext=(1e3, 1e3 * direct / 1.12),
             arrowprops=dict(arrowstyle='<|-|>', mutation_scale=4.5, lw=0.6, color=C['ink2'], shrinkA=0, shrinkB=0))
ax2.text(760, np.sqrt(build * 1e3 * direct), r'%.0f$\times$' % ratio, fontsize=5.8, color=C['black'], ha='right', va='center')
ax2.set_xscale('log'); ax2.set_yscale('log')
ax2.set_xlim(0.85, 1.25e3); ax2.set_ylim(YLO, 3e4)
ax2.set_xticks([1, 10, 100, 1000]); ax2.set_xticklabels(['1', '10', '100', '1000'])
ax2.set_yticks([1, 10, 100, 1000, 10000])
ax2.minorticks_off()
ax2.set_xlabel('number of queries', labelpad=1.5)
ax2.set_ylabel('cumulative time (s)', labelpad=1.5)
fs.style(ax2)
from matplotlib.legend_handler import HandlerTuple
h2 = [(Patch(facecolor=C['salmon'], alpha=0.40, edgecolor='none'),
       Line2D([], [], color=cD, ls=lD, marker=mD, ms=3.0, lw=1.0)),
      Line2D([], [], color=cR, ls=lR, marker=mR, ms=3.4, lw=1.1)]
lg2 = ax2.legend(h2, ['Direct-CVaR', 'RiskIndex'], handler_map={tuple: HandlerTuple(ndivide=1)},
                 loc='upper left', frameon=True, framealpha=0.95, edgecolor='#cccccc', fancybox=False,
                 fontsize=5.3, handlelength=2.0, handleheight=0.9, borderpad=0.35, labelspacing=0.22, borderaxespad=0.3)
lg2.get_frame().set_linewidth(0.4)
fs.sub(ax2, '(b) Amortization on DBLP', y=-0.37)

fs.save(fig, 'fig1_intro', OUT)
