#!/bin/bash
# Aggregate the raw results and produce every table and figure of the paper (vector PDF + PNG, LaTeX table bodies).
# Usage: scripts/make_figures.sh [paper_dir]     (paper_dir: figures and table bodies are written into the paper)
set -e
cd "$(dirname "$0")/.."
export PYTHONPATH=experiments
PAPER=${1:-}

# ---- Fig. 1(a): exhaustive equal-mean / different-tail pair on a tiny graph, re-evaluated on fresh worlds
bin/riskindex pair graph=er:48:8:1 model=IC p=0.2 k=3 R=4000 alpha=0.1 tol=0.01 seed=7 2>/dev/null > results/fig1_pair.txt
bin/riskindex dist graph=er:48:8:1 model=IC p=0.2 sets="$(cat results/fig1_pair.txt)" R=20000 seed=909 > results/fig1_dist.csv

# ---- Fig. 5 (estimation fidelity, results/fidelity) and Fig. 6 (sensitivity of the GRR index, results/sens, sens_v2)
python3 experiments/analyze_fidelity.py > results/fidelity/summary.txt
python3 experiments/fig_fidelity.py                                                                    # figures/fig3_fidelity.pdf
(cd experiments && FIGDIR=../figures RESDIR=../results python3 fig_sens_v2.py)                         # figures/fig5_sens.pdf
python3 experiments/analyze_alloc_v2.py > results/alloc_v2/summary.txt                                 # budget split M vs b (Sec. V-G)
python3 experiments/analyze_cap_v2.py   > results/alloc_v2/cap_summary.txt                             # DBLP IC memory caps (Sec. V-G)

# ---- Fig. 2: schematic overview of the index (TikZ; needs XeLaTeX)
xelatex -interaction=nonstopmode -output-directory=figures experiments/fig2_arch.tex >/dev/null && rm -f figures/fig2_arch.aux figures/fig2_arch.log

# ---- Baseline and ablation tables that the final analysis joins: Tier B baselines (results/quality) paired with the
# index without the forward layer (results/quality_v2); Tier C aggregates of the index without the layer (results/scale_v2)
python3 experiments/merge_v2.py         > results/quality_v2/merge_log.txt
python3 experiments/analyze_scale_v2.py > results/scale_v2/summary.txt
# ---- Everything else: Tables II-IV, Figs. 1, 3, 4 and every number of Section V, from the final runs (results/*_ad)
# joined with the two ablation runs (results/*_v2: w/o forward layer; results/*_lg: w/o need gate).
scripts/ad_analyze.sh $PAPER
if [ -n "$PAPER" ]; then
  cp figures/fig3_fidelity.pdf figures/fig5_sens.pdf figures/fig2_arch.pdf $PAPER/figures/
fi
echo FIGURES_DONE
