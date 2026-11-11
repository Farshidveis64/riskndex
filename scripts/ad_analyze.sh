#!/bin/bash
# All analyses, tables and figures of the paper, from the final runs (scripts/ad_runner.py, THETA=5).
# Usage: scripts/ad_analyze.sh [paper_dir]
set -e
cd "$(dirname "$0")/.."
PAPER=${1:-}
mkdir -p results/tables_ad
python3 experiments/merge_ad.py | tee results/quality_ad/merge_log.txt
python3 experiments/headline_ad.py > results/quality_ad/headline.txt
python3 experiments/make_quality_table_ad.py results/quality_ad/merged.csv results/tables_ad/quality.tex > /dev/null
[ -s results/oracle_ad/oracle_ad.csv ] && python3 experiments/analyze_oracle_ad.py > results/oracle_ad/analysis.txt
if ls results/scale_ad/scale_hyb_?.csv > /dev/null 2>&1; then
  cat $(ls results/scale_ad/scale_hyb_?.csv) > results/scale_ad/scale_hyb.csv
  cat $(ls results/scale_ad/scale_?.csv) > results/scale_ad/scale.csv
  python3 experiments/make_tables_ad.py > results/tables_ad/make_tables_log.txt
fi
[ -s results/headroom_ad/headroom.csv ] && [ -s results/sweep_ad/sweep.csv ] && python3 experiments/ad_large_text.py > results/ad_large_text.txt
python3 experiments/ad_theta_table.py > results/ad_dev/theta_table_log.txt                     # need-gate threshold, Sec. V-E (development data)
[ -s results/ad_dev/tierA_gate_diag.log ] && python3 experiments/ad_tierA_gate.py > results/ad_dev/tierA_gate_analysis.txt
if [ -n "$PAPER" ]; then
  cd experiments
  RESDIR=../results FIGDIR=$PAPER/figures python3 fig_quality_ad.py > ../results/quality_ad/fig_quality_log.txt
  [ -s ../results/headroom_ad/headroom.csv ] && FIGDIR=$PAPER/figures python3 fig_large_ad.py > ../results/ad_fig_large_log.txt
  [ -s ../results/tables_ad/scale_agg_ad.csv ] && python3 fig_intro_ad.py ../results $PAPER/figures > ../results/ad_fig_intro_log.txt
  cd ..; python3 experiments/insert_tables.py $PAPER/main.tex results
fi
[ -s results/quality_sgg/quality_a.csv ] && python3 experiments/sgg_control.py > /dev/null                     # Sec. V-E control (results/quality_sgg/analysis.txt)
echo AD_ANALYZE_DONE
