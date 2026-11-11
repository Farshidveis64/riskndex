#!/bin/bash
# Reproduce every table and figure of the paper from scratch (single core each; the stages are independent and
# may run concurrently on separate cores).  See README.md for what each stage produces.
set -e
cd "$(dirname "$0")/.."
scripts/build.sh
scripts/download_data.sh
# 1. GRR sketch alone: estimation fidelity (Fig. 5) and parameter sensitivity of the GRR index (Fig. 6, Sec. V-G)
scripts/run_fidelity.sh                       # -> results/fidelity   (~30 min)
scripts/run_sensitivity.sh                    # -> results/sens       (~2 h)   Wiki-Vote, email
scripts/run_sensitivity_epinions.sh           #                       (~1.5 h) Epinions, appends
scripts/run_sensitivity_v2.sh                 # -> results/sens_v2    (~3 h)   same sweeps for the index with RIS chain + validation
scripts/run_alloc_v2.sh; scripts/run_dblp_cap_v2.sh   # -> results/alloc_v2 budget split M vs b; DBLP IC memory caps
# 2a. Per-query baselines (RIS, GRR mean chain, Direct-CVaR, SG, exhaustive optimum), Tiers A/B
scripts/run_oracle.sh                         # -> results/oracle     (~2.5 h)
scripts/run_quality.sh; scripts/run_quality_part2.sh   # -> results/quality (~10 h; part 2 = Epinions, Slashdot, DBLP)
# 2b. Ablation "w/o forward layer" (GRR index with RIS chain and validation sketch, and its variants), same sketches and worlds
scripts/run_quality_v2.sh small; scripts/run_quality_v2.sh dblp; scripts/run_dblp_inddirect.sh   # -> results/quality_v2 (~12 h)
scripts/run_oracle_v2.sh                      # -> results/oracle_v2  (~2.5 h)
scripts/run_scale_v2.sh; scripts/run_scale_v2_lj.sh   # -> results/scale_v2 (~2 h; LiveJournal needs ~2 GB RAM)
# 3. Ablation "w/o need gate" (forward layer with the cost gate only; also the forward-simulated references)
scripts/lg_runner.py A; scripts/lg_runner.py B        # -> results/quality_lg, headroom_lg, oracle_lg, sweep_lg, scale_lg (~10 h)
# 4. The index of the paper (forward layer gated by need and cost, theta = 5)
THETA=5 scripts/ad_runner.py A; THETA=5 scripts/ad_runner.py B; THETA=5 scripts/ad_runner.py C   # -> results/*_ad (~8 h)
scripts/run_quality_sgg.sh a; scripts/run_quality_sgg.sh b   # -> results/quality_sgg (~25 min): matched-grid control of the ablation "w/o tail selection"
# 5. Analysis, tables, figures
scripts/make_figures.sh
