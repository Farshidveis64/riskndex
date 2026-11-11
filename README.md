# RiskIndex: Fast Tail-Risk Influence Maximization via Reusable Graph Indexing

Reproducibility artifact for the ICDE 2027 submission. Every number in every table and figure of the paper is
produced by the code in this directory from the public SNAP graphs and deterministic seeds; `AUDIT.md` maps each
claim of the paper to the raw result file and the script that computes it.

## The index in one paragraph

RiskIndex answers tail-risk influence queries (k, alpha) from a structure built once per graph. A grouped
reverse-reachability (GRR) sketch draws b RR sets per sampled world for M worlds (`src/grr`), so that a seed set has a
scenario profile instead of a single hit rate. One lazy-greedy chain per threshold of a family of truncated spreads
(`src/risk_basis`), plus the classical RIS chain on an independent RR sketch, provides the candidates; the sorted
variance-matched profile of every chain prefix is materialized on a second, independent validation sketch and stored
as prefix sums (`src/index`), so a query is O(H) array reads and returns the seed set with an out-of-sample
certificate, a resolution ratio rho_hat and a standard error. The index is tail-resolution adaptive: for the leading
budgets k <= k_f whose RIS-prefix tail the sketch cannot resolve (need gate R_tail(k) = b * C_0.05 / n < theta) and
whose spreads are small enough to simulate (cost gate sigma_hat <= sigma_max), a forward layer (`src/forward`) stores
profiles of exact per-world spreads from W forward-simulated worlds and serves those budgets. R_tail is read from the
variance-matched RIS-chain profile stored on the validation sketch (`tail_ratios` in `src/main.cpp`, `cvar_from_prefix`
at alpha = 0.05), the same quantity the paper defines.

## Method labels in the raw result files

| Label | Meaning in the paper |
|---|---|
| `RiskIndex-AD` | RiskIndex (the index of the paper): GRR index + forward layer gated by need and cost (`fl=1 fltheta=5`) |
| `RiskIndex-LG` | ablation "w/o need gate": forward layer with the cost gate only (`fl=1`, `fltheta` unset) |
| `RiskIndex-HV` | ablation "w/o forward layer": GRR index with RIS chain and validation sketch (`hyb=1 val=1`) |
| `RiskIndex-H`, `RiskIndex-V`, `RiskIndex-HV-raw`, `RiskIndex` | ablations w/o validation, w/o RIS chain, w/o variance matching, w/o both (Table III rows) |
| `RIS-IndRR` | classical RIS (IMM selection step on the independent RR sketch); reference of all gains |
| `Mean-GRR` | the same greedy on the GRR build sketch (mean chain) |
| `DirectCVaR-GRR`, `DirectCVaR-GRR-VM`, `DirectCVaR-IndRR` | per-query CVaR greedy on the build sketch / on the variance-matched objective / on an independent RR sketch |
| `SG`, `SG-cached` | sequential greedy of Zhou & Tokekar, per query / chains cached |
| `Oracle-CVaR`, `Oracle-Mean`, `Oracle-basis` | exhaustive optima (Tier A, Sec. V-C) and the two forward-simulated references of Sec. V-D |

The three index rows come from three runs over identical sketches and evaluation worlds (the seeds derive every
random stream), so they are paired with each other and with every baseline; `merge_ad.py`, `merge_v2.py` and
`make_tables_ad.py` check the shared rows for exact equality (0 mismatches, see `AUDIT.md`).

## Layout

```
src/                  single-threaded C++17 implementation (header-only library + driver)
  common.hpp          hashing, RNG, empirical lower-tail CVaR, variance-matched profile (Prop. IV.1)
  graph.hpp           CSR graph, SNAP loader, synthetic generators (ER/BA/WS/SBM), padding for the size sweep
  diffusion.hpp       IC/WC/LT live-edge worlds (scenario-consistent lazy edge states), RR traversal
  grr/grr.hpp         grouped RR sketch (Step 1): M worlds x b targets, inverted lists, profiles; grouped=false gives ordinary RR
  risk_basis/         truncated-spread objective, marginal gains (Eq. gain), lazy-greedy chains (Step 2)
  index/index.hpp     threshold grid, RIS chain, validation sketch, prefix sums (+ squares, rho_hat, mean), O(H) query (Steps 3-4),
                      certificate / raw certificate / standard error, rho-guard, budget query, NoMat ablation
  forward/forward.hpp forward layer (Step 3b): screening, forward chains, exact profiles, need and cost gates
  baselines/          Direct-CVaR greedy (valid pruning bounds, raw and variance-matched), mean greedy, SG and SG-cached
  eval.hpp            fresh-world evaluation (common random numbers), exhaustive oracle (n <= 64)
  main.cpp            driver: quality | oracle | headroom | scale | fidelity | selftest | info | pair | dist
  tools/toy_fig3.cpp  small worked example of the library API (not used by the paper)
scripts/              build, data download, one script per experiment, the resumable runners, make_figures.sh
  dev/, lg_dev/, ad_dev/   development logs (dev seeds 100-102): selection rule and RIS chain + validation sketch;
                      forward layer and cost gate; need gate with the pre-registered rule for theta
experiments/          analysis (pandas/scipy) and figure generators (matplotlib); fig_toy_grr.tex, fig2_arch.tex (TikZ)
configs/parameters.md every parameter of every experiment
results/              raw outputs of every run (see "Results" below) and the generated tables
figures/              generated figures
AUDIT.md              claim in the paper -> result file -> script
```

## Results directories

| Directory | Content | Produced by |
|---|---|---|
| `results/quality_ad`, `oracle_ad`, `headroom_ad`, `sweep_ad`, `scale_ad` | final runs of the index of the paper (Tier B, Tier A, forward-simulated references, size sweep, Tier C) | `THETA=5 scripts/ad_runner.py A\|B\|C` |
| `results/quality_lg`, `oracle_lg`, `headroom_lg`, `sweep_lg`, `scale_lg` | ablation w/o need gate (same jobs, `fltheta` unset) | `scripts/lg_runner.py A\|B` (= `scripts/run_*_lg.sh`) |
| `results/quality_v2`, `oracle_v2`, `scale_v2`, `sens_v2`, `alloc_v2` | ablation w/o forward layer and its variants; H/M/b sweeps and price of grouping (Fig. 7); budget split and DBLP memory caps (Sec. V-G) | `scripts/run_*_v2.sh`, `run_dblp_inddirect.sh`, `run_dblp_cap_v2.sh` |
| `results/quality`, `oracle`, `scale`, `sens` | per-query baselines, exhaustive optima and the sketch-only sweeps, same sketches and worlds | `scripts/run_quality.sh`, `run_quality_part2.sh`, `run_oracle.sh`, `run_scale.sh`, `run_sensitivity*.sh` |
| `results/quality_sgg` | matched-grid control for the ablation row "w/o tail selection" (SG rule on the index's own grid, six graphs) | `scripts/run_quality_sgg.sh a\|b`, `experiments/sgg_control.py` |
| `results/fidelity` | estimation fidelity of GRR vs. independent RR (Fig. 6) | `scripts/run_fidelity.sh` |
| `results/dev`, `lg_dev`, `ad_dev` | development runs behind the frozen configuration; `ad_dev` also holds the Tier A gate diagnostic of the reported configuration (analysis only) | `scripts/dev/*`, `lg_dev/*`, `ad_dev/*` |
| `results/tables_ad` | generated LaTeX bodies of Tables II-IV and the Tier C aggregates | `scripts/ad_analyze.sh` |

## Requirements

- A C++17 compiler (tested with GCC 13.3); no external C++ libraries (header-only code, single-threaded).
- Python 3 (tested with 3.13) and the packages in `requirements.txt` (`pip install -r requirements.txt`) for analysis and figures.
- TeX Live (pdflatex, bibtex, xelatex) only to rebuild Fig. 2 and Fig. 3 and the paper.
- All reported runs fit in 8 GB of RAM (2-vCPU VM, see the environment below).

Quick check (about one minute):

```bash
scripts/build.sh && bin/riskindex selftest     # prints three "0 failures" lines and exits 0
scripts/ad_analyze.sh                          # recomputes all tables and headline numbers from results/ (about 30 s)
```

## Reproduce

```bash
scripts/build.sh                 # -> bin/riskindex (g++ -std=c++17 -O3)
scripts/download_data.sh         # SNAP graphs into datasets/raw
scripts/run_all.sh               # all runs (stages 1-4, about 60 h on one core; stages are independent)
scripts/make_figures.sh [paper]  # analysis, tables, figures; with a paper directory, writes figures/ and table bodies into main.tex
bin/riskindex selftest           # lazy-vs-exhaustive greedy, certificate recovery, query/query_score/budget_query agreement (GRR and routed); exits 0 on success
```

`scripts/run_all.sh` lists the stages: (1) sketch-only experiments (fidelity, sensitivity, budget split), (2) per-query
baselines and the index without the forward layer, (3) the ablation without the need gate, which also produces the
forward-simulated references, (4) the final runs with `THETA=5`. `scripts/ad_status.sh` shows the state of the
resumable runners. `scripts/ad_analyze.sh` alone redoes the analysis from existing raw results (about 30 s).

Environment used for the reported numbers: Intel Xeon @ 2.80 GHz (2 vCPU VM, 8 GB RAM), Linux, GCC 13.3.0
(`-O3 -march=native`), Python 3.13, numpy, scipy, pandas, matplotlib 3.11, TeX Live (pdflatex, bibtex, xelatex).
All timings are single-threaded wall-clock on one core.

## Driver usage

```
bin/riskindex quality graph=snap:PATH:directed|undirected name=NAME model=WC|LT|IC [p=..] \
    M=1024 brule=100 bmax=65536 H=8 Kmax=100 ks=1,2,3,5,10,25,50,100 alphas=0.05,0.1,0.2,0.5,1.0 seeds=0,1,2 R=3000 \
    hyb=1 val=1 rhomin=0.5 fl=1 fltheta=5 [flk=5 flsmax=1000 flWsel=1000 flWval=1000 flPscreen=2000 flP=60 flT=7 flalphas=0.05,0.1,0.2] \
    [ind=1 direct=1 inddirect=1 nomat=1 rawidx=1 directvm=1 sg=1 sgc=1] out=... sysout=...
      # .hcert file: tag,graph,model,cfg,method,k,alpha,h,cert,rho,mc_cert,se,se_w,raw,guarded
bin/riskindex oracle   graph=sbm:48:4:0.5:0.02:1 model=IC p=0.2 k=3 M=4096 b=64 H=8 Kmax=3 lite=1 hyb=1 val=1 rhomin=0.5 fl=1 fltheta=5 ...
bin/riskindex headroom graph=... model=... k=5 ...        # forward-simulated references of Sec. V-D (screened pool, 3000 fresh worlds)
bin/riskindex scale    graph=... M=512 brule=100 bmax=16384 H=8 Kmax=100 hyb=1 val=1 fl=1 fltheta=5 Q=1000 R=2000 ...
bin/riskindex fidelity graph=... M=256 bs=1,4,16,64,256,1024,4096 nsets=100 k=5 R=20000 ...
bin/riskindex info     graph=...                          # n, m, degrees, spectral radius, p_c = 1/lambda_1
bin/riskindex pair     graph=er:48:8:1 model=IC p=0.2 k=3 R=4000 alpha=0.1 tol=0.01   # Fig. 1(a) pair search
bin/riskindex dist     graph=... model=... sets="a b c;d e f" R=20000               # spread distributions of given sets
```

Graph specs: `snap:path:directed|undirected`, `pad:path:directed|undirected:factor` (padding with isolated nodes,
size sweep), `ba:n:m:seed`, `er:n:avgdeg:seed`, `ws:n:k:beta:seed`, `sbm:n:blocks:p_in:p_out:seed`. Every run is
deterministic in its seeds; selection (sketch) and evaluation streams are disjoint. `fl=1` without `fltheta` builds
the layer with the cost gate only (ablation w/o need gate); `rhomin` is the rho-guard.

## Mapping paper <-> code

| Paper | Code / script |
|---|---|
| Eq. (6), scenario estimate | `GRRSketch::profile`, `CoverState::profile` |
| Prop. IV.1 (variance identity, matched profile), rho_hat | `vm_shrink` in `src/common.hpp` |
| Eq. (7), Lemma IV.5, Eq. (8) | `CoverState::gain_tau` |
| Algorithm 1, Steps 1-3 (sketches, chains, materialization) | `GRRSketch::build`, `build_chain`, `ThresholdGrid::make`, `RiskIndex::build`, `RiskIndex::rematerialize` |
| RIS chain (tau = inf) on the independent RR sketch | `mean_greedy` on `GRRSketch::build(..., grouped=false)`, `materialize_sequence` |
| Step 3b, forward layer and its gates; Corollary IV.9 | `ForwardLayer::build` / `ForwardLayer::query`, `tail_ratios` and `prefix_spreads` in `src/main.cpp` |
| Routed query / budget query (Step 4: N_k in {M, W}, guard on the GRR directory only) | `adaptive_query`, `adaptive_query_score`, `adaptive_budget_query` in `src/forward/forward.hpp` |
| Step 4, Eqs. (10), (11); rho-guard; standard error | `RiskIndex::select` (shared by `query`, `query_score`, `budget_query`), `cvar_from_prefix`, `tail_se_from_prefix`, `tail_se_within` |
| Corollary IV.7 (raw certificate) | `RiskIndex::raw_cert` |
| Budget query | `RiskIndex::budget_query` |
| Resolution rule for b | `auto_b` in `src/main.cpp` |
| Direct-CVaR with valid pruning bounds (raw / VM) | `direct_cvar_greedy` |
| SG, SG-cached | `seq_greedy`, `SeqGreedyCache` |
| Fresh-world evaluation / exhaustive oracle | `evaluate_sets`, `Oracle` in `src/eval.hpp` |
| Fig. 1 | `experiments/fig_intro_ad.py` (`results/fig1_dist.csv`, `results/tables_ad/scale_agg_ad.csv`, `results/scale_v2/scale_queries.csv`) |
| Fig. 2 (toy example of tail blindness; hand-built, 12 RR sets) | `experiments/fig_toy_grr.tex` (pdfLaTeX) -> `fig_toy_grr.pdf` |
| Fig. 3 (overview; schematic, instances illustrative) | `experiments/fig2_arch.tex` (XeLaTeX) -> `fig2_arch.pdf` |
| Fig. 4 (decision quality) | `experiments/fig_quality_ad.py` -> `fig4_quality.pdf` |
| Fig. 5 (small budgets on large graphs) | `experiments/fig_large_ad.py` -> `fig_large.pdf`; numbers: `experiments/ad_large_text.py` |
| Fig. 6 (estimation fidelity) | `experiments/fig_fidelity.py` -> `fig3_fidelity.pdf` |
| Fig. 7 (sensitivity of the GRR index) | `experiments/fig_sens_v2.py` -> `fig5_sens.pdf` |
| Table II (Tier B) | `experiments/make_quality_table_ad.py` -> `results/tables_ad/quality.tex` |
| Table III (Tier C), Table IV (ablation) | `experiments/make_tables_ad.py` -> `results/tables_ad/{systems,ablation}.tex` |
| Need-gate threshold paragraph (Sec. V-E, development data) | `experiments/ad_theta_table.py` -> `results/tables_ad/theta.tex`, `results/ad_dev/theta_table_log.txt` |
| Tier A gate behaviour (Sec. V-C: 29% served, -0.26 points on served cases) | `scripts/ad_dev/tierA_gate_diag.sh` -> `results/ad_dev/tierA_gate_diag.log`; `experiments/ad_tierA_gate.py` -> `results/ad_dev/tierA_gate_analysis.txt` |
| Section V numbers | `experiments/headline_ad.py`, `analyze_oracle_ad.py`, `ad_large_text.py`, `analyze_alloc_v2.py`, `analyze_cap_v2.py`; all claims: `AUDIT.md` |
| Inserting the generated tables into the paper | `experiments/insert_tables.py <main.tex> results` |
| Everything above | `scripts/make_figures.sh <paper_dir>` |

## Citation and license

Please cite the paper if you use this code (`CITATION.cff`):

```bibtex
@inproceedings{veisi2027riskindex,
  author    = {Farshid Veisi and Kamala Berahmand and Sona Taheri and Mahnaz Moradi and Mahdi Jalili},
  title     = {{RiskIndex}: Fast Tail-Risk Influence Maximization via Reusable Graph Indexing},
  booktitle = {Submitted to the 43rd IEEE International Conference on Data Engineering (ICDE)},
  year      = {2027}
}
```

The code is released under the MIT License (`LICENSE`). The SNAP graphs are not redistributed; `scripts/download_data.sh`
fetches them from https://snap.stanford.edu/data.
