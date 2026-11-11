# Forward layer: development log (stage 2 of 3; stage 3 = the need gate, scripts/ad_dev/README.md)

Goal: find out why the small-graph gains of the GRR index (RiskIndex-HV: +14.1% on email, Wiki-Vote, Facebook and ca-HepPh vs.
+1.1% on Epinions, Slashdot and DBLP at k <= 5) do not carry over to large graphs, and keep them where they exist.
All development runs use sketch seed 100 (101, 102 where stated) and world streams disjoint from the reported
runs: selection 555000+seed, evaluation 777000+seed (headroom), forward layer seed*7919+505 / +707, Tier B
evaluation stream 31337+seed with dev seeds only.

## 1. Headroom diagnosis (run_headroom.sh, run_headroom_youtube.sh; experiments/lg_headroom.py)

`riskindex headroom` builds RiskIndex-HV exactly as in Tier B and compares it on 3000 fresh worlds with two
forward-simulated references over a screened pool: the 2000 best-covered nodes are screened by their exact
single-seed mean and CVaR on 1000 selection worlds, and the pool is the top-60 by each statistic, the RIS prefix and
the index's own seeds. Oracle-CVaR is a greedy on the exact empirical CVaR. Oracle-basis grows exact truncated-spread
chains and picks the best candidate, the RIS prefix and RiskIndex-HV included, by exact CVaR on the selection worlds.
A first version that used only the top-100 nodes by mean as the pool missed the CVaR-optimal seeds (it lost to the
index on ca-HepPh LT). It was replaced before any large-graph conclusion was drawn.

Outcome (seed 100; gain over RIS, k in {1,2,3,5}, alpha <= 0.2):
- Small graphs: HV 8.96%, Oracle-CVaR 12.25%, Oracle-basis 12.48%, best 13.03%. The index leaves headroom at k = 1-3
  on ca-HepPh: k = 1 gives +125% for HV vs. +253% for the oracle.
- Epinions, Slashdot and DBLP: the best reference gains <= 1.4% except DBLP LT (4.9%; HV 4.2%). The tail headroom of
  real large social graphs is small. Their top seeds are reliable: tail ratio 0.3-0.6, or, under IC-critical, the
  risk is world-level and no seed choice hedges it.

## 2. Forward layer (src/forward/forward.hpp, `fl=1`; called "exact small-budget layer" in this log)

For k <= kf, exact forward chains (mean, CVaR_alpha, truncated spreads) are grown over the screened pool, and the
profiles of all k-prefixes, the GRR chains and the RIS chain included, are materialized on independent forward
validation worlds. Selection is by exact CVaR. Defaults: kf = 5, Wsel = Wval = 1000, Pscreen = 2000, P = 60, T = 7,
alpha grid {0.05, 0.1, 0.2}.

PRE-REGISTERED DECISION RULE (fixed before the development runs of the layer on all graphs): adopt the layer with
these defaults if, on the development seeds,
 (i)  the mean gain over RIS at k <= 5 (alpha <= 0.2) on the four small graphs exceeds RiskIndex-HV by >= 1 point;
 (ii) on Epinions, Slashdot and DBLP it is not below RiskIndex-HV by more than 0.2 points;
 (iii) on ca-HepPh replicated 20 and 100 times (run_sweep_dev.sh) it keeps at least half of its 1-copy gain.
Otherwise the layer is reported as a negative result. k >= 10 is unchanged by construction (kf = 5).

### Cost cap (decided on a development cost measurement, before any evaluation run)
On YouTube WC (seed 100, Tier C configuration) the layer did not finish within 5 minutes. There the expected
spread of the RIS 1-prefix is ~2*10^4 nodes, and forward simulation costs proportionally to the spreads. The layer
therefore serves only the budgets k <= kf whose RIS prefix has an estimated expected spread (from the RR sketch) of
at most sigma_max = 1000 nodes (`flsmax`). Budgets above the cap use RiskIndex-HV unchanged.
- Applies to: every small graph, every IC-critical setting, and DBLP WC/LT.
- Does not apply to: Epinions and Slashdot under WC/LT (RIS 1-prefix spreads of 1400-3000 nodes) and the Tier C
  social graphs under WC.
- Effect on the development results: Epinions LT drops from +0.64% to HV's +0.31%, and Slashdot WC from +1.20% to
  HV's 0.00%.
- Criteria (i)-(iii) are unaffected.

### Decision (development runs complete)
Development results (gain over RIS, alpha <= 0.2; experiments/lg_dev.py):
- Small graphs (seeds 100-102): k <= 5 HV 10.97% -> LG 13.79% (+2.8 points); k >= 10 unchanged.
- Large graphs (seed 100): k <= 5 HV 1.31% -> LG 2.21%. LG >= HV in every cell, e.g. DBLP IC 6.51% -> 10.63%.
- Size sweep (ca-HepPh padded 1-100x, seed 100): LG keeps 66.6-67.0% (LT) and 74.4-74.9% (WC) at every size. HV
  gets 39.8-48.6% (LT) and 46.2-63.0% (WC), lowest at 1.2M nodes in both models.
Criteria (i)-(iii) are met, so the forward layer is ADOPTED with the pre-registered defaults and the cost cap
sigma_max = 1000. At this stage the index with the (ungated) layer is RiskIndex-LG, which the paper reports as the ablation "w/o need gate"; RiskIndex-HV appears as
"w/o forward layer". After the decision, the candidate sets shared by several objectives in a greedy pass are
evaluated once (identical results, verified on Wiki-Vote IC; about 1.6x faster).
