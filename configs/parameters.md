# Parameter configurations used in the paper

| Experiment | Script | M | b | H | Kmax | k grid | alpha grid | seeds | eval worlds R |
|---|---|---|---|---|---|---|---|---|---|
| Tier A oracle (n<=64) | run_oracle.sh (baselines), run_oracle_v2.sh, ad_runner.py (index) | 4096 | 64 | 8 | k | 3 (4 on 4 families, IC 0.2) | .05,.1,.2,.5,1 | 5 | 5000 select / 5000 eval |
| RQ1 fidelity, b sweep | run_fidelity.sh | 256 | 1..4096 | -- | -- | random sets, k=5 (10 Epinions) | .05,.1,.2 | 5 (3 Epinions) | 20000 (5000 Epinions) |
| RQ1 fidelity, M sweep | run_fidelity.sh | 64..2048 | 1024 | -- | -- | k=5 | .05,.1,.2 | 5 | 20000 |
| Tier B quality | run_quality.sh, run_quality_part2.sh (baselines), run_quality_v2.sh, ad_runner.py (index) | 1024 | rule c=100, memory cap 2^16 | 8 | 100 | 1,2,3,5,10,25,50,100 | .05,.1,.2,.5,1 | 5 (3 for >70K nodes) | 5000 / 3000 (Epinions, Slashdot) / 2000 (DBLP) |
| Budget allocation (fixed M*b) | run_alloc_v2.sh | 512 / 1024 / 2048 | rule c=100 / c=50 / c=25 | 8 | 100 | 1,2,3,5,10,25,50,100 | same | 3 | 5000 |
| Sensitivity H | run_sensitivity.sh | 512 | rule | 1,2,4,8,16,32 | 100 | per graph | same | 5 (3 Epinions) | 5000 / 3000 |
| Sensitivity M | run_sensitivity.sh | 64..1024 | rule | 8 | 100 | per graph | same | 5 (3) | 5000 / 3000 |
| Sensitivity b | run_sensitivity.sh | 512 | 64..16384 | 8 | 100 | per graph | same | 5 (3) | 5000 / 3000 |
| Price of grouping | run_sensitivity.sh | 128..4096 | rule | 8 | 100 | per graph | same | 5 (3) | 5000 / 3000 |
| Kmax | run_sensitivity.sh | 512 | rule | 8 | 25..400 | 10 | .1 | 3 | 100 |
| Tier C scale | run_scale_v2.sh, run_scale_v2_lj.sh, ad_runner.py | 512 | rule, cap 2^14 | 8 | 100 | 10,25,50,100 (LJ: 10,50,100) | .05,.1,.2,.5,1 (LJ: .05,.2,1) | 3 (LJ: 1) | 2000 (LJ: 1000) |
| Matched-grid SG control (ablation) | run_quality_sgg.sh | 1024 | rule, cap 2^16 | 8 | 100 | as Tier B | as Tier B | as Tier B, six graphs (no DBLP) | as Tier B |

Diffusion: WC (p_uv = 1/indeg v), LT (w_uv = 1/indeg v), IC with p = 1/lambda_1(A) per graph
(email 0.01622, wiki 0.02215, facebook 0.00616, hepph 0.00408, epinions 0.00939, slashdot 0.00808, dblp 0.00863);
Tier A IC p in {0.1, 0.2}.
Resolution rule: b = next power of two >= 100 * n / sigma_pilot, sigma_pilot = mean spread of the 10 highest
out-degree nodes over 100 cascades on a separate random stream; clamped to [64, 2^16] in Tier B (the cap binds only on DBLP IC) and to [64, 2^14] in Tier C.

Index contents: per prefix, prefix sums P_0..P_M and prefix sums of squares Q_0..Q_M of the sorted
variance-matched profile (doubles), the resolution ratio rho_hat, the profile mean, and the seed id:
H*Kmax*(16(M+1)+20) bytes.  Query output: certificate, raw certificate (recovered from rho_hat and the mean),
tail standard error sqrt((Q_r/r-(P_r/r)^2)/r) with r=ceil(alpha*M), rho_hat, chain id, seeds.
Selection rule (shared by query, query_score and budget_query): argmax of the certificate followed by the rho-guard (`rhomin=0.5`: a non-mean entry with rho_hat < 0.5 is replaced by the mean-chain prefix). The noise-margin rule `z=...` is kept as an option; it gave no
measurable improvement on the development seeds 100-102 and is not used in the paper).
Development protocol: all configuration choices (M/b allocation, selection rule, RIS chain, validation sketch, forward layer and its gates) were made on sketch seeds 100-102
with evaluation streams disjoint from those of the reported runs (seeds 0-4), on email, Wiki-Vote, Facebook and ca-HepPh only.
Baseline SG [Zhou & Tokekar]: uniform grid of H'=H thresholds on (0, Gamma], Gamma = largest per-world spread of the
mean-greedy k-set on the sketch; selection by surrogate value; SG-cached builds the H' chains once to Kmax.
Direct-CVaR-VM: valid pruning bound DeltaCVaR_vm <= (mean - CVaR_vm) + DeltaMean (the raw 1/alpha bound is not valid
under variance matching); `directvm_legacy=1` reproduces the earlier (invalid) bound for the self-test.

Self-test: `bin/riskindex selftest [inst=6 M=96 b=24]` checks that lazy Direct-CVaR greedy with the valid bounds (raw and
variance-matched objective) selects the same objective value as exhaustive greedy without pruning on small synthetic graphs,
that prefix-sum CVaR equals the direct computation, and that the raw certificate recovered from (rho, mean, C_vm) equals the
raw empirical CVaR (and is flagged invalid when rho = 0).
Tier C (run_scale_v2.sh, ad_runner.py) keeps cap 2^14 at M=512.
Random streams: the build sketch R_B and the RR sketch R_I are built from the same sketch seed, so they share the sequence of
target nodes; their worlds differ (one world per group in R_B, one fresh world per RR set in R_I).  The validation sketch,
the forward selection and validation worlds and the evaluation worlds use separate streams.

Forward layer (`fl=1`; src/forward/forward.hpp; Step 3b of the paper). All defaults were fixed before the reported runs
(cost gate and layer parameters: scripts/lg_dev/README.md; need gate and theta: pre-registered rule in scripts/ad_dev/README.md):
  kf = 5 (`flk`); a budget k <= kf is served by the layer only if the leading budgets 1..k all pass both gates:
    need gate  R_tail(k) = b * C_hat_0.05(RIS k-prefix on the validation sketch) / n < theta = 5 (`fltheta`;
               unset = no need gate, the ablation "w/o need gate", label RiskIndex-LG);
    cost gate  RR-sketch spread estimate of the RIS k-prefix <= sigma_max = 1000 (`flsmax`);
  screening: the Ps = 2000 (`flPscreen`) best-covered nodes of the RR sketch, by exact single-seed mean and CVaR_alpha
  (alpha in A = {0.05, 0.1, 0.2}, `flalphas`) on W = 1000 forward selection worlds (`flWsel`, stream seed*7919+505);
  pool = top-P (P = 60, `flP`) per statistic plus the first kf seeds of every index chain;
  forward chains: greedy on the exact mean, the exact empirical CVaR_alpha (alpha in A) and the exact truncated spread
  F_tau for T = 7 (`flT`) thresholds, geometric between the 5% quantile of the RIS 1-prefix and the 20% quantile of the
  RIS kf-prefix on the selection worlds;
  candidates: all forward chains and the index chains (RIS chain included), truncated to kf, identical sets once;
  validation: sorted exact profiles of every k-prefix on W = 1000 fresh worlds (`flWval`, stream seed*7919+707), stored as
  prefix sums and prefix sums of squares; query: argmax of the exact empirical CVaR (ties: RIS chain, then larger mean);
  standard error: between-world only; rho_hat = 1 by convention and no rho-guard for forward-served budgets; candidates per
  served budget H_F = H + 1 + |A| + T_F = 19 (T_F = `flT`); memory 16 H_F k_f' (W+1) bytes; routing (k_f') fixed at build.
  R_tail(k) is read from the variance-matched profile of the RIS k-prefix stored on the validation sketch (Step 3), i.e.
  the same certificate C_hat_0.05 that a query at alpha = 0.05 would return for that prefix.
  theta was chosen from {5, 10, 20, 50} on development seeds 100-102 by the rule in scripts/ad_dev/README.md:
  Tier A dev within 0.05 points of the index without the layer at alpha = 0.05 (theta = 5: -0.03; 10: -0.13; 20: -0.18),
  then keep >= 90% of the ungated layer's Tier B dev gain (theta = 5 keeps 88%, reported as a limitation).
Final runs: THETA=5 scripts/ad_runner.py (Tier B as above with fl=1 fltheta=5; Tier A oracle M=4096 b=64 fl=1 fltheta=5;
Tier C scale M=512 cap 2^14 fl=1 fltheta=5; headroom and size sweep as in scripts/run_headroom_lg.sh, run_sweep_lg.sh).
