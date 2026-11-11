Development runs used to choose the frozen configuration (sketch seeds 100-102, evaluation streams disjoint
from the reported seeds 0-4; email, Wiki-Vote, Facebook, ca-HepPh only).
  run_variants.sh   A: M=512 rule-b (+ margin rules z=0.5,1,2) | B: M=1024 b/2 | C: M=1024 b/2 split 0.5 | D: M=2048 b/4
  run_variants2.sh  E: M=1024 rule-b (adopted) | F: M=2048 b/2
Analysis: experiments/an.py (per-method vs mean chain), experiments/cross.py (fixed-reference cross-variant),
experiments/cert.py (certificate calibration, SE coverage).
Outcome (alpha<=0.2, 1440 common cases): k<=5 gain A 14.1 / B 12.9 / D 11.2 / E 16.0 / F 14.1 %; k>=10 -0.28 / +0.26 / +0.72 / +0.40 / +0.84 %.
The noise-margin selection rule changed results by <0.2 points and was dropped; cross-fitting (C) calibrated k>=10
certificates but hurt selection and was dropped.
  run_underres.sh   rho-guard / margin ablation on artificially under-resolved sketches (b=128, 512) and rule-b,
                    wiki/facebook/hepph x IC/LT/WC, dev seeds 100-102.  alpha<=0.2, gain vs mean chain (%):
                      b=128:  k<=5  0.94 / rho0.3 0.96 / rho0.5 0.81 / rho0.7 0.62 ; k>=10 -1.19 / -0.25 / -0.18 / -0.10
                      b=512:  k<=5 11.91 / 11.91 / 11.91 / 10.61            ; k>=10 -0.73 / -0.39 / -0.22 / -0.17
                      rule:   k<=5 18.58 (all equal)                         ; k>=10 -0.25 (all equal)
                    -> rho_min = 0.5 adopted (removes most under-resolution loss at k>=10, <0.05 points at k<=5, inert at rule-b).
                    Noise-margin z=1,2: no improvement at k>=10, 0.5-1 point loss at k<=5 -> dropped.

## Hybrid mean chain and selection margin (run_hybrid_dev.sh; RIS chain and validation sketch)

Variants (all on the same GRR sketch, rho-guard rho_min = 0.5 throughout):
- RiskIndex: the submitted index (tau = inf chain = lazy greedy on the GRR sketch).
- RiskIndex-H: tau = inf chain = RIS greedy on an independent RR sketch of the same size (Mb sets), materialized
  on the GRR sketch; plain argmax selection.
- RiskIndex-H-z<z>: RiskIndex-H with the noise margin: a threshold chain replaces the RIS chain only if its
  certificate exceeds the RIS chain's by more than z * sqrt(se_w(chain)^2 + se_w(RIS)^2), where se_w includes the
  within-world noise of the tail scenarios (tail_se_within); z in {0.5, 1, 2, 3}.
- RiskIndex-H-zb<z>: the same margin with the between-world standard error only (z in {1, 2}).

Sketches: dev seeds 100-102 (evaluation stream 31337+seed, R = 3000), email, Wiki-Vote, Facebook, ca-HepPh,
WC/LT/IC-critical; deployed sketch (M = 1024, rule b, cap 2^16) and under-resolved sketches (b = 512, b = 128).

Decision rule, fixed before the runs were inspected: the reference is RIS (ind. RR). Among z in {0, 0.5, 1, 2, 3}
(z = 0 is RiskIndex-H), keep the variants whose mean gain in realized CVaR (alpha <= 0.2) at k >= 10 is at least
-0.1% on the deployed and on the b = 512 sketches; among these choose the largest mean gain at k <= 5 on the deployed
sketches; differences below 0.1 points are ties, broken towards the smaller z. b = 128 is reported as a stress test.
The chosen rule is then frozen for all reported runs (seeds 0-4).

Outcome, first round (results/dev/hyb_*.csv, experiments/hyb_dev.py; gain over RIS ind., alpha <= 0.2):
no RiskIndex-H variant met the k >= 10 criterion (deployed sketch: plain -0.69%, z=1 -0.53%, z=3 -0.17%;
b=512: -0.60 / -0.41 / -0.13%).  Cause: the threshold chains are evaluated on the worlds that selected them
(optimistic, in-sample certificates), the RIS chain on worlds that did not, so the argmax prefers the
in-sample chains at large k, where the RIS chain is truly better (price of grouping).

Second round (run_hv_dev.sh): all chains re-materialized on an independent validation GRR sketch
(RiskIndex-HV), so that every certificate is out-of-sample; RiskIndex-V = validation with the GRR mean chain.
Same decision rule applied to {HV, HV-z0.5, HV-z1, HV-z2} (results/dev/hv_summary.txt):
  deployed sketch  k<=5: RiskIndex 11.14 | HV 10.97 | z0.5 10.87 | z1 10.78 | z2 9.85
                   k>=10: RiskIndex -0.77 | HV -0.05 | z0.5 -0.01 | z1 +0.02 | z2 +0.02
  b = 512          k<=5: 7.35 | 8.68 | 8.39 | 7.96 | 6.77 ;  k>=10: -0.72 | +0.05 | +0.07 | +0.08 | +0.06
  b = 128          k<=5: -0.50 | +0.06 | +0.09 | +0.05 | +0.03 ;  k>=10: -0.79 | -0.03 | 0.00 | 0.00 | 0.00
All four HV variants are feasible; the largest k<=5 gain is HV (z = 0).  DECISION (frozen for all reported runs):
RIS chain for tau = inf + validation sketch, plain argmax with the rho-guard rho_min = 0.5 (no margin).

rho-guard of the final index (run_guard_v2_dev.sh; experiments/analyze_guard_v2.py), alpha <= 0.2, gain over RIS:
  b = 128 (rho 0.57): rho_min 0: k<=5 +0.14 / k>=10 -0.43 ; rho_min 0.5: +0.06 / -0.03 (guard fires on 28.5%)
  b = 512 (rho 0.76): rho_min 0: +8.68 / +0.01 ; rho_min 0.5: +8.68 / +0.05 (fires on 6.0%)
  -> the frozen rho_min = 0.5 is kept; it never fires on the reported Tier B sketches.
