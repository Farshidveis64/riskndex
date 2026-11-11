# Need gate of the forward layer: development log (stage 3 of 3)

Goal: keep the small-budget gains of the (ungated) forward-simulation layer where the GRR sketch cannot
resolve the tail, and pay its cost (build time; 0.2 points in Tier A) only there.  Diagnosis (on the reported runs of
the ungated-layer runs, used as a diagnostic only): rho_hat of the selected entries is >= 0.80 in every Tier B cell at
k <= 5, so rho_hat cannot gate the layer; the tail-resolution ratio

    R_tail(k) = b * C_hat_0.05(RIS k-prefix on the validation sketch) / n

(the number of GRR resolution units n/b in the mean of the worst 5% scenarios of the RIS answer) separates the cells
where the layer helps (R_tail < 10: 165 of the 173 gain points) from those where it does not (mean +0.13 points).

## Gate (src/forward/forward.hpp, `fltheta`)
The layer is built for the leading budgets k <= kf = 5 that pass both gates:
  need gate  R_tail(k) < theta        (new)
  cost gate  sigma_hat(RIS k-prefix) <= sigma_max = 1000   (as in the ungated layer, unchanged)
Everything else (screening, chains, W, validation, selection) is unchanged.

## PRE-REGISTERED RULE for theta (written before the Tier A development runs were analysed)
Candidates: theta in {5, 10, 20, 50}.  Development data: Tier B dev seeds 100-102 (results/lg_dev; the gated answer is
LG where R_tail(k) < theta, HV otherwise, an approximation of the rebuilt layer) and Tier A dev seeds 100-102
(scripts/ad_dev/run_oracle_dev.sh, exact gated layer).
Adopt the smallest theta such that
 (i)  Tier A dev: the gated index is within 0.05 points of RiskIndex-HV in mean % of the optimum at alpha = 0.05, and
 (ii) Tier B dev: the gated index keeps at least 90% of the ungated layer's gain over RiskIndex-HV, pooled over all
      seven graphs (k <= 5, alpha <= 0.2, gain over RIS).
If no candidate satisfies both, adopt the smallest theta satisfying (i) and report (ii) as a limitation; if none
satisfies (i), the need gate is reported as a negative result and the ungated layer stands.
After the decision, theta is frozen and all reported runs (seeds 0-4 / 0-2, identical sketches and evaluation worlds as
the ungated-layer and no-layer ablations) are run once.

## Development results and decision
Tier B dev (results/lg_dev, seeds 100-102; gain over RIS, k <= 5, alpha <= 0.2, pooled over the seven graphs; HV 9.03%,
LG 11.47%): gated 11.18% (theta 5, 88% of the layer's gain kept), 11.23% (10, 90%), 11.37% (20, 96%), 11.45% (50, 99%).
Fraction of (graph, model, seed, k <= 5) cells with the layer: 19%, 33%, 49%, 83%.
Tier A dev (seeds 100-102, 720 cases per alpha; mean % of the optimum at alpha = 0.05; HV 99.65): gated 99.63 (theta 5,
-0.03), 99.52 (10, -0.13), 99.47 (20, -0.18); see results/ad_dev/ for theta 50.
No candidate satisfies both (i) and (ii): theta = 5 satisfies (i) but keeps 88% < 90% of the layer's gain; theta = 10
keeps 90% but loses 0.13 points in Tier A.  Per the rule, theta = 5 is ADOPTED and (ii) is reported as a limitation:
the gate gives up about 12% of the layer's development gain (the cells with 5 <= R_tail < 10, e.g. ca-HepPh WC k = 2
and Wiki-Vote IC k = 1) to keep the layer off where the sketch resolves the tail.  theta is frozen; the reported runs
(scripts/ad_runner.py, THETA=5) follow.
