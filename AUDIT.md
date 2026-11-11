# Audit: every claim in the paper -> raw results -> code

All numbers are recomputed from the raw outputs by `scripts/ad_analyze.sh`. Method labels in the raw files:
`RiskIndex-AD` = the index of the paper (GRR index + forward layer gated by need and cost, `fltheta=5`);
`RiskIndex-LG` = the same layer without the need gate (`fltheta` unset; ablation "w/o need gate", `results/*_lg`);
`RiskIndex-HV` = the GRR index without the forward layer (ablation "w/o forward layer", `results/*_v2`).
Reproduction checks: the RiskIndex-HV, Mean-GRR and RIS rows of the final runs equal those of the ablation runs exactly
(Tier B 3480 rows each, Tier A 6000 rows each, 0 mismatches); k >= 10 answers equal RiskIndex-HV and RiskIndex-LG
(348/348 ties at alpha = 0.05).

Design decisions were frozen before the reported runs; the development logs are `scripts/dev/README.md` (RIS chain,
validation sketch, selection rule), `scripts/lg_dev/README.md` (forward layer, cost gate) and `scripts/ad_dev/README.md`
(need gate; pre-registered rule for theta; theta = 5 adopted, 88% < 90% of the layer's development gain kept, reported
as a limitation).

| # | Claim (Section) | Value | Source / script |
|---|---|---|---|
| 1 | Tier B, k<=5, alpha=.05: RiskIndex vs RIS | +13.6% (108/215/25 of 348, p=1.2e-15) | `results/quality_ad/headline.txt` |
| 2 | vs GRR mean chain / SG / Direct-CVaR | +17.4% (p=3.4e-19) / +11.1% (93/214/41, p=1.5e-6) / +3.7% (80/225/43, p=0.038) | same |
| 3 | w/o forward layer vs Direct-CVaR | +1.6%, p=0.63 (w/o forward layer) | `results/quality_v2/headline_v2.txt` |
| 4 | layer serves 17% of k<=5 queries (ungated: 85%); ungated +0.5 points more (p=0.07) | 0.170 / 0.848; AD +13.63 vs LG +14.11, paired -0.32%, p=0.068 | `headline.txt`, `results/quality_ad/certificates.csv` (`exact`) |
| 5 | dispersion quartiles (Fig. 4d) | Q1 +25.1; Q2-Q4 +0.3, +0.1, -0.0 | `results/quality_ad/fig_quality_log.txt` |
| 6 | four small / three large graphs (k<=5, alpha=.05) | +19.0% / +1.7% (LG +19.6 / +2.0; HV +14.1 / +1.1) | `headline.txt` |
| 7 | ca-HepPh WC k=1; ratios | 20.2 vs RIS 5.3 (3.8x), Mean 5.0, SG 6.7, Direct 17.0; LT: 3.85x RIS, 4.35x SG, 2.09x Direct | `headline.txt` |
| 8 | k>=10 | -0.03% vs RIS (53/212/83, p=0.021); identical to HV and LG | `headline.txt` |
| 9 | DBLP IC k<=5 (alpha=.05) | +7.2% (HV +1.6%); R_tail(k) <= 2.8 for k <= 5 | `headline.txt`; R_tail from hcert mc_cert (`ad_gate_dev.py` formula) |
| 10 | Table II | generated | `results/tables_ad/quality.tex` (Table III = Tier C `systems.tex`, Table IV = ablation `ablation.tex`) |
| 11 | Tier A (alpha=.05) | AD 99.64; HV 99.72; LG 99.50; Direct 99.47; SG 99.51; RIS 98.03; mean-opt 98.03; mean chain 97.79; Direct (ind.) 97.80 | `results/oracle_ad/analysis.txt` |
| 12 | Tier A other levels, regret, exact optimum | 99.67-99.85%; mean regret 0.028; exact optimum 51.1% | same |
| 13 | Tier A tests | vs Direct 1368/3545/1087 p=3.2e-15; vs SG 1341/3661/998 p=1.7e-23; vs mean-based p<=4.7e-29; vs HV -0.08 points p=5.8e-7; vs LG +0.11 p=4.7e-13 | same |
| 14 | mean-opt != CVaR-opt instances | 69.6%; AD 99.74% vs 97.17%; max regret 1.10 vs 4.37 | same |
| 15 | Ablation (six graphs, alpha=.05) | k<=5: HV +10.9, LG +15.1, AD +14.7; DBLP: +2.8, +5.3, +4.6; k>=10 unchanged | `results/tables_ad/ablation.tex`, inline check |
| 16 | layer build time, Tier B total | AD 424 s vs LG 1698 s (4.0x); per cell: DBLP WC 0 vs 207 s, DBLP LT 0 vs 47 s, DBLP IC 103 vs 101 s | `quality_*_sys.csv` fl lines (both versions) |
| 17 | layer built in 9 of 21 settings, 1-3 budgets except DBLP IC | kf > 0 (mean over seeds) in dblp IC, email IC/LT, facebook IC, hepph IC/LT/WC, slashdot IC, wiki IC | `results/ad_large_text.txt` (Tier B forward layer table) |
| 18 | Certificates (alpha<=.2) | median -0.15/-0.11/-0.03%; 96.5% within 10%; raw median -3.45%; k>=10 +0.2%; k<=5 GRR-served -1.3%, layer-served +0.15%; coverage +-2se_w 79.7% (k>=10 85.7%), between-world only 62.2%, in-sample 38.8%; median se_w/cert 1.28% | `headline.txt`, `certificates.csv` |
| 19 | V-D headroom (seed 0; k<=5, alpha<=.2) | small: AD 12.19, HV 9.85, LG 12.89, refs 11.97/12.64; large: AD 3.09, HV 1.99, LG 3.26, refs 3.07/3.08; DBLP IC: refs 19.43/19.48, HV 9.84, AD 19.73 | `results/ad_large_text.txt` |
| 20 | Epinions/Slashdot | best reference per cell <= 1.01%; tail ratio 0.33-0.61; R_tail(1) >= 5.8 | same; `lg_large_text.py` R_tail table |
| 21 | size sweep (seeds 0-2) | HV WC 42.2-52.1, LT 41.5-44.5 -> 32.3 at 1.2M; AD WC 64.4-67.7, LT 63.5-64.9; layer 2.6-10.6 s up to 240K, 35.8-37.4 s at 1.2M | same; `results/sweep_ad/` |
| 22 | Tier C | layer built only on DBLP (WC: kf=1, 108.5 s mean; LT: 24.9 s) and email LT (0.2 s); R_tail(1) on DBLP WC = 16384*63.2/317080 = 3.3; forward lookup P50 0.12-0.21 us, 0.3 MB | `results/scale_ad/scale_hyb.csv`, `results/tables_ad/scale_agg_ad.csv` |
| 23 | break-even with the layer | real graphs except YouTube: 1.2-15.6 (index), 1.5-21.9 (+sketch); YouTube 67.3/99.2; BA 3.8-60.3 (+sketch); DBLP 10.3/11.0 | `scale_agg_ad.csv`, Table IV |
| 24 | Fig. 1(b) | DBLP build 150.9 s (10.3 + 32.1 + 108.5), Q* = 11.0, 91x at 1000 queries | `results/ad_fig_intro_log.txt` |
| 25 | V-G: DBLP IC memory caps 2^16 -> 2^17 raise the GRR index's k<=5 gain only from 1.6 to 1.7% | cap65536 1.60, cap131072 1.72 (alpha=.05) | `experiments/analyze_cap_v2.py` |
| 26 | Step 3b: rho_hat >= 0.8 in every Tier B setting at k<=5 | min 0.80 (DBLP IC) | `results/quality_v2/certificates.csv` (RiskIndex-HV rows) |


| 27 | Tier A gate behaviour (V-B): layer serves 29% of cases (88% of IC p=0.1), max R_tail(k) of served cases < 5 at b/n = 4/3; on served cases AD is 0.26 points below the GRR index (57/203/86), elsewhere identical; overall -0.08 | 28.8%; 0.884; -0.256 (alpha=.05); unserved diff 0.000; overall -0.074 | `scripts/ad_dev/tierA_gate_diag.sh` (reproduces `oracle_ad.csv` exactly: 54000 rows, 0 mismatches), `experiments/ad_tierA_gate.py` -> `results/ad_dev/tierA_gate_analysis.txt` |
| 28 | Need-gate threshold (V-C, development data): theta = 5/10/20/50 serve 19/33/49/83% of Tier B dev cells, gain +11.2/+11.2/+11.4/+11.5% (no layer +9.0, ungated +11.5); Tier A dev 99.62/99.52/99.47/99.46% (no layer 99.65, ungated 99.46) | as stated | `experiments/ad_theta_table.py` -> `results/tables_ad/theta.tex`; `results/ad_dev/oracle_dev_theta{5,10,20,50,1e300}.csv`, `results/lg_dev/lgdev_*.csv(.hcert)` |
| 29 | Online cost O(H) for GRR-served and O(H_F) for forward-served budgets, H_F <= H+1+|A|+T = 19 | stored candidates per served budget: H = 8 index chains (RIS included) + 1 + |A| = 3 + T = 7 forward chains = 19 exactly (identical sets are evaluated once during chain growth) | `ForwardLayer::build` (`src/forward/forward.hpp`), `ForwardLayer::query` |
| 30 | R_tail is computed from the variance-matched RIS-chain profile stored on R_V | `tail_ratios` reads `chains[mean_h()].prefix_sum` (sorted variance-matched profile, Prop. IV.1) and applies `cvar_from_prefix(., 0.05)` | `src/main.cpp`, `src/index/index.hpp` |
| 31 | query(), query_score() and budget_query() use one selection (argmax, noise margin, rho-guard) within the GRR directory | self-test: 240 agreement checks, 0 failures; reported numbers unchanged (all experiments call `query()`) | `RiskIndex::select`, `cmd_selftest` in `src/main.cpp` |

| 32 | Routed interfaces: k_f' = length of the longest prefix of budgets passing both gates (loop stops at the first failure); a budget k <= k_f' is served from the forward directory (W worlds, rho_hat = 1, no rho-guard), others from the GRR directory (M worlds, guard); query, query_score and budget query share the routing; budget query O(Kmax max(H, H_F)) | self-test: 160 routed checks (30 forward-served), 0 failures; routing fixed at build (`ForwardLayer::kf`) | `adaptive_query`, `adaptive_query_score`, `adaptive_budget_query` in `src/forward/forward.hpp`; `cmd_selftest` |
| 33 | theta rule with fallback (smallest theta within 0.05 Tier A points that keeps 90% of the Tier B gain; otherwise the smallest theta meeting the Tier A criterion, shortfall reported); theta = 5 by the fallback, 88% | rule and decision as written before the Tier A development analysis | `scripts/ad_dev/README.md` (sections "PRE-REGISTERED RULE" and "DECISION") |
| 34 | Forward layer memory: 16 H_F k_f' (W+1) bytes, at most 1.5 MB in our runs | 19 x 5 x 1001 x 16 = 1.52 MB; measured 0.3-1.5 MB (Table III column) | `ForwardLayer::bytes`, `results/scale_ad/*_sys.csv` |

| 35 | Tier B k <= 5 gain is concentrated: median 0, ca-HepPh +67%, +2.5% without it (70/200/18 of 288, p = 2e-9); per-graph gains Wiki +2.9, Facebook +6.1, DBLP +4.6 (significant), none on Epinions/Slashdot, email -0.3 | recomputed from `results/quality_ad/merged.csv` (alpha = 0.05, k <= 5, paired on graph/model/seed/k) | Sec. V-B |
| 36 | Margins over Direct-CVaR and SG at k <= 5 rest on ca-HepPh: +0.6% (p = 0.82) and +0.4% (p = 0.11) without it; Direct-CVaR leads in 7 of 21 rows of Table II | same file; Table II row counts | Sec. V-B |
| 37 | Forward layer changes answers only on ca-HepPh (+8.9% over RiskIndex-HV) and DBLP (+1.8%); 1 of 252 cases elsewhere | same file, RiskIndex-AD vs RiskIndex-HV | Sec. V-E |
| 38 | Matched-grid control for "w/o tail selection": SG rule on the index's grid gains +10.4% at k <= 5 and loses 12.1% at k >= 10 (index: +10.9%, 0.0%) | `scripts/run_quality_sgg.sh`, `experiments/sgg_control.py`, `results/quality_sgg/analysis.txt`; reproduction of RIS, mean chain, RiskIndex-HV: 0 mismatches | Sec. V-E |
| 39 | Size sweep: GRR index 42-52% (WC), 32-45% (LT), seed s.d. up to 14 points; RiskIndex 63-68%, s.d. below 5 | `results/ad_large_text.txt` | Sec. V-D, Fig. 5b |
| 40 | Index size 13,136,000 bytes (M = 1024) and 6,582,400 bytes (M = 512) | `*_sys.csv` byte columns | Sec. IV-A (index size), Sec. V-F |

Layout: body ends on page 12; acknowledgment (AI-generated content disclosure) and references follow, 14 pages in total; fonts embedded, no overfull boxes; author block present (ICDE 2027 is single-blind).