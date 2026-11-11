#!/bin/bash
# Tier A, index without the forward layer (RiskIndex-HV: RIS chain + validation sketch + rho-guard 0.5) and the
# submitted index with the rho-guard (RiskIndex-G), on the identical sketches and evaluation worlds as
# run_oracle.sh.  lite=1 skips the exhaustive search and the baselines; the exhaustive optimum and all baselines
# are joined from results/oracle/oracle.csv (same selection worlds RA, seed 101, and evaluation worlds RB, seed 202).
set -e
cd "$(dirname "$0")/.."
B=${RIBIN:-bin/riskindex}; OUT=results/oracle_v2/oracle_v2.csv; mkdir -p results/oracle_v2; rm -f $OUT
SEEDS=0,1,2,3,4
run() { $B oracle graph=$1 model=$2 p=$3 k=$4 alphas=0.05,0.1,0.2,0.5,1.0 M=4096 b=64 H=8 Kmax=$4 seeds=$SEEDS RA=5000 RB=5000 lite=1 hyb=1 val=1 rhomin=0.5 out=$OUT tag=$5; }
for inst in 1 2 3 4 5; do
  for spec in "er:48:4:$inst" "er:48:8:$inst" "ba:48:2:$inst" "ba:48:4:$inst" "sbm:48:4:0.5:0.02:$inst" "sbm:48:4:0.3:0.08:$inst" "ws:48:4:0.2:$inst" "ws:48:6:0.5:$inst" "er:64:4:$inst" "ba:64:2:$inst" "sbm:64:4:0.4:0.02:$inst"; do
    for p in 0.1 0.2; do run $spec IC $p 3 "IC$p"; done
    run $spec WC 0 3 WC
    run $spec LT 0 3 LT
  done
  for spec in "er:48:4:$inst" "ba:48:2:$inst" "sbm:48:4:0.5:0.02:$inst" "ws:48:4:0.2:$inst"; do
    run $spec IC 0.2 4 IC0.2
  done
done
echo ORACLE_V2_DONE
