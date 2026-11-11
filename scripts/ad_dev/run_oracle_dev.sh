#!/bin/bash
# Need-gate development, Tier A runs: sketch seeds 100-102 (disjoint from the reported seeds 0-4), one run per
# candidate threshold theta of the need gate; RiskIndex-HV (no layer) and RiskIndex-LG (gated layer) per run.
# Usage: run_oracle_dev.sh <theta>   -> results/ad_dev/oracle_dev_theta<theta>.csv
set -e
cd "$(dirname "$0")/../.."
TH=$1; B=bin/riskindex; OUT=results/ad_dev/oracle_dev_theta$TH.csv; mkdir -p results/ad_dev; rm -f $OUT
run() { $B oracle graph=$1 model=$2 p=$3 k=$4 alphas=0.05,0.1,0.2,0.5,1.0 M=4096 b=64 H=8 Kmax=$4 seeds=100,101,102 RA=5000 RB=5000 lite=1 hyb=1 val=1 rhomin=0.5 fl=1 fltheta=$TH out=$OUT tag=$5; }
for inst in 1 2 3 4 5; do
  for spec in "er:48:4:$inst" "er:48:8:$inst" "ba:48:2:$inst" "ba:48:4:$inst" "sbm:48:4:0.5:0.02:$inst" "sbm:48:4:0.3:0.08:$inst" "ws:48:4:0.2:$inst" "ws:48:6:0.5:$inst" "er:64:4:$inst" "ba:64:2:$inst" "sbm:64:4:0.4:0.02:$inst"; do
    for p in 0.1 0.2; do run $spec IC $p 3 "IC$p"; done
    run $spec WC 0 3 WC; run $spec LT 0 3 LT
  done
  for spec in "er:48:4:$inst" "ba:48:2:$inst" "sbm:48:4:0.5:0.02:$inst" "ws:48:4:0.2:$inst"; do run $spec IC 0.2 4 IC0.2; done
done
echo ORACLE_DEV_${TH}_DONE
