#!/bin/bash
# Analysis only (no parameter choice): rerun the Tier A configuration of the paper with the gate diagnostic enabled and
# record, per instance and seed, the tail-resolution ratios R_tail(k) of the RIS prefixes and the budgets served by the
# forward layer (stderr lines "[flgate] ...").  Answers must reproduce results/oracle_ad/oracle_ad.csv exactly.
set -e
cd "$(dirname "$0")/../.."
B=${RIBIN:-bin/riskindex}; OUT=results/ad_dev/tierA_gate_diag; mkdir -p results/ad_dev; rm -f $OUT.csv $OUT.log
run() { $B oracle graph=$1 model=$2 p=$3 k=$4 alphas=0.05,0.1,0.2,0.5,1.0 M=4096 b=64 H=8 Kmax=$4 seeds=0,1,2,3,4 RA=5000 RB=5000 lite=1 hyb=1 val=1 rhomin=0.5 fl=1 fltheta=5 out=$OUT.csv tag=$5 2>&1 >/dev/null | grep "^\[flgate\]" >> $OUT.log || true; }
for inst in 1 2 3 4 5; do
  for spec in er:48:4:$inst er:48:8:$inst ba:48:2:$inst ba:48:4:$inst sbm:48:4:0.5:0.02:$inst sbm:48:4:0.3:0.08:$inst ws:48:4:0.2:$inst ws:48:6:0.5:$inst er:64:4:$inst ba:64:2:$inst sbm:64:4:0.4:0.02:$inst; do
    for p in 0.1 0.2; do run $spec IC $p 3 IC$p; done
    run $spec WC 0 3 WC; run $spec LT 0 3 LT
  done
  for spec in er:48:4:$inst ba:48:2:$inst sbm:48:4:0.5:0.02:$inst ws:48:4:0.2:$inst; do run $spec IC 0.2 4 IC0.2; done
done
echo TIERA_GATE_DIAG_DONE
