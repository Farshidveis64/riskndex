#!/bin/bash
# Forward-simulated references: tail headroom at k <= 5 on the Tier B graphs (seed 0; selection worlds 555000+seed,
# evaluation worlds 777000+seed).  RIS, RiskIndex-HV, RiskIndex-LG and the forward-simulated pool references
# (Oracle-CVaR, Oracle-basis) on the same 3000 evaluation worlds.
set -e
cd "$(dirname "$0")/.."
B=${RIBIN:-bin/riskindex}; D=datasets/raw; mkdir -p results/headroom_lg; OUT=results/headroom_lg/headroom.csv; rm -f $OUT
run() { echo "=== $1 $4"; $B headroom graph=snap:$2:$3 name=$1 model=$4 p=$5 M=1024 brule=100 bmax=65536 seed=${SEED:-0} P=60 Pscreen=2000 Wsel=1000 Weval=3000 fl=1 out=$OUT 2>&1 | grep -v auto_b; }
for g in "email email-Eu-core.txt directed 0.01622" "wiki wiki-Vote.txt directed 0.02215" "facebook facebook_combined.txt undirected 0.00616" "hepph ca-HepPh.txt undirected 0.00408" "epinions soc-Epinions1.txt directed 0.00939" "slashdot soc-Slashdot0902.txt directed 0.00808" "dblp com-dblp.ungraph.txt undirected 0.00863"; do
  set -- $g; for mo in WC LT IC; do p=0; [ $mo = IC ] && p=$4; run $1 $D/$2 $3 $mo $p; done; done
echo HEADROOM_LG_DONE
