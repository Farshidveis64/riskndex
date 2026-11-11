#!/bin/bash
# LargeGraph diagnosis (development seed 100 only; evaluation streams 555000+seed / 777000+seed are disjoint from
# the reported runs). Usage: run_headroom.sh small|large
set -e
cd "$(dirname "$0")/../.."
B=${RIBIN:-bin/riskindex}; D=datasets/raw; mkdir -p results/lg_dev
PART=$1; OUT=results/lg_dev/headroom_$PART.csv; rm -f $OUT
run() { $B headroom graph=snap:$2:$3 name=$1 model=$4 p=$5 M=1024 brule=100 bmax=65536 seed=100 P=60 Pscreen=2000 Wsel=1000 Weval=3000 out=$OUT 2>&1 | grep -v auto_b; }
if [ $PART = small ]; then
for g in "email email-Eu-core.txt directed 0.01622" "wiki wiki-Vote.txt directed 0.02215" "facebook facebook_combined.txt undirected 0.00616" "hepph ca-HepPh.txt undirected 0.00408"; do
  set -- $g; for mo in WC LT IC; do p=0; [ $mo = IC ] && p=$4; run $1 $D/$2 $3 $mo $p; done; done
else
for g in "epinions soc-Epinions1.txt directed 0.00939" "slashdot soc-Slashdot0902.txt directed 0.00808" "dblp com-dblp.ungraph.txt undirected 0.00863"; do
  set -- $g; for mo in WC LT IC; do p=0; [ $mo = IC ] && p=$4; run $1 $D/$2 $3 $mo $p; done; done
fi
echo HEADROOM_${PART}_DONE
