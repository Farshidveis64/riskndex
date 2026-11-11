#!/bin/bash
# LargeGraph development: RiskIndex-LG vs RiskIndex-HV vs RIS on the Tier B graphs with development seeds.
# Usage: run_lg_dev.sh small|large
set -e
cd "$(dirname "$0")/../.."
B=${RIBIN:-bin/riskindex}; D=datasets/raw; PART=$1; OUT=results/lg_dev/lgdev_$PART.csv; SYS=results/lg_dev/lgdev_${PART}_sys.csv
rm -f $OUT $SYS $OUT.hcert $OUT.cert
V2="hyb=1 val=1 hvz=0 rhomin=0.5 sq=1 rawidx=0 direct=0 sg=0 sgc=0 inddirect=0 nomat=0 fl=1"
run() { echo "=== $1 $4"; $B quality graph=snap:$2:$3 name=$1 model=$4 p=$5 M=1024 brule=100 bmax=65536 H=8 Kmax=100 ks=1,2,3,5,10 alphas=0.05,0.1,0.2,0.5,1.0 seeds=$6 R=$7 $V2 out=$OUT sysout=$SYS tag=dev 2>&1 | grep -v "auto_b\|^\[hyb\]\|^\[val\]"; }
if [ $PART = small ]; then
for g in "email email-Eu-core.txt directed 0.01622" "wiki wiki-Vote.txt directed 0.02215" "facebook facebook_combined.txt undirected 0.00616" "hepph ca-HepPh.txt undirected 0.00408"; do
  set -- $g; for mo in WC LT IC; do p=0; [ $mo = IC ] && p=$4; run $1 $D/$2 $3 $mo $p 100,101,102 3000; done; done
else
for g in "epinions soc-Epinions1.txt directed 0.00939" "slashdot soc-Slashdot0902.txt directed 0.00808" "dblp com-dblp.ungraph.txt undirected 0.00863"; do
  set -- $g; for mo in WC LT IC; do p=0; [ $mo = IC ] && p=$4; run $1 $D/$2 $3 $mo $p 100 2000; done; done
fi
echo LGDEV_${PART}_DONE
