#!/bin/bash
# Tier B, index without the forward layer (RiskIndex-HV: RIS chain for tau = inf, all chains materialized on a
# validation GRR sketch, plain argmax with the rho-guard 0.5) on the identical sketches (same seeds) and
# evaluation worlds as the submitted Tier B runs, so all baselines of results/quality/quality_merged.csv stay paired.
# Records: RiskIndex (submitted index; reproduction check), RiskIndex-H (RIS chain, no validation),
# RiskIndex-V (validation, GRR mean chain), RiskIndex-HV (final), RiskIndex-HV-z1 (margin ablation),
# Mean-GRR and RIS-IndRR (reproduction checks).
# Usage: run_quality_v2.sh [small|dblp|all]
set -e
cd "$(dirname "$0")/.."
B=${RIBIN:-bin/riskindex}; D=datasets/raw; mkdir -p results/quality_v2
PART=${1:-all}; OUT=results/quality_v2/quality_$PART.csv; SYS=results/quality_v2/quality_${PART}_sys.csv
rm -f $OUT $SYS $OUT.cert $OUT.hcert
ALPHAS=0.05,0.1,0.2,0.5,1.0; KS=1,2,3,5,10,25,50,100
V2="hyb=1 val=1 hvz=0 hvzs=1 hvraw=1 rhomin=0.5 sq=1 rawidx=0 direct=0 sg=0 sgc=0 inddirect=0 nomat=0"
run() { echo "=== $1 $4"; $B quality graph=snap:$2:$3 name=$1 model=$4 p=$5 M=1024 brule=100 bmax=65536 H=8 Kmax=100 ks=$KS alphas=$ALPHAS seeds=$6 R=$7 $V2 out=$OUT sysout=$SYS tag=tierB 2>&1 | grep -v auto_b; }
if [ $PART != dblp ]; then
for mo in WC LT IC; do p=0; [ $mo = IC ] && p=0.01622; run email $D/email-Eu-core.txt directed $mo $p 0,1,2,3,4 5000; done
for mo in WC LT IC; do p=0; [ $mo = IC ] && p=0.02215; run wiki $D/wiki-Vote.txt directed $mo $p 0,1,2,3,4 5000; done
for mo in WC LT IC; do p=0; [ $mo = IC ] && p=0.00616; run facebook $D/facebook_combined.txt undirected $mo $p 0,1,2,3,4 5000; done
for mo in WC LT IC; do p=0; [ $mo = IC ] && p=0.00408; run hepph $D/ca-HepPh.txt undirected $mo $p 0,1,2,3,4 5000; done
for mo in WC LT IC; do p=0; [ $mo = IC ] && p=0.00939; run epinions $D/soc-Epinions1.txt directed $mo $p 0,1,2 3000; done
for mo in WC LT IC; do p=0; [ $mo = IC ] && p=0.00808; run slashdot $D/soc-Slashdot0902.txt directed $mo $p 0,1,2 3000; done
fi
if [ $PART != small ]; then
for mo in WC LT IC; do p=0; [ $mo = IC ] && p=0.00863; run dblp $D/com-dblp.ungraph.txt undirected $mo $p 0,1,2 2000; done
fi
echo QUALITY_V2_${PART}_DONE
