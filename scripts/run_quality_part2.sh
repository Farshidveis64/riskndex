#!/bin/bash
# Tier B, part 2 (appends): Epinions, Slashdot (three seeds, R=3000) and DBLP (three seeds, R=2000).
# Same frozen configuration as part 1.  Direct-CVaR-VM is not run (cost); Direct-CVaR (ind. RR) is not run on DBLP.
set -e
cd "$(dirname "$0")/.."
B=bin/riskindex; D=datasets/raw
OUT=results/quality/quality.csv; SYS=results/quality/quality_sys.csv
ALPHAS=0.05,0.1,0.2,0.5,1.0; KS=1,2,3,5,10,25,50,100
run() { echo "=== $1 $4"; $B quality graph=snap:$2:$3 name=$1 model=$4 p=$5 M=1024 brule=100 bmax=65536 H=8 Kmax=100 sq=1 rhomin=0.5 ks=$KS alphas=$ALPHAS seeds=0,1,2 R=$6 inddirect=$7 rawidx=1 sg=1 sgc=1 out=$OUT sysout=$SYS tag=tierB; }
for mo in WC LT IC; do p=0; [ $mo = IC ] && p=0.00939; run epinions $D/soc-Epinions1.txt directed $mo $p 3000 1; done
for mo in WC LT IC; do p=0; [ $mo = IC ] && p=0.00808; run slashdot $D/soc-Slashdot0902.txt directed $mo $p 3000 1; done
for mo in WC LT IC; do p=0; [ $mo = IC ] && p=0.00863; run dblp $D/com-dblp.ungraph.txt undirected $mo $p 2000 0; done
echo QUALITY_PART2_DONE
