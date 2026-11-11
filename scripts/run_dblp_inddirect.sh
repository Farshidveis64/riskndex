#!/bin/bash
# Fills the Direct-CVaR (ind. RR) cells of DBLP in Table II (alpha = 0.05, the level of the table; three seeds).
# Same sketches (seeds) and evaluation worlds as the Tier B runs.
set -e
cd "$(dirname "$0")/.."
B=${RIBIN:-bin/riskindex}; D=datasets/raw; mkdir -p results/quality_v2
OUT=results/quality_v2/dblp_inddirect.csv; SYS=results/quality_v2/dblp_inddirect_sys.csv; rm -f $OUT $SYS $OUT.cert
for mo in WC LT IC; do p=0; [ $mo = IC ] && p=0.00863
  $B quality graph=snap:$D/com-dblp.ungraph.txt:undirected name=dblp model=$mo p=$p M=1024 brule=100 bmax=65536 H=8 Kmax=100 sq=1 rhomin=0.5 \
     ks=1,2,3,5,10,25,50,100 alphas=0.05 seeds=0,1,2 R=2000 rawidx=0 direct=0 inddirect=2 sg=0 sgc=0 nomat=0 out=$OUT sysout=$SYS tag=tierB 2>&1 | grep -v auto_b
done
echo DBLP_INDDIRECT_DONE
