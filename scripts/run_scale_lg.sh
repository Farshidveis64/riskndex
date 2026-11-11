#!/bin/bash
# Tier C, ablation w/o need gate: build cost and small-budget lookup latency of the ungated forward layer (fl=1) on the Tier C
# graphs, with the same sketches and seeds as results/scale_v2.  Direct-CVaR costs and quality at k >= 10 are
# unchanged and are taken from results/scale_v2 (the forward layer serves only k <= 5).
set -e
cd "$(dirname "$0")/.."
B=${RIBIN:-bin/riskindex}; D=datasets/raw; mkdir -p results/scale_lg
OUT=results/scale_lg/scale.csv; HO=results/scale_lg/scale_hyb.csv; rm -f $OUT $HO
KS=10,25,50,100; ALPHAS=0.05,0.1,0.2,0.5,1.0
LG="hyb=1 val=1 rhomin=0.5 sq=1 fl=1 direct=0"
for n in 10000 100000 1000000; do
  $B scale graph=ba:$n:5:1 name=ba$n model=WC M=512 brule=100 bmax=16384 H=8 Kmax=100 $LG ks=$KS alphas=$ALPHAS seeds=0,1,2 Q=1000 R=0 out=$OUT qout=/dev/null hout=$HO tag=ba
done
for g in "email $D/email-Eu-core.txt directed" "wiki $D/wiki-Vote.txt directed" "epinions $D/soc-Epinions1.txt directed" "slashdot $D/soc-Slashdot0902.txt directed" "dblp $D/com-dblp.ungraph.txt undirected" "youtube $D/com-youtube.ungraph.txt undirected"; do
  set -- $g
  for mo in WC LT; do
    $B scale graph=snap:$2:$3 name=$1 model=$mo M=512 brule=100 bmax=16384 H=8 Kmax=100 $LG ks=$KS alphas=$ALPHAS seeds=0,1,2 Q=1000 R=0 out=$OUT qout=/dev/null hout=$HO tag=real
  done
done
[ -f $D/soc-LiveJournal1.txt ] || gunzip -k -f $D/soc-LiveJournal1.txt.gz
$B scale graph=snap:$D/soc-LiveJournal1.txt:directed name=livejournal model=WC M=512 brule=100 bmax=16384 H=8 Kmax=100 $LG ks=10,50,100 alphas=0.05,0.2,1.0 seeds=0 Q=1000 R=0 out=$OUT qout=/dev/null hout=$HO tag=real
echo SCALE_LG_DONE
