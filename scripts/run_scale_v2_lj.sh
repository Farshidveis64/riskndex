#!/bin/bash
# Tier C, index without the forward layer: LiveJournal (single seed, reduced direct grid), appended to results/scale_v2.
set -e
cd "$(dirname "$0")/.."
B=${RIBIN:-bin/riskindex}; D=datasets/raw
OUT=results/scale_v2/scale.csv; QO=results/scale_v2/scale_queries.csv; EO=results/scale_v2/scale_quality.csv; HO=results/scale_v2/scale_hyb.csv
[ -f $D/soc-LiveJournal1.txt ] || gunzip -k -f $D/soc-LiveJournal1.txt.gz
$B scale graph=snap:$D/soc-LiveJournal1.txt:directed name=livejournal model=WC M=512 brule=100 bmax=16384 H=8 Kmax=100 hyb=1 val=1 rhomin=0.5 sq=1 ks=10,50,100 alphas=0.05,0.2,1.0 seeds=0 Q=1000 R=1000 out=$OUT qout=$QO eout=$EO hout=$HO tag=real
echo SCALE_V2_LJ_DONE
