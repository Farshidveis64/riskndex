#!/bin/bash
# Tier C, index without the forward layer (independent RR sketch + RIS chain for tau = inf, built and
# freed first; GRR build sketch and H-1 threshold chains; validation GRR sketch on which all chains are
# materialized, then freed; rho-guard 0.5; prefix sums of squares stored).  The RR sketch, the RIS chain and the
# validation sketch are charged to the index.  Quality on R = 2000 fresh worlds (1000 on LiveJournal).
set -e
cd "$(dirname "$0")/.."
B=${RIBIN:-bin/riskindex}; D=datasets/raw; mkdir -p results/scale_v2
OUT=results/scale_v2/scale.csv; QO=results/scale_v2/scale_queries.csv; EO=results/scale_v2/scale_quality.csv; HO=results/scale_v2/scale_hyb.csv
rm -f $OUT $QO $EO $HO
KS=10,25,50,100; ALPHAS=0.05,0.1,0.2,0.5,1.0
V2="hyb=1 val=1 rhomin=0.5 sq=1"
for n in 10000 30000 100000 300000 1000000; do
  $B scale graph=ba:$n:5:1 name=ba$n model=WC M=512 brule=100 bmax=16384 H=8 Kmax=100 $V2 ks=$KS alphas=$ALPHAS seeds=0,1,2 Q=1000 R=2000 out=$OUT qout=$QO eout=$EO hout=$HO tag=ba
done
for g in "email $D/email-Eu-core.txt directed" "wiki $D/wiki-Vote.txt directed" "epinions $D/soc-Epinions1.txt directed" "slashdot $D/soc-Slashdot0902.txt directed" "dblp $D/com-dblp.ungraph.txt undirected" "youtube $D/com-youtube.ungraph.txt undirected"; do
  set -- $g
  for mo in WC LT; do
    $B scale graph=snap:$2:$3 name=$1 model=$mo M=512 brule=100 bmax=16384 H=8 Kmax=100 $V2 ks=$KS alphas=$ALPHAS seeds=0,1,2 Q=1000 R=2000 out=$OUT qout=$QO eout=$EO hout=$HO tag=real
  done
done
gunzip -k -f $D/soc-LiveJournal1.txt.gz 2>/dev/null || true
$B scale graph=snap:$D/soc-LiveJournal1.txt:directed name=livejournal model=WC M=512 brule=100 bmax=16384 H=8 Kmax=100 $V2 ks=10,50,100 alphas=0.05,0.2,1.0 seeds=0 Q=1000 R=1000 out=$OUT qout=$QO eout=$EO hout=$HO tag=real
echo SCALE_V2_DONE
