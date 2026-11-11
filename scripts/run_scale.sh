#!/bin/bash
# Tier C: offline build cost, index size, online latency, workload, break-even (M=512, b by the rule, prefix sums of squares stored).
set -e
cd "$(dirname "$0")/.."
B=bin/riskindex; D=datasets/raw; mkdir -p results/scale
OUT=results/scale/scale.csv; QO=results/scale/scale_queries.csv; EO=results/scale/scale_quality.csv
rm -f $OUT $QO $EO
KS=10,25,50,100; ALPHAS=0.05,0.1,0.2,0.5,1.0
# synthetic BA sweep (avg degree ~10), WC, fixed M=256, b by the resolution rule
for n in 10000 30000 100000 300000 1000000; do
  $B scale graph=ba:$n:5:1 name=ba$n model=WC M=512 brule=100 bmax=16384 H=8 Kmax=100 sq=1 ks=$KS alphas=$ALPHAS seeds=0,1,2 Q=1000 R=500 out=$OUT qout=$QO eout=$EO tag=ba
done
# real graphs, WC and LT
for g in "email $D/email-Eu-core.txt directed" "wiki $D/wiki-Vote.txt directed" "epinions $D/soc-Epinions1.txt directed" "slashdot $D/soc-Slashdot0902.txt directed" "dblp $D/com-dblp.ungraph.txt undirected" "youtube $D/com-youtube.ungraph.txt undirected"; do
  set -- $g
  for mo in WC LT; do
    $B scale graph=snap:$2:$3 name=$1 model=$mo M=512 brule=100 bmax=16384 H=8 Kmax=100 sq=1 ks=$KS alphas=$ALPHAS seeds=0,1,2 Q=1000 R=500 out=$OUT qout=$QO eout=$EO tag=real
  done
done
# LiveJournal (4.8M nodes, 69M edges): single seed, reduced direct grid
gunzip -k -f $D/soc-LiveJournal1.txt.gz 2>/dev/null || true
$B scale graph=snap:$D/soc-LiveJournal1.txt:directed name=livejournal model=WC M=512 brule=100 bmax=16384 H=8 Kmax=100 sq=1 ks=10,50,100 alphas=0.05,0.2,1.0 seeds=0 Q=1000 R=200 out=$OUT qout=$QO eout=$EO tag=real
# Orkut (3.1M nodes, 234M directed edges): optional, needs ~6 GB RAM
if [ "${WITH_ORKUT:-0}" = "1" ]; then
  gunzip -k -f $D/com-orkut.ungraph.txt.gz 2>/dev/null || true
  $B scale graph=snap:$D/com-orkut.ungraph.txt:undirected name=orkut model=WC M=512 brule=100 bmax=16384 H=8 Kmax=100 sq=1 ks=10,50,100 alphas=0.05,0.2,1.0 seeds=0 Q=1000 R=0 out=$OUT qout=$QO eout=$EO tag=real
fi
echo SCALE_DONE
