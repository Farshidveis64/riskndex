#!/bin/bash
# LargeGraph diagnosis on YouTube (development seed 100).
set -e
cd "$(dirname "$0")/../.."
B=${RIBIN:-bin/riskindex}; D=datasets/raw; OUT=results/lg_dev/headroom_youtube.csv; rm -f $OUT
for mo in WC LT IC; do p=0; [ $mo = IC ] && p=0.00475
  $B headroom graph=snap:$D/com-youtube.ungraph.txt:undirected name=youtube model=$mo p=$p M=1024 brule=100 bmax=65536 seed=100 P=60 Pscreen=2000 Wsel=1000 Weval=2000 out=$OUT 2>&1 | grep -v "^auto_b" ; done
echo HEADROOM_youtube_DONE
