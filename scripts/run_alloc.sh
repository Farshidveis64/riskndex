#!/bin/bash
# Budget allocation: worlds M vs targets b at a FIXED sketch budget M*b, on the four small graphs (three seeds).
#   allocA: M=512,  b = rule (c=100)       -- the v2 configuration
#   allocB: M=1024, b = rule/2 (c=50)
#   allocD: M=2048, b = rule/4 (c=25)
# The frozen default (M=1024, b = rule) in results/quality/ has twice this budget.
set -e
cd "$(dirname "$0")/.."
B=bin/riskindex; D=datasets/raw; mkdir -p results/alloc
OUT=results/alloc/alloc.csv; SYS=results/alloc/alloc_sys.csv; rm -f $OUT $SYS $OUT.cert
ALPHAS=0.05,0.1,0.2,0.5,1.0; KS=1,2,3,5,10,25,50,100
COMMON="H=8 Kmax=100 sq=1 ks=$KS alphas=$ALPHAS seeds=0,1,2 R=5000 rawidx=1 nomat=0 direct=0 inddirect=0"
GRAPHS=("email $D/email-Eu-core.txt directed 0.01622" "wiki $D/wiki-Vote.txt directed 0.02215"
        "facebook $D/facebook_combined.txt undirected 0.00616" "hepph $D/ca-HepPh.txt undirected 0.00408")
for g in "${GRAPHS[@]}"; do set -- $g
  for mo in WC LT IC; do p=0; [ $mo = IC ] && p=$4
    $B quality graph=snap:$2:$3 name=$1 model=$mo p=$p M=512  brule=100 bmax=16384 bscale=1.0  $COMMON out=$OUT sysout=$SYS tag=allocA
    $B quality graph=snap:$2:$3 name=$1 model=$mo p=$p M=1024 brule=100 bmax=16384 bscale=0.5  $COMMON out=$OUT sysout=$SYS tag=allocB
    $B quality graph=snap:$2:$3 name=$1 model=$mo p=$p M=2048 brule=100 bmax=16384 bscale=0.25 $COMMON out=$OUT sysout=$SYS tag=allocD
  done
done
echo ALLOC_DONE
