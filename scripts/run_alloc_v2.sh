#!/bin/bash
# Budget allocation (index without the forward layer): worlds M vs targets b at a fixed sketch budget M*b, with the frozen
# RiskIndex; same sketches, seeds and evaluation worlds as run_alloc.sh.
set -e
cd "$(dirname "$0")/.."
B=${RIBIN:-bin/riskindex}; D=datasets/raw; mkdir -p results/alloc_v2
OUT=results/alloc_v2/alloc.csv; SYS=results/alloc_v2/alloc_sys.csv; rm -f $OUT $SYS $OUT.cert $OUT.hcert
ALPHAS=0.05,0.1,0.2,0.5,1.0; KS=1,2,3,5,10,25,50,100
COMMON="H=8 Kmax=100 sq=1 ks=$KS alphas=$ALPHAS seeds=0,1,2 R=5000 rawidx=0 nomat=0 direct=0 inddirect=0 hyb=1 val=1 rhomin=0.5"
GRAPHS=("email $D/email-Eu-core.txt directed 0.01622" "wiki $D/wiki-Vote.txt directed 0.02215"
        "facebook $D/facebook_combined.txt undirected 0.00616" "hepph $D/ca-HepPh.txt undirected 0.00408")
for g in "${GRAPHS[@]}"; do set -- $g
  for mo in WC LT IC; do p=0; [ $mo = IC ] && p=$4
    $B quality graph=snap:$2:$3 name=$1 model=$mo p=$p M=512  brule=100 bmax=16384 bscale=1.0  $COMMON out=$OUT sysout=$SYS tag=allocA 2>&1 | grep -v auto_b
    $B quality graph=snap:$2:$3 name=$1 model=$mo p=$p M=1024 brule=100 bmax=16384 bscale=0.5  $COMMON out=$OUT sysout=$SYS tag=allocB 2>&1 | grep -v auto_b
    $B quality graph=snap:$2:$3 name=$1 model=$mo p=$p M=2048 brule=100 bmax=16384 bscale=0.25 $COMMON out=$OUT sysout=$SYS tag=allocD 2>&1 | grep -v auto_b
  done
done
echo ALLOC_V2_DONE
