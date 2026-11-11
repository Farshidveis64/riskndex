#!/bin/bash
# Sensitivity of the GRR index (no forward layer): H, M and b sweeps and the price-of-grouping sweep
# (RIS chain + validation sketch + rho-guard 0.5).  Identical sketches (seeds) and evaluation worlds as
# run_sensitivity.sh, so Direct-CVaR (GRR) of the H sweep -- which does not depend on H -- is joined from
# results/sens/sens.csv instead of being recomputed.
set -e
cd "$(dirname "$0")/.."
B=${RIBIN:-bin/riskindex}; D=datasets/raw; mkdir -p results/sens_v2
OUT=results/sens_v2/sens.csv; SYS=results/sens_v2/sens_sys.csv; rm -f $OUT $SYS $OUT.cert $OUT.hcert
ALPHAS=0.05,0.1,0.2,0.5,1.0; S=0,1,2
V2="ind=1 hyb=1 val=1 rhomin=0.5 sq=1 direct=0 inddirect=0 nomat=0 rawidx=0 sg=0 sgc=0"
for g in "wiki $D/wiki-Vote.txt directed 0.02215 WC" "wiki $D/wiki-Vote.txt directed 0.02215 IC" "email $D/email-Eu-core.txt directed 0.01622 WC" "email $D/email-Eu-core.txt directed 0.01622 IC" "epinions $D/soc-Epinions1.txt directed 0.00939 WC"; do
  set -- $g; name=$1; path=$2; dir=$3; pc=$4; mo=$5
  case $name in epinions) KS=10,25,50,100; R=2000;; *) KS=2,5,10,25,50,100; R=4000;; esac
  p=0; [ $mo = IC ] && p=$pc
  $B quality graph=snap:$path:$dir name=$name model=$mo p=$p M=512 brule=100 bmax=16384 Hs=1,2,4,8,16,32 Kmax=100 ks=$KS alphas=$ALPHAS seeds=$S R=$R $V2 out=$OUT sysout=$SYS tag=Hsweep 2>&1 | grep -v auto_b
  $B quality graph=snap:$path:$dir name=$name model=$mo p=$p Ms=64,128,256,512,1024 brule=100 bmax=16384 H=8 Kmax=100 ks=$KS alphas=$ALPHAS seeds=$S R=$R $V2 out=$OUT sysout=$SYS tag=Msweep 2>&1 | grep -v auto_b
  $B quality graph=snap:$path:$dir name=$name model=$mo p=$p M=512 bs=64,256,1024,4096,16384 H=8 Kmax=100 ks=$KS alphas=$ALPHAS seeds=$S R=$R $V2 out=$OUT sysout=$SYS tag=bsweep 2>&1 | grep -v auto_b
  $B quality graph=snap:$path:$dir name=$name model=$mo p=$p Ms=128,256,512,1024,2048,4096 brule=100 bmax=16384 H=8 Kmax=100 ks=$KS alphas=$ALPHAS seeds=$S R=$R $V2 out=$OUT sysout=$SYS tag=Mprice 2>&1 | grep -v auto_b
done
echo SENS_V2_DONE
