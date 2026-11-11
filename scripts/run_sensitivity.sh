#!/bin/bash
# RQ3/RQ6: threshold grid H, worlds M, targets-per-world b, Kmax, price of grouping.
set -e
cd "$(dirname "$0")/.."
B=bin/riskindex; D=datasets/raw; mkdir -p results/sens
OUT=results/sens/sens.csv; SYS=results/sens/sens_sys.csv; rm -f $OUT $SYS $OUT.cert
ALPHAS=0.05,0.1,0.2,0.5,1.0; S=0,1,2
for g in "wiki $D/wiki-Vote.txt directed 0.02215 WC" "wiki $D/wiki-Vote.txt directed 0.02215 IC" "email $D/email-Eu-core.txt directed 0.01622 WC" "email $D/email-Eu-core.txt directed 0.01622 IC" "epinions $D/soc-Epinions1.txt directed 0.00939 WC"; do
  set -- $g; name=$1; path=$2; dir=$3; pc=$4; mo=$5
  case $name in epinions) KS=10,25,50,100; R=2000;; *) KS=2,5,10,25,50,100; R=4000;; esac
  p=0; [ $mo = IC ] && p=$pc
  # H sweep (b fixed by rule, M=512), compared against Direct-CVaR on the same sketch
  $B quality graph=snap:$path:$dir name=$name model=$mo p=$p M=512 brule=100 bmax=16384 Hs=1,2,4,8,16,32 Kmax=100 ks=$KS alphas=$ALPHAS seeds=$S R=$R ind=0 direct=1 nomat=0 out=$OUT sysout=$SYS tag=Hsweep
  # M sweep (b fixed by rule)
  $B quality graph=snap:$path:$dir name=$name model=$mo p=$p Ms=64,128,256,512,1024 brule=100 bmax=16384 H=8 Kmax=100 ks=$KS alphas=$ALPHAS seeds=$S R=$R ind=0 direct=0 nomat=0 out=$OUT sysout=$SYS tag=Msweep
  # b sweep (M=512)
  $B quality graph=snap:$path:$dir name=$name model=$mo p=$p M=512 bs=64,256,1024,4096,16384 H=8 Kmax=100 ks=$KS alphas=$ALPHAS seeds=$S R=$R ind=0 direct=0 nomat=0 out=$OUT sysout=$SYS tag=bsweep
  # price of grouping: mean-chain quality vs independent RIS at equal budget, as M grows
  $B quality graph=snap:$path:$dir name=$name model=$mo p=$p Ms=128,256,512,1024,2048,4096 brule=100 bmax=16384 H=8 Kmax=100 ks=$KS alphas=$ALPHAS seeds=$S R=$R ind=1 inddirect=0 direct=0 nomat=0 out=$OUT sysout=$SYS tag=Mprice
  # Kmax sweep (index size / build cost only)
  $B quality graph=snap:$path:$dir name=$name model=$mo p=$p M=512 brule=100 bmax=16384 H=8 Kmaxs=25,50,100,200,400 ks=10 alphas=0.1 seeds=0,1,2 R=100 ind=0 direct=0 nomat=0 out=$OUT sysout=$SYS tag=Ksweep
done
echo SENS_DONE
