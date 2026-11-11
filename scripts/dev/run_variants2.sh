#!/bin/bash
cd $(dirname "$0")/../..; B=bin/riskindex; D=datasets/raw; O=results/dev/var
rm -f ${O}_E* ${O}_F*
GR=("email $D/email-Eu-core.txt directed 0.01622" "wiki $D/wiki-Vote.txt directed 0.02215" "facebook $D/facebook_combined.txt undirected 0.00616" "hepph $D/ca-HepPh.txt undirected 0.00408")
COMMON="H=8 Kmax=100 ks=1,2,3,5,10,25,50,100 alphas=0.05,0.1,0.2,0.5,1.0 seeds=100,101,102 R=3000 rawidx=1 nomat=0 direct=0 inddirect=0"
for g in "${GR[@]}"; do set -- $g; for mo in WC LT IC; do p=0; [ $mo = IC ] && p=$4;
  $B quality graph=snap:$2:$3 name=$1 model=$mo p=$p M=1024 brule=100 bmax=16384 $COMMON out=${O}_E.csv sysout=${O}_E_sys.csv tag=E 2>&1 | grep -v auto_b
  $B quality graph=snap:$2:$3 name=$1 model=$mo p=$p M=2048 brule=100 bmax=16384 bscale=0.5 $COMMON out=${O}_F.csv sysout=${O}_F_sys.csv tag=F 2>&1 | grep -v auto_b
done; done; echo VARIANTS2_DONE
