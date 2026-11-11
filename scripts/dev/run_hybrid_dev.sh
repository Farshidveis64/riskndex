#!/bin/bash
# Development runs for the hybrid mean chain and the selection margin (dev sketch seeds 100-102, disjoint
# from the reported seeds 0-4; evaluation stream 31337+seed is therefore also disjoint).  Four smallest graphs only.
# Settings: the deployed sketch (M=1024, rule b, cap 2^16) and two under-resolved sketches (b=512, b=128).
cd $(dirname "$0")/../..; B=bin/riskindex; D=datasets/raw; O=results/dev/hyb
rm -f ${O}_*.csv*
GR=("email $D/email-Eu-core.txt directed 0.01622" "wiki $D/wiki-Vote.txt directed 0.02215" "facebook $D/facebook_combined.txt undirected 0.00616" "hepph $D/ca-HepPh.txt undirected 0.00408")
COMMON="H=8 Kmax=100 sq=1 rhomin=0.5 ks=1,2,3,5,10,25,50,100 alphas=0.05,0.1,0.2,0.5,1.0 seeds=100,101,102 R=3000 rawidx=0 nomat=0 direct=0 inddirect=0 sg=0 sgc=0 hyb=1 hzs=0.5,1,2,3 hzbs=1,2"
for g in "${GR[@]}"; do set -- $g; for mo in WC LT IC; do p=0; [ $mo = IC ] && p=$4;
  $B quality graph=snap:$2:$3 name=$1 model=$mo p=$p M=1024 brule=100 bmax=65536 $COMMON out=${O}_rule.csv sysout=${O}_rule_sys.csv tag=rule 2>&1 | grep -v auto_b
  $B quality graph=snap:$2:$3 name=$1 model=$mo p=$p M=1024 b=512 $COMMON out=${O}_b512.csv sysout=${O}_b512_sys.csv tag=b512 2>&1 | grep -v auto_b
  $B quality graph=snap:$2:$3 name=$1 model=$mo p=$p M=1024 b=128 $COMMON out=${O}_b128.csv sysout=${O}_b128_sys.csv tag=b128 2>&1 | grep -v auto_b
done; done; echo HYBDEV_DONE
