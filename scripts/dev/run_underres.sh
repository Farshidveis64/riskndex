#!/bin/bash
cd $(dirname "$0")/../..; B=bin/riskindex; D=datasets/raw; O=results/dev/ur
rm -f $O*.csv*
GR=("wiki $D/wiki-Vote.txt directed 0.02215" "facebook $D/facebook_combined.txt undirected 0.00616" "hepph $D/ca-HepPh.txt undirected 0.00408")
COMMON="H=8 Kmax=100 sq=1 ks=1,2,3,5,10,25,50,100 alphas=0.05,0.1,0.2,0.5,1.0 seeds=100,101,102 R=3000 rawidx=0 nomat=0 direct=0 inddirect=0 rhomins=0.3,0.5,0.7 zs=1,2 sgg=1 sgc=1"
for g in "${GR[@]}"; do set -- $g; for mo in IC LT WC; do p=0; [ $mo = IC ] && p=$4;
  for b in 128 512; do $B quality graph=snap:$2:$3 name=$1 model=$mo p=$p M=1024 b=$b $COMMON out=${O}_b$b.csv sysout=${O}_b${b}_sys.csv tag=b$b 2>&1 | grep -v auto_b; done
  $B quality graph=snap:$2:$3 name=$1 model=$mo p=$p M=1024 brule=100 bmax=16384 $COMMON out=${O}_rule.csv sysout=${O}_rule_sys.csv tag=rule 2>&1 | grep -v auto_b
done; done; echo UR_DONE
