#!/bin/bash
# Development check of the rho-guard for the final index (RIS chain + validation): rho_min = 0 vs 0.5 on the
# under-resolved dev sketches (b = 128, 512) and the deployed sketch; dev seeds 100-102, four small graphs.
cd $(dirname "$0")/../..; B=${RIBIN:-bin/riskindex}; D=datasets/raw; O=results/dev/guard
rm -f ${O}_*.csv*
GR=("email $D/email-Eu-core.txt directed 0.01622" "wiki $D/wiki-Vote.txt directed 0.02215" "facebook $D/facebook_combined.txt undirected 0.00616" "hepph $D/ca-HepPh.txt undirected 0.00408")
for rm in 0 0.5; do
COMMON="H=8 Kmax=100 sq=1 rhomin=$rm ks=1,2,3,5,10,25,50,100 alphas=0.05,0.1,0.2,0.5,1.0 seeds=100,101,102 R=3000 rawidx=0 nomat=0 direct=0 inddirect=0 sg=0 sgc=0 hyb=1 val=1"
for g in "${GR[@]}"; do set -- $g; for mo in WC LT IC; do p=0; [ $mo = IC ] && p=$4;
  for b in 128 512; do $B quality graph=snap:$2:$3 name=$1 model=$mo p=$p M=1024 b=$b $COMMON out=${O}_b${b}_r$rm.csv sysout=${O}_sys.csv tag=b$b 2>&1 | grep -v auto_b; done
done; done; done; echo GUARD_DEV_DONE
