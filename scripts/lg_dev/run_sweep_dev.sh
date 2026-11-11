#!/bin/bash
# LargeGraph development: size sweep at fixed influence structure (ca-HepPh padded with isolated nodes to 1, 5, 20
# and 100 times its size), development seed 100.  RIS, RiskIndex-HV (no forward layer) and RiskIndex-LG (ungated layer).
set -e
cd "$(dirname "$0")/../.."
B=${RIBIN:-bin/riskindex}; D=datasets/raw; OUT=results/lg_dev/sweep_dev.csv; SYS=results/lg_dev/sweep_dev_sys.csv
V2="hyb=1 val=1 hvz=0 rhomin=0.5 sq=1 rawidx=0 direct=0 sg=0 sgc=0 inddirect=0 nomat=0 fl=1"
for f in ${FACTORS:-1 5 20 100}; do for mo in ${MODELS:-WC LT}; do
  echo "=== hepph pad x$f $mo"
  $B quality graph=pad:$D/ca-HepPh.txt:undirected:$f name=hepph-pad$f model=$mo p=0.00408 M=1024 brule=100 bmax=65536 H=8 Kmax=100 ks=1,2,3,5,10 alphas=0.05,0.1,0.2,0.5,1.0 seeds=${SEEDS:-100} R=3000 $V2 out=$OUT sysout=$SYS tag=dev 2>&1 | grep -v "^\[hyb\]\|^\[val\]"
done; done
echo SWEEP_DEV_DONE
