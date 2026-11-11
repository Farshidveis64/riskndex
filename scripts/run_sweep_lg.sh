#!/bin/bash
# Size sweep at fixed influence structure (ablation w/o need gate).  ca-HepPh padded with isolated nodes to 1, 5, 20 and
# 100 times its size (12K to 1.2M nodes); RIS, RiskIndex-HV (without the forward layer) and RiskIndex-LG (ungated layer) on the same
# sketches and evaluation worlds (stream 31337+seed), seeds 0-2, WC and LT, M = 1024, rule b with cap 2^16.
set -e
cd "$(dirname "$0")/.."
B=${RIBIN:-bin/riskindex}; D=datasets/raw; mkdir -p results/sweep_lg
OUT=results/sweep_lg/sweep.csv; SYS=results/sweep_lg/sweep_sys.csv; rm -f $OUT $SYS $OUT.cert $OUT.hcert
LG="hyb=1 val=1 hvz=0 rhomin=0.5 sq=1 fl=1 rawidx=0 direct=0 sg=0 sgc=0 inddirect=0 nomat=0"
for f in 1 5 20 100; do for mo in WC LT; do
  echo "=== ca-HepPh x$f $mo"
  $B quality graph=pad:$D/ca-HepPh.txt:undirected:$f name=hepph-pad$f model=$mo M=1024 brule=100 bmax=65536 H=8 Kmax=100 ks=1,2,3,5,10 alphas=0.05,0.1,0.2,0.5,1.0 seeds=0,1,2 R=3000 $LG out=$OUT sysout=$SYS tag=sweep 2>&1 | grep -v "auto_b\|^\[hyb\]\|^\[val\]"
done; done
echo SWEEP_LG_DONE
