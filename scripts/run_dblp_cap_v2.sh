#!/bin/bash
# DBLP IC-critical (the under-resolved Tier B cell): the frozen RiskIndex at the memory caps b = 2^14 and 2^17
# (2^16 is the Tier B configuration, results/quality_v2).  Same seeds and evaluation worlds as Tier B.
set -e
cd "$(dirname "$0")/.."
B=${RIBIN:-bin/riskindex}; D=datasets/raw; mkdir -p results/alloc_v2
OUT=results/alloc_v2/dblp_cap.csv; SYS=results/alloc_v2/dblp_cap_sys.csv; rm -f $OUT $SYS $OUT.cert $OUT.hcert
for b in 16384 131072; do
  $B quality graph=snap:$D/com-dblp.ungraph.txt:undirected name=dblp model=IC p=0.00863 M=1024 b=$b H=8 Kmax=100 sq=1 rhomin=0.5 \
     ks=1,2,3,5,10,25,50,100 alphas=0.05,0.1,0.2,0.5,1.0 seeds=0,1,2 R=2000 hyb=1 val=1 hvzs=1 rawidx=0 direct=0 sg=0 sgc=0 inddirect=0 nomat=0 \
     out=$OUT sysout=$SYS tag=cap$b 2>&1 | grep -v auto_b
done
echo DBLP_CAP_V2_DONE
