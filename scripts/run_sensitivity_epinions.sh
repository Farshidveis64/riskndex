#!/bin/bash
# Epinions part of the sensitivity study (same settings as run_sensitivity.sh; appends).
set -e
cd "$(dirname "$0")/.."
B=bin/riskindex; D=datasets/raw
OUT=results/sens/sens.csv; SYS=results/sens/sens_sys.csv
ALPHAS=0.05,0.1,0.2,0.5,1.0; S=0,1,2; KS=10,25,50,100; R=2000
G="graph=snap:$D/soc-Epinions1.txt:directed name=epinions model=WC p=0"
$B quality $G M=512 brule=100 bmax=16384 Hs=1,2,4,8,16,32 Kmax=100 ks=$KS alphas=$ALPHAS seeds=$S R=$R ind=0 direct=1 nomat=0 out=$OUT sysout=$SYS tag=Hsweep
$B quality $G Ms=64,128,256,512,1024 brule=100 bmax=16384 H=8 Kmax=100 ks=$KS alphas=$ALPHAS seeds=$S R=$R ind=0 direct=0 nomat=0 out=$OUT sysout=$SYS tag=Msweep
$B quality $G M=512 bs=64,256,1024,4096,16384 H=8 Kmax=100 ks=$KS alphas=$ALPHAS seeds=$S R=$R ind=0 direct=0 nomat=0 out=$OUT sysout=$SYS tag=bsweep
$B quality $G Ms=128,256,512,1024,2048,4096 brule=100 bmax=16384 H=8 Kmax=100 ks=$KS alphas=$ALPHAS seeds=$S R=$R ind=1 inddirect=0 direct=0 nomat=0 out=$OUT sysout=$SYS tag=Mprice
$B quality $G M=512 brule=100 bmax=16384 H=8 Kmaxs=25,50,100,200,400 ks=10 alphas=0.1 seeds=0,1,2 R=100 ind=0 direct=0 nomat=0 out=$OUT sysout=$SYS tag=Ksweep
echo SENS_EPINIONS_DONE
