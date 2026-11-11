#!/bin/bash
# RQ1: estimation fidelity, GRR vs independent RR at equal total budget M*b.
set -e
cd "$(dirname "$0")/.."
B=bin/riskindex; OUT=results/fidelity/fidelity.csv; mkdir -p results/fidelity; rm -f $OUT
D=datasets/raw
# b sweep at fixed M = 256 (budget M*b grows with b); 100 random seed sets of size k; truth on R fresh worlds
$B fidelity graph=snap:$D/email-Eu-core.txt:directed name=email model=WC fb=2048 M=256 bs=1,4,16,64,256,1024,4096 seeds=0,1,2,3,4 nsets=100 k=5 R=20000 out=$OUT tag=bsweep
$B fidelity graph=snap:$D/email-Eu-core.txt:directed name=email model=IC fb=2048 p=0.01622 M=256 bs=1,4,16,64,256,1024,4096 seeds=0,1,2,3,4 nsets=100 k=5 R=20000 out=$OUT tag=bsweep
$B fidelity graph=snap:$D/wiki-Vote.txt:directed name=wiki model=WC fb=4096 M=256 bs=1,4,16,64,256,1024,4096 seeds=0,1,2,3,4 nsets=100 k=5 R=20000 out=$OUT tag=bsweep
$B fidelity graph=snap:$D/wiki-Vote.txt:directed name=wiki model=IC fb=4096 p=0.02215 M=256 bs=1,4,16,64,256,1024,4096 seeds=0,1,2,3,4 nsets=100 k=5 R=20000 out=$OUT tag=bsweep
$B fidelity graph=snap:$D/wiki-Vote.txt:directed name=wiki model=LT fb=4096 M=256 bs=1,4,16,64,256,1024,4096 seeds=0,1,2,3,4 nsets=100 k=5 R=20000 out=$OUT tag=bsweep
$B fidelity graph=snap:$D/soc-Epinions1.txt:directed name=epinions model=WC fb=4096 M=256 bs=16,64,256,1024,4096 seeds=0,1,2 nsets=100 k=10 R=5000 out=$OUT tag=bsweep
# M sweep at fixed b = 1024 (between-world resolution)
$B fidelity graph=snap:$D/wiki-Vote.txt:directed name=wiki model=WC fb=4096 M=64 bs=1024 seeds=0,1,2,3,4 nsets=100 k=5 R=20000 out=$OUT tag=Msweep
$B fidelity graph=snap:$D/wiki-Vote.txt:directed name=wiki model=WC fb=4096 M=128 bs=1024 seeds=0,1,2,3,4 nsets=100 k=5 R=20000 out=$OUT tag=Msweep
$B fidelity graph=snap:$D/wiki-Vote.txt:directed name=wiki model=WC fb=4096 M=512 bs=1024 seeds=0,1,2,3,4 nsets=100 k=5 R=20000 out=$OUT tag=Msweep
$B fidelity graph=snap:$D/wiki-Vote.txt:directed name=wiki model=WC fb=4096 M=1024 bs=1024 seeds=0,1,2,3,4 nsets=100 k=5 R=20000 out=$OUT tag=Msweep
$B fidelity graph=snap:$D/wiki-Vote.txt:directed name=wiki model=WC fb=4096 M=2048 bs=1024 seeds=0,1,2,3,4 nsets=100 k=5 R=20000 out=$OUT tag=Msweep
echo FIDELITY_DONE
