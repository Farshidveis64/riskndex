#!/bin/bash
# Tier B, part 1: email-Eu-core, Wiki-Vote, Facebook, ca-HepPh (five sketch seeds).
# Frozen configuration: M=1024, b by the resolution rule (c=100, memory cap 2^16), H=8, Kmax=100, rho-guard rho_min=0.5,
# variance-matched profiles with prefix sums of squares (sq=1).  Fresh-world evaluation on R=5000 worlds.
# Baselines: RiskIndex-raw, Direct-CVaR (GRR), Direct-CVaR (ind. RR), Direct-CVaR-VM (valid bound; not on ca-HepPh,
# where it needs minutes per query), SG and SG-cached (Zhou & Tokekar), RIS-Mean (GRR and ind. RR).
set -e
cd "$(dirname "$0")/.."
B=bin/riskindex; D=datasets/raw; mkdir -p results/quality
OUT=results/quality/quality.csv; SYS=results/quality/quality_sys.csv
rm -f $OUT $SYS $OUT.cert
ALPHAS=0.05,0.1,0.2,0.5,1.0; KS=1,2,3,5,10,25,50,100
GRAPHS=("email $D/email-Eu-core.txt directed 0.01622" "wiki $D/wiki-Vote.txt directed 0.02215"
        "facebook $D/facebook_combined.txt undirected 0.00616" "hepph $D/ca-HepPh.txt undirected 0.00408")
for g in "${GRAPHS[@]}"; do set -- $g; name=$1
  DV=1; [ $name = hepph ] && DV=0
  for mo in WC LT IC; do p=0; [ $mo = IC ] && p=$4
    echo "=== $name $mo p=$p"
    $B quality graph=snap:$2:$3 name=$name model=$mo p=$p M=1024 brule=100 bmax=65536 H=8 Kmax=100 sq=1 rhomin=0.5 ks=$KS alphas=$ALPHAS seeds=0,1,2,3,4 R=5000 \
       rawidx=1 directvm=$DV sg=1 sgc=1 inddirect=1 out=$OUT sysout=$SYS tag=tierB
  done
done
echo QUALITY_PART1_DONE
