#!/bin/bash
# Download the public SNAP graphs used in the paper (official URLs).
set -e
cd "$(dirname "$0")/.."
mkdir -p datasets/raw && cd datasets/raw
get() { [ -f "$2" ] || curl -sS -f -L -o "$2" "$1"; [ -f "${2%.gz}" ] || gunzip -k "$2"; }
get https://snap.stanford.edu/data/email-Eu-core.txt.gz email-Eu-core.txt.gz
get https://snap.stanford.edu/data/wiki-Vote.txt.gz wiki-Vote.txt.gz
get https://snap.stanford.edu/data/facebook_combined.txt.gz facebook_combined.txt.gz
get https://snap.stanford.edu/data/ca-HepPh.txt.gz ca-HepPh.txt.gz
get https://snap.stanford.edu/data/soc-Epinions1.txt.gz soc-Epinions1.txt.gz
get https://snap.stanford.edu/data/soc-Slashdot0902.txt.gz soc-Slashdot0902.txt.gz
get https://snap.stanford.edu/data/bigdata/communities/com-dblp.ungraph.txt.gz com-dblp.ungraph.txt.gz
get https://snap.stanford.edu/data/bigdata/communities/com-youtube.ungraph.txt.gz com-youtube.ungraph.txt.gz
get https://snap.stanford.edu/data/soc-LiveJournal1.txt.gz soc-LiveJournal1.txt.gz
echo DATA_DONE
