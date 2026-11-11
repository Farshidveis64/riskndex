#!/bin/bash
cd "$(dirname "$0")/.."
date
for q in A B C; do
  n=$(wc -l < results/ad_jobs/$q.done 2>/dev/null || echo 0)
  L=results/ad_jobs/runner_$q.log; [ -f results/ad_jobs_$q.log ] && L=results/ad_jobs_$q.log
  last=$(tail -1 $L 2>/dev/null)
  if pgrep -f "ad_runner.py $q" > /dev/null; then st=running; elif grep -q ALL_DONE $L 2>/dev/null; then st=finished; else st=waiting/STOPPED; fi
  echo "queue $q: $st, $n done; $last"
done
