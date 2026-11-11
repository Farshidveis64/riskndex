#!/bin/bash
# two queues on the two cores
cd "$(dirname "$0")/../.."
if [ "$1" = 1 ]; then scripts/lg_dev/run_headroom.sh large > results/lg_dev/headroom_large.log 2>&1
else scripts/lg_dev/run_headroom.sh small > results/lg_dev/headroom_small.log 2>&1; scripts/lg_dev/run_headroom_youtube.sh > results/lg_dev/headroom_youtube.log 2>&1; fi
