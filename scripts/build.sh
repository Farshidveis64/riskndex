#!/bin/bash
# Build the single-threaded C++17 driver.
set -e
cd "$(dirname "$0")/.."
mkdir -p bin
g++ -std=c++17 -O3 -march=native -Wall -Wno-unused-result -o bin/riskindex src/main.cpp
g++ -std=c++17 -O2 -Wall -Wno-unused-result -Wno-unused-function -o bin/toy_fig3 src/tools/toy_fig3.cpp
echo "built bin/riskindex and bin/toy_fig3 ($(g++ --version | head -1))"
