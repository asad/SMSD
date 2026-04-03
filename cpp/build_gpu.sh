#!/bin/bash
set -e
cd /Users/asad/tool/SMSD/cpp && rm -rf build && \
cmake -B build -DSMSD_BUILD_METAL=OFF -DCMAKE_BUILD_TYPE=Release && \
cmake --build build --config Release -j 8 && \
cd build && ctest --output-on-failure -V && \
cd /Users/asad/tool/SMSD && mvn clean test && \
cd benchmarks && python3 benchmark_python.py && \
echo "ALL DONE"

