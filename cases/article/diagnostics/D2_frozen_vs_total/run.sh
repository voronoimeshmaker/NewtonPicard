#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target diag_D2 -j 4
python3 support/run_suite.py D2
