#!/usr/bin/env bash
set -euo pipefail
case_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(cd -- "$case_dir/../../.." && pwd)"
cmake -S "$project_dir" -B "$project_dir/build" -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build "$project_dir/build" --target fadeev_newton_case --parallel 4
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
"$project_dir/build/cases/fadeev/fadeev_newton_case" "$case_dir/config.yaml" "$case_dir/results"
