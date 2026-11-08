#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/../../.."
export FC=frtpx Q_COMPILE_FLAGS='-O3 -Kpic -Cpp -Nalloc_assign' Q_OMP_FLAGS='-Kopenmp -Nfjomplib' Q_BLAS_FLAGS='-SSL2'
python3 tools/hsex/build_q_openmp.py "${1:?provide a new out-of-source build directory}"
