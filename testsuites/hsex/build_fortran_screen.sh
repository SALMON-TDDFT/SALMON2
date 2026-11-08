#!/bin/bash
# Standalone verification uses production screening and communication sources.
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
out=${1:?provide a new out-of-source build directory}
if test -e "$out"; then echo 'Build directory must be new' >&2; exit 1; fi
python3 "$root/tools/check_forbidden_eigensolvers.py" --source-root "$root/src"
python3 "$root/tools/check_forbidden_eigensolvers.py" --source-root "$root/tools/hsex/native"
mkdir -p "$out"
out=$(cd "$out" && pwd)
fc=${FC:-mpif90}
flags=${SCREEN_FFLAGS:--O2 -fopenmp -cpp -fallow-argument-mismatch -ffree-line-length-none}
blas=${SCREEN_BLAS:--framework Accelerate}
native_flags=${SCREEN_NATIVE_FLAGS--DHSEX_ACCELERATE}
cd "$out"
# Compiler flag variables are intentionally word-split, as in build tooling.
$fc $flags -c "$root/src/misc/nvtx_wrapper.f90" -o nvtx.o
$fc $flags -c "$root/src/parallel/communication.f90" -o comm.o
$fc $flags $native_flags -c "$root/tools/hsex/native/q_response.f90" -o q.o
$fc $flags -c "$root/src/xc/hsex_screen_math.f90" -o math.o
$fc $flags -c "$root/src/xc/hsex_sternheimer.f90" -o sternheimer.o
$fc $flags -c "$root/src/xc/hsex_screen_driver.f90" -o driver.o
for name in driver head pipeline refine; do
  $fc $flags -c "$root/testsuites/hsex/fortran_screen_${name}.f90" -o "test_${name}.o"
done
$fc $flags -c "$root/testsuites/hsex/fortran_sternheimer_static.f90" -o test_static.o
python3 "$root/tools/check_forbidden_eigensolvers.py" --object-tree "$out"
$fc $flags test_driver.o math.o $blas -o assembly
$fc $flags test_head.o math.o $blas -o head
$fc $flags test_pipeline.o driver.o sternheimer.o math.o q.o comm.o nvtx.o $blas -o pipeline

$fc $flags test_refine.o q.o $blas -o refine
$fc $flags test_static.o sternheimer.o math.o $blas -o static_sternheimer
