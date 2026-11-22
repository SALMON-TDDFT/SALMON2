#!/bin/bash
# Native focused tests; no production computation is orchestrated in Python.
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
out=${1:?provide a new out-of-source directory}
test ! -e "$out" || { echo 'Build directory must be new' >&2; exit 1; }
python3 "$root/tools/check_forbidden_eigensolvers.py" --source-root "$root/src"
mkdir -p "$out"
cd "$out"
fc=${FC:-mpif90}
flags=${SCREEN_FFLAGS:--O2 -g -cpp -fopenmp -fallow-argument-mismatch -fcheck=all -fbacktrace}
blas=${SCREEN_BLAS:--framework Accelerate}
$fc $flags -c "$root/src/xc/hsex_screen_math.f90" -o math.o
$fc $flags -c "$root/testsuites/hsex/fortran_grid_translation.f90" -o test_grid_translation.o
$fc $flags -c "$root/testsuites/hsex/fortran_q2_head.f90" -o test_q2_head.o
$fc $flags -c "$root/src/xc/hsex_projected_derivative.f90" -o projected_derivative.o
$fc $flags -c "$root/testsuites/hsex/fortran_projected_derivative.f90" -o test_projected_derivative.o
$fc $flags -c "$root/src/xc/hsex_frozen_window.f90" -o window.o
$fc $flags -c "$root/src/xc/hsex_frozen_response.f90" -o frozen.o
$fc $flags -c "$root/src/xc/hsex_frozen_kinetic.f90" -o frozen_kinetic.o
$fc $flags -c "$root/testsuites/hsex/fortran_frozen_kinetic.f90" -o test_frozen_kinetic.o
$fc $flags -c "$root/src/xc/hsex_longwave.f90" -o longwave.o
$fc $flags -c "$root/src/xc/hsex_coherent_optical.f90" -o coherent_optical.o
$fc $flags -c "$root/testsuites/hsex/fortran_coherent_optical.f90" -o test_coherent_optical.o
$fc $flags -c "$root/src/xc/hsex_frozen_limit.f90" -o frozen_limit.o
$fc $flags -c "$root/testsuites/hsex/fortran_frozen_limit.f90" -o test_limit.o
$fc $flags -c "$root/testsuites/hsex/fortran_frozen_window.f90" -o test_window.o
$fc $flags -c "$root/testsuites/hsex/fortran_frozen_response.f90" -o test.o
$fc $flags -c "$root/testsuites/hsex/fortran_frozen_time.f90" -o time.o
$fc $flags -c "$root/src/misc/nvtx_wrapper.f90" -o nvtx.o
$fc $flags -c "$root/src/parallel/communication.f90" -o comm.o
$fc $flags -c "$root/src/xc/hsex_impulse_resolvent.f90" -o impulse.o
$fc $flags -c "$root/src/xc/hsex_equilibrium_q.f90" -o equilibrium.o
$fc $flags -c "$root/testsuites/hsex/fortran_equilibrium_q.f90" -o test_equilibrium.o
$fc $flags -c "$root/testsuites/hsex/fortran_impulse_resolvent.f90" -o test_impulse.o
$fc $flags -c "$root/testsuites/hsex/fortran_impulse_replay.f90" -o replay_impulse.o
$fc $flags -c "$root/src/xc/hsex_frozen_q.f90" -o frozen_q.o
$fc $flags -c "$root/src/xc/hsex_density_gamma.f90" -o density_gamma.o
$fc $flags -c "$root/testsuites/hsex/fortran_density_gamma_mpi.f90" -o test_density_gamma.o
$fc $flags -c "$root/src/xc/hsex_frozen_full.f90" -o frozen_full.o
${CC:-cc} -c "$root/src/xc/hsex_memory.c" -o memory.o
$fc $flags -c "$root/testsuites/hsex/fortran_frozen_q_mpi.f90" -o test_q.o
$fc $flags -c "$root/testsuites/hsex/fortran_density_gamma_full.f90" -o test_density_full.o
$fc $flags -c "$root/testsuites/hsex/fortran_frozen_full_mpi.f90" -o test_full.o
python3 "$root/tools/check_forbidden_eigensolvers.py" --object-tree "$out"
$fc $flags math.o window.o frozen.o test_window.o $blas -o test_frozen_window
$fc $flags math.o frozen.o window.o test.o $blas -o test_frozen_response
$fc $flags math.o frozen_kinetic.o test_frozen_kinetic.o $blas -o test_frozen_kinetic
$fc $flags math.o test_grid_translation.o $blas -o test_grid_translation
$fc $flags math.o test_q2_head.o $blas -o test_q2_head
$fc $flags math.o comm.o nvtx.o impulse.o test_impulse.o $blas -o test_impulse_resolvent
$fc $flags math.o comm.o nvtx.o impulse.o replay_impulse.o $blas -o replay_impulse_operator
$fc $flags math.o frozen.o window.o comm.o nvtx.o equilibrium.o test_equilibrium.o $blas -o test_equilibrium_q
$fc $flags math.o projected_derivative.o test_projected_derivative.o $blas -o test_projected_derivative
$fc $flags math.o frozen.o window.o longwave.o coherent_optical.o test_coherent_optical.o $blas -o test_coherent_optical
$fc $flags math.o frozen.o window.o time.o $blas -o test_frozen_time
$fc $flags math.o frozen.o window.o frozen_q.o comm.o nvtx.o memory.o test_q.o $blas -o test_frozen_q
$fc $flags math.o frozen.o window.o longwave.o frozen_limit.o test_limit.o $blas -o test_frozen_limit
$fc $flags math.o frozen.o window.o longwave.o frozen_limit.o frozen_full.o density_gamma.o comm.o nvtx.o test_full.o memory.o $blas -o test_frozen_full
$fc $flags math.o density_gamma.o comm.o nvtx.o test_density_gamma.o $blas -o test_density_gamma
$fc $flags math.o frozen.o window.o longwave.o frozen_limit.o frozen_full.o density_gamma.o comm.o nvtx.o test_density_full.o memory.o $blas -o test_density_full
shasum -a 256 test_density_full > validated.sha256
shasum -a 256 test_frozen_window >> validated.sha256
shasum -a 256 test_frozen_kinetic >> validated.sha256
shasum -a 256 test_density_gamma >> validated.sha256
shasum -a 256 test_frozen_response >> validated.sha256
shasum -a 256 test_grid_translation >> validated.sha256
shasum -a 256 test_q2_head >> validated.sha256
shasum -a 256 test_impulse_resolvent >> validated.sha256
shasum -a 256 test_equilibrium_q >> validated.sha256
shasum -a 256 replay_impulse_operator >> validated.sha256
shasum -a 256 test_frozen_time >> validated.sha256
shasum -a 256 test_frozen_q >> validated.sha256
shasum -a 256 test_frozen_full >> validated.sha256
shasum -a 256 test_frozen_limit >> validated.sha256
shasum -a 256 test_coherent_optical >> validated.sha256
shasum -a 256 test_projected_derivative >> validated.sha256
shasum -a 256 -c validated.sha256
./test_frozen_window
./test_frozen_kinetic
./test_frozen_response
./test_grid_translation
./test_q2_head
./test_frozen_time
./test_frozen_limit
./test_coherent_optical
./test_projected_derivative
