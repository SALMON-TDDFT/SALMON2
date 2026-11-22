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
fftw_flags=${SCREEN_FFTW_FLAGS:-}
fftw_libs=${SCREEN_FFTW_LIBS:--lfftw3}
if test -z "$fftw_flags" && test -f /opt/homebrew/opt/fftw/include/fftw3.f03; then
  fftw_flags=-I/opt/homebrew/opt/fftw/include
  fftw_libs="-L/opt/homebrew/opt/fftw/lib -lfftw3"
fi
# Compiler flag variables are intentionally word-split, as in build tooling.
$fc $flags -c "$root/src/common/hsex_update_schedule.f90" -o update_schedule.o
$fc $flags -c "$root/src/misc/nvtx_wrapper.f90" -o nvtx.o
$fc $flags -c "$root/src/parallel/communication.f90" -o comm.o
$fc $flags $native_flags -c "$root/tools/hsex/native/q_response.f90" -o q.o
$fc $flags -c "$root/src/xc/hsex_screen_math.f90" -o math.o
$fc $flags -c "$root/src/xc/hsex_frozen_window.f90" -o window.o
$fc $flags -c "$root/src/xc/hsex_frozen_response.f90" -o frozen.o
$fc $flags -c "$root/src/xc/hsex_frozen_kinetic.f90" -o frozen_kinetic.o
$fc $flags -c "$root/testsuites/hsex/fortran_frozen_kinetic.f90" -o test_frozen_kinetic.o
$fc $flags -c "$root/src/xc/hsex_impulse_resolvent.f90" -o impulse.o
$fc $flags -c "$root/src/xc/hsex_equilibrium_q.f90" -o equilibrium.o
$fc $flags -c "$root/src/xc/hsex_projected_derivative.f90" -o projected_derivative.o
$fc $flags -c "$root/src/xc/hsex_longwave.f90" -o longwave.o
$fc $flags -c "$root/src/xc/hsex_coherent_optical.f90" -o coherent_optical.o
$fc $flags -c "$root/src/xc/hsex_frozen_q.f90" -o frozen_q.o
$fc $flags -c "$root/src/xc/hsex_frozen_limit.f90" -o frozen_limit.o
$fc $flags -c "$root/src/xc/hsex_density_gamma.f90" -o density_gamma.o
$fc $flags -c "$root/testsuites/hsex/fortran_density_gamma_mpi.f90" -o test_density_gamma.o
$fc $flags -c "$root/src/xc/hsex_frozen_full.f90" -o frozen_full.o
$fc $flags -c "$root/src/xc/hsex_projection.f90" -o projection.o
$fc $flags -c "$root/src/xc/hsex_pdep.f90" -o pdep.o
$fc $flags -c "$root/src/xc/hsex_sternheimer.f90" -o sternheimer.o
$fc $flags -c "$root/src/xc/hsex_response_space.f90" -o response_space.o
$fc $flags -c "$root/src/xc/hsex_response_symmetry.f90" -o response_symmetry.o
$fc $flags -c "$root/src/xc/hsex_reduced_q.f90" -o reduced_q.o
$fc $flags -c "$root/src/xc/hsex_pdep_q.f90" -o pdep_q.o
${CC:-cc} -c "$root/src/xc/hsex_memory.c" -o memory.o
$fc $flags -c "$root/src/xc/hsex_pdep_cache.f90" -o pdep_cache.o
$fc $flags $fftw_flags -c "$root/src/xc/hsex_gamma_pdep.f90" -o gamma_pdep.o
$fc $flags -c "$root/src/xc/hsex_screen_driver.f90" -o driver.o
for name in driver head pipeline refine; do
  $fc $flags -c "$root/testsuites/hsex/fortran_screen_${name}.f90" -o "test_${name}.o"
done
$fc $flags -c "$root/testsuites/hsex/fortran_projection_capture.f90" -o test_projection.o
$fc $flags -c "$root/testsuites/hsex/fortran_pdep.f90" -o test_pdep.o
$fc $flags -c "$root/testsuites/hsex/fortran_pdep_q_mpi.f90" -o test_pdep_q.o
$fc $flags -c "$root/testsuites/hsex/fortran_q_batch.f90" -o test_q_batch.o
$fc $flags -c "$root/testsuites/hsex/fortran_screen_forms.f90" -o test_forms.o
$fc $flags -c "$root/testsuites/hsex/fortran_sternheimer_static.f90" -o test_static.o
$fc $flags -c "$root/testsuites/hsex/fortran_longwave.f90" -o test_longwave.o
$fc $flags -c "$root/testsuites/hsex/fortran_gamma_optical.f90" -o test_gamma_optical.o
$fc $flags -c "$root/testsuites/hsex/fortran_gamma_pdep_mpi.f90" -o test_gamma_pdep.o
$fc $flags -c "$root/testsuites/hsex/fortran_dielectric_inverse.f90" -o test_inverse.o
$fc $flags -c "$root/testsuites/hsex/fortran_response_space.f90" -o test_response_space.o
$fc $flags -c "$root/testsuites/hsex/fortran_response_star_mpi.f90" -o test_response_star.o
$fc $flags -c "$root/testsuites/hsex/fortran_response_symmetry.f90" -o test_response_symmetry.o
$fc $flags -c "$root/testsuites/hsex/fortran_reduced_q_mpi.f90" -o test_reduced_q.o
$fc $flags -c "$root/testsuites/hsex/fortran_update_schedule.f90" -o test_update_schedule.o
$fc $flags -c "$root/testsuites/hsex/fortran_frozen_window.f90" -o test_window.o
$fc $flags -c "$root/testsuites/hsex/fortran_frozen_response.f90" -o test_frozen_response.o
$fc $flags -c "$root/testsuites/hsex/fortran_frozen_time.f90" -o test_frozen_time.o
$fc $flags -c "$root/testsuites/hsex/fortran_frozen_q_mpi.f90" -o test_frozen_q.o
$fc $flags -c "$root/testsuites/hsex/fortran_density_gamma_full.f90" -o test_density_full.o
$fc $flags -c "$root/testsuites/hsex/fortran_frozen_full_mpi.f90" -o test_frozen_full.o
$fc $flags -c "$root/testsuites/hsex/fortran_frozen_limit.f90" -o test_frozen_limit.o
$fc $flags -c "$root/testsuites/hsex/fortran_coherent_optical.f90" -o test_coherent_optical.o
$fc $flags -c "$root/testsuites/hsex/fortran_projected_derivative.f90" -o test_projected_derivative.o
$fc $flags -c "$root/testsuites/hsex/fortran_grid_translation.f90" -o test_grid_translation.o
$fc $flags -c "$root/testsuites/hsex/fortran_q2_head.f90" -o test_q2_head.o
$fc $flags -c "$root/testsuites/hsex/fortran_impulse_resolvent.f90" -o test_impulse.o
$fc $flags -c "$root/testsuites/hsex/fortran_equilibrium_q.f90" -o test_equilibrium.o
python3 "$root/tools/check_forbidden_eigensolvers.py" --object-tree "$out"
$fc $flags math.o window.o frozen.o test_window.o $blas -o frozen_window
$fc $flags test_driver.o math.o $blas -o assembly
$fc $flags test_head.o math.o $blas -o head
$fc $flags test_pipeline.o driver.o longwave.o sternheimer.o response_space.o response_symmetry.o math.o q.o comm.o nvtx.o $blas -o pipeline

$fc $flags test_refine.o q.o $blas -o refine
$fc $flags test_projection.o projection.o math.o $blas -o projection_capture
$fc $flags test_q_batch.o q.o $blas -o q_batch
$fc $flags test_forms.o math.o $blas -o forms
$fc $flags test_static.o sternheimer.o math.o $blas -o static_sternheimer
$fc $flags test_longwave.o longwave.o math.o $blas -o longwave
$fc $flags test_inverse.o math.o $blas -o dielectric_inverse

$fc $flags test_pdep.o pdep.o math.o $blas -o pdep

$fc $flags test_pdep_q.o pdep_q.o pdep.o sternheimer.o math.o comm.o nvtx.o memory.o $blas -o pdep_q

$fc $flags test_gamma_optical.o longwave.o math.o $blas -o gamma_optical

$fc $flags test_gamma_pdep.o gamma_pdep.o pdep_cache.o pdep.o longwave.o sternheimer.o math.o comm.o nvtx.o $blas $fftw_libs -o gamma_pdep

$fc $flags test_response_space.o response_space.o math.o $blas -o response_space

$fc $flags test_reduced_q.o reduced_q.o response_symmetry.o response_space.o sternheimer.o math.o comm.o nvtx.o $blas -o reduced_q

$fc $flags test_response_symmetry.o response_symmetry.o math.o $blas -o response_symmetry

$fc $flags test_response_star.o reduced_q.o response_symmetry.o response_space.o sternheimer.o math.o comm.o nvtx.o $blas -o response_star

$fc $flags test_update_schedule.o update_schedule.o -o update_schedule
$fc $flags test_frozen_response.o frozen.o window.o math.o $blas -o frozen_response
$fc $flags math.o frozen_kinetic.o test_frozen_kinetic.o $blas -o frozen_kinetic
$fc $flags test_frozen_time.o frozen.o window.o math.o $blas -o frozen_time
$fc $flags test_frozen_q.o frozen_q.o frozen.o window.o math.o comm.o nvtx.o memory.o $blas -o frozen_q
$fc $flags test_frozen_limit.o frozen_limit.o longwave.o frozen.o window.o math.o $blas -o frozen_limit
$fc $flags test_coherent_optical.o coherent_optical.o longwave.o frozen.o window.o math.o $blas -o coherent_optical
$fc $flags test_projected_derivative.o projected_derivative.o math.o $blas -o projected_derivative
$fc $flags test_grid_translation.o math.o $blas -o grid_translation
$fc $flags test_q2_head.o math.o $blas -o q2_head
$fc $flags test_impulse.o impulse.o math.o comm.o nvtx.o $blas -o impulse_resolvent
$fc $flags test_equilibrium.o equilibrium.o frozen.o window.o math.o comm.o nvtx.o $blas -o equilibrium_q
$fc $flags test_frozen_full.o frozen_full.o density_gamma.o frozen_limit.o longwave.o frozen.o window.o math.o comm.o nvtx.o memory.o $blas -o frozen_full

$fc $flags test_density_gamma.o density_gamma.o math.o comm.o nvtx.o $blas -o density_gamma

$fc $flags test_density_full.o frozen_full.o density_gamma.o frozen_limit.o longwave.o frozen.o window.o math.o comm.o nvtx.o memory.o $blas -o density_full
