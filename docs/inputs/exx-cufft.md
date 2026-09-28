# Experimental cuFFT backend for compact exact exchange

The default remains `exx_local_backend='cpu'`. The optional `cufft` backend
moves the existing compact-support convolution batches to an NVIDIA GPU, in
double-complex precision. It uses the same discrete kernel/filter and selected
pairs as the CPU implementation. Pair selection, spatial MPI communication,
localization and the rest of the hybrid calculation remain on the CPU.

This backend is experimental. The current development machine has no NVHPC
compiler or NVIDIA device: CPU stub behavior and the CPU callback dispatch can
be tested here, but GPU compilation, numerical parity, and performance require
validation on the target GPU system. No GPU speedup or GPU parity is claimed
from those CPU checks.

## Build and input

Use an independent build directory and the NVIDIA HPC Fortran compiler, together
with a compatible MPI installation and the normal hybrid dependencies (FFTW,
Libxc and BLAS/LAPACK). Add your site's usual dependency options to, for example:

```sh
cmake -S . -B build-exx-cufft \
  -DCMAKE_Fortran_COMPILER=nvfortran -DCMAKE_C_COMPILER=nvc \
  -DUSE_SCALAPACK=ON -DUSE_MPI=ON -DUSE_HSE=ON -DUSE_EXX_CUFFT=ON -DUSE_OPENACC=OFF
cmake --build build-exx-cufft
```

`USE_EXX_CUFFT` defaults to `OFF`. It requires `USE_HSE=ON` and NVIDIA HPC
Fortran. Keep the generic `USE_OPENACC` path disabled: this is a separate compact
exchange backend, not a switch to the generic accelerated SALMON calculation.

Add to `&functional` in a supported native hybrid calculation:

```fortran
 exx_local_fft='auto'
 exx_local_backend='cufft'
 exx_gpu_batch_size=8
 exx_mlwf_norm_fraction=.999d0
 exx_mlwf_radius=0d0
```

`exx_local_backend` is case-insensitive and accepts `cpu` or `cufft`.
`exx_gpu_batch_size` is a positive integer, default 8, bounding each spatial
worker's selected target columns per backend call. A larger batch increases
GPU work storage; it does not change the retained-norm target or pair selection.
The `.999` support fraction is an example, not a certified spectral accuracy.
Convergence with respect to support and all other numerical settings remains
necessary for a production calculation.

The GPU option requires native hybrid orbitals at unshifted Gamma,
`num_kgrid=1,1,1`, `nproc_k=1`, automatic local FFT selection, zero fixed radius,
and an adaptive retained-norm fraction strictly between zero and one. Static
SCF and fixed-ion native response/pulse RT are admitted; LCFO RT, moving-ion
calculations, and multi-k GPU exchange are outside this backend's scope.
Other native hybrid input restrictions still apply. Selecting `cufft` in a
build without the backend is an error. A nonempty GPU call also errors when no
NVIDIA device is available; it does not silently pass the GPU test on the CPU.

## Current implementation limits

Each nonempty compact batch creates and destroys a cuFFT plan and transfers
its filter and work arrays. Persistent device plans/filter caching are not part
of this first backend, so transfer and planning overhead may dominate small
batches. Performance has not been measured on a GPU.

Only compact local convolutions use this backend. Broad-support/full-grid FFT
fallbacks, ACE construction/application, and rVV10 remain CPU work. A rank with
no selected local pairs does not call cuFFT; selecting the backend alone does
not prove that a particular step executed a GPU convolution. Check the existing
`EXX_ADAPTIVE` local/global pair counters to identify compact work versus
full-grid fallback. Those counters are workload diagnostics, not GPU timings.

Ordinary backend failures are propagated across the participating ranks.
Unrecoverable OpenACC runtime failures use the registered fatal-error handler;
when MPI is active, it aborts the MPI job so peer ranks do not hang in a
collective waiting for the failed rank. This is an explicit failure path, not a
CPU retry or successful GPU validation.

## MPI rank to GPU assignment

Each rank uses its OpenACC-selected NVIDIA device. SALMON does not assign
local MPI ranks to GPUs automatically. Arrange a device binding for every rank
that performs compact exchange. If each rank sees exactly one GPU through
`CUDA_VISIBLE_DEVICES`, set `ACC_DEVICE_NUM=0` for that rank. Otherwise select a
valid local visible device number separately for each rank.

For example, on a Slurm site configured to expose one bound GPU per task:

```sh
srun --ntasks=2 --gpus-per-task=1 --gpu-bind=single:1 \
  env ACC_DEVICE_NUM=0 ./build-exx-cufft/salmon < inputfile
```

Confirm the site's device visibility and MPI binding policy before using that
example. Setting the same `ACC_DEVICE_NUM` for multiple ranks with the same
visible GPU list makes them share a device; this is not automatic distribution.
The host exchanges compact arrays through MPI; GPU-aware MPI is not required
by this backend.

## Verification

CPU stub validation needs GNU Fortran (`GNU_FC` may override `gfortran`):

```sh
python3 testsuites/unit_pbeh_rvv10/test_cufft.py -v
```

This explicitly skips GPU parity. It checks that disabled nonempty work fails,
empty work is a no-op, and invalid shapes/indices are rejected. To additionally
fault-inject the host fatal callback on one MPI rank while another waits in a
barrier (without needing a GPU), set `SALMON_TEST_MPIEXEC=/path/to/mpiexec` when
running the same test. This checks job termination, not registration with the
real OpenACC runtime.

On an NVIDIA system, request the actual device test:

```sh
FFTW_ROOT=/path/to/fftw NVFC=nvfortran \
  python3 testsuites/unit_pbeh_rvv10/test_cufft.py --gpu -v
```

Alternatively set `SALMON_TEST_CUFFT_GPU=1`. An explicitly requested GPU test
fails if the compiler, FFTW dependency or functioning device is absent. FFTW
flags may instead be supplied through both `FFTW_FFLAGS` and `FFTW_LIBS`, or
resolved through `pkg-config fftw3`. `CUFFT_TEST_FFLAGS` accepts additional
NVHPC flags, such as the site's GPU architecture setting. Configure runtime
library paths as required by the site's FFTW installation.

The device fixture compares the production cuFFT batch with the production
FFTW scalar convolution on unequal padded dimensions `3 x 5 x 8`, irregular
periodic support indices, complex sources and targets, a zero column, multiple
batch sizes/tails, repeated calls, and an empty batch. It compares normalized
complex actions at a tolerance of `2e-12`.

The independent CPU MPI callback test exercises the selected-pair batching and
scatter path without requiring CUDA:

```sh
python3 testsuites/unit_hse_wannier/test_spatial_local.py \
  --build /path/to/cpu-hybrid-build --ranks 1 2 4
```

Its callback uses real FFTW transforms and checks skipped/zero targets, partial
batches, conserved pair counters, and collective propagation of backend errors.
It validates CPU dispatch and MPI handling; it does not validate cuFFT itself.
