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

The compact plan owns resident GPU data for one exchange-action call. The
source tile is gathered once per source orbital, before the target-batch loop.
Support indices, the Fourier kernel and the source remain on the GPU; prepare
uploads them only when their contents change. Each batch uploads its target
columns and downloads its action columns. Pair-density generation, forward FFT,
kernel multiplication, inverse FFT and action assembly stay on the GPU.

The cuFFT plan and work buffers are reused when padded dimensions, support size,
batch capacity and selected device match. Tail batches zero unused columns and
use the same plan. A geometry/capacity change rebuilds the resources. Explicit
release at the end of the exchange action frees them; this cache does not yet
persist across separate Hamiltonian applications or RT steps. CPU-side MPI
collectives still require per-batch target/result transfers. Performance has
not been measured on a GPU.

Only compact local convolutions use this backend. Broad-support/full-grid FFT
fallbacks, ACE construction/application, and rVV10 remain CPU work. A rank with
no owned target columns does not prepare the GPU backend; selecting the backend alone does
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
python3 developer_tests/653_functional/test_cufft.py -v
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
  python3 developer_tests/653_functional/test_cufft.py --gpu -v
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
complex actions at a tolerance of `2e-12`. The stateful fixture also checks
operation counters: repeated prepare leaves uploads unchanged, source/filter/
index changes upload only that data, and tail batches reuse the plan. Capacity
and geometry changes, repeated release, and owner finalization are covered.

The independent CPU MPI callback test exercises the selected-pair batching and
scatter path without requiring CUDA:

```sh
python3 developer_tests/651_hybrid_exchange/wannier/test_spatial_local.py \
  --build /path/to/cpu-hybrid-build --ranks 1 2 4
```

Its callback and stateful oracle use real FFTW transforms and check skipped/zero
targets, partial batches, conserved pair counters, one prepare per source,
changed-source results, and collective propagation of prepare/apply errors.
A negative release failure on one rank is tested during reinitialization; the
pair-screen fixture also tests release failures at the exchange-action boundary.
It validates CPU dispatch and MPI handling; it does not validate cuFFT itself.

## Full k-mesh hybrid route

The separate `exx_kpoint_backend='cufft'` now accelerates the distributed
k-mesh density convolution. It does not relax the compact Gamma backend's
restrictions. See [Si CPU/GPU tests](../../benchmarks/benchmark_si_kpoint_cufft/README.md)
for supported inputs, resident data lifetime and MIYABI-G jobs. This route also
requires actual NVHPC/GPU validation before numerical or performance claims.

The common `exx_k_exchange` density-tile engine handles HSE06, PBE0, PBEh(40)
and PBEh(40)+rVV10 with the same MPI layout. The screened or spherical Coulomb
kernel and exchange fraction select the functional. Rectangular orthorhombic
real-space grids and full uniform k meshes are supported; the k backend requires
more than one k point and fully occupied spin pairs. Fractional/empty-state
sources, DC, localized support and Wannier snapshots retain their specialized
CPU routes and cannot select this cuFFT backend. rVV10 is a separate correlation
calculation, not part of the exchange accelerator.

Run `python3 developer_tests/651_hybrid_exchange/k_exchange/run.py` for MPI1/2/3/4 comparisons
against the retained full-support Wannier reference. `661_unit_kpoint_cufft/run.py`
also exercises rectangular and global-Coulomb convolution; add `--gpu` only on
an NVHPC/NVIDIA machine. GNU oracle/stub tests do not validate the GPU binary.
