# Experimental cuFFT backend for compact exact exchange

> 2026-10-06：以下の`benchmarks`は当時のローカル性能測定です。マージ対象から除外し、測定コードはブランチ外へ保存しています。

> 2026-10-06：本記録中の`developer_tests`は当時のローカル開発検証です。GitHub配布から除外しました。通常の回帰試験は`testsuites`を使用します。

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
SCF and fixed-ion native response/pulse RT are admitted; moving-ion
calculations and multi-k GPU exchange are outside this backend's scope.
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

## Validation scope

The standard CPU Si regression does not validate this experimental GPU backend.
Actual device compilation and CPU/GPU numerical comparison are required for the
selected compiler, GPU and MPI layout before using it for production results.

## Full k-mesh hybrid route

The separate `exx_kpoint_backend='cufft'` now accelerates the distributed
k-mesh density convolution. It does not relax the compact Gamma backend's
restrictions. This route requires actual NVHPC/GPU validation before numerical or performance claims.

The common `exx_k_exchange` density-tile engine handles HSE06, PBE0, PBEh(40)
and PBEh(40)+rVV10 with the same MPI layout. The screened or spherical Coulomb
kernel and exchange fraction select the functional. Rectangular orthorhombic
real-space grids and full uniform k meshes are supported; the k backend requires
more than one k point and fully occupied spin pairs. Fractional/empty-state
sources, DC, localized support and Wannier snapshots retain their specialized
CPU routes and cannot select this cuFFT backend. rVV10 is a separate correlation
calculation, not part of the exchange accelerator.

Historical standalone GPU/oracle probes are not part of this proposal.
The submitted Si regression is CPU-only and does not validate the GPU binary.
