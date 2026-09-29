# Multi-k cuFFT Implementation Plan

**Goal:** Add an optional GPU path for the existing distributed multi-k exchange convolution, and prepare a CPU/GPU Si 8³ timing comparison for MIYABI-G.

**Architecture:** Preserve MPI density-tile transposes and CPU source/target BLAS. A resident backend uploads fixed kernel/geometry/layout once, then packs received density tiles, transforms the k mesh, multiplies the discrete kernel, inverse transforms and unpacks entirely on device. CPU remains the reference/default. No multi-k MLWF approximation is introduced.

**Tech Stack:** Fortran, OpenACC, cuFFT, FFTW, MPI, NVHPC, PBS.

## Tasks

1. Add a neutral owning multi-k backend interface in `src/xc/exx_k_backend.f90`; implement `src/xc/exx_k_cufft.f90` with CPU unavailable stub, shape/overflow validation, persistent plans/data, explicit release and fatal MPI handling. Test direct convolution parity and buffer reuse on optional GPU tests.
2. Connect the backend at the distributed convolution boundary without changing density construction, MPI layout, normalization or action BLAS. Validate failures collectively before each subsequent MPI transpose. CPU oracle injection tests must cover tails, uneven partitions and one-rank failures.
3. Expose a distinct multi-k backend namelist control, with lowercase conversion, broadcast, logging and admission checks. Start from the selected Si functional; do not silently route unsupported functional/layout combinations to CPU.
4. Prepare paired CPU/GPU Si inputs, a small-k correctness check, and an 8³ fixed-work timing job. Keep grid, SCF iterations, exchange parameters, MPI and threads matched; record setup/iteration time and separate host/GPU memory. Document restrictions and MIYABI module/build steps.
5. Build MPI and serial CPU variants, run relevant hybrid regressions and injected backend tests, review changes. Real GPU tests remain opt-in and unverified locally; never claim actual MIYABI timings from CPU/stub execution.

The user has approved adding the multi-k GPU route for CPU/GPU comparison. The user selected HSE06. PBE0/PBEh multi-k Wannier exchange is outside this first extension.

## Verification

- MPI and non-MPI CPU builds passed.
- CPU oracle at the GPU boundary agrees with the original k-mesh exchange for MPI 1/2/3/4, shuffled/shifted k order, uneven k distribution, tail and idle rows, changed sources, zero targets and repeated actions.
- One-rank negative prepare/apply failure and mixed CPU/GPU backend selection are rejected collectively before the next variable-size transpose.
- Independent Python exchange reference passed; CPU input/conventional hybrid RT regressions: 12 passed. Numbered CTest 422–437: 48 stages passed.
- CPU disabled-backend stub tests passed. A small Si 2³-k/6³-grid/2-SCF CPU pilot verified input generation and timing, energy, rank-memory output capture; it is not an 8³ benchmark result.
- GPU direct tests include odd 3³ k transforms, padding, changed constants and reuse counters. Default GPU tests explicitly skip; an explicit GPU request fails because nvfortran is absent locally. GPU compilation, MPI GPU execution and MIYABI timings remain unverified.
- Independent static review completed. Backend-selection handshake and explicit multi-node test placement were corrected during review.
