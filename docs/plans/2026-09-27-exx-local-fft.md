# Local periodic EXX convolution implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Turn compact MLWF source support into reduced exchange FFT volume,
while preserving the existing discrete periodic HSE/PBEh operator exactly.

**Architecture:** The user approved proceeding toward local exchange shared by
optical response and DC-MD, retaining extended-state support. First implement a
standalone local convolution engine and connect the existing Wannier/ACE static
adapter. The inverse FFT of the global multiplier defines the discrete periodic
kernel. Restrict its displacement table to a compact cyclic bounding box and
zero-pad each axis to at least 2*m-1. A local FFT then computes the exact same
potential at the compact source points. Fall back to the existing global FFT
when the padded volume is not smaller. No extra density/pair magnitude cutoff.

**Tech Stack:** Fortran, FFTW, Python unittest, existing MPI adapter.

The current stage is a shared exchange-engine prerequisite, not an activation of
PBEh optical propagation or DC-MD. Their force, moving-basis and restart work
remains necessary. The first local engine caches one box shape and executes
local pairs serially; full-domain pairs retain existing OpenMP/batching. Record
FFT volume and avoid claiming overall linear scaling from FFT-volume reduction.

### Task 1: standalone convolution
- Create `src/xc/exx_local_fft.f90`, `developer_tests/651_hybrid_exchange/wannier/local_fft_probe.f90`
  and `test_local_fft.py`.
- First fail tests for a missing module. Cover wrapped boxes in three axes,
  singleton/full supports, complex pair densities, rectangular meshes, arbitrary
  periodic discrete multipliers, zero values, changing box sizes and origins.
- Compare against an independent direct sum with an explicit inverse DFT kernel.
- Build minimal cyclic intervals by excluding each axis's largest empty gap.
  Embed the kernel displacement table in a padded box; retain the periodic
  global-kernel zero mode. Cache kernel and plans; clean every FFTW resource.

### Task 2: Wannier adapter
- Modify `src/xc/hse_wannier.f90` to prepare a local box for each translated
  nonzero source. Compute and scatter only source-support potentials; leave
  source/target weighting and ACE unchanged. Skip exactly zero pair densities.
- Add opt-out via `exx_local_fft='auto'/'off'` in `&functional` for regression
  and performance comparison. Report local/global pair counts and FFT volumes.
- Register build dependency in `src/xc/CMakeLists.txt` and standalone probes.
- Add whole-operator parity for HSE/PBEh, multi-k translations, sparse and
  extended sources, complex targets, Hermiticity and zero/tiny pair densities.

### Task 3: integration and evidence
- Run native auto/off SCF parity at fixed finite support and DC MPI checks;
  keep legacy/full-support regressions and HSE-disabled build working.
- Add a bounded large-grid/local-support benchmark with measured FFT-volume
  reduction and wall time; do not extrapolate a production scaling claim.
- Update the input contract and plan ledger. Request independent review before
  committing. Preserve all existing finite-radius MD/convergence restrictions.

## Execution ledger

- Implemented the standalone module, exact cyclic-box embedding and fallback.
  A missing-module test was observed first; independent inverse-DFT/direct-sum
  tests now pass, including complex kernels and 1e-100 nonzero densities.
- Connected the native Wannier/ACE adapter with auto/off selection and pair-grid
  accounting. Full-domain worker FFT buffers are allocated only on fallback;
  local plans run serially outside OpenMP. One local box shape is cached.
- Whole-operator HSE/PBEh tests pass for Gamma and multiple k points with compact
  and extended sources. Legacy 10-case Wannier and exact-pair probes pass.
- Native SCF auto/off parity converges at step 30 on the larger hydrogen fixture;
  local pair-grid volume is 54,000 versus 262,144. Measured results are stored in
  docs/results/pbeh40-rvv10/local-fft.json.
- HSE-enabled and HSE-disabled builds succeed. Existing MPI/input/functional
  regressions pass; a further native-local-path regression was added afterwards.
- Independent review found no blocker in normalization, wrapping, lifetimes or
  threading. Its evidence suggestions were addressed: system_clock wall timing,
  changed-origin/cache-reuse tests, and tiny nonzero local densities.
- Remaining scope: local-batch/OpenMP optimization, distributed source storage,
  pair reduction for localized reference states, PBEh optical adapters, and
  moving-fragment/force consistency for DC-MD. No production scaling claim.
