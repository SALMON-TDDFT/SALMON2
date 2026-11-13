# TDDFT connection validation

## Scope and ownership

The native adapter supports periodic, unpolarized, fixed-ion PZ with middlepoint
propagation, CPU execution and one MPI rank (OpenMP configurable). It holds one
persistent SALMON process per state; it is not a reentrant shared-library handle.
The input bundle is snapshotted and hashed; each state owns independent input
copies plus all native wavefunctions, potential histories and module globals.

`worker_configure` is called only when SALMONPY_TDDFT_WORKER=1; ordinary input
and update paths retain their previous behavior. Existing initialization and
propagation routines execute in place, including predictor-corrector and
nonlocal-phase conventions. Public A means total Ac=A/c, and returned J is the
SALMON electronic matter-current density. Electronic charge current is -J.

## Local acceptance

macOS arm64, GNU Fortran 15, Accelerate, OMP_NUM_THREADS=1:

- Before/after ordinary current and energy files are byte-identical, with
  predictor-corrector disabled and enabled.
- Interactive and batch Python currents match ordinary SALMON at rtol=1e-12,
  atol=1e-14. Interleaved independent states reproduce the same accepted sequence.
- Invalid A does not advance state. Killing a worker between steps raises an
  error, reconstructs the previous accepted state and permits an identical retry.
- Nonzero initial A with all three components varying matches the original
  ordinary SALMON current at the same tolerance.
- Closing the initial loader leaves already-created states usable.
- Unsafe output stems and external k-grid paths are rejected before initialization.
- Python suite: 34 tests. Full source/object prohibition audits pass.

Before executable SHA256:
`32f9123cd400aab16750e6ee4c16253c26f257a99824d25dd957379088f345f7`

After executable SHA256:
`08d7e7e127d81f8dd35886f2d607135db2264e33c69066bdf6876136e83628a4`

Reproduce with audited before/after executables:

```
PYTHONPATH=. python tools/test_tddft_entrypoint.py BEFORE AFTER NEW_DIRECTORY
SALMONPY_TEST_PC=y PYTHONPATH=. python tools/test_tddft_entrypoint.py BEFORE AFTER OTHER_DIRECTORY
```

The fixture uses two analytically normalized orbitals and a repository-supplied
Si pseudopotential on an 8×8×8 periodic grid. It is intentionally **not a converged
GS** and does not establish material accuracy. No external GS dataset was collected
or ground-state solver run. At this initial synthetic-test stage, a real GS regression was pending;
the subsequently authorized local calculation is recorded below. binary compatibility is inherited from SALMON's
existing loader and is not claimed across arbitrary compilers/architectures.

AETRS/HSE, other XC functionals, MPI decomposition, GPUs, longitudinal internal
feedback, moving ions, spin-orbit, +U and excited-state checkpoint loading are not
part of this validation. No new unitarity/adaptive-step claim is made.


## CI acceptance

[Native TDDFT run 43363609508](https://github.com/SALMON-TDDFT/SALMON2/actions/runs/43363609508)
passed on Ubuntu 22.04 with GNU Fortran/OpenBLAS at commit `514e3030`.
It builds and audits both the pinned pre-worker baseline (`f419cb79`) and the
new solver, then runs the full regression with predictor-corrector off/on.
Both report ordinary output byte equality and passing live-state/recovery tests.

[Wheel run 43363121675](https://github.com/SALMON-TDDFT/SALMON2/actions/runs/43363121675)
passed at Python/native-LD-SBE package commit `e837ebac`: 4 platforms × 5 CPython
versions × 34 tests = 680 passing installed-package tests, plus the sdist build.
The full TDDFT executable is **not** part of these wheels. Native TDDFT was
executed on Linux CI and local macOS arm64, not Windows or Intel macOS.

Linux integration setup exposed two portability defects: raw LAPACK vendor flags
were placed before objects (fixed by CMake library discovery), and an early stdio
include disabled nftw feature declarations (fixed with POSIX.1-2008 feature macros
before system headers, also validated by the local macOS build).


## Authorized local GS-to-response run

On 2026-11-13 the user requested a full local Mac calculation. A new Si GS was
computed from the repository exercise geometry/pseudopotential with PZ,
12×12×12 real-space grid, 2×2×2 k points and 32 states. No existing user GS was
modified. OpenBLAS was selected after a minimal ZDOTC diagnostic isolated a
GNU Fortran/Accelerate complex-return incompatibility; no solver math changed.

- SCF converged at iteration 37, residual 6.0476828e-10 (threshold 1e-9).
- TDDFT: 1000 steps, dt=0.05 a.u., 1.209442 fs; direct/native versus Python
  maximum absolute J difference 5.014435047745458e-19 a.u.
- SBE from actual GS text exports, including nonlocal correction: 1000 steps;
  maximum absolute J difference 3.6930636500287495e-17 a.u.
- Two TDDFT states driven by opposite A: 128 interleaved steps, each agrees
  with its separate trajectory; neither mutates the other's current.
- LD: 501 samples, batch/step maximum J difference exactly zero.
- FDTD: ellipsoid and box, 200 steps each; direct/wrapped observations and
  material-grid output files byte-identical.

Native execution: OpenMP 2 threads, OpenBLAS 1 thread. Executable SHA256:
`17aec5e4bb4abf2b2a4031e1a098b88d3abc5776c9c6f9f2dab91f4486c2b785`.

Inputs, outputs, GS bundle, executable, JSON summary and comparison plot are in:
`/Users/otobetoshihito/.codex/visualizations/2026/11/13/01a1d275-1ce2-7121-b8f4-cfd9dbc2d897/salmonpy-local-mac-20261113`.

These are functional/regression results. Grid, k-point and physical-observable
convergence were not assessed; this short pulse is not a material benchmark.

Environment note: absolute paths in this report identify the historical validation host. They are not portable setup defaults or required user paths.
