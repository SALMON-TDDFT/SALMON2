# TDCDFT implementation validation — 2026-09-23

Repository: SALMON-TDDFT/SALMON2, branch `TDCDFT`, base `ae7d8953` (`origin/develop-2.0.0`).
This is an implementation/smoke-test record, not a converged Si exciton calculation.

## Implemented

Opt-in PZ + head-only LRC and normalized Proca-type field equation; midpoint propagation;
consistent endpoint kinetic/nonlocal current; initial half step; separate xc output;
versioned shared-file RT restart with current history; parameter/compatibility validation.
The classical Ac_tot history and legacy rtdata.bin layout remain unchanged.
GPU and other unsupported physics/propagation paths are rejected explicitly.

The new model is documented in `samples/exercise_si_tdcdft/README.md`.
SALMON-DOCS is a separate repository: `docs/salmon-docs-tdcdft.patch` contains the matching
change to `source/input_keyword_list.rst`, made against a separate checkout.
No remote repository was modified.

## Build environment

macOS ARM64, Homebrew gfortran/GCC 15, OpenBLAS, CMake; serial and Open MPI builds.
An unmodified base executable was built before source changes for baseline comparison.
Legacy source compatibility needed `-fallow-argument-mismatch` (serial communication stubs)
and `-include stdio.h` (posix.c uses snprintf without the corresponding include).
These were build flags, not unrelated source edits.

Serial configuration used:

```sh
cmake -S . -B /private/tmp/salmon-tdcdft-work \
  -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-15 \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/gfortran \
  -DCMAKE_BUILD_TYPE=Release -DUSE_MPI=OFF \
  -DCMAKE_Fortran_FLAGS=-fallow-argument-mismatch \
  '-DCMAKE_C_FLAGS=-include stdio.h' \
  -DBLA_VENDOR=OpenBLAS -DCMAKE_PREFIX_PATH=/opt/homebrew/opt/openblas
cmake --build /private/tmp/salmon-tdcdft-work -j 8
```

MPI used `USE_MPI=ON`, `mpifort`, and `OMPI_CC=/opt/homebrew/bin/gcc-15` for both
configuration and build. The latter avoids the Apple Clang wrapper rejecting OpenMP.
Local MPI execution needed permission for process-communication sockets outside the sandbox.

## Tests and evidence

- Unmodified baseline: tests 111 (Si GS) and 112 (Si response), six CTest steps passed.
- Modified serial and four-process MPI: tests 111, 112, 142, analytic field test;
  ten CTest steps passed per build.
- Analytic field update: constant drive, zero coupling, smooth sinusoidal drive,
  damped oscillator with second-order timestep convergence.
- Si integration regression: disabled output matches unmodified executable;
  enabled alpha=0 impulse output is identical to disabled; LRC and Proca feedback are nonzero.
- Restart at step 90 of 180 reproduces current and the full response spectrum.
  Missing/incompatible restart state is guarded by the reader; changed-parameter rejection tested.
- Weak impulse 0.001 vs 0.0005: normalized-current relative difference
  approximately 8.11e-5. dt=0.16 vs 0.08: relative current difference approximately 3.54e-4.
  These short, coarse-grid checks demonstrate numerical consistency, not physical convergence.
- Laser alpha=0 vs alpha=1e-10: continuous current limit verified.
- Negative/NaN coefficients, unknown mode, unsupported propagator and changed restart coefficient rejected.
- Full integration regression also passed with two MPI processes.
- Analysis tests: damped-oscillator peak and amplitude scaling, pump subtraction and time-origin shift,
  rejection of inconsistent time grids/zero impulse/empty peak intervals (three Python tests).
- Python Fourier transform matches SALMON's complex epsilon on the Si test within rtol=1e-10, atol=1e-9.
- Six input templates exercised with shortened 180-step runs; pump/probe timing shortened to include the
  probe. Differential current vanishes before the probe, is nonzero afterward, and gives finite spectra.
- Independent read-only code review found a missing initial half step and a laser zero-coupling
  discontinuity; both fixed and covered by tests. Re-review found no further definite correctness issue.

CTest fixture setup originally omitted the execution fixture used by verification; adding the new
GS-dependent test exposed out-of-order verification. The shared helper now registers that fixture,
so the new RT test actually depends on completed GS execution. No reference tolerances were weakened.

Representative commands (paths can be changed):

```sh
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 ctest --test-dir BUILD \
  -R '111_bulk_Si|112_bulk_Si|142_Si|unit_tdcdft' --output-on-failure
python3 testsuites/142_Si_tdcdft/test_update.py
python3 testsuites/142_Si_tdcdft/test_regression.py EXE BASELINE_GS_DIR BASELINE_RT_DATA
# MPI variant appends: /path/to/mpiexec -n 2
python3 samples/exercise_si_tdcdft/test_analyze.py
```

The legacy verification scripts require a `python` command with NumPy; locally this was supplied
by a temporary shim to python3. The new field-unit test only needs Python and gfortran.

## Limits and next scientific work

No fitted Si Proca parameter set, converged absorption curve, E1/E2 shift, bleaching, or exciton
binding energy has been established. The sample gamma=0.01 a.u. is a numerical demonstration.
The small 12³ real-space / 2³ k grid and 180-step tests cannot support such claims.
Full sample-duration runs and systematic k/grid/dt/time scans remain to be performed.

For lasers, enabled alpha=0 is the comparison baseline: the new path uses endpoint nonlocal phases,
whereas the legacy path uses midpoint phases in that part of the current. Impulse comparisons are
identical because the classical A is constant. The model-field energy is not added to total-energy
output. Zero-field model stability does not guarantee stability of the coupled electronic system.

Before interpreting pump-induced peaks, verify a no-pump delayed probe against equilibrium response,
probe-amplitude convergence, pump-only subtraction, observation-window convergence and coefficient
sensitivity. Fixed alpha does not describe dynamically evolving screening. Test against literature
or BSE with consistent quasiparticle-gap treatment before drawing quantitative conclusions.


## Follow-up: fixed 4³ k grid

The user removed the k-point convergence requirement. The examples now use 4³ k points,
and completed 12000-step ALDA/LRC runs and their initial analysis are recorded in
`docs/results/si-k4/README.md`. The earlier 2³-grid entries above remain the implementation-test record.
