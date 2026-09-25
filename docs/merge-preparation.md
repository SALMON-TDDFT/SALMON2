# HSE merge preparation — 2026-09-25

## Scope and base

Candidate branch: `hse-merge-prep`, based on upstream `develop-2.0.0` at
`ea309eed`. Extracted from TDCDFT at `f946904d` plus the pending HSE laser,
screening-input/unit conversion and default-propagator changes. The original
TDCDFT worktree is preserved. No TDCDFT runtime, input or ELF implementation is
included, and no large trajectory/spectrum archives are imported.

Includes native screened exchange, ACE, distributed Bloch/symmetry processing,
OpenMP/BLAS/FFT optimizations, short-pulse support, restart guards, user examples,
and bounded native numerical/integration tests. PT-CN and full-action Taylor
remain explicit developer paths. Ordinary HSE RT needs no propagation namelist.

The newer upstream DC-LCFO input validation is retained. A generic RT fix uses
iteration relative to the loaded checkpoint for wavefunction buffer alternation.
The HSE metadata writer now removes stale metadata on non-HSE checkpoint writes.

## Verification actually run on the extracted tree

- HSE ON / MPI ON: complete Release build, GNU Fortran 15, OpenMPI,
  FFTW, Libxc and OpenBLAS on Apple Silicon.
- HSE OFF / MPI OFF: complete serial CPU Release build.
- Eight native numerical tests passed, including MPI decomposition, arbitrary
  targets, symmetry reconstruction, ACE, developer CN and Libxc screening.
- Thirteen bounded input/restart cases passed: omitted versus explicit Taylor,
  atomic versus A_eV_fs screening input/defaults, custom omega, invalid settings,
  same-physics restart, changed-omega rejection and serial PZ odd-step restart.
  The PZ restarted endpoint current and total energy exactly match the uninterrupted
  two-step run. Stale HSE metadata removal passed.
- Full-k short Acos2 + impulse-probe: 16 steps completed. Both direct pulse restart
  and reopening its checkpoint as an impulse were rejected as intended.
- Standard CTest cases 420 and 421: all six prep/run/verify stages passed on MPI4.
  These use a one-iteration GS producer and two-step RT consumer for integration
  coverage; they are deliberately not convergence or physical-spectrum tests.
- `git diff --check` passed. Read-only review of both pending HSE changes and
  extracted source found no remaining concrete source correctness blocker.

The pre-extraction branch's broader 83-test reference suite also passed, but
that result is not substituted for tests on the extracted tree.

Build reproduction (adjust dependency paths for your platform):

```sh
cmake -S . -B /tmp/salmon-hse-build -DUSE_HSE=ON -DUSE_MPI=ON \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=mpifort \
  -DCMAKE_PREFIX_PATH='/path/openblas;/path/fftw;/path/libxc'
cmake --build /tmp/salmon-hse-build -j 8
OMP_NUM_THREADS=1 ctest --test-dir /tmp/salmon-hse-build \
  -R '(420_bulk_Si_hse_gs|421_bulk_Si_hse_rt)' --output-on-failure
SALMON_TEST_MPI=1 OMP_NUM_THREADS=1 python3 -m unittest discover \
  -s samples/hse_mlwf_reference -p 'test_native*.py'
```

This Mac build additionally used existing toolchain workarounds
`CMAKE_C_FLAGS=-include stdio.h` and
`CMAKE_Fortran_FLAGS=-fallow-argument-mismatch`, with gcc-15 as C compiler.
Compiler warnings in legacy MPI argument interfaces remain; no new warning-free
or cross-platform claim is made.

The input/restart harness is `samples/hse_native/check_input_smoke.py --help`;
it accepts an existing compatible converged Si8 GS and retains every input/log.
The short pulse harness is `samples/hse_mlwf_reference/validate_native_pulse.py`.

## Remaining before final integration

- Maintainer review of the extracted HSE API, allowed scope and numerical method.
- Transfer the input documentation in `docs/inputs/hse.md` into SALMON-DOCS,
  as required by CODING_RULES.md; that separate repository is not modified here.
- Platform CI beyond this Mac, especially target compiler/MPI/BLAS combinations.
  Fugaku performance, accelerators and unsupported physical modes are not certified.
- Full upstream test-suite execution is not claimed. The new integration cases,
  native tests and representative HSE-disabled restart checks are the current evidence.

No merge, force push or deletion of the research branch has been performed.
