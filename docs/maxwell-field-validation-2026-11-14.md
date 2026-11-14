# Field components (j→A): validation against SALMON, 2026-11-14

Code: `salmonpy/maxwell.py`. Reference: pinned upstream SALMON `3322d05` with the
SALMONpy runtime patch (`runtime/source.json`), built with GNU Fortran on
x86-64 Linux, OpenBLAS, serial, `-fallow-argument-mismatch`, one OpenMP thread.

## Direct comparison with SALMON

`tools/test_longitudinal_salmon.py SALMON NEW_DIR --trans-longi MODE` propagates
the same synthetic normalized two-orbital Si fixture (8³ grid, one k-point, PZ,
middlepoint, dt = 0.02, 400 steps, Acos2 pulse along z) twice:

1. ordinary SALMON, with the field computed inside SALMON;
2. the SALMONpy TDDFT worker with `trans_longi='tr'`, coupled through
   `salmonpy.maxwell` (`couple` or an explicit loop for the grid check).

| Mode | SALMON input | Largest relative difference |
|---|---|---|
| `lo` | `trans_longi='lo'` | Ac_tot 8.5e-16, matter current 1.7e-15 |
| `2d` | `trans_longi='2d'`, film 10 bohr, ε = (1, 2.25) | Ac_tot 3.1e-15, matter current 1.7e-15 |
| `ms1d` | `multi_scale_maxwell_tddft`, `fdtddim='1d'`, 1 point, h = 5, vacuum (300, 50), Mur | A on all 351 cells at every step 8.5e-16, current 2.3e-15 |

The differences are at the 16-digit precision of SALMON's text output. The
synthetic orbital carries a nonzero current, so the induced field dominates the
total field in `lo` and `ms1d`; the comparison is therefore sensitive to the
feedback sign, the time alignment and the initial conditions.

Notes found while setting up the comparison:

- SALMON prints `rt.data` rows for itt = 1..nt (t = dt first); align by time.
- `film_thickness` and `epsilon_em` belong to `&maxwell`, not `&emfield`.
- The multiscale input checker requires the legacy `nxvacl_m`/`nxvacr_m`
  (≥ 2; the error message reads "must not larger than 2"). Then
  nxvac = (|nxvacl_m|+1, |nxvacr_m|).
- In 1D, x is the propagation axis; an x-polarized incident wave does not
  propagate. `Multiscale1D` rejects it.

## Unit tests (`salmonpy/tests/test_maxwell.py`)

Analytic checks: leapfrog plasma-oscillation period, thin-film sheet
transmission −iω/(−iω+γ), dispersion-free 1D vacuum propagation at c·dt = h,
one-cell 1D layer versus the thin-film model; bit-identity of `Transverse`
coupling with `crystal.respond` (bundled SBE kernel) and of `Longitudinal`
coupling with SALMON's recurrence; atomic failure, capacity, closing of
matter states and output archives.

## Charge-current contract and E/H timing (same day, later)

The `crystal` contract now returns the electric (charge) current density j.
`SBEState`/`TDDFTState` keep SALMON's matter current J internally
(`matter_current`) and declare `timing` with output j. `TDDFTResponseBackend`
no longer negates. `crystal.respond` and `maxwell.couple` resolve declarations
through `salmonpy.timing.matter_timing` (crystal default plus a warning when a
backend has none). `NativeEHField.timing` declares j and A at offset 1/2;
`CoupledStepper` gathers matter declarations over MPI, shifts them by the
electronic clock offset 1/2 and checks them before the first step.

Environment additionally used: OpenMPI 4 and mpi4py from Ubuntu/PyPI, Libxc 6
(Ubuntu), oversubscribed on 2 cores, run as root (`OMPI_ALLOW_RUN_AS_ROOT`).

| Check | Result |
|---|---|
| `salmonpy/tests` (unit) | 133 run, all pass (6 skipped: external runtime) |
| `tools/test_tddft_entrypoint.py` | worker A→J matches baseline; batch, interleaved, rejection, recovery pass |
| `tools/test_sbe_entrypoint.py` | 16 full-application configurations; Python currents agree |
| `tools/run_coupled_mpi_tests.py` (embedded MPI library) | all 13 cases and the feedback/rank-parity check pass |
| `mpiexec -n 3 tools/test_tddft_distributed_mpi.py` | four-point TDDFT charge current, repeats, serial equality pass |
| `tools/run_parallel_mpi_tests.py` (host-MPI bridge) | custom, LD and sample pass |
| `tools/test_longitudinal_salmon.py` lo / 2d / ms1d | unchanged: 8.5e-16 / 3.1e-15 / 8.5e-16 (A) |

Not run: `tools/test_tddft_outputs.py` (needs a converged GS bundle) and
`tools/test_sbe_native_regression.py` (regenerates the fixture from the
pre-extraction source; only its sign handling changed).

Build issue found and fixed: `src/rt/embedded_tddft.f90` in
`runtime/patches/salmonpy-runtime.patch` had three lines longer than 132
characters, which GNU Fortran on Linux rejects (`-Werror=line-truncation`).
They were split with continuation lines only; with spaces, newlines and `&`
removed the file is token-identical. The patch was regenerated from a clean
checkout of the pinned commit (the method reproduces the old patch
byte-for-byte), `runtime/source.json` now records patch SHA256
`0f9faa9db8365cad1ecb457e5672d43c8618b440a6eaac63b3ee782229f4ae7f`
(previously `14b82f79...44f7`). `tools/test_runtime_preparation.py` passes; the
embedded MPI library builds without `-ffree-line-length-none`; the coupled MPI
suite passes again and its 64-step sample history is bit-identical to the
library built from the old patch.
