# Scalar pseudopotential grid ownership (2026-11-15)

Local Apple Silicon/GNU Fortran/OpenBLAS/Libxc/FFTW checks, synthetic small systems and licensed bundled Si KY input. No CI/network run was started. Build parallelism at most 2, numerical tests 1–2 threads, aggregate worker budget at most 4. Outputs/builds remain outside the repository.

## Implementation and contract

- `s_pp_grid_config` is copied into `ppg` on construction. Legacy calls capture once; explicit scalar calls bypass shared input. Periodicity/Poisson/FFT/diagnostic/SOC/+U selection remains fixed.
- Scalar grid setup reads atom count/species from `system`, species count from PP arrays. Cache, interpolation and projector communication-table helpers no longer read global atom/species counts.
- NLCC placement uses system species and the loaded PP periodicity. The periodic nonlocal k-phase update dispatches using owned SOC/+U settings.
- Position updates retain the original deallocate/rebuild sequence. Atom mapping/count, cell/grid geometry, radial layout and source config must stay fixed. Mismatches are rejected before updating grid arrays.
- Loaded radial numerical tables are immutable by caller contract. Checks cover configuration/layout, not full table contents; there is no per-step full numerical-table copy/hash.

Explicit independence is scalar, without +U. Legacy SOC/+U specialized modules, communicator lifetimes, global timers, and other whole-solver inputs remain outside this extraction. Full TDDFT's common grid/k/XC guard remains. No general native failure atomicity, concurrent calls, MPI isolation or physical convergence is claimed.

## Focused regression

The original cache test failed with `STOP 1` against the previous audited build: a two-atom object's cache was allocated using poisoned global `natom=9`. This failure preceded production edits.

The final probe covers one- and two-atom owners with distinct species mapping/local channels, poisoned atom counts/species/boundaries/FFT/SOC/+U globals, copied config, moved-position update versus fresh creation, and nonzero synthetic NLCC A/B/A comparison. Periodic phase values are checked against the analytic negative-exponential expression. Four subprocess cases require the expected diagnostics for changed grid geometry, radial layout, atom mapping and radial source.

This focused probe exercises isolated scalar grid setup plus periodic phase calculation, not the complete poisoned-global periodic FFT path. Full periodic solver preservation is assessed separately with native fast-HSE regressions. A test-fixture unallocated global kion was caught with bounds checking and fixed; the test runner retains `-fcheck=all` on the probe.

## Final receipts

| Check | Result |
| --- | --- |
| Focused native grid/phase/NLCC probe, OMP 2 | Passed |
| Four incompatible-update diagnostic cases | Passed |
| Full native context interleaving, OMP 2 | 8 cases, bit-identical |
| Original fast-HSE periodic driver, OMP 1 | 32 cases, current/energy output difference 0 |
| Python source tests | 305 passed |
| Source distribution → wheel | Passed |
| Fresh installed-wheel tests | 305 run, 1 optional plotting skip, no failures |
| Installed prepared source | All 24 replacement files byte-identical |
| Installed Japanese/English/AI guides | Byte-identical to checkout |
| Source and object forbidden-eigensolver audit | Passed before link |

Executable SHA256: `ddf428e129903c73fb06001de2b62a41953ce3f06feafc6e5e754c0b405d1451`.

Replacement manifest SHA256: `774c5a99c4c0c3f45b99336a22bd8ded9c1c6a0776724ec5350203306fd5106d`.

Local receipts: `/tmp/salmonpy-ppgrid-final-build-20261115.json`, `/tmp/salmonpy-ppgrid-final-{probe,contexts,default}-20261115/result.json`. Wheel: `/tmp/salmonpy-ppgrid-final-dist-20261115`; installed source: `/tmp/salmonpy-ppgrid-installed-source-20261115`. These paths record this machine’s validation, not installation requirements.
