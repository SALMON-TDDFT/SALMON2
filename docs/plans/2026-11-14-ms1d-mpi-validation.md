# Multiscale1D × two TDDFT groups: validation

Date: 2026-11-14. Local macOS 26.7.2 arm64, Python 3.14.5, host OpenMPI,
GNU Fortran 15, OpenBLAS; one OpenMP/BLAS thread per process. These are local
connection tests, not a converged Si optical calculation or a scaling study.

## Implementation and scope

`CoupledStepper.from_multiscale(context, bindings, field)` reuses the outer
coupling schedule and `ResponseSession`. One Python field rank owns the existing
Multiscale1D state. Two material groups each own an independent native TDDFT state
created from the same immutable GS input. Plan point order maps to cells 1..N.

Material callbacks use `step(A)`; legacy `advance(A_old, A_new)` callbacks are
wrapped with accepted A history. Physical charge current is not negated again.
Multiscale1D uses clock offset 0, native E/H keeps 1/2. Initial feedback stays
zero, separately from the observed initial current. Nonzero material A(0) is
rejected before material creation, with the point/value and vacuum/delay remedy.
No cross-group rollback or unrestricted native reentrancy is claimed.

## Provenance

The runtime was freshly prepared outside the repository from
`runtime/source.json`, not copied into this checkout. Shared LD/SBE kernels were
verified byte-identical and the forbidden-eigensolver checks passed.

- Upstream: `3322d0561c0cad5180679136b3eddbe8b909487c`
- Patch SHA256: `0f9faa9db8365cad1ecb457e5672d43c8618b440a6eaac63b3ee782229f4ae7f`
- SALMON executable SHA256: `12b91fa5d1e18035ebd01673824459d6aad257f8c00c3b9cdefc784579ba022e`
- Embedded library SHA256: `782fc8552062e900574d3c898957e1710a71cdd204270d0fc8e4eb01908d783d`

## Numerical results

Synthetic normalized Si orbitals (8³ grid, two bands, one k point, PZ), 400 steps,
dt=0.02 a.u.; two material cells, spacing 5 bohr, vacuum (300,50), incident pulse
specified in `tools/test_longitudinal_salmon.py`. Initial orbitals are not
self-consistently converged. Native SALMON runs with two macro ranks.

| Comparison | Result |
|---|---:|
| Whole-grid A vs native SALMON, max relative error | 8.48275e-16 |
| Point 1 j vs native SALMON, max relative error | 2.64602e-15 |
| Point 2 j vs native SALMON, max relative error | 2.35997e-15 |
| 3 vs 5 total ranks, A/j/grid A max absolute difference | 0 (exact) |
| Serial independent-worker comparison | passed |
| Same-GS points have distinct currents | max difference 5.3703e-6 |

SALMON field file n records A[n] and feedback j[n-1]. Comparisons preserve this
placement; the first two feedback records are zero. Histories preserve physical
initial current independently. All-grid comparison includes vacuum and boundary
propagation, not just material sites. A sign flip, one-step shift and vacuum-only
perturbation were checked to fail the reference comparison.

## Communication measurement

Per-rank accumulated wall times divided by 400; the table takes the maximum
rank average for each phase (not the average of per-step maxima).

| Total ranks | A distribution | j collection | Material phase including waits |
|---|---:|---:|---:|
| 3 (1 field + 1 + 1) | 0.798 µs/step | 8.006 µs/step | 152.531 µs/step |
| 5 (1 field + 2 + 2) | 0.887 µs/step | 8.650 µs/step | 143.106 µs/step |

Each A/j logical payload is 48 bytes/step. This is not network traffic volume.
A timing covers leader and subgroup broadcasts. j timing covers leader gather,
point ordering and field-group distribution. Collective wait is included; native
work, subgroup reductions and error agreement are in the material phase. Other
validation and Python field work are outside these three measurements. No extra
barrier is inserted. One Mac and two tiny material points cannot establish
multi-node performance. Clock failure records null times and measurement errors;
a forced clock failure produced bit-identical physical histories.

## Reproduction and regression

Use fresh output directories outside the checkout:

```sh
PYTHONPATH=. OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  python tools/run_ms1d_mpi_tests.py BUILD_JSON PREPARED_SOURCE NEW_OUTPUT
```

This bounded runner creates the two-point native/serial reference, executes
three- and five-rank coupling, checks rank parity, and executes MPI contract
regressions. The latter cover mixed step/legacy callbacks, nonconsecutive IDs,
initial-A rejection before callbacks, timing mismatch, capacity and backend
failures, cleanup, and measurement failure isolation. JSON and NPZ evidence is
saved beneath the output directory. Linux runtime CI now invokes this runner;
remote CI has not been executed for these unpushed changes.

Additional checks passed:

- Source unit suite: 140 tests, six unavailable-native tests skipped.
- Source distribution → wheel build and installed native hash verification.
- Installed wheel unit suite: all 140 tests, no skips.
- Existing E/H/LD/TDDFT MPI regression: all 13 cases, exact field-rank parity.
- New user sample: 400 steps with five ranks, from checkout and installed wheel.
- English/Japanese manuals and AI JSON snapshot synchronization checks.

Local full evidence for this run is in `/tmp/salmonpy-ms1d-acceptance`,
`/tmp/salmonpy-ms1d-eh-final`, and `/tmp/salmonpy-ms1d-wheel-example` (temporary,
not distributable fixtures or runtime locations to assume on another machine).
