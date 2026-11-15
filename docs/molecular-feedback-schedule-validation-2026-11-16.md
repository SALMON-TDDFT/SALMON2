# Molecular feedback schedule regression — 2026-11-16

## Scope and skill use

Used the pinned materials-simulation-skills benchmark-and-mms-planner and
convergence-study dt_refinement helpers. Planner explicitly uses `general`,
not a claimed SALMON/Maxwell-specific family. The planner output is a plan,
not validation. No production spatial/material convergence campaign was run.

Evidence: `/tmp/salmonpy-feedback-vv-plan-20261116.json`,
`/tmp/salmonpy-feedback-oscillator-20261116.json` and
`/tmp/salmonpy-feedback-oscillator-order-20261116.json`.
Free-oscillator four-level refinement yields observed orders 1.99903 and
1.99976. This applies only to the fixture midpoint oscillator, not coupled TDDFT.
The standalone schedule test has two passing checks, including a reversed-work
sign negative control.

## Native results

Apple Silicon/macOS, Python 3.14. Two independent native vacuum field instances
in separate processes; four MPI ranks total, OMP/BLAS one thread each. Grid 40³,
spacing 2 bohr; unit-mass/unit-charge oscillator omega=2, deposition volume 8
bohr³. Environment has illumination, reference has none; one molecule only.
Explicit first-interval hold. Forcing uses E_n, deposit is j_(n+1/2), then field
advances to n+1. No future field is used before its required current exists.

| dt (au) | Steps | Selected E vs illumination-only max | Oscillator work residual max | Held-drive work | Temporal work discrepancy |
|---|---:|---:|---:|---:|---:|
| 0.003 | 64 | 2.304e-19 | 1.105e-24 | 1.869e-11 | -3.598e-11 |
| 0.0015 | 128 | 3.253e-19 | 2.262e-24 | 1.492e-11 | -1.748e-11 |
| 0.00075 | 256 | 5.421e-19 | 8.306e-24 | 1.301e-11 | -8.633e-12 |

Self fields are nonzero; cancellation is not a zero-source test. Selected-field
agreement passes rtol=1e-10/atol=1e-14, and oscillator work residual passes 1e-12.
A production coupled-energy balance is **not established**: the temporal work discrepancy is
large relative to absorbed work. It roughly halves when dt halves, identifying
first-order coupling error despite the second-order isolated oscillator.
Do not promote this held-field schedule as a production energy-conserving
molecular loop or attach the isolated TDDFT midpoint worker by relabelling time.
A centered/predictor-corrector schedule and its state restoration requirements
remain design work. Vacuum radiation reaction is removed by this policy.

Native electric work here is evaluated from sampled fields and deposited dipole
increments; no full-volume Maxwell energy or boundary flux was measured.
There is no material-environment scattering validation or multi-molecule loop.

Final runs: `/tmp/salmonpy-feedback-reviewed-{external,feedback}-level{0,1,2}-20261116`
with adjacent logs. Each report includes library/script/native-input hashes,
oscillator parameters, policy and baseline report/history hashes. NPZ stores
E/p/j clocks and separate total/reference/drive work arrays. Read-only review
confirmed causal indexing and work bookkeeping; its provenance request was fixed
and all six final native runs repeated. No GitHub Actions or push.

Final package checks: source suite 435 tests/one skip; sdist-to-wheel build,
installed suite 435 tests/two skips, no failures. The field-selection sample
also runs with installed Python `-I`. Dependencies were reused offline; this is
not a fresh dependency-resolution test. Logs:
`/tmp/salmonpy-feedback-source-tests-20261116.log`,
`/tmp/salmonpy-feedback-build-20261116.log`,
`/tmp/salmonpy-feedback-installed-tests-20261116.log`.
