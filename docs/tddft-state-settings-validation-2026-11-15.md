# Fixed timestep per TDDFT context — 2026-11-15

This extends the full-context interleaving milestone with an owned, fixed timestep.
The internal Fortran rt_open accepts optional delta_t; omitted values capture the
input dt. Invalid nonpositive/nonfinite values are rejected before allocation or
output. Each s_rt owns time_step; initialization, Taylor coefficients, PT-CN,
field differences and output times use it. No numerical active-context/global-dt
switching is used. Fast HSE arithmetic and update order are retained.

Scope: local Apple Silicon, GNU Fortran, serial CPU/FFTE Poisson, fixed ions,
periodic scalar HSE, shared grid/k/XC and other configuration. HSE FFT resources
remain independently owned. No public Python API, concurrent calls, adaptive dt,
FDTD subcycling, heterogeneous grid/XC or general native failure atomicity is
claimed. Module XC scratch and code-optimization tables are the next extraction.

## Evidence

- Red: the new delta_t test failed compilation against the previous adapter,
  which had no delta_t argument (`/tmp/salmonpy-dt-red-20261115/build.log`).
- Eight cases each with OMP=1 and OMP=2: k meshes 2x1x1 and 3x1x1; Taylor/ACE,
  learned Taylor, full/direct Taylor, learned PT-CN; twelve steps each.
- State A dt=0.001, B dt=0.0015. Separate-process references use their own input dt.
  Mixed runs deliberately poison the old global dt. Current/energy/output times,
  orbitals, density, potentials and HSE histories are bit-identical to those refs.
- Reversed order and close/reopen match independent runs; invalid dt/A and guarded
  setting changes are rejected. The original native Acos2 full-driver comparison
  also passes at rtol=1e-12, atol=1e-14, including the returned charge-current sign.
- Mutation evidence: an external adapter copy ignoring selected_dt fails with
  `native output differs: b_rt.data` in `/tmp/salmonpy-context-dt-mutation.log`.
  Neither the validated engine nor the repository was modified for that mutation.
- Source and all application objects passed mandatory eigensolver audits before
  linking. The source inventory and completed-file SHA256 identities were updated.

- OMP1 receipt: `/tmp/salmonpy-context-dt-test1-20261115/result.json`; SHA256 `892e7c1782a99fc1e1e55fb00c6762f3e8db5484fe572a22132ae27bb2ff11f3`.
- OMP2 receipt: `/tmp/salmonpy-context-dt-test2-20261115/result.json`; SHA256 `ed50540430854722dc30295509add50269af4eb2cde44a5cbc3a634aedfd80fd`.
- Engine: `/tmp/salmonpy-context-dt-build-20261115/salmon`; SHA256 `d0580133a34f583ebeeb9072e9de977fc3702b56e8c8f28b9c3570e6b6674621`.

All fixtures are synthetic. This is implementation regression, not a converged
material or physical-accuracy validation. Output is outside the repository; no
real GS was collected. GitHub Actions was not invoked. Changes are uncommitted
at the time of this report.

## Legacy path and packaging

The unchanged default-construction path passed 32 full-driver cases (k=2/3/6/8,
dt=0.001/0.002, four propagation/factor modes) against the immutable fast-PR
reference: formatted current and energy differences are zero. Receipt:
`/tmp/salmonpy-context-dt-default-20261115/result.json`.

Python source suite: 305 passed. Sdist-to-wheel build passed, and a new isolated
wheel installation ran 305 tests with one optional plotting skip and no failure.
Installed runtime preparation also passed its source inventory/audit, including
the new timestep-owned completed files. Logs:
`/tmp/salmonpy-context-dt-python.log`,
`/tmp/salmonpy-context-dt-packaging.log`,
`/tmp/salmonpy-context-dt-installed-tests.log`,
`/tmp/salmonpy-context-dt-wheel-prepare.log`.

Read-only independent review found no remaining global-dt use in the supported
fixed-ion scalar-HSE step path and no blocker within this experimental scope.
