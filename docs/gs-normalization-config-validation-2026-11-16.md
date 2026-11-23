# GS electron-count normalization isolation — 2026-11-16

## Implementation

`s_gs_run_config%normalization` captures `nelec` and `nelec_spin(2)` as values.
`gs_run` copies them before initialization. Density convergence uses these owned
values, or legacy input capture if no GS config is owned. The selection remains
nelec for nspin=1, positive spin-count sum for nspin=2 otherwise nelec; the DC
state's elec_num_tot still overrides the result. Spin summation uses an int64
intermediate to avoid integer overflow. All other operations/reductions/order
remain unchanged. Only non-DC rho_dne checks a positive divisor at the boundary;
invalid values return status1 before activation.

This is normalization of convergence, not an electron-count/charge/occupation API.
Material initialization and occupation-setting inputs remain shared. No new public
Python interface or concurrent GS claim is made.

## Evidence

Local Apple Silicon/macOS, GNU Fortran/OpenBLAS/static Libxc, serial CPU fast HSE,
-O2 -fallow-argument-mismatch (no unlimited-line-length option), build jobs2,
GS OMP/BLAS1. No GitHub Actions or network download.

- Red probe failed against previous mixing build: missing normalization config
  and resolver, `/tmp/salmonpy-gs-normalization-red-20261116/build.log`.
- Build/source/object audits passed:
  `/tmp/salmonpy-gs-normalization-build-20261116.json`.
  Original fast-HSE reference: `/tmp/salmonpy-gs-metric-reference-20261116.json`.
- Five convergence metrics, 14 cases (three density metrics with pre-SCF0/1/2/3;
  two potential metrics with pre-SCF disabled): two sequential states' printed
  residuals, gs_info/eigen/k outputs and every restart payload match the original
  executable exactly. Results: `/tmp/salmonpy-gs-normalization-{rho_dne,norm_rho,norm_rho_dng,norm_pot,norm_pot_dng}-20261116/result.json`.
- Four additional rho_dne cases double only B's normalization count. Electronic
  outputs and every restart byte still match the unchanged original input;
  B's printed residual is half the original, A is unchanged. Compare the exact
  decimal values with a bound derived solely from independent E-format rounding:
  half one printed quantum for B plus one quarter quantum for the original.
  Initial fixed-relative log tolerance was too strict for the printed digits;
  no solver modification or payload tolerance relaxation was made.
  Final records: `/tmp/salmonpy-gs-normalization-double-final-20261116/result.json`.
- Retained A includes its normalization config in byte traces, unchanged across
  caller-value mutation and B's run/close. Pure resolver checks cover unpolarized,
  positive spin sum and zero-spin fallback; nonpositive divisor rejected before
  activation. Spin selection checks do not constitute a polarized SCF run.
- Source Python suite: 435 tests, 1 skip (`/tmp/salmonpy-gs-normalization-python-20261116.log`).

- Explicit sdist-to-wheel build passed: `/tmp/salmonpy-gs-normalization-dist-20261116`.
  Installed wheel suite: 435 tests, 2 skips (`/tmp/salmonpy-gs-normalization-installed-tests-20261116.log`).
  Local dependencies were reused; no online dependency resolution test. Packaged
  replacement sources, inventories, native probe, regression tool and AI manual
  match the checkout bytes.
- Fast-HSE TDDFT interleaving: 8 cases matched the original, 2/3 k points,
  Taylor4/PTCN/factor modes, different dt/OpenMP/capacity, 12 steps:
  `/tmp/salmonpy-gs-normalization-tddft-20261116`.
- Japanese/English runtime and embedded AI guides synchronized; diff checks clean.

Both references' hash sets and residual traces remain in result records; temporary
wavefunction/density payloads are deleted after comparison. Same tiny Si diagnostic
fixture (8³ grid, 2 k points, loose stopping threshold), not physical convergence.

## Remaining limits

MPI/DC/polarized/spinor/DFT-MD/geometry-optimization/GPU/other-platform runs were
not newly performed. DC priority and legacy fallback received code review only.
Electronic occupations, gap reporting and other material initialization still
read legacy electron-number inputs; these have not been silently rerouted.
Molecular work remains paused. Native errors outside the explicit validation may
still terminate the host; this is not a fully independent/concurrent GS API.
