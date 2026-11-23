# Explicit cell resolution — 2026-11-15

s_cell_config and resolve_cell extract cell-vector/length/spacing precedence into
a value-configured routine. Explicit lengths precede diagonal vectors, then
points times spacing; tilted vectors retain original norm arithmetic. Explicit
spacing precedes length/point spacing. The compatibility caller preserves legacy
global al normalization for downstream dependencies. No whole-context common
cell/grid/k/XC guard was relaxed; k/pseudo/XC/HSE input extraction remains.

The resolver validates finite positive resolved lengths/spacings and rejects
zero scaled determinant before committing its value-only result. This does not
certify near-singular conditioning. Other downstream lattice checks remain.

## Verification

- Red probe failed against previous build on missing configuration/resolver:
  `/tmp/salmonpy-cell-red-20261115/build.log`.
- Final native probe: explicit-length precedence, diagonal-vector fallback,
  point-times-spacing fallback, tilted geometry and independent results pass.
  NaN, invalid points and dependent-vector cells reject without changing the
  prior result. The singular-cell case was added after independent review and
  its scaled determinant fix was re-reviewed without remaining blocker.
- Full-context regression:8 synthetic cases, OMP2, differing dt/capacities/OpenMP,
  separate processes/reversed order/close-reopen; current/energy/state traces
  bit-identical to isolated references. The cell/grid remains common.
- Legacy full driver:32 cases vs immutable fast-PR original; formatted current
  and energy differences zero. Original native Acos2 comparison also passes the
  existing rtol1e-12/atol1e-14 contract.
- Source/all-object eigensolver audits pass. The unvalidated development build
  was incrementally rebuilt for the singularity fix before issuing its first
  build receipt; all replacement bytes match fresh preparation. Earlier validated
  reference builds were untouched. New executable native lines stay within132.
- Python source suite305 passed. Sdist-to-wheel and fresh installed tests pass
  (305 tests, one optional plotting skip). Installed preparation passes audits;
  all replacement/probe/tool bytes match the checkout. Japanese/English/AI guides
  synchronized. No Actions launched; uncommitted at report creation.

## Limits

Local Apple Silicon synthetic fixed-ion scalar-HSE regression. Lower-level tilted
cell resolution is not heterogeneous or general tilted-cell TDDFT validation.
No MPI/GPU/concurrent host/AVX-SVE/force-microAc/physical convergence claim.

## Receipts

- `/tmp/salmonpy-cell-probe-20261115/result.json`; SHA256 `9a164b278fcc3d21b8afb560834eb47fef879fa03611507c602b7b9abf80b507`.
- `/tmp/salmonpy-cell-context-20261115/result.json`; SHA256 `e8a5e6a98d81e19eab763572eb9c9349d8de5740a34133a322e9d87fb0d5f60d`.
- `/tmp/salmonpy-cell-default-20261115/result.json`; SHA256 `614185ef60c5179efa41ce4790dfc2cca634bb609198fba7878342ace86bcf80`.
- `/tmp/salmonpy-cell-build-20261115.json`; SHA256 `c0367ce1c0986d2abbcc810833fcf8f66fc3616561043ed4190f73d8fcdf4f87`.

Logs: `/tmp/salmonpy-cell-python.log`, `/tmp/salmonpy-cell-packaging.log`,
`/tmp/salmonpy-cell-installed-tests.log`, `/tmp/salmonpy-cell-wheel-prepare.log`.
