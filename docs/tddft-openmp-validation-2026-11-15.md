# Owned OpenMP selection — 2026-11-15

The explicit fast-HSE profile stores stencil/current/force OpenMP selection
flags in each s_stencil. Numerical consumers (Hamiltonian, normal and microAc
stencils, current and force) receive the owning stencil; no bare shared flag
consumers remain. Legacy input selection and the 32^3 auto threshold are kept.
The internal rt_open optional openmp_mode accepts auto/orbital/rgrid and rejects
other values before allocation/output. Omission preserves individual inputs.
This does not set process thread counts or support concurrent host calls.

## Local verification

- Red native probe failed against the preceding immutable build because the
  owning flag fields/configuration routine did not exist:
  `/tmp/salmonpy-openmp-red-20261115/build.log`.
- Owned workspace probe with OMP=2 passes explicit orbital/rgrid selection,
  automatic threshold, independent owners, and the previous XC/table checks.
- Eight full-step cases with OMP=2: k=2/3 and four propagation/factor modes;
  A uses orbital splitting/dt=0.001, B real-space splitting/dt=0.0015. Separate
  processes, interleaving, reverse order, close/reopen and old-input poisoning
  pass. Current/energy/full diagnostic state traces are bit-identical to each
  state's isolated reference. Invalid constructor selectors are rejected.
- The unchanged original fast-HSE Acos2 driver comparison passes the existing
  rtol=1e-12/atol=1e-14 contract and charge-current sign test.
- Legacy default path: 32 full-driver cases against the original fast PR;
  formatted current and energy maximum differences are zero.
- Mandatory source/object eigensolver audits pass before linking. Completed-file
  inventory hashes match freshly prepared sources. New numerical lines do not
  exceed 132 columns.
- Source Python suite: 305 passed. Sdist-to-wheel build passed. Fresh installed
  wheel: 305 tests, one optional plotting skip, no failures. Installed preparation
  and source byte equality pass. Japanese/English guides and AI JSON synchronized.
- Read-only independent review found no scoped ownership/OMP/default-path blocker.

## Limits

Synthetic fixed-ion scalar-HSE regression on local Apple Silicon only. Force,
microAc and AVX/SVE optimized branches are wired but not numerically validated by
these current tests. MPI, GPU, host concurrency and physical convergence are not
claimed. Common grid/k/XC and other global input restrictions remain. No Actions
were launched; changes are uncommitted at report creation.

## Receipts

- `/tmp/salmonpy-owned-omp-context-20261115/result.json`; SHA256 `2b968f6c90dafce6ebd1a6246f0344ec69bbe180e4e953b6c71e6f8f3f0893bc`.
- `/tmp/salmonpy-owned-omp-default-20261115/result.json`; SHA256 `8d11737cc2f8f026cbc0a8525d63232a66bc5b1b8044ecc1457b4240b74befc9`.
- `/tmp/salmonpy-owned-omp-probe-20261115/result.json`; SHA256 `a362c4f9cdd3d0b6ed2f320d19acce1265a11e4ed7238f1f364626c1bf79a0a6`.
- `/tmp/salmonpy-owned-omp-build-20261115.json`; SHA256 `526437ea12ac03e28e7c93b43902b10e1f74f2a807d9c769574421114122580d`.

Logs: `/tmp/salmonpy-owned-omp-python.log`,
`/tmp/salmonpy-owned-omp-packaging.log`,
`/tmp/salmonpy-owned-omp-installed-tests.log`,
`/tmp/salmonpy-owned-omp-wheel-prepare.log`.
