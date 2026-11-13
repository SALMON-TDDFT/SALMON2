# Material-ID-independent FDTD validation, 2026-11-13

Implemented in the external-runtime compatibility patch. E edges use their
lower-coordinate material when both endpoints are nonvacuum. Vacuum boundaries
keep the legacy coefficients. H uses fixed spatial priority, with complete
material face/edge/corner exchange during initialization. LD impulse initialization
only contributes the owning material. Numerical time-update order and the shared
LD kernel are unchanged. This is staircase assignment, not dispersive averaging.
Runtime interface is now `salmonpy-runtime/2`; old multi-material binaries are
rejected by the Python adapter. The new sample and paired material guides explain
arbitrary per-material pole counts/coefficients and supported scope.

## Local evidence

Apple Silicon macOS, GNU Fortran 15, OpenBLAS; serial runtime includes Libxc 5.2.3.
All prepared-source/shared-kernel and forbidden-eigensolver object audits passed.

- Original finite two-layer fixture: 64x40x40, dt=0.003 a.u., 500 steps,
  one process/thread. Permuting material IDs previously changed rear Ez by
  relative L2 0.0275428472 (2.75%); after the fix all four E/H observation arrays
  are bitwise equal (maximum absolute difference 0).
- Original A-only and B-only fixtures: all eight observation files byte-identical
  to the old runtime. This is tested singleton compatibility, not a guarantee for
  all shapes, periodic boundaries or parallel decompositions.
- `test_fdtd_material_interfaces.py`: 19 comparisons pass. Three orientations of
  unequal-pole LD/LD, LD/dielectric and magnetic dielectric pairs; nonoverlapping
  shape order; duplicate identical material IDs; singleton old/new pulse and LD
  impulse; two-material impulse ID permutations. All comparisons use exact arrays.
- `test_fdtd_material_mpi.py`: magnetic interfaces placed at a rank corner;
  1x1x1 versus 1x2x2, 2x1x2 and 2x2x1 are bitwise equal. The pre-halo implementation
  fails this same regression (max field difference 8.32709e-6), demonstrating the
  test detects missing diagonal communication. No periodic-boundary claim.
- Python checkout: 91 tests pass. Fresh sdist -> wheel -> repaired wheel ->
  isolated installation: 91 tests pass. The installed two-layer sample executes.
- Relocatable v2 runtime packaging and loading pass. Native TDDFT relocation test
  reports ordinary output byte equality and A-to-J/interleaved-state/recovery pass.
- CI configuration now includes serial interface tests on all runtime-bundle
  platforms and Linux native MPI corner tests. These new CI changes have not yet
  been pushed/executed; local results must not be called cross-platform CI results.

The old runtime used for singleton comparison is the previously audited v1
relocated runtime. The new relocated binary SHA256 is
`63d6aff963e11e1000401ce208e4eaf3416b6fbe54a9a647c5d9ae3e6bfbb72a`.
Patch SHA256 is recorded in runtime/source.json, with the pinned upstream revision.

## Limits and next interface

The half-timestep check still changes rear Ez by about 1.24%; this short finite
fixture is not converged infinite-film R/T. ID independence does not establish
subpixel accuracy, rotational invariance or production convergence. Public
materials remain isotropic LD and constant dielectric models. E/H-TDDFT mixing
and the two-level Python/Fortran MPI interface are designed separately in
`plans/2026-11-13-two-level-mpi-interface.md`, not implemented by this patch.

Local artifacts (historical paths, not setup defaults):
`/Users/otobetoshihito/.codex/visualizations/2026/11/13/01a1d275-1ce2-7121-b8f4-cfd9dbc2d897/ld-two-layers-v2`.
Includes all original-fixture inputs/outputs, `verification.json`, native/Python
logs and `interface-invariance.png/.pdf`. Regression tools in this repository
allow rerunning with caller-supplied audited binaries and new output directories.
