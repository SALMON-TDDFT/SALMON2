# XC and stencil workspace ownership — 2026-11-15

The explicit fast-HSE completed-source profile now owns XC 3D scratch and PZ
flat work arrays in each s_dft_system. A standalone calc_xc call can pass its own
PZ workspace; omission uses call-local allocatable storage. The functional and
Libxc handles remain read-only, and distinct output/workspace components avoid
parent/component argument aliasing. Shape/spin resize frees only that owner.

Modulo tables and optimized-callable flag belong to each s_stencil. Both RT and
DFT initialization and all located Hamiltonian/stencil consumers route that
object explicitly. No active-owner switching is introduced. Owning contexts
release their allocatable workspace/table components on destruction.

## Local verification

- The new native probe first failed against the previous build because owned
  workspace types and explicit stencil arguments did not exist. Evidence:
  `/tmp/salmonpy-workspace-red-20261115/build.log`.
- Final workspace probe passes with OMP=1 and OMP=2: PZ unpolarized and spin
  potential equality against original built-in routines, distinct shapes,
  resizing by shape/spin, modulo-table contents/bounds, freeing B and continuing A.
- Full native TDDFT: eight cases with OMP=2, k=2/3, four Taylor/ACE/direct/PT-CN
  modes, dt A=0.001/B=0.0015. Interleaving, reverse order and close/reopen agree
  bit-for-bit with separate-process references for current, energy, orbitals,
  density, potentials and HSE histories. The original-driver Acos2 comparisons
  pass at rtol=1e-12, atol=1e-14. Global-dt poisoning remains enabled in mixed runs.
- Default full-driver path: 32 cases, k=2/3/6/8 and dt=0.001/0.002, against the
  immutable fast-PR original. Formatted current and energy differences are zero.
- Sources and all application objects passed mandatory eigensolver audits before
  linking. Source inventories and replacement SHA256 identities were updated.
  New executable source lines do not introduce lines exceeding 132 columns.
- Python source suite: 305 passed. Sdist-to-wheel build passed. Fresh wheel
  installation: 305 tests, one optional plotting skip, no failures. Installed
  preparation passed inventory/audits; all prepared replacement bytes and the
  installed probe/tool match the checkout.
- Independent read-only review found no ownership, aliasing or OpenMP blocker
  within the stated scope. The requested spin-PZ numerical test was then added
  and passes. No GitHub Actions were run.

## Limits

Synthetic fixtures only, no collected real GS or physical convergence claim.
This is **not** heterogeneous full TDDFT: common grid/k/XC guards remain, as input
and OpenMP scheduling globals still exist. Native concurrent calls, MPI, GPU,
AVX/SVE tuned-kernel dispatch, general failure atomicity and RSS limits remain
outside this validation. Changes are uncommitted at report creation.

## Receipts

- workspace_omp1: `/tmp/salmonpy-workspace-final-probe1-20261115/result.json`, SHA256 `404373e91608bd07a0995db2b4733005c96e75ba7e8539eaf2a5285974e9893a`.
- workspace_omp2: `/tmp/salmonpy-workspace-final-probe2-20261115/result.json`, SHA256 `d4990386134fb36c81edbd56fa3d899693ec0f7bf0c3a577385935b4555d0632`.
- full_context_omp2: `/tmp/salmonpy-workspace-final-context2-20261115/result.json`, SHA256 `fcf22cb7d3703032916056a2ab20ebca8e0106d598b003680791f5e617fff397`.
- legacy_full_driver: `/tmp/salmonpy-workspace-final-default-20261115/result.json`, SHA256 `8978f0bdef9a439f5cc04e28098c49db5ef4d82b61c40042254d1f9de64d77ef`.
- Engine: `/tmp/salmonpy-workspace-final-build-20261115/salmon`, SHA256 `2d1811cbfeb344486a6b10cb6b445126e41ac901450b75708a7c99d2356f5932`.

Supporting logs: `/tmp/salmonpy-workspace-python.log`,
`/tmp/salmonpy-workspace-final-packaging.log`,
`/tmp/salmonpy-workspace-installed-tests.log`,
`/tmp/salmonpy-workspace-wheel-prepare.log`.
