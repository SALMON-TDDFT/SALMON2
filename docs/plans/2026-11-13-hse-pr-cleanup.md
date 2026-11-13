# HSE PR cleanup implementation plan

**Goal:** Prepare the 2026-11-06 HSE/SOC proposal for upstream review without the developmental giant-Gamma routes.

**Approved architecture:** Retain ordinary HSE06/PBE0/PBEh, dense-k MPI+OpenMP, SOC/NLCC, occupied-support ACE, ordinary Gamma and restart. Remove giant-Gamma SR group/wave FFT and source-support ACE; preserve shared full-support spatial/DC kernels. Reject removed modes explicitly; never silently change an approximation. Keep learned dense-k history distinct from the excluded Gamma source adapter. No remote jobs or fixed executables are modified.

**Base:** codex/hse-soc-pr-20261106 at 5926312f66a19ee6afdd0a6cbd68ac99039933de. Old hse-merge-prep at 483058cf restored separately at /private/tmp/salmon-hse-merge-prep; do not build it (prohibited reference in lcfo_diag_chefsi_complex.f90).

**Tech stack:** Fortran, CMake, Python prohibition checker, GNU Fortran, MPI, FFTW/Libxc/BLAS.

## Tasks

1. Restore source/configure/incremental/application-object prohibition gates from the pre-cleanup implementation. Verify configure and incremental rejection using a synthetic fixture; never execute PZHEEV.
2. Remove SR group/wave implementations and source-support ACE construction. Preserve the canonical Gamma Hermitian pair backend and ordinary full-support spatial exchange. Remove unsupported input switches from the namelist and source support from the accepted choices. Verify removed modes reject, defaults/ordinary Gamma remain valid, and no removed module dependencies remain.
3. Update input/review documentation, list missing features versus old hse-merge-prep, and remove excluded-route examples. Keep historical source and evidence in their original branches.
4. Run the workspace source checker before any build. Build serial and MPI CPU/HSE configurations, audit application objects, and run bounded synthetic numerical checks for the retained exchange paths. Do not collect GS data. Report unavailable platform checks. Run diff whitespace checks and inspect the final diff.

## Inventory before changes

Old PR source versus selected base: 75 added, 7 removed and 58 modified files under src (path-based comparison; renames count as add/remove). The selected base includes EXX functional unification, global hybrids/rVV10, distributed ACE/MLWF, conventional/DC-to-real-space RT, fractional occupations, dense-k pair/orbit FFT reuse, worker services/layout, SOC/noncollinear/NLCC, factor-history restart and source-level compiler fixes. These are absent or substantially older in the September PR. The selected base also includes upstream updates through 376635eb and later scalar/SOC worker optimizations. The November 12 tip only preserves standalone Si512 candidates; those are not linked runtime modules.

## Known preparation defect

Commit 05359b43 removed cmakefiles/pzheev_guard.cmake, its src/CMakeLists.txt hooks and the bundled checker. The selected base passes the workspace source audit but lacks enforced configure/incremental gates. Restore them before building. Historical validation prose claiming active guards must be corrected against actual source.

## Completion

All four preparation tasks are complete. The ordinary source/ACE path and
full-support Gamma/DC exchange remain; dedicated SR group/wave/source-ACE
implementations and the unlinked standalone Si512 snapshots are omitted.
Three CPU application builds, nine input fixtures, eleven spatial
oracle/rejection cases and four build-guard fixture outcomes passed.
Scope, reproduction, independent-review findings and validation limits are
recorded in `docs/hse-pr-scope-20261113.md`. The separate SALMON-DOCS working
copy received a narrowly scoped proposal note; its existing edits were preserved.
The preparation stays on its own local branch; no remote job, fixed binary,
upstream branch, PR publication or remote push is changed.
