# SCF exchange-control owner routing validation

Local Mac verification on 2026-11-15, following commit `7b5fcb95`. No CI triggered.

## Change

GS SCF now accesses adaptive readiness, exchange freezing and pending support changes through the `s_dft_system`-bound HSE resource owner. Optional-field setter and individual queries preserve the existing flag storage in `s_exx_rt_state`. Unbound systems resolve to the legacy owner. Localization validation accepts the system and checks the same owner. Collective reduction, support-event reset and temporary freeze/restore positions remain unchanged.

This is internal Fortran routing, not a public Python API. It does not create independent full GS contexts, remove remaining adaptive configuration globals/DC material swapping, or permit concurrent calls.

## Evidence

- The new native test initially failed to compile against the preceding audited build because the owner-aware contract was absent (`/tmp/salmonpy-scf-owner-red-20261115/build.log`).
- The native resource test now passes: flags on A do not affect B; changing the unbound legacy owner does not affect B; changing one optional flag preserves others; localization validation selects the bound owner while the legacy owner would reject. Existing native refresh/action comparisons remain included.
- Serial and MPI HSE builds completed with source/pre-link/all-object forbidden-eigensolver audits. Build records: `/tmp/salmonpy-scf-owner-{serial,mpi}-build-20261115.json`.
- Four tiny serial GS cases match the immutable original fast-HSE executable: `/tmp/salmonpy-scf-owner-gs-20261115/result.json`.
- Two DC fragments with two and four total MPI ranks, four pre-SCF settings per layout: all eight cases match original phase/energy traces and fragment scalar files. Records: `/tmp/salmonpy-scf-owner-dc2-20261115/result.json`, `/tmp/salmonpy-scf-owner-dc4-20261115/result.json`.
- Eight independent TDDFT context cases, differing timestep/OpenMP/capacity, remain identical: `/tmp/salmonpy-scf-owner-contexts-20261115/result.json`.
- Source Python suite: 305 passed. sdist → wheel → fresh environment: 305 tests, one optional plotting skip, no failures. Logs: `/tmp/salmonpy-scf-owner-python-tests.log`, `/tmp/salmonpy-scf-owner-installed-tests.log`.
- Installed runtime preparation replacement bytes and all four runtime/AI guides match the tested checkout. Distribution: `/tmp/salmonpy-scf-owner-dist-20261115`; installed environment: `/tmp/salmonpy-scf-owner-installed-20261115`.

These are same-platform implementation regressions with small synthetic/coarse inputs, not physical convergence, complete adaptive-localization accuracy or cross-platform verification. The DC fixtures use full-support exchange; independent flag/localization routing is exercised by the native probe. Build jobs and execution threads were limited to four in aggregate. Historical audited sources/binaries were not modified.
