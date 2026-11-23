# Fast-PR HSE SCF k-only repartition

Extend optional native k-repartition protocol to version2 in HSE-enabled builds;
ordinary-only builds retain version1. Preserve default exact-layout behavior.

Map dense ACE/scf-damped/initial/midpoint factors and local source/action/work
arrays by global representative k index in Fortran, keeping replicated occupation,
phase/mixing/history data. Slice factor/PT-CN per-k histories when present.
Rebind metric communicator IDs, never copy FFT/backend handles. Clear only the
rank-specific geometric phase cache; new-layout kernel/map setup occurs at the
next natural exchange refresh, not by replaying SCF or rebuilding accepted ACE.
Packed/distributed ACE is outside this k-only dense profile. Preserve diagnostics
conservatively when scalar ACE condition estimates cover several old k points.

## Validation completed

The user explicitly authorized source/wheel HSE validation GS, density and SCF
restart datasets on2026-11-20. Artifacts were saved under canonical calc2026;
cloud synchronization was not independently confirmed.

- Source/object prohibition audits and ordinary/HSE builds passed.
- Compiled synthetic Fortran probe covered dense ACE/cache/history merge/slice,
  packed-factor rejection and malformed local-cache shape rejection.
- Full-k and reduced-k HSE1↔2 plus PBE pre-stage1↔2 fresh-process restart passed.
  Restored cut arrays exact, executed bootstrap iterations0, then5 new iterations
  match source-initialized continuous8 at rtol1e-9/atol1e-10.
- Installed wheel HSE1→2 reduced and2→1 full passed; ordinary PBE2→1 passed.
- Ordinary/HSE fixed-previous-runtime controls: SCF arrays and subsequent16-step
  RT current exact; sibling state and memory-snapshot continuation preserved.
- Source and installed suites637 tests each:636pass,1skip. No tests added in count.
- First parity invocation lacked required MPI launch, then an unconstrained
  OMP18 invocation failed an exact info comparison. These failed logs are retained;
  the established OMP1/BLAS1 controls passed without changing tolerances.
- Generated serializers, guide JSON synchronization and diff whitespace checked.

Evidence: calc2026/salmonpy-scf-fast-hse-k-record-20261120.
No CI dispatch. Orbital/spatial/packed/distributed ACE, changed-runtime migration
and physical convergence certification remain outside this validated scope.
