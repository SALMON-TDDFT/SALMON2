# Archive and custom-backend SCF initial guesses — design draft

Status: approved by project lead: direct common-array handoff. Existing retained-GS implementation
was committed/pushed as8f3a7663. This document does not claim array import exists.

## Goal and choices

Reuse the established wavefunction/occupation-only initial-guess semantics for
GS/2 archives and user-written electronic backends, without a retained SALMON
source owner or intermediate GS record files.

Recommended: one explicit immutable array-seed contract with archive and custom
provider adapters; Fortran copies into a fresh SCF owner. Alternative: reuse the
existing archive-to-native converter then load a retained owner. That avoids a
new array ABI but requires intermediate files, explicit converter tools and an
extra owner. A file-only custom-provider protocol would also unnecessarily force
user code to serialize data. Do not silently use either as fallback.

## Proposed ownership and interface

- Keep `session.run(...,initial_guess=live_state)` unchanged.
- Add an initial-guess data object containing declared source model/sampling,
  wavefunction/occupation arrays, axis order, local/global bounds and provenance.
  Public constructor and archive factory produce the same contract. Custom
  backends may provide that object without inheritance or SALMON internal types.
- Source arrays are copied or read through owned immutable storage. The adapter
  must not retain mutable user buffers or a live custom-backend lifetime.
- Archive adapter verifies existing GS/2 metadata/chunk hashes and reads only
  required k/orbital chunks where feasible. No hidden network/build/converter or
  job launch. Explicit pseudopotential bytes/references remain required.
- Source metadata is checked against target model; permit atomic displacement,
  preserve current cell/grid/k/spin/electron/band/XC/potential compatibility.
  Compare explicit k vectors and weights, never infer a matching mesh from count.
  No interpolation or implicit real/complex truncation. No convergence claim is
  fabricated from an imported seed.
- Define one typed C/Fortran import entry; document array order, extents, ownership,
  communicator scope and rejection semantics. Python schedules and packs data;
  numerical checks/copy/density/potential calculations belong to Fortran/C.
- Fortran imports values into a fresh candidate, recalculates density/potentials,
  and starts fresh mixing/HSE resources. Temporary seed arrays are released after
  initialization. Source/siblings remain unchanged on success or supported input
  rejection; do not promise arbitrary OOM/process/MPI failure recovery.
- Use only owning-group setup agreements, not per-step application-wide checks.
  Keep large MPI seeds rank-local where possible; no mandatory full-state gather
  or a full-wavefunction broadcast to every rank.

## Verification and documentation

1. Synthetic archive/custom-provider seeds exercise ownership, layouts, complex
   values, k weights, metadata mismatch and NaN/Inf rejection before mutation.
2. Compare the array route against the established retained-owner seed route on
   PZ and fast-PR HSE, including MPI k/spatial layouts. Interleave receivers,
   mutate/release provider buffers and continue; results must remain identical.
3. Keep existing cold-start/retained-seed regression results unchanged. Test
   source and installed wheel; update Japanese/English guide, AI task JSON and
   checklist. Add a short custom-provider example.
4. Never collect real wavefunction/density/restart data without applicable user
   authorization. Synthetic archives may test disk I/O; actual numerical values
   can be compared in RAM. Scalar verification output follows existing grants.
5. Source/object eigensolver audits, scalar IEEE inquiries and Git-only inventory
   build are required for any changed compiled source. A changed source tree also
   needs refreshed binary-bound NLCC evidence before advancing its capability
   identity; do not transfer the8f3a7663 runtime receipt to a new binary.
6. Required CI remains integrity/packaging/synthetic/API only. Do not dispatch or
   add material-observable reproduction CI.

## Boundaries

No new restriction on unrelated custom backends, no automatic fallback to old
formats, no automatic acceptance of an unconverged GS. Archive import is an
initial estimate, not complete mixing-history SCF restoration. Broader source
models can be represented in the data contract, while a particular solver reports
only the combinations it can initialize correctly.
