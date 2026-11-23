# Fast-PR HSE initialization child communicator ownership

## Implemented scope

A current explicit fast-PR runtime exposes context_mpi_protocol1. Every rt_open
owns communicators created through the communication wrapper during its synchronous
initialization scope. Its close releases persistent requests, XC/exchange/FFT
owners and then those children. Shared-model siblings, attachment-owned owner0
children and the host communicator remain live. Monotonic attachment-local owner
IDs distinguish reused registry slots; explicit comm_free_group clears handle and
owner tag to avoid later double release. This is not arbitrary heterogeneous-input
reentrancy or ownership of every MPI resource created anywhere in the solver.

The common density symmetry map no longer borrows the most recently initialized
context's retireable r communicator for k-only groups: a size-one r group uses
stable MPI_COMM_SELF. Its singleton sum arithmetic is unchanged. Other distributed
r groups retain the supplied communicator and are not newly supported for shared
HSE. Complete-file HSE overrides and their inventories are updated; original
source/provenance and fixed binaries are retained.

## Verification

Local macOS26 Apple Silicon, Python3.14, Open MPI; OMP/BLAS1; Mac64GiB/18 cores.
Fresh bundled-source preparation and new external build; source/object forbidden-
eigensolver audits pass. No GitHub Actions dispatched.

- New MPI resource probe first rejects the previous binary's missing capability.
- Taylor learned/PT-CN learned × full/reduced-k × internal MPI1/2:8 cases,32
  sibling create/step/close cycles each, against an independent solver group.
  One state tracks10 children; two track20; closing the sibling returns to10;
  closing all returns to0. Slots peak at20 and do not grow with cycles. Idle model
  recreation starts from original GS; host communicator remains usable. Survivor
  current/energy matches independent control exactly in every case.
- Existing shared/isolated/full/reduced/direct-resume matrix: Taylor strict and
  learned, PT-CN off and learned × MPI1/2 × full/reduced-k. Shared/independent and
  resumed traces remain exact; full/reduced stays within unchanged tolerances.
  Symmetry-breaking fields still reject without changing the accepted state.
- Normal PBE/LD/SBE/shared full-k HSE/FDTD mixture,6 ranks,128 steps: continuous
  fields/currents match the preceding runtime baseline exactly. Save at63, fresh
  direct resume with a handled late-candidate failure preserves all originals and
  output streams; retry continues exactly with65 new steps and no replay.
- Source and isolated installed wheel:603 Python tests each, one skip each.
  sdist→wheel/native hashes pass. Installed-source preparation verifies the new
  complete overrides. Installed reduced-k learned PT-CN MPI2 resource probe runs
  64 cycles with the same bounded counts and exact independent control.
- Paired guides, AI JSON and latest checklist synchronized. Read-only ownership
  review found no blocker. Initial declaration-order compiler error and temporary
  unsynchronized-guide test failure are retained with successful corrected logs.

Evidence retained outside Git at developer outputs salmonpy-context-mpi-20261118.
No calibrated process-memory claim, physical-convergence certificate or other-
platform guarantee follows from tracked communicator counts. Standard TDDFT,
unscoped/common groups, broader GS/different-input ownership, I/O/timers and
initialization/allocation/MPI failure recovery remain separate backlog items.
