# Group5: retain MPI child ownership during RT propagation

Initialization already assigns child communicators to a state token. Add an
internal resume operation and scope standard embedded step / fast-PR rt_step
with that existing token. Children created lazily during propagation then retire
with the state, rather than waiting for attachment teardown. Restore the inactive
selector on every returned path, collectively rejecting conflicting activation.
Standalone calls remain unchanged. No numerical payload, codec, Python API or
MPI split change; main-thread sequential and one attachment/process remain.

Extend the existing native directory/I/O fixture: late children for two owners,
nested-scope rejection, independent release, attachment-owned child retained and
host MPI survives detach. Compare ordinary and full/reduced learned HSE native
currents/energies/history with old controls and memory rollback. Existing GS
references only, no GS copies/SCF/CI. Source/object audit before all builds,
source/wheel and JP/EN/AI/list synchronization.

## Acceptance

Native late-child ownership/reused slots/nested rejection/host MPI test passes
in source and installed tooling. Old communication module rejects compilation
of the new fixture. PBE/rgrid2 and full/reduced learned HSE MPI2/group currents
(and HSE energy/history) exactly match prior independent controls; rollback and
sibling timing remain exact. Final old HSE controls use hse-timers-build, not
an intermediate MPI-owner build. Source/installed638:637pass1skip. Source/object
audits passed, JP/EN/AI/list updated. No GS copies/SCF/CI.
Evidence: calc2026/salmonpy-rt-mpi-owner-record-20261120; uncommitted.
