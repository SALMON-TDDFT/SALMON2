# Native protected mixed steps and append-output rollback

Continue the approved common protected-step interface. Implement explicit native
same-lifetime transaction protocol1 for standard TDDFT, native E/H and fast-PR
HSE. Modes0 capture,1 rollback+close owned RT streams,2 retire/flush,3 reopen.
Backups retain complete existing numerical payloads/XC/HSE histories in memory;
never clone FFT/MPI owners or replay accepted propagation. Preserve original
wavefunction/overlap flags rather than disk-LOAD normalization. Snapshot slots
are separate from staged checkpoint candidates; reject concurrent ownership.

RT output units stay owned by Fortran. A shared full-source helper flushes,
closes and reopens only the two audited persistent RT streams on their root.
Python root journals regular *.data file lengths/inodes in the owned output
folder. Rollback closes writers, truncates appended attempt records and reopens
at the accepted end before retry. Simple FDTD observation writers are already
closed per call. Reject auxiliary spatial/projection/TTM/frequency/restart output
profiles before capture. No concurrent readers/writers, durable crash-atomic
publication or output-device failure recovery is promised; tentative appends may
be externally visible until the collective vote. Unknown native failure poisons
owners; never infer rollback from old protocols.

Expose native_step_targets(coupled, custom_targets=()) with exact supported owner
types and explicit overrides for user backends. Built-in SBE endpoint wrappers
use their existing handle snapshots and restore wrapper latches/clocks. Native
Python A/E/current/clock refs are retained. New opt-in --protect-steps in ordinary
mixed example uses common CoupledStepper.step_protected; default path unchanged.

Acceptance: MPI field2 + SBE1 + TDDFT2 and field2 + SBE1 + HSE2 Taylor/PTCN,
plus actual mixed TDDFT+HSE if fixture compatibility allows. First/warm/final
one-rank deep/late failures preserve every native checkpoint byte, wrapper aliases
and ordinary files; retry trajectories/output exact to uninterrupted and old
normal builds. No replay. Explicit new external builds, manifests/overlay bases,
source/wheel tests and JP/EN/AI sync; persistent local evidence, no Actions,
OMP/BLAS1 and warning-only48GiB on user64GB/18-core host.
