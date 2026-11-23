# Isolated TDDFT: measured electron number

Prior molecular feedback committed/pushed as b7657d4a to the independent SALMONpy
branch. This follow-up is uncommitted. No manual heavy CI dispatch.

## Implementation

E2 returns electronic dipole change3, total charge dipole3 and occupied electron
number1. Fortran reuses rt%rIe(itt), already integrated by subdip, and worker_ne
at initialization. No new density integration, reduction or physical update.
Scalar IEEE finite check precedes output. Python exposes frame.electron_number
and state.electron_number; Fortran convergence kernel computes loss fractions.
Target1e-5 is diagnostic and does not stop a finite calculation.

Evidence revision2 and isolated replay checkpoint revision2 require seven
observables. Arrays are E0/observables0/observables/E_history/observable_history.
No old-protocol/checkpoint fallback. Snapshot and replay comparisons include
measured electron number. NaN/Inf replies are rejected as numerical failures.

## Actual C2H2 / LD / FDTD

Reuse the GS from the explicitly authorized C2H2 campaign; no new SCF.
New runtime built from a Git-tracked-files-only export with source/configure/
incremental/object prohibition guards. PZ, fixed ions,24×24×32 at0.35Å.
Outer MPI4, one field rank plus bookkeeping rank per field,64 steps,dt=.003au.
E2/batch reference dipole difference remains3.85523175341039e-12e·bohr.
Independent E-driven states, memory/disk replay and worker recovery pass.

Maximum electron loss9.947598300641403e-14, relative loss9.947598300641637e-15.
Maximum relative gain4.227729277772696e-14. Loss target1e-5 passes.
Inject rejection at interval30, restore and retry: every retained array including
measured electrons exactly equals continuous execution. Original E/p/j/drive/time
arrays exactly equal the earlier six-dipole-component version. Installed wheel
LD retry also equals source continuous execution, including electron number.

## Source identity and limits

Source inventory updated for main_tddft.f90/tddft_worker.f90. GS-free prescribed
XC and core-force tests rerun using a freshly audited Libxc-enabled build:
MPI1/2 × OMP1/2, PZ/Libxc-PZ/PBE, all pass. XC discrepancy up to2.854e-10,
core-force discrepancy up to5.951e-8. Only then advance source capability identity.
This is not a new whole-application NLCC force receipt for either new executable.

Source/wheel suites684 tests each, one optional skip each. Paired guides/AI catalog
synchronized. No quantum-work/total-Maxwell-energy certificate, optical convergence,
multinode or other-platform claim. Electronic restoration remains history replay.
All campaign evidence stays under canonical calc2026/salmonpy-c2h2-feedback-20261122.
