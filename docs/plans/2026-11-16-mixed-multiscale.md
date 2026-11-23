# Mixed LD / SBE / TDDFT multiscale target

User priority: one multiscale calculation may mix LD, SBE, TDDFT and user material backends. Continue input isolation where it removes a measured integration obstacle; do not make complete Fortran reentrancy a prerequisite for the existing process/MPI routes. Molecular work remains paused.

- Material ID identifies an immutable material definition; point ID identifies its independent mutable state. Different points may share initial data but never propagation history. Report both IDs with spatial point ordering.
- Crystalline SBE/TDDFT use A/c endpoints and physical charge current j. No E-driven crystal API. LD uses E; endpoint demonstration derives E=-(A_new-A_old)/dt. This explicit schedule is distinct from native E/H's semi-implicit coefficients/Q update.
- Start with 1D and one material/point per electronic MPI group. Internal TDDFT decomposition stays in that group. Large LD regions remain many local states on field ranks, not one MPI group per LD cell in production.
- One point initially has one response model. Do not add LD and electronic j at the same point without a separately specified physical decomposition (double-counting risk).
- Keep field solver, material preparation/response and execution topology separate. User backends retain the same contract.

## Acceptance sequence

1. LD+SBE on one Multiscale1D field: nonzero response, exact independent serial/MPI schedule comparison, timing and point order. Implemented sample: mpi_mixed_multiscale.py (3 ranks).
2. Add an actual TDDFT point (4 ranks for managed serial worker), compare serial/MPI currents and field. Validated for the Si starter PZ GS over100 steps; all A/j values match an independent serial schedule exactly. Real GS/runtime must be supplied explicitly; no downloads or builds on import.
3. Repeat with internally parallel native TDDFT groups; verify communicator ownership, run IDs and communication measurements.
4. Native MPI E/H FDTD with local multi-material LD plus remote SBE/TDDFT: preserve LD native update order and interface deposition, compare pure-LD and pure-electronic limits before mixed layers. One-field-rank demonstration does not establish this topology.
5. Persist material/point/state identity, snapshots and whole-session restart, with explicit capability refusal for unsupported backends.

Numerical regression is not convergence of a real composite. Keep 4GiB/4-core local budget, explicit resource reporting and no routine remote CI.

## Native acceptance update

Steps3/4 now have a bounded native fixture: one E/H rank with two field-local LD materials, one SBE rank and one/two TDDFT orbital ranks. Native/serial material replay,1/2-rank full-history parity,LD vacuum contrast and pure-LD batch parity pass. See docs/native-mixed-material-validation-2026-11-16.md. Multiple field ranks were subsequently verified in the update below. Large numbers of LD cells and whole-session restart remain open.

## Distributed-field acceptance update

Two field ranks with SBE1+TDDFT1 pass all x/y/z split comparisons against field1+SBE1+TDDFT1. Both LD materials span both y/z field blocks. Four-core campaign does not test simultaneous field2+SBE1+TDDFT2 (five ranks). Scaling,aggregate memory measurement and whole-session restart remain open. See docs/distributed-field-mixed-validation-2026-11-16.md.

## Usability follow-up after topology checks

Topology demonstrations do not yet construct an entire native mixed run from user material IDs,layer geometry and rank budgets. Reuse the existing configuration/identity helpers for that next integration; avoid a duplicate input schema. Preserve explicit LD E versus crystalline A/c timing. Track the known prefilled-incident longitudinal validation gap before broadening admission. Whole-session restart and recoverable collective errors remain separate milestones. The active checklist contains the current prioritized audit.

## Layer entry and repeated-material update

The input entry and prefilled-incident validation gap above are addressed in01d43cc1; see2026-11-16-native-layer-assembly.md. Next acceptance step5 remains open: complete coupled-field/LD/electronic checkpoint and restart is not implemented.

Repeated Si and repeated SBE layouts now each run within four ranks, sharing material definitions/initial data while assigning distinct point IDs, MPI groups and evolving states. Every point matches an independent serial replay exactly over128 steps. This verifies process/group isolation, not multiple native TDDFT owners in one process. Evidence: docs/repeated-mixed-material-validation-2026-11-16.md.
