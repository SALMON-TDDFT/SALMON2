# Mixed LD / SBE / TDDFT multiscale target

User priority: one multiscale calculation may mix LD, SBE, TDDFT and user material backends. Continue input isolation where it removes a measured integration obstacle; do not make complete Fortran reentrancy a prerequisite for the existing process/MPI routes. Molecular work remains paused.

- Material ID identifies an immutable material definition; point ID identifies its independent mutable state. Different points may share initial data but never propagation history. Report both IDs with spatial point ordering.
- Crystalline SBE/TDDFT use A/c endpoints and physical charge current j. No E-driven crystal API. LD uses E; endpoint demonstration derives E=-(A_new-A_old)/dt. This explicit schedule is distinct from native E/H's semi-implicit coefficients/Q update.
- Start with 1D and one material/point per electronic MPI group. Internal TDDFT decomposition stays in that group. Large LD regions remain many local states on field ranks, not one MPI group per LD cell in production.
- One point initially has one response model. Do not add LD and electronic j at the same point without a separately specified physical decomposition (double-counting risk).
- Keep field solver, material preparation/response and execution topology separate. User backends retain the same contract.

## Acceptance sequence

1. LD+SBE on one Multiscale1D field: nonzero response, exact independent serial/MPI schedule comparison, timing and point order. Implemented sample: mpi_mixed_multiscale.py (3 ranks).
2. Add an actual TDDFT point (4 ranks for managed serial worker), compare serial/MPI currents and field. Sample options are prepared, but three-way execution is not yet validated. Real GS/runtime must be supplied explicitly; no downloads or builds on import.
3. Repeat with internally parallel native TDDFT groups; verify communicator ownership, run IDs and communication measurements.
4. Native MPI E/H FDTD with local multi-material LD plus remote SBE/TDDFT: preserve LD native update order and interface deposition, compare pure-LD and pure-electronic limits before mixed layers. One-field-rank demonstration does not establish this topology.
5. Persist material/point/state identity, snapshots and whole-session restart, with explicit capability refusal for unsupported backends.

Numerical regression is not convergence of a real composite. Keep 4GiB/4-core local budget, explicit resource reporting and no routine remote CI.
