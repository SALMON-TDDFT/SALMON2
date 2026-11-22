# Frozen response implementation plan

> **For Claude:** Execute task-by-task with failing tests before implementation.

**Goal:** Implement and validate a native Fortran frozen-H finite-q response without instantaneous eigenstate projection.

**Architecture:** Represent the commutator source as low-rank factors. Build independent block Krylov spaces with left/right real-space H actions, solve the small projected Sylvester resolvent for all shifts, and measure the full equation residual from low-rank factors. No real-grid density matrix or real-grid Hamiltonian is constructed by this solver. Small projected eigensystems are internal numerical reductions, not a conduction-state expansion for main RT.

**Tech Stack:** Fortran, BLAS/LAPACK, existing finite-value helpers, SALMON k-point MPI adapter. Existing jobs and fixed binaries remain unchanged.

1. Add `testsuites/hsex/fortran_frozen_response.f90` and an audited standalone builder. Verify missing-module failure. Test stationary/coherent sources against an explicit small Liouvillian, shifts, uniform-potential response, occupancy rotation, capacity rejection and scalar NaN/Infinity rejection.
2. Add `src/xc/hsex_frozen_response.f90`: commutator factors, two-sided orthonormal spaces, H-action caching, projected multi-shift solve, true residual and diagonal-vertex contraction. Verify red-to-green. Distinguish success, rank cap and invalid/callback failure. Use no mutable module state.
3. Add independent impulse wavefunction propagation and small Hartree-RPA references. Check chi0 agreement and measure the difference between snapshot-RPA W and a self-consistent evolving background for coherent states. Add frozen-background commutator diagnostic.
4. Register focused native tests and document their scope. Run production-source/object prohibition audits before new executables; include all linked application objects. Check finite inputs and invalid-value rejection.
5. Add an opt-in Si finite-q diagnostic, reusing real-space H actions and occupied k-point exchange. First compare equilibrium chi0*x at matched eta against a reference, then a coherent snapshot. Retain MPI distribution and report basis sizes, true residuals, time and memory. Do not publish a new production W or change optical head/wings as part of this diagnostic.
6. Record results under canonical calc2026, no GS payload collection. Report numerical accuracy and cost, including cases where the response space grows too large; do not infer speedup from initial source rank.

Acceptance is staged: successful small-model verification does not establish a production screened interaction for arbitrary coherent pump states. Full-q W updates follow only after finite-q and static-W assumptions are validated.
