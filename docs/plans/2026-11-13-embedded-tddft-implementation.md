# Embedded one-point TDDFT implementation plan

Goal: advance one A-driven electronic state on each supplied Python MPI solver
communicator, preserving the native initialization, odd/even orbital swap and
current evaluation. Then connect field phases using an explicitly documented
scheme; update the paired manuals only after numerical validation.

Architecture: extract initialization/advance/finalization from main_tddft into
an explicit native state shared by the CLI and embedded entry. Each process may
own one electronic state initially; distinct solver groups own distinct state.
Python owns MPI lifetime, native code duplicates only the supplied group.
Use an explicitly built host-MPI library with source/object provenance. No nested
launcher, import-time compilation, or claim of arbitrary multi-state reentrancy.

1. Add a native contract test: same prepared GS and endpoint A history on isolated
   groups of one and two ranks; compare every current to the serial worker. Test
   independent point histories and preservation of Python MPI after close.
2. Extract native RT state and step without changing operation order; retain the
   CLI path as a numerical reference. Add host library init/step/close ABI and
   explicit build provenance. Audit global communicator and ownership behavior.
3. Add Python collective provider/backend integration; restrict one point per
   group and one native state per process until lifetime/isolation tests pass.
4. Run source/object guards, CLI parity, internal-MPI tests and installed-wheel
   tests; document measured restrictions and sample commands.
5. Field coupling: inspect native multiscale temporal/spatial transfer before
   implementing E/H phases. The current pinned main_ms uses Weyl A propagation;
   its centered E/H observables are not the staggered standalone E/H integrator.
   Do not claim bitwise equivalence between these distinct discretizations.

## Completed checkpoint

Implemented the restricted native context, explicit build tool, collective Python
provider/backend and E/H field/response stepper. Internal orbital, real-space and
k-point decompositions are checked; the E/H map uses native collocation and its
transpose with an explicit half-step clock. The installed-wheel regression and
paired documentation are complete for this scope. See
`docs/native-eh-tddft-validation-2026-11-13.md` for evidence and limitations.
General multi-state reentrancy, restart and production convergence/scaling remain
outside the one-point first implementation.
