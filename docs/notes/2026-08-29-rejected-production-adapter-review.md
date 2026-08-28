# Rejected production-adapter experiment

The uncommitted Task 7a.6 adapter experiment was preserved long enough to run
the following focused checks before the design was reconsidered:

- production adapter MPI fixture: PASS at 1, 2, 4, and 8 ranks;
- continuation SCF MPI fixture: PASS at 1, 2, 4, and 8 ranks;
- continuation acceptance, SIPG operator, production face traces, and divided
  controls: PASS;
- `git diff --check`: PASS.

The experiment was not committed because review found:

- a rank-local volume-kernel failure could deadlock later collectives;
- Hermiticity and symmetry were accepted as cached booleans and numerical
  residual values did not independently gate publication;
- final refresh restored older density and trace fields before reacceptance;
- fingerprints depended on rank decomposition and did not prove unique global
  row ownership;
- `state%projector` held `C^dagger S C`, not the basis-space occupied map
  `C_occ C_occ^dagger S`;
- the adapter payload did not prove complete face coverage, effective
  WF/PW-selection provenance, or one uniform lambda;
- the adapter fixture manually called phases and therefore did not exercise
  the complete solver, rollback, nontrivial metric, or rank-local failure.

The replacement design is
`docs/plans/2026-08-29-dg-concrete-continuation-solver-design.md`.  The
experimental source and fixture were intentionally removed before restarting
Task 7a.6 under TDD; this note preserves their test evidence and rejection
rationale.
