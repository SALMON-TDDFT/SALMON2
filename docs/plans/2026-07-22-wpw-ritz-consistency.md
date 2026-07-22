# WPW Ritz Consistency Diagnostic Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Locate any inconsistency between the reduced Ritz update and the next direct fixed-H operator application.

**Architecture:** Arm diagnostics after updates 31/95/159, retain the complete rank-local post-update residual vectors, and compare them state-by-state with direct H/S residuals at iterations 32/96/160. Also compute matched per-state reduced-coordinate residuals from the solved reduced matrices. Keep all solver decisions unchanged.

**Tech Stack:** Fortran, MPI, Python source contracts, CMake.

---

### Task 1: Add deterministic Ritz-boundary diagnostics

**Files:**
- Modify: `tests/dg/test_dg_wpw_matrix_free_scf.py`
- Modify: `src/gs/dc/dg_wpw_matrix_free_scf.f90`

1. Add RED source-contract assertions requiring post-update `HZ*C`/`SZ*C`
   residual construction, `R_red=H_red*C-S_red*C*lambda`, next-iteration
   same-state residual-vector comparison, collective finite-value validation,
   `Q_new^H S Q_new-I`, and `[DG-WPW-RITZ-CONSISTENCY]` output.
2. Assert the existing `if(residual<tol.and.orth<tol)return` predicate precedes
   diagnostics and remains unchanged.
3. Run `python3 tests/dg/test_dg_wpw_matrix_free_scf.py`; expect failure because
   the diagnostic is absent.
4. Add local post-update H/S images, complete pending residual vectors, per-state
   reduced residual norms, and a rank-uniform pending flag.
5. Define comparison iterations as 32/96/160; arm pending post-update state only
   at 31/95/159, so each comparison refers to the same Ritz vectors and carried
   eigenvalues.
6. Before entering or skipping any diagnostic `gram`, collectively reduce local
   validity.  Store and clear the pending transaction identically on all ranks.
7. Use production `gram` to summarize occupied/extra post-update physical
   residuals, reduced-coordinate residuals, and the same-state vector defect
   `R_direct-R_post`, including worst state IDs. Define each relative defect as
   `||delta R||/max(||R_direct||,||R_post||)`, map `(0,0)` to zero, and use the
   first state on ties. Also compute and log `max|Q_new^H S Q_new-I|` with the
   production global Gram inside the same collective transaction.
8. Return a distinct nonzero `info` only for non-finite or structurally invalid
   diagnostic data; do not gate on defect magnitude.
9. State in the source contract that diagnostics run only after the unchanged
   convergence test fails and add no H/S/preconditioner application.
10. Run the matrix-free test; expect PASS.
11. Commit the focused implementation and test.

### Task 2: Verify, review, and run the physical comparison

**Files:**
- Modify: `docs/plans/2026-07-22-wpw-ritz-consistency.md`
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

1. Run:
   - `python3 tests/dg/test_dg_wpw_matrix_free_scf.py`
   - `python3 tests/dg/check_dg_wpw_fixed_h_relaxation.py`
   - `python3 tests/dg/test_wpw_generalized_algebra.py`
   - `python3 tests/dg/run_dg_wpw_core_wannier_seed_mpi.py`
   - `cmake --build build-mpi-eigenexa -j 8`
   - `git diff --check`
2. Request code review for numerical meaning, MPI collectives, unchanged solver
   state, and selected-iteration pairing. Resolve all Critical/Important findings.
3. Commit only intended source/test changes.
4. Create a fresh B=6 run from Task 16 with
   `yn_dg_wpw_preconditioner='n'` and default retained search history.
5. Run 8 MPI ranks through the unchanged 160-inner fixed-H limit.
6. Compare the matched per-state reduced-coordinate residual, post-update
   physical occupied/extra residual, next-direct residual, same-state vector
   boundary defect, post-update metric orthogonality, and subsequent
   re-diagonalized residual at direct comparison iterations 32/96/160 (paired
   with post updates 31/95/159).
7. Record the evidence and next hypothesis in this plan and the session handoff.
8. Run `git diff --check` and commit only the two result documents.

## Result (2026-07-23)

Implemented in `f9fe8b4`.  The two-rank long-spectrum fixture executes all three
transactions and verifies near-zero boundary defects.  Focused tests, the full
MPI/EigenExa build, and two review rounds reported no blocking issue.

The fresh physical run is
`stage2d_wpw_runs/20260723_task18_ritz_consistency_b6/run.log`.  It uses the
Task 16 B=6 input with no preconditioner and retained search history, reproduces
the residual history, and reaches the unchanged `info=40` limit.  Matched values
are:

| post→direct | reduced occupied / extra | physical occupied / extra | relative boundary defect occupied / extra | metric orthogonality |
|---:|---:|---:|---:|---:|
| 31→32 | `3.21E-08 / 7.37E-08` | `1.97E-03 / 4.04E-03` | `3.29E-12 / 9.71E-13` | `1.62E-11` |
| 95→96 | `3.72E-09 / 9.71E-09` | `3.57E-04 / 1.03E-03` | `1.29E-11 / 2.21E-11` | `2.34E-12` |
| 159→160 | `1.99E-09 / 1.48E-09` | `2.56E-04 / 2.09E-04` | `3.90E-11 / 1.02E-10` | `1.98E-13` |

The post-update `HZ*C`/`SZ*C` residual and the following direct H/S residual
agree to at least about `1E-10` relatively, and the updated Ritz vectors remain
metric orthonormal.  The reduced eigenproblem is also solved much more accurately
than the physical Ritz residual.  Therefore stale operator images, an inconsistent
Ritz update, and the update-to-next-application recurrence are not the plateau
source.  The large residual lies outside the current Rayleigh--Ritz search span;
the remaining diagnosis should address why retained residual correction
directions reduce it so slowly (especially metric-aware preconditioning), not
change the occupied-W representation, cutoff, tolerance, or publication gates.
