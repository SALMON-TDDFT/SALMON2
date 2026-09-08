# Hybrid Fixed-Frame Preconditioner Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Supply the approved gauge-covariant, energy-shifted preconditioner needed before Task 8 production wiring.

**Architecture:** Store distributed rows of the unitary map Q from an immutable physical reference into the current uncompressed fragment coordinates. Build reference H/S diagonals once per operator epoch without a Hamiltonian eigensolve, then apply Q D(epsilon)^-1 Q^dagger to each residual using its current Rayleigh value. Keep the existing residual-only fixture callback compatible and require exactly one explicit callback.

**Tech Stack:** Fortran 2008, MPI, existing bounded LOBPCG fixture, CMake, Python test runners.

---

Use the existing `wpw-s-orthogonal-complement` worktree only. Preserve all dirty
changes and logs. No Si64 calculation, new worktree or production activation.
The design is the 2026-09-03 amendment in
`2026-09-02-hybrid-divided-one-shot-lcfo-design.md`.

### Task 8a: Fixed-frame preconditioner and shifted bounded-update callback

**Files:**

- Create: `src/gs/dc/dg_hybrid_fragment_preconditioner.f90`
- Modify: `src/gs/dc/dg_hybrid_fragment_subspace.f90`
- Modify: `src/gs/dc/CMakeLists.txt` (only the new source line)
- Create: `tests/dg/test_dg_hybrid_fragment_preconditioner_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_fragment_preconditioner_mpi.py`
- Modify: `tests/dg/test_dg_hybrid_fragment_subspace_mpi.f90`

**Step 1: RED fixed-frame action and provenance tests**

Provide `s_dg_hybrid_preconditioner_key` containing fragment ID, generation,
operator epoch, basis/metric/operator/reference fingerprints. Provide a cache
with private components. Expose:

```fortran
call prepare_dg_hybrid_fragment_preconditioner(comm,n,row_ids,q_rows,h_rows,s_rows,&
  key,tolerance,cache,fingerprint,ok,message)
call apply_dg_hybrid_fragment_preconditioner(comm,row_ids,key,cache,shifts,residual,&
  output,fingerprint,ok,message)
```

The matrix rows have extent `(size(row_ids),n)`; shifts have one entry per
residual column. Q spans the complete local fragment basis, never a selected
occupied subset or the terminal complete-system compressed catalog. Output is
allocatable and remains unallocated on failure.

Test diagonal physical H/S against the analytic regularized answer, then
rotate H/S, Q and r with a complex unitary mixing WF columns while leaving PW
columns fixed. Assert identical physical action for unequal per-state shifts.
Reject nonunitary Q, duplicate/missing rows, nonfinite matrices/residuals/shifts,
nonpositive reference metric norms, rank-disagreeing controls/keys, stale epoch,
operator/basis/metric/reference identity and changed row layout. Cover zero-row
ranks and state-count growth. Reject a failed rebuild transactionally without
destroying an earlier valid cache. Check denominator signs and exact-zero floor.

Run `python3 tests/dg/run_dg_hybrid_fragment_preconditioner_mpi.py`.
Expected: RED because the module/API is absent.

**Step 2: Implement the distributed frame action**

Validate all shapes/ownership/key fields collectively, then stream reference
columns to compute `h_a=q_a^dagger H q_a`, `s_a=q_a^dagger S q_a` and certify
Q unitarity. Validate Hermitian input rows without gathering an n-by-n matrix
on one rank. Keep distributed Q rows and replicated O(n) reference diagonals.
Bind their integrity, layout and complete key into an exact receipt.

At application reduce `Q^dagger r`, divide by regularized denominators and
multiply local Q rows. Define `scale_j=max(1,maxabs(h)+abs(epsilon_j)*maxabs(s))`,
`roundoff_j=64*epsilon_machine*n*scale_j`,
`floor_j=max(tolerance*scale_j,roundoff_j)`. Resolvable negative denominators
retain a negative floor; unresolved zeros use positive floor. Detect overflow
and nonfinite results collectively, restoring floating-point trap settings.
No reference-frame eigensolve or default identity recovery is permitted.

**Step 3: RED current-shift delivery**

Append optional `apply_shifted_preconditioner(input,rayleigh_values,output,ok)`
to `advance_dg_hybrid_fragment_subspace`. Make the original residual-only
callback optional without changing old positional calls. Require exactly one
of the two callbacks; a failed shifted callback must never fall back to the
old callback. The test records independently computed current Rayleigh values
at each step and verifies that the supplied shifts change when X changes.
An invariant in-span problem must still perform its small Ritz rotation.

Run `python3 tests/dg/run_dg_hybrid_fragment_subspace_mpi.py` and observe RED
for the missing shifted callback API before modifying production code.

**Step 4: Implement and link shifted callbacks**

At the existing preconditioner invocation pass the `values` calculated from
the current accepted X. Collectively agree callback selection before any
callback execution. Retain successful cap/rollback semantics, no greater trial
dimension, and the strict solver's unchanged contract. Link the new module.
Exercise the real frame preconditioner through a three-step bounded update
in both original and rotated coordinates and compare physical projectors,
residuals and common-mu density. Include an occupied-plus-guard subset smaller
than the fragment basis, with a nontrivial metric and nonidentity action.

**Step 5: GREEN, review and checkpoint**

Run fresh:

```text
python3 tests/dg/run_dg_hybrid_fragment_preconditioner_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_subspace_mpi.py
python3 tests/dg/run_dg_hybrid_block_cg_mpi.py
python3 tests/dg/run_dg_hybrid_divided_operator_mpi.py
python3 tests/dg/run_dc_fragment_occupation_mpi.py
cmake --build build-hybrid-commit -j2
git diff --check
```

Require PASS at configured 1/2/4/8 ranks and independent review without
Critical/Important issues. Save new verification logs without overwriting old
ones. Commit only these Task 8a files/hunks and the approved plan amendments;
exclude the unrelated canonical-PP CMake line. Report this prerequisite as a
checkpoint, not completion of Task 8. Resume production wiring in the parent
plan after review.
