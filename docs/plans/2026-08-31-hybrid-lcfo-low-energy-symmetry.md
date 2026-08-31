# Hybrid LCFO Low-Energy Symmetry Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the fatal all-WF+PW closure gate with post-LCFO symmetry acceptance of the occupied 128 states and the degeneracy-complete low-energy target beginning at 384 states.

**Architecture:** Preserve the existing full retained-basis representation as a finite diagnostic projected action.  Add a focused low-energy symmetry module that selects a degeneracy-complete eigenvalue window and measures the induced action on occupied and target eigenspaces.  The continuation solves enough empty states to certify the 384-state boundary, accepts only the final occupied/target projectors and density, and records the full-basis defect without using it as a stop condition.

**Tech Stack:** Fortran 2008, MPI, distributed row-owned generalized eigensystems, LAPACK/ScaLAPACK, Python test runners.

---

### Task 1: Specify the degeneracy-complete low-energy symmetry contract

**Files:**
- Create: `src/gs/dc/dg_hybrid_low_energy_symmetry.f90`
- Create: `tests/dg/test_dg_hybrid_low_energy_symmetry_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_low_energy_symmetry_mpi.py`

**Step 1: Write the failing MPI test**

Define the public routines:

```fortran
select_dg_hybrid_symmetry_target(eigenvalues,requested_rank,tolerance,target_rank,ok,message)
evaluate_dg_hybrid_low_energy_symmetry(comm,metric,representation,coefficients,eigenvalues,&
  occupied_rank,target_rank,tolerance,occupied_defect,target_defect,energy_defect,ok,message)
```

Use a three-dimensional fixture with `S=I`, orthonormal eigenvectors, and a
projected symmetry action `diag(1,1,0.5)`.  Require the first two states to pass
even though the complete three-dimensional action is not unitary.  Change the
action so state two leaks into state three and require target rejection.

Use eigenvalues `[-1,0,0,1]` with requested rank two and require selection of
rank three.  Require rejection if the supplied spectrum ends inside the
degenerate cluster, so production cannot silently split an unresolved
multiplet.

**Step 2: Run RED**

Run:

```text
python3 tests/dg/run_dg_hybrid_low_energy_symmetry_mpi.py
```

Expected: compilation fails because `dg_hybrid_low_energy_symmetry` does not
exist.

**Step 3: Implement the minimal module**

For target selection, compare adjacent eigenvalues using
`tolerance*max(1,abs(e_i),abs(e_j))`, extend through the complete boundary
cluster, and require one resolved state above the selected cluster unless the
entire basis was solved.

For each operation form the induced target action

```text
D = C_target^H S R C_target
```

and measure `D^H D-I`.  Measure the occupied block independently and measure
`D^H diag(e) D-diag(e)` for energy covariance.  Validate dimensions,
orthonormality, finite values, MPI agreement, and report defects separately;
do not inspect or require unitarity of unused high-energy directions.

**Step 4: Run GREEN**

Run the new runner on 1, 2, 4, and 8 ranks.  Expected: PASS.

**Step 5: Commit**

Commit only the new module and focused test:

```text
test(dg): define low-energy LCFO symmetry acceptance
```

### Task 2: Downgrade complete retained-basis closure to a diagnostic

**Files:**
- Modify: `src/gs/dc/dg_hybrid_retained_basis_symmetry.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_hybrid_retained_basis_symmetry.py`
- Modify: `tests/dg/check_dg_hybrid_continuation_route.py`

**Step 1: Write RED route assertions**

Require the representation builder to return both a finite
`full_basis_closure_defect` and structural success.  Require `main_dft` to emit

```text
[HYBRID-RETAINED-BASIS-SYMMETRY] closed=<0|1> defect=<value>
```

and forbid `error stop` based only on `closed=0`.  Continue to require fatal
handling for invalid dimensions, incomplete metric rows, failed metric solve,
collective failure, or non-finite defects.

**Step 2: Run RED**

Run both static contract checks.  Expected: FAIL because closure and structural
validity still share one callback flag.

**Step 3: Separate structural validity from closure**

Add `closure_defect_arg` and `closure_ok_arg` outputs.  Set the existing
callback success true after a finite representation and defect are produced,
even when `closure_ok_arg` is false.  In `main_dft`, print the receipt once on
rank zero and pass the finite projected representation into continuation.

**Step 4: Run GREEN and focused protection**

Run:

```text
python3 tests/dg/check_dg_hybrid_retained_basis_symmetry.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
python3 tests/dg/run_dg_hybrid_low_energy_symmetry_mpi.py
```

Expected: PASS.

**Step 5: Commit**

```text
fix(dg): defer retained-basis symmetry acceptance to LCFO
```

### Task 3: Solve and certify the 384-state low-energy window

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_hybrid_continuation_controller.f90`
- Modify: `tests/dg/test_dg_hybrid_continuation_controller_mpi.f90`
- Modify: `tests/dg/check_dg_hybrid_continuation_route.py`

**Step 1: Write failing solve-window tests**

Extend `dg_hybrid_continuation_state_count` with a requested symmetry target.
Require an occupied count of 128, target 384, and basis count 1560 to request
at least 385 eigenpairs, retaining the extra eigenvalue as proof that the
boundary cluster is complete.  Cover clipping at the full basis and invalid
targets smaller than the occupied window.

Require the production route to pass the existing `ntarget` into the concrete
continuation and to keep all solver coefficients/eigenvalues needed by the
degeneracy-complete target analysis.  Density reconstruction must continue to
use only the occupied columns and physical occupations.

**Step 2: Run RED**

Run the controller MPI runner and continuation route check.  Expected: FAIL on
the old occupied-only solve count and missing target analysis.

**Step 3: Implement target solve and post-solve checks**

Request `max(occupied-gap count, ntarget+1)` eigenpairs, bounded by the retained
basis size.  After every expensive candidate check, select the
degeneracy-complete target and call
`evaluate_dg_hybrid_low_energy_symmetry` for occupied and target spaces.
Replace the full-Hamiltonian covariance acceptance gate with the target energy
covariance defect.  Keep the existing occupied-projector covariance and add a
direct mapped-core density defect using the authoritative pencil maps.

Fail a stage only on occupied, target, target-energy, or density defects above
`dg_ow_symmetry_tolerance`.  Record requested/extended target ranks and all
four defects in the final acceptance receipt.  Do not symmetrize coefficients,
Hamiltonian, occupations, or density.

**Step 4: Run GREEN**

Run:

```text
python3 tests/dg/run_dg_hybrid_continuation_controller_mpi.py
python3 tests/dg/run_dg_hybrid_low_energy_symmetry_mpi.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
python3 tests/dg/run_dg_hybrid_continuation_acceptance_mpi.py
```

Expected: PASS on every configured rank count.

**Step 5: Commit**

```text
feat(dg): certify low-energy LCFO symmetry
```

### Task 4: Protect the physical post-LCFO acceptance route

**Files:**
- Modify: `tests/dg/check_dg_hybrid_retained_basis_symmetry.py`
- Modify: `tests/dg/check_dg_hybrid_continuation_route.py`
- Modify: `tests/dg/run_dg_hybrid_si64_continuation_rt.py`

**Step 1: Add RED acceptance assertions**

Require the Si64 runner to reject a log that lacks the retained-basis
diagnostic, requested/extended target ranks, occupied defect, target defect,
target-energy defect, or density defect.  Require it to reject a run where a
post-LCFO physical defect exceeds tolerance, while permitting
`closed=0` for the full retained basis.

**Step 2: Run RED**

Run the static route checks.  Expected: FAIL because the receipts are not yet
protected.

**Step 3: Update receipt validation**

Add the exact parser and finite/tolerance checks without weakening existing
GS-to-RT provenance, energy, checkpoint, or zero-field RT assertions.

**Step 4: Run GREEN and protected regression**

Run every focused runner referenced by the parent handoff plan, the new
low-energy runner, all continuation runners, and:

```text
cmake --build build-hybrid-commit -j2
```

Expected: all PASS and successful build.

**Step 5: Commit**

```text
test(dg): protect low-energy symmetry handoff
```

### Task 5: Perform one fresh Si64 GS-to-RT acceptance calculation

**Files:**
- Preserve evidence under `verification-si64-dg-continuation-20260830/` without staging it.

**Step 1: Run one fresh eight-rank calculation**

Use a new output directory, `OMP_NUM_THREADS=1`, no timeout, and the existing
capture hook so any later failure retains the fixed variational matrices.  Do
not start another conventional DC calculation while this run is active.

**Step 2: Verify the GS receipts**

Require full-basis diagnostic publication, successful target extension,
post-LCFO occupied/target/energy/density defects within tolerance, final
Hybrid GS acceptance, and checkpoint publication.

**Step 3: Verify RT**

Require the runner's existing exact-checkpoint provenance, zero-field RT
stability, norm, charge, energy, and symmetry assertions.

**Step 4: Preserve evidence and report**

Do not stage the large output directory.  Report whether the 384-state window
was extended and identify the worst symmetry operation and defect.

**Step 5: Resume the parent handoff plan**

Complete the remaining protected verification and checkpoint requirements in
`docs/plans/2026-08-30-hybrid-dc-symmetry-handoff.md`.
