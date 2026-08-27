# WF+PW Fragment-Compatible Symmetry Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Guarantee a valid divided WF+PW production catalog when none of the supplied nonidentity affine generators preserves the DC fragment partition.

**Architecture:** Prepend a distributed identity spatial map and identity reciprocal rotation before testing generator compatibility. Retain the identity and every supplied operation that passes the existing whole-fragment collective contract, then build the reciprocal and windowed-PW catalogs from those aligned selected columns.

**Tech Stack:** Fortran 2008, MPI, standalone Python MPI runners, SALMON source-contract tests.

---

### Task 1: Specify identity fallback and implement it

**Files:**
- Modify: `tests/dg/test_dg_hybrid_production_pw_basis_mpi.f90`
- Modify: `src/common/dg_hybrid_production_pw_basis.f90`

**Step 1: Write the failing MPI case**

Pass a production catalog exactly one nonidentity spatial permutation that splits both synthetic fragments. Its reciprocal rotation may be identity, but the spatial operation must be `[1,3,2,4]`. Require catalog construction to succeed because the production adapter supplies identity itself.

**Step 2: Run RED**

Run:

```bash
python3 tests/dg/run_dg_hybrid_production_pw_basis_mpi.py --build-dir build-hybrid-commit --mpi-ranks 1 2 4 8
```

Expected: FAIL with `no symmetry operation maps whole fragments`.

**Step 3: Add the minimal identity column**

In `build_dg_hybrid_production_pw_basis`, construct candidate arrays with one leading column:

```fortran
candidate_row_action(:,1)=[(i,i=1,global_point_count)]
candidate_reciprocal_rotation(:,:,1)=0d0
do i=1,3
  candidate_reciprocal_rotation(i,i,1)=1d0
enddo
candidate_row_action(:,2:)=normalized_row_action
candidate_reciprocal_rotation(:,:,2:)=reciprocal_rotation
```

Run the existing per-column whole-fragment selector over the candidate arrays. Keep all subsequent permutation, reciprocal-closure, and covariance checks unchanged.

**Step 4: Run GREEN and regress window distribution**

Run:

```bash
python3 tests/dg/run_dg_hybrid_production_pw_basis_mpi.py --build-dir build-hybrid-commit --mpi-ranks 1 2 4 8
python3 tests/dg/run_dg_hybrid_window_distribution_mpi.py --build-dir build-hybrid-commit --mpi-ranks 1 2 4 8
```

Expected: both runners PASS for every requested rank count.

**Step 5: Build and run route contracts**

Run:

```bash
cmake --build build-hybrid-commit -j 8
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_self_consistent_route.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: build succeeds and all contracts PASS.

### Task 2: Validate Si64 and commit the Task 7 fix

**Files:**
- Modify: `src/common/dg_hybrid_production_pw_basis.f90`
- Modify: `src/common/dg_hybrid_window_distribution.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_hybrid_divided_dc_controls.py`
- Modify: `tests/dg/test_dg_hybrid_production_pw_basis_mpi.f90`
- Create: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in`
- Create: `tests/dg/run_dg_hybrid_si64_divided_lcfo.py`

**Step 1: Preserve the previous result and run Si64**

Move the existing result directory to a descriptive failed-result name. Run without a time cutoff:

```bash
OMP_NUM_THREADS=1 python3 tests/dg/run_dg_hybrid_si64_divided_lcfo.py --mpi-ranks 8
```

Expected: normal completion, accepted Wannier90/WF+PW basis, converged divided SCF, exactly one final LCFO solve, and an occupied checkpoint.

**Step 2: Inspect fresh evidence**

Require return code zero, eight MPI ranks, one OpenMP thread, one final LCFO receipt, finite residuals, and a checkpoint. Do not add a post-LCFO density SCF or convergence check.

**Step 3: Stage only Task 7 hunks**

Use `git add -p` for every pre-dirty modified file. Add only the two new Task 7 files. Then run:

```bash
git diff --cached --check
git diff --cached
```

Reject unrelated hunks and preserve every user change.

**Step 4: Commit**

```bash
git commit -m "test(dg): validate divided WF+PW LCFO on Si64"
```
