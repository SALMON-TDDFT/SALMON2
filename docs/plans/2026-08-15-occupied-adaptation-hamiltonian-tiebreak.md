# Hamiltonian-resolved occupied adaptation implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Resolve a group-averaged occupied-projector rank boundary with a covariant physical Hamiltonian, while preserving exact rank, symmetry, MPI determinism, and bounded memory.

**Architecture:** Extend the distributed group-average eigensolver with an optional occupied-frame Hamiltonian. The projector remains the primary discriminator; only its boundary-degenerate eigenspace is diagonalized with the repeated-orbit Hamiltonian. Production constructs the Hamiltonian in the translation-adapted occupied frame from the smoothly composed LCFO eigenfunctions and their eigenvalues, then retains all existing affine closure and density gates.

**Tech Stack:** Fortran 2008, MPI, EigenExa, LAPACK, Python source contracts, Si8/MPI8 monitored integration.

---

### Task 1: Lock the boundary-block Hamiltonian algebra

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_eigenexa_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_eigenexa_mpi.py`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`

**Step 1: Write the failing GREEN fixture**

Add an EigenExa case whose averaged-projector spectrum has a two-dimensional block crossing the requested boundary. Supply a Hermitian occupied-frame Hamiltonian which separates that block. Assert success at the requested rank, equality with a dense reference projector, finite secondary edges and gap, and the same fingerprint on MPI 1/2/4/8.

**Step 2: Run RED**

```bash
SALMON_GP5_OVERLAY=/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/build-mpi-eigenexa-wannier-lib \
python3 tests/dg/run_dg_overlapping_wannier_eigenexa_mpi.py
```

Expected: FAIL because the API still rejects the primary split block.

**Step 3: Add the optional Hamiltonian contract**

Extend the group-average routine and cocycle wrapper with optional `occupied_hamiltonian`, secondary selected/rejected edges and gap, and primary boundary dimension. Require Hermiticity, finiteness, raw replicated-payload agreement, and matching `noccupied x noccupied` shape before shape-dependent collectives.

**Step 4: Implement the boundary projection**

Collectively determine the primary boundary block from adjacent projector gaps. Gather only its orbit-Gram eigenvectors. With projector eigenvalues `lambda` and repeated occupied Hamiltonian `E`, form

```text
H_ab = sqrt(lambda_a*lambda_b) * v_a^H E_orbit v_b
```

where `E_orbit` is block diagonal with one occupied Hamiltonian per point operation. Hermitize, diagonalize, measure residuals, and reject if the secondary requested boundary remains tolerance-degenerate. Rotate primary orbit eigenvectors by the secondary eigenvectors before reconstruction.

**Step 5: Run GREEN and adverse cases**

Add non-Hermitian, nonfinite, rank-disagreeing, and still-degenerate Hamiltonian cases. Every failure must be collective and leave outputs deallocated.

**Step 6: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_construction.f90 \
  tests/dg/test_dg_overlapping_wannier_eigenexa_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_eigenexa_mpi.py
git commit -m "feat(dg): resolve occupied projector boundary with Hamiltonian"
```

### Task 2: Add checked memory and provenance receipts

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_eigenexa_mpi.f90`

**Step 1: Write receipt REDs**

Require nonzero Hamiltonian provenance, primary boundary dimension, secondary edges/gap, eigensystem residual, and conservative peak workspace. Add unsafe extent and quantization-bound rejection cases.

**Step 2: Run RED**

Run the EigenExa MPI fixture and confirm the missing receipts fail.

**Step 3: Implement checked accounting**

Before allocation, checked-add the boundary eigenvector slab, secondary dense matrix/eigenvectors/eigenvalues, one gathered orbit vector, and reconstruction vector. Reduce rank-local receipt validity before returning, reduce peak bytes with MPI MAX, and bind Hamiltonian raw payload plus secondary spectrum into the returned fingerprint.

**Step 4: Verify and commit**

Run MPI 1/2/4/8 and `git diff --check`, then commit Task 2 files.

### Task 3: Construct the production occupied Hamiltonian

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Write production-route RED**

Require production to construct a Hamiltonian from actual smoothly composed occupied LCFO rows and `lcfo_retained_eigenvalues(1:nstate)`, transform it into the translation-adapted frame, and pass it to the cocycle-aware adaptation. Forbid treating a bare eigenvalue diagonal as already expressed in the adapted frame.

**Step 2: Run RED**

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: FAIL because production passes no Hamiltonian.

**Step 3: Stream the overlap and form the small operator**

Before releasing `lcfo_occupied_core`, compute

```text
Q_ni = <lcfo_n | translation_adapted_i>
H_ij = sum_n conj(Q_ni) epsilon_n Q_nj
```

using weighted local overlaps plus one `nstate*nstate` MPI reduction. Check the MPI count, safe magnitude, finiteness, Hermiticity, allocation consensus, and a frame/eigenvalue-bound fingerprint. This is a small `nstate x nstate` matrix, never an `Ngrid x Nstate` duplicate.

**Step 4: Connect and publish diagnostics**

Pass the operator to the cocycle-aware group average. Print primary block dimension, secondary selected/rejected edges and gap, and bind them into checkpoint provenance.

**Step 5: Verify and commit**

Run construction, EigenExa, route, and checkpoint MPI 1/2/4/8 tests plus `git diff --check`; commit Task 3 files.

### Task 4: Re-run monitored Si8 and classify the next boundary

**Files:**
- Modify only if a genuine new defect is reproduced.
- Store evidence under a fresh `/tmp` directory, not the repository.

**Step 1: Build committed full-feature source**

Configure MPI, EigenExa, Wannier90, SPGLIB, and GNU 15 `-fallow-invalid-boz`. Build EigenExa serially, then SALMON with `-j2`.

**Step 2: Run monitored Si8/MPI8**

```bash
python3 tests/dg/run_si8_overlapping_wannier_memory.py \
  /tmp/salmon-wpw-full/salmon /tmp/si8-hamiltonian-adaptation-YYYYMMDD \
  --interval 10
```

Keep eight MPI ranks, one OpenMP/BLAS thread per rank, `nice 15`, the 3 GiB/rank ceiling, and the 8 GiB available-memory floor.

**Step 3: Acceptance**

Require the run to pass the former `0.3178877` boundary with a finite positive secondary gap, unchanged rank 16, bounded density/closure receipts, and stable memory. Record the next diagnostic or Wannier90 result.

**Step 4: Final verification**

Freshly run:

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
SALMON_GP5_OVERLAY=/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/build-mpi-eigenexa-wannier-lib \
python3 tests/dg/run_dg_overlapping_wannier_eigenexa_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
git diff --check
```

Do not claim completion unless every applicable command exits zero and the Si8 log demonstrates that the previous boundary was crossed.

