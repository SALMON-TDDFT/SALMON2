# Global-Covariant Fragment-Local Wannier Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Define one full-system symmetry-covariant Wannier gauge and distribute only its core/buffer restrictions so arbitrary fragment boundaries preserve exact instantaneous symmetry through GS, V3, and Exp coefficient RT.

**Architecture:** Full-system affine operations act on global grid IDs and may split one source fragment among several owners.  A distributed construction builds and verifies a fixed-rank symmetry-closed occupied-plus-seed subspace before fragment-local localization; the same global representation and phase convention are then used for localization and `S/H/X/V` publication.  No full-system real-space KS/Wannier gather or fragment site-symmetry prerequisite is permitted.

**Tech Stack:** Fortran 2008, MPI, LAPACK/ScaLAPACK, EigenExa, spglib, complex generalized-overlap algebra, Python source contracts and MPI fixture runners.

---

### Task 1: Freeze pointwise affine ownership and phase semantics

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_symmetry.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_fragment_symmetry_mpi.f90`
- Modify: `tests/dg/check_dg_fragment_symmetry_production.py`

**Step 1: Write the RED fixtures**

Add a synthetic periodic grid where an inversion center lies outside each fragment and one source fragment is split between two target ranks.  Require a deterministic map from every source global grid ID to `(owner_rank, owner_local_index, lattice_wrap)`.  Add a complex Bloch-phase case requiring the product phase for `g*h` to equal the sequential `h` then `g` phase.

**Step 2: Verify RED**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_fragment_symmetry_mpi.py
```

Expected: FAIL because the pointwise affine map does not publish lattice wraps/phases and the production helper is incomplete.

**Step 3: Implement the pointwise map**

Add a pure helper that consumes global grid dimensions, rank-local physical IDs, integer rotation, fractional translation, and tolerance.  Convert each operation to an exactly commensurate integer grid action, map every local point, resolve the mapped physical ID against the collective owner table, and return owner/local-index/wrap arrays.  Permit different target owners for different points of one source fragment.

Generalize `assemble_dg_distributed_basis_symmetry_overlap` to fetch mapped values owner by owner in bounded batches and multiply by the returned periodic phase.  Reject missing/duplicate owners, noncommensurate operations, nonfinite phases, and inconsistent collective shapes.

**Step 4: Focused verification**

Run both MPI fixtures on 1/2/4/8 ranks, `python3 tests/dg/check_dg_fragment_symmetry_production.py`, and `git diff --check`.

**Step 5: Specification and quality reviews**

Confirm that no test assumes a fragment permutation or an internal symmetry center.  Review integer overflow, negative modulo, operation composition order, phase convention, MPI failure consistency, and bounded allocation.  Resolve every Critical/Important finding and rerun Step 4.

**Step 6: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_construction.f90 \
  src/gs/dc/dg_overlapping_wannier_symmetry.f90 tests/dg
git commit -m "feat(dg): map full-system symmetry pointwise across fragments"
```

### Task 2: Build a fixed-rank global symmetry-closed subspace

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_types.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write the RED fixture**

Construct fragment-local occupied vectors plus atomic seeds whose unsymmetrized retained space loses inversion.  Require a wished-for `build_dg_global_covariant_retained_subspace` to return exactly the requested target rank, contain the occupied projector, close under inversion and a noncommuting second operation, and give identical fingerprints on 1/2/4/8 ranks.  Add an insufficient-rank case that must fail rather than drop an operation.

**Step 2: Verify RED**

Run the construction MPI fixture and confirm failure for the missing API.

**Step 3: Implement distributed closure**

Generate pointwise symmetry images of occupied states and projection seeds only on rank-local core/buffer storage.  Assemble their distributed Gram, occupied-projector, localization, and representation matrices.  Use rank-revealing Hermitian algebra to select exactly `target_rank` vectors while preserving the occupied projector and complete group-invariant blocks.  Tie degeneracies deterministically using the ordered seed manifest and global operation order.

Do not allocate `N_grid_global x N_candidate`, `N_grid_global x N_W`, or a full-system KS array on one rank.  Fail if the smallest symmetry-closed space containing the occupied states exceeds the requested target rank.

**Step 4: Focused verification**

Run construction, metric, projection, point-group, and route fixtures on 1/2/4/8 ranks plus `git diff --check`.

**Step 5: Specification and quality reviews**

Review occupied inclusion, degeneracy handling, rank thresholds, deterministic ordering, collective memory scaling, and absence of symmetry averaging of the ionic structure.  Resolve all Critical/Important findings and rerun.

**Step 6: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_construction.f90 \
  src/gs/dc/dg_overlapping_wannier_types.f90 tests/dg
git commit -m "feat(dg): select a fixed-rank symmetry-closed Wannier space"
```

### Task 3: Carry one global representation through fragment-local localization

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_localization.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_symmetry.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_localization_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_fragment_symmetry_mpi.f90`

**Step 1: Write the RED fixtures**

Use a dense global representation whose inversion splits every fragment.  Require localization to reduce periodic spread while its transformation commutes with every `D(g)`.  Verify metric unitarity, group closure, center-orbit/phase consistency, and identical results on 1/2/4/8 ranks.  Add a lower-spread unconstrained transform that must be rejected because it breaks inversion.

**Step 2: Verify RED**

Run localization and fragment-symmetry MPI fixtures.  Expected: FAIL because localization currently consumes only a fragment site-stabilizer representation.

**Step 3: Implement the minimum change**

Pass the global representation/product table from construction into localization.  Build Riemannian gradients in the commutant by a Reynolds projection over the complete global group, transport the conjugate direction in the same commutant, and apply one common gauge to every distributed restriction.  Validate symmetry after every accepted Armijo step and roll back collectively on violation.

Remove fragment-site symmetry as a prerequisite for localization.  The instantaneous full-system group remains the sole symmetry authority.

**Step 4: Focused verification and reviews**

Run localization and symmetry fixtures on 1/2/4/8 ranks with bounds checking and floating-point traps.  Review conjugation/composition conventions, commutant projection, deterministic reductions, rollback, and convergence.  Resolve all Critical/Important findings and rerun.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_overlapping_wannier_localization.f90 \
  src/gs/dc/dg_overlapping_wannier_symmetry.f90 tests/dg
git commit -m "fix(dg): localize one global-covariant Wannier gauge"
```

### Task 4: Integrate global-covariant construction into OW GS and V3

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_checkpoint.f90`
- Modify: `tests/dg/check_dg_fragment_symmetry_production.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `tests/dg/test_dg_overlapping_wannier_checkpoint_mpi.f90`

**Step 1: Write RED production contracts**

Require the full-system catalog and pointwise action before retained-subspace selection, the global representation before localization, and measured post-localization representation before SCF/operator publication.  Require V3 evidence for group order, inversion presence, subspace leakage, raw/final unitarity, closure, and pre/post projection covariance.  Forbid fragment-permutation or site-stabilizer promotion as the full-system gate.

**Step 2: Verify RED**

Run the production contracts and checkpoint fixture; confirm the expected failures.

**Step 3: Implement production wiring**

Replace representative-fragment replication and checkpoint-time coset reconstruction with the global-covariant retained-subspace result.  Store only local core/buffer restrictions and the dense global `D(g)`.  Pass the rank-consistent exact closure into SCF and V3.  Measure the final distributed action and allow projection only for residuals already within the configured correction bound.

Invalidate older checkpoints whose provenance does not bind the full-system operation list, pointwise ownership/phase fingerprint, and global representation fingerprint.

**Step 4: Focused verification**

Run route, checkpoint, construction, localization, symmetry, SCF, operator, observable, and coefficient-RT fixtures on 1/2/4/8 ranks; run the obsolete-route removal checker and `git diff --check`.

**Step 5: Specification and quality reviews**

Confirm the retained path is only OW GS → V3 → generalized-eigen Exp RT, with normal DC LCFO/EigenExa untouched.  Review restart compatibility, publication ordering, correction limits, failure paths, and memory accounting.  Resolve all Critical/Important findings and rerun.

**Step 6: Commit**

```bash
git add src/gs/main_dft.f90 src/gs/dc/dg_overlapping_wannier_checkpoint.f90 tests/dg
git commit -m "fix(dg): publish globally covariant fragment-local Wanniers"
```

### Task 5: Genuine Si64 inversion and clean-first prerequisite verification

**Files:**
- Modify only as required by verified failures
- Record evidence using the existing focused Si64 output convention

**Step 1: Run clean-first parent-prerequisite build**

Create a clean source archive from committed `HEAD`, overlay only the reviewed task diff, and configure Release with `FC=mpifort`, MPI, ScaLAPACK, EigenExa, and spglib enabled, Wannier90 disabled, and the required EigenExa BOZ compatibility flag.  Build `eigenexa-project-build -j1` before the full `-j4` build.

**Step 2: Run genuine Si64 GS**

Run the ideal undisplaced Si64 case on 8 MPI ranks.  Require localization convergence, OW-SCF convergence, occupied-subspace inclusion, full-system inversion in the catalog, subspace leakage within tolerance, metric-unitary and closed global representation, covariance of `S/H/X/V`, `global_inversion_promoted=T`, and successful V3 publication.

Do not accept loose-SCF results as final physical evidence.  The earlier loose run may be used only to diagnose failures.

**Step 3: Audit memory and scaling**

Require reported storage to remain local-grid/buffer plus dense Wannier-space matrices.  Reject any full-system real-space KS/Wannier gather or point-group-factor expansion of the retained rank.

**Step 4: Reviews and rebuild**

Perform specification and code-quality reviews of the genuine log and diff.  Resolve every Critical/Important finding.  After any source change, repeat the clean-first build and focused Si64 run.

**Step 5: Commit**

Commit only reviewed fixes or evidence-contract changes with a focused message.

### Task 6: Polarization-derived LR/HHG, final verification, commit, and push

**Files:**
- Modify: `tests/dg/run_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/check_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/test_si64_harmonic_morphology.py`
- Create: semi-log HHG figure under the focused result directory

**Step 1: Run accepted coefficient RT**

Use only the accepted V3 checkpoint and generalized-eigenvalue Exp coefficient RT.  Run field-off, weak impulse/linear response, and a long laser pulse.  Use the largest `dt` that passes a smaller-step comparison.  Record polarization as the primary observable and current as secondary evidence.

**Step 2: Compute spectra from polarization**

Apply the documented window to polarization, compute LR/HHG spectra from it, and produce a semi-log harmonic-order plot.  Do not use current as the primary spectrum and do not introduce lattice displacement or phonons into the ideal inversion test.

**Step 3: Apply physical gates**

Require field-off stationarity, cubic-axis agreement, small transverse polarization, a clear third-harmonic peak, and suppression of even harmonics.  Classify each even order explicitly as a peak or dip relative to neighboring bins.  Interpret amplitudes qualitatively because Si64 is too small for quantitative convergence.

**Step 4: Final focused and clean-first verification**

Run all retained-route contracts and MPI fixtures, the genuine GS and RT/HHG checkers, `git diff --check`, and one final clean-first parent-prerequisite overlay build.  Perform final specification and code-quality reviews and resolve all Critical/Important findings.

**Step 5: Commit and push**

Commit the verified LR/HHG runner and evidence.  Verify both remotes, push `codex/wpw-s-orthogonal-complement` to `origin`, and push the identical commit to the configured upstream branch.  Report commit IDs, exact verification commands, inversion and harmonic diagnostics, and the absolute path to the semi-log HHG figure.
