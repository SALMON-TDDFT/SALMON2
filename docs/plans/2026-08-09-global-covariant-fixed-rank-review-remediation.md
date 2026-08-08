# Global-Covariant Fixed-Rank Review Remediation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Resolve every Critical/Important review finding while retaining exact full-system atomic symmetry and a fixed rank-384 fragment-local Wannier space.

**Architecture:** Build and verify a complete distributed symmetry-closed candidate, then select exact invariant blocks with measured metric, metric-orthonormal occupied coefficients, and a physical localizer.  Use the same point action and phase convention through construction, localization, publication, and ideal-Si64 evidence.

**Tech Stack:** Fortran 2008, MPI, LAPACK/ScaLAPACK, EigenExa, spglib, Python source contracts and physical-evidence checkers.

---

### Task 1: Repair candidate closure semantics

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

1. Add RED cases for a non-first identity operation, mandatory processing of every required seed, minimum-rank orbit crossing, measured projector leakage, metric unitarity, and group closure.
2. Run the construction fixture and confirm the intended failures.
3. Find the identity from both sides of the product table; reject missing or duplicate identities.
4. Prevent minimum-rank exit until all required seeds have been processed and publish a valid required retained rank.
5. Carry point-action phases into each streamed symmetry image.
6. Measure and reject candidate leakage before representation synchronization.
7. Run construction on 1/2/4/8 ranks and `git diff --check`.

### Task 2: Select exact rank with measured algebra

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

1. Add RED cases for nonidentity measured metric, metric-orthonormal occupied coefficients, a physical seed/position localizer, exact target rank, and rejection of a cut irreducible block.
2. Build the occupied Gram inverse square root and verify occupied inclusion.
3. Assemble the physical localizer from distributed projection overlap and periodic moments, then symmetry-average it.
4. Select complete blocks at rank 384, remeasure the selected metric and representation, and gate leakage/unitarity/closure.
5. Run construction, metric, projection, point-group, symmetry, and route fixtures on 1/2/4/8 ranks.

### Task 3: Make localization objective and convergence identical

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_localization.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_localization_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Add a RED case where the sparse pair graph is empty but the complete periodic-spread gradient is nonzero.
2. Remove graph emptiness as a convergence gate; convergence is based only on the symmetry-projected complete gradient.
3. Retain bounded point blocks and verify analytic gradients against real/imaginary finite differences.
4. Run localization on 1/2/4/8 ranks with bounds and floating-point traps.

### Task 4: Align ideal-Si64 and HHG evidence

**Files:**
- Modify: `tests/dg/run_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/check_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/test_si64_harmonic_morphology.py`
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_gs.in`

1. Add RED source/evidence contracts that forbid fragment site-symmetry prerequisites and displaced-structure acceptance.
2. Require full-system identity/group order/inversion, candidate and selected leakage, metric unitarity, closure, and covariance.
3. Require explicit peak/dip/slope output for H2/H4 plus even-to-H3 suppression.
4. Recalibrate localization tolerance only after the repaired strict-SCF path produces a convergence history.
5. Run morphology, route, obsolete-route, and ideal-evidence contracts.

### Task 5: Genuine verification, reviews, commit, and push

1. Perform specification and code-quality reviews; resolve all Critical/Important findings.
2. Commit the reviewed repair.
3. Build a clean committed-HEAD archive with EigenExa `-j1`, then full `-j4`.
4. Run strict undisplaced genuine Si64 GS on 8 MPI ranks through V3.
5. Run field-off, LR, and long-pulse Exp coefficient RT; compute spectra from polarization and generate the semi-log HHG figure.
6. Run final focused 1/2/4/8 verification and repeat both reviews.
7. Commit evidence, push the identical head to origin and upstream, and report diagnostics and figure path.
