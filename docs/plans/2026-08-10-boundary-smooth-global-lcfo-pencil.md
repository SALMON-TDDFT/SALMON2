# Boundary-Smooth Global LCFO Pencil Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Use DC-LCFO only as local expansion data and derive the physical 128-state occupied space from a partition-of-unity stitched, full-system-symmetrized generalized eigenproblem.

**Architecture:** Assemble distributed `H/S/rho` tiles from overlapping core-plus-buffer functions with smooth partition weights and their kinetic derivatives.  Average matrix tiles under the full atomic affine group through a deterministic generating set, solve `H C = S C epsilon` with EigenExa, and combine complete occupied blocks with the 256 complete `s+p` channels.

**Tech Stack:** Fortran 2008, MPI/OpenMP, ScaLAPACK/EigenExa, spglib, Wannier90 3.1, Python contracts, CMake clean overlays.

---

## Mandatory discipline

Each task must show a genuine RED, run focused MPI 1/2/4/8 verification,
receive specification and code-quality reviews, resolve every Critical and
Important finding, and pass a clean-first committed-parent prerequisite
overlay.  Preserve normal DC LCFO+EigenExa and generalized-eigenvalue Exp-only
V3 RT.  Never use raw LCFO occupied states as the accepted physical projector.

### Task GP1: Construct a globally normalized smooth partition of unity

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

Write RED fixtures for arbitrary fragment layouts, unequal boundary coverage,
periodic wrap, off-fragment symmetry centers, zero-row ranks, and malformed
coverage.  Require pointwise partition sum one, nonnegative bounded weights,
paired-face value/first-difference continuity, deterministic MPI results, and
bounded tiles.  Implement from global physical IDs and core/buffer geometry;
fragment rank must not affect weights.  Verify/review/overlay and commit GP1.

### Task GP2: Assemble stitched overlap and density matrix tiles

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_metric.f90`
- Modify: `src/gs/main_dft.f90`
- Modify focused metric and route fixtures

Write RED dense references for `S` and `rho` from overlapping basis values,
including duplicate buffer coverage and arbitrary ownership.  Require Hermitian
tiles, electron-count preservation, positive-definite retained overlap,
process-count invariance, and measured memory independent of affine order.
Assemble only row-owned tiles using `sqrt(w_f)` weighted LCFO/local-projector
values.  Reject missing coverage or rank loss.  Verify/review/overlay and
commit GP2.

### Task GP3: Assemble the boundary-correct Hamiltonian

**Files:**
- Create or modify the smallest overlapping-Wannier operator module
- Modify: `src/gs/main_dft.f90`
- Add focused MPI operator fixtures

Write RED references that distinguish the correct kinetic form using
derivatives of `sqrt(w_f) phi_fa` from the incorrect value-only weighting.
Cover local potential, nonlocal pseudopotential, periodic faces, and buffer
cutoff.  Require Hermiticity, face-flux balance, finite energy, and bounded
row-tile storage.  Reuse normal DC operators without modifying their route.
Verify/review/overlay and commit GP3.

### Task GP4: Symmetrize distributed `H/S/rho` under the full affine group

**Files:**
- Modify overlapping-Wannier symmetry/construction modules
- Modify: `src/gs/main_dft.f90`
- Add focused symmetry MPI fixtures

Write RED matrices containing both nonsymmetric error and a symmetric
fragment artifact.  Require generator-streamed group averaging to remove the
former, retain and independently report the latter, match full dense averaging,
and prove full-group covariance from the product/cocycle catalog.  No
`Naffine*Nbasis**2` storage is allowed.  Publish before/after commutators and
workspace.  Verify/review/overlay and commit GP4.

### Task GP5: Solve and select the symmetry-complete occupied space

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_solver.f90`
- Modify: `src/gs/main_dft.f90`
- Add focused generalized-eigenproblem MPI fixtures

Write RED for a small stitched pencil with known eigenpairs, overlap
conditioning, degeneracy at and away from the occupied boundary, non-real
Gamma contamination, and residual failure.  Solve distributed `H C=S C e`
with EigenExa and select exactly 128 states without splitting a complete
cluster/irrep.  Require electron count, `S` orthonormality, residual, density
reconstruction, and deterministic gauge receipts.  Verify/review/overlay and
commit GP5.

### Task GP6: Form the 384-state Wannier seed and publish V3 evidence

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_checkpoint.f90`
- Modify checkpoint, DMN, W90, and route fixtures

Write RED requiring the new 128 occupied states plus all 256 complete `s+p`
channels, exact final rank 384, strict fixed-center/full-affine closure, density
and face-jump receipts, pencil fingerprints, nonzero measured workspace, and
restart validation.  Run constrained Wannier90 and center closure without any
fallback.  Verify on MPI 1/2/4/8, review, clean overlay, and commit GP6.

### Task GP7: Genuine Si64 GS/V3 acceptance

Run ideal undisplaced Si64 on 8 MPI ranks and one OpenMP thread using the clean
full-feature overlay.  Require normal DC convergence, partition continuity,
bounded boundary energy/density jumps, 1536/32/48 affine provenance, an
inversion-containing fixed-center group, stitched pencil rank and residual,
128+256=384, converged MLWFs, center closure, accepted V3, and restart reuse.
Record exact evidence and qualitative small-cell limitations.  Resolve all
Critical/Important findings and commit.

### Task GP8: Polarization-primary LR and long-pulse HHG

Resume generalized-eigenvalue Exp coefficient RT only after GP7.  Use the
accepted undisplaced checkpoint, a long pulse, and the justified large time
step.  Compute spectra from polarization, retain current as secondary, and
produce the requested semilog HHG figure with explicit H2/H4 peak-versus-dip,
background slope, and even/odd suppression analysis.  Run focused MPI 1/2/4/8,
reviews, and the final clean overlay; commit results.

### Task GP9: Final publication

Run all retained route contracts, focused suites, Si64 GS/restart/LR/HHG, and
`git diff --check`.  Confirm normal DC LCFO+EigenExa and Exp-only V3 RT remain
the only accepted routes.  Resolve final review findings, commit, then push
`codex/wpw-s-orthogonal-complement` to both `origin` and `upstream`.
