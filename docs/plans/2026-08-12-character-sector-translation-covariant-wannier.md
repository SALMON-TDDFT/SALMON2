# Character-Sector Translation-Covariant Wannier Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Construct an exact translation-covariant localized Wannier gauge by aligning finite-translation character sectors to the Wannier90 localization reference.

**Architecture:** Canonicalize the finite abelian translation group and its characters, stream the LCFO representation into equal-multiplicity character sectors, align their internal gauges by polar factors, then inverse-transform characters into exact real-space translation orbits. Validate point-cogroup covariance with the existing cocycle without retaining the full affine representation.

**Tech Stack:** Fortran 2008, MPI/OpenMP, BLAS/LAPACK, EigenExa, Spglib, Wannier90, Python test runners.

---

### Task 1: Canonical finite-abelian character table

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: RED**

Add `Z2`, `Z2 x Z2`, and `Z4` fixtures with permuted operation numbering.
Require canonical identity, inverses, minimal generators, element words,
character values, orthogonality, multiplication, conjugate pairing, and a
rank-independent fingerprint. Reject zero order, nonassociative, nonabelian,
duplicate/nonfaithful, missing inverse, integer-overflow, and nonfinite phase
contracts.

Run `python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py` and
require failure because the character-table API is absent.

**Step 2: GREEN**

Implement `build_dg_finite_abelian_character_table` using only `O(|T|)` integer
metadata and `O(|T|^2)` character output explicitly requested by the caller.
Canonicalize independently of input numbering.  Use wide-integer checked
arithmetic and strict total-order keys.

**Step 3: Verify and review**

Run MPI 1/2/4/8, route, obsolete-route, and `git diff --check`. Review group
orientation, canonicalization, overflow, allocation failure, and adverse
tables. Resolve all Critical/Important findings with new REDs.

**Step 4: Overlay and commit**

Run the committed-parent full-feature clean-first overlay, rerun focused tests,
and commit Task 1 only.

### Task 2: Stream translation generators and split character sectors

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_eigenexa_mpi.py`

**Step 1: RED**

Build complex commuting generator representations with repeated character
multiplicity, arbitrary input gauge, and distributed row ownership. Require
the known sector spectrum and equal multiplicity. Reject noncommuting,
nonunitary, wrong-order, rank-losing, split-cluster, and nonfinite generators.

**Step 2: GREEN**

Stream only minimal generator overlaps. Refine their distributed eigenspaces
with EigenExa, producing one sector at a time. Publish identity, unitarity,
commutator, order, multiplicity, Gamma-pairing, fingerprint, and peak-byte
receipts. Do not retain all projectors or all translation operations.

**Step 3: Verify, review, overlay, commit**

Run construction and EigenExa MPI 1/2/4/8. Review eigenspace orientation,
degeneracies, memory scaling, MPI agreement, and normal DC isolation. Resolve
Critical/Important issues, clean-first overlay, and commit.

### Task 3: Align sector gauges to Wannier90

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_w90_mpi.py`

**Step 1: RED**

Use repeated-center, complex-gauge sector fixtures. Require invariance under
sector-frame rotations and Wannier numbering, the closest polar alignment,
complete degenerate-cluster handling, and conjugate-sector Gamma reality.
Reject singular localization links and split clusters.

**Step 2: GREEN**

Project the Wannier90 reference into each sector and align its internal frame
to a deterministic reference using localization link matrices and unitary polar
factors. Store one sector workspace and small internal matrices only.

**Step 3: Verify, review, overlay, commit**

Run W90 MPI 1/2/4/8 and bundled noncommutative DMN tests. Review phase/gauge
invariance, SVD thresholds, spread cost, complex conjugation, and memory.
Resolve Critical/Important issues, clean-first overlay, and commit.

### Task 4: Inverse character transform and exact translation orbits

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_construction_mpi.py`

**Step 1: RED**

Require inverse transforms for `Z2 x Z2` and `Z4` to produce exact translated
orbits with repeated centers/internal channels. Check density, rank,
orthogonality, known translation permutation, real Gamma output, and input
gauge invariance. Corrupt one character phase and require rejection.

**Step 2: GREEN**

Apply the normalized inverse character transform channel by channel. Apply the
unitary result identically to values and gradients. Validate translation
identity/unitarity/products and measured periodic center orbits.

**Step 3: Verify, review, overlay, commit**

Run construction, physical matrices, checkpoint, and Exp-only RT on MPI
1/2/4/8. Review transform normalization, action direction, buffer gradients,
density invariance, and fingerprints. Resolve all Critical/Important issues,
clean-first overlay, and commit.

### Task 5: Point-cogroup cocycle and production integration

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: relevant MPI fixtures

**Step 1: RED**

Require production order: full-affine LCFO proof, fixed-center DMN, Wannier90,
character sectors, sector alignment, inverse transform, point-cogroup cocycle,
center redistribution, V3. Forbid all-translation/all-affine dense storage and
tolerance relaxation.

**Step 2: GREEN**

Integrate the transform and validate how each of 48 point representatives
permutes characters and acts on internal channels. Use the existing translation
cocycle for products. Propagate the gauge through all operators and provenance.

**Step 3: Full focused verification and reviews**

Run fragment symmetry, construction, EigenExa, checkpoint, W90, route,
obsolete-route, normal DC isolation, DMN, and Exp RT on MPI 1/2/4/8. Review
physics, arbitrary fragments/materials, cocycle orientation, memory, and error
collectives. Resolve every Critical/Important finding.

**Step 4: Clean overlay and commit**

Run the full-feature clean-first committed-parent overlay and commit.

### Task 6: Genuine ideal Si64 through V3

Run the undisplaced Si64 case with MPI 8/OMP 1. Require rank 384, 32 characters
of multiplicity 12, fixed-center DMN, Wannier90, exact translation orbits,
48-representative cocycle closure, density/electron preservation, bounded
spread increase and memory, fragment redistribution, and V3 publication.

For any genuine failure, add a minimal RED, fix, rerun focused verification,
perform both reviews, and clean-first overlay before committing. Do not relax
symmetry tolerances. Only after V3 passes, resume polarization-primary linear
response and long-pulse Exp HHG with semi-log spectra and even-order peak/dip
classification.
