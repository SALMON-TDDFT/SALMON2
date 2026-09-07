# Wannier90 Block-Monomial SAWF Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Feed Wannier90 a symmetry representation that moves complete Wannier-center blocks rather than merely preserving the total retained subspace.

**Architecture:** First add a fail-closed analyzer for block-monomial generator actions in a center-labelled trial frame. Integrate it before DMN writing and only emit the measured representation when its leakage, internal unitarity, and covariance gates pass. Keep all generator matrices streamed.

**Tech Stack:** Fortran 2008, MPI, Wannier90 SAWF DMN, existing DG MPI fixtures.

---

### Task 1: Block-monomial analyzer

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

1. Add a repeated-center synthetic RED requiring one center permutation and one nontrivial internal unitary block.
2. Add leakage and ambiguous-center REDs.
3. Implement streamed center-block assignment, leakage, and internal-unitarity receipts.
4. Verify construction MPI 1/2/4/8.

### Task 2: Pre-Wannier trial-frame feasibility gate

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

1. Compute periodic centers/confidence for the trial frame before DMN writing.
2. Analyze each affine generator one at a time.
3. Publish worst leakage, ambiguity, and internal-unitarity receipts.
4. Reject before Wannier90 if the frame is not block-monomial.

### Task 3: DMN production binding

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/lcfo_wannier_sawf_dmn.f90`
- Modify: `tests/dg/check_sawf_dmn_format.py`

1. Write the measured block-monomial action as `D_wann`.
2. Retain the exact `D_band A=A D_wann` validator.
3. Add a mismatched trial-frame/representation RED.
4. Verify focused SAWF and route tests.

### Task 4: Si64 end-to-end verification

1. Build the lightweight executable.
2. Run Si64 with 8 MPI ranks through Wannier90 and the immediate center gate.
3. Require generator covariance and full center-orbit closure.
4. Run focused MPI 1/2/4/8, route, and `git diff --check`.
