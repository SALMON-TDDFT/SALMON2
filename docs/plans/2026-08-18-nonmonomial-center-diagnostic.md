# Nonmonomial Wannier-Center Diagnostic Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Allow physically valid dense point-cogroup representations to proceed while preserving center-orbit diagnostics and all physical symmetry gates.

**Architecture:** Keep the existing center matcher and distributed point-gauge diagnostic. Convert only the successful-unitary nonmonomial outcome from fatal to diagnostic; keep diagnostic contract/unitarity failures fatal and leave downstream owner assignment and physical validation unchanged.

**Tech Stack:** Fortran 2008, MPI, Python route contracts, existing construction MPI fixtures.

---

### Task 1: Add the route RED

1. Update `tests/dg/check_dg_overlapping_wannier_route.py` to require the center diagnostic before owner assignment.
2. Reject the unconditional `error stop 'localized Wannier center orbit failed'`.
3. Require diagnostic failure itself to remain fatal.
4. Run the route checker and observe the expected failure.

### Task 2: Make nonmonomial center action diagnostic

1. In `src/gs/main_dft.f90`, retain the center mismatch message and distributed diagnostic call.
2. If the diagnostic fails, print its reason and stop.
3. If it succeeds, print both the mismatch and dense-unitary diagnostic, then continue.
4. Do not change basis, tolerance, centers, or symmetry receipts.
5. Run route and construction MPI tests.

### Task 3: Verify and run Si64

1. Run W90 MPI 1/2/4/8, construction MPI 1/2/4/8, route, build, and diff checks.
2. Rerun Si64 MPI 8 / OMP 1.
3. Require progression past center ownership and into stitched overlap/Hamiltonian/observable gates.
4. Record the first genuine downstream physical failure without weakening it.

