# Wannier Center Closure Diagnostics Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Emit actionable center-orbit mismatch and periodic-moment diagnostics while preserving the existing rejection gate.

**Architecture:** Extend the center validator with an optional magnitude payload and enrich only its failure reporting. Production passes its existing measured magnitudes and emits their range before the gate.

**Tech Stack:** Fortran 2008, MPI fixture, Python route checker, CMake.

---

### Task 1: Add the diagnostic RED

**Files:**
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

Require the broken center-orbit fixture to report operation, source, nearest residual, and moment magnitude. Run the MPI fixture and confirm failure.

### Task 2: Add diagnostic reporting

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`

Add the optional magnitude contract, calculate the nearest periodic mismatch on failure, format the detailed error, pass production magnitudes, and print their range before validation.

### Task 3: Verify and commit

Run construction MPI 1/2/4/8, route checker, Release build, and `git diff --check`. Commit only the diagnostic implementation, fixture, and plan documents.

