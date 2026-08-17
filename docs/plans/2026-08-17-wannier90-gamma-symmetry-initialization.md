# Wannier90 Gamma Symmetry Initialization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Preserve complex symmetry representations when initializing a Gamma-only symmetry-adapted Wannier90 calculation.

**Architecture:** Patch Wannier90's projection-path selection so `lsitesymmetry` forces the existing complex SVD and initial-U symmetrization path. Keep the real Gamma fast path for calculations without site symmetry.

**Tech Stack:** CMake script-mode patching, Fortran, MPI, Wannier90, Python regression checks.

---

### Task 1: Add the failing branch regression

Modify `tests/dg/check_sawf_dmn_format.py` to require the patched condition
`.not. gamma_only .or. lsitesymmetry`. Run it and observe failure.

### Task 2: Patch projection-path selection

Extend `cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake`
to replace the original Gamma branch condition. Re-run the regression and
inspect the patched `overlap.F90`.

### Task 3: Verify and run Si64

Rebuild the Wannier90 library and SALMON, run W90 MPI 1/2/4/8, DMN and route
tests, then rerun Si64 with 8 MPI ranks. Preserve diagnostics if the covariance
gate still fails.

