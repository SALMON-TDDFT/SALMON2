# Wannier90 Symmetry-Constrained Monotonic Backtracking Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the unreliable parabolic line search in the site-symmetry route with bounded-memory monotonic backtracking based on the actual spread.

**Architecture:** Patch Wannier90's existing line-search block and reuse its saved U/M state for every retry.  Keep upstream behavior for non-symmetry and explicit fixed-step calculations, and retain the last accepted step across iterations.

**Tech Stack:** CMake source patching, Fortran 2008, Wannier90 3.1, Python source-route tests.

---

### Task 1: Add failing route contracts

**Files:**
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

Require a site-symmetry-only retry loop, state restoration before retry,
spread-based acceptance, machine-precision lower bound, accepted-step reuse,
and CG reset after failure.  Require removal of temporary slope and matrix-form
diagnostics.  Run the test and observe the intended failure.

### Task 2: Implement minimal backtracking

**Files:**
- Modify: `cmakefiles/Builder/patches/apply_wannier90_generator_symmetry.cmake`

Patch the line-search block to save U/M once, restore before each trial, halve
the retained step on spread increase, accept the first non-increasing spread,
or restore the original state and take zero step at the precision floor.  Do
not allocate another U/M-sized array.  Run the route test to green.

### Task 3: Rebuild and verify generated source

Apply the patch to the bundled Wannier90 source, rebuild its executable/library
and SALMON, then inspect the generated retry/rollback order.

### Task 4: Short Si64 validation

Restart from the established 934.321851 Å² checkpoint for a short segment.
Verify every accepted spread is non-increasing, record step attenuation,
runtime and memory, and stop before a long production run.

