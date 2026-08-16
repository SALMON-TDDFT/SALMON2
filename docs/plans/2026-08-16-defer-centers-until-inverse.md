# Defer Periodic Centers Until Character Inversion Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Avoid projecting a translation-changing position operator into one character sector and measure centers only after the complete inverse character transform.

**Architecture:** Preserve W90/LCFO anchoring and the existing streamed per-character pipeline. Delete only the production single-sector position tuple/joint canonicalization block and its dead workspace/provenance fields.

**Tech Stack:** Fortran 2008, MPI, EigenExa, Wannier90, Python route contracts.

---

### Task 1: Add the production-route RED

Modify `tests/dg/check_dg_overlapping_wannier_route.py` to reject production
calls to `build_dg_sector_periodic_position_tuple` and
`jointly_canonicalize_dg_sector_periodic_position_gauge`, while preserving the
post-inverse center measurement and affine orbit validation. Run and confirm
failure, then commit.

### Task 2: Remove the invalid production stage

Modify `src/gs/main_dft.f90` to allocate/fill spatial IDs, materialize the
anchored reference sector, and continue directly to translation action
preparation. Remove dead point-overlap, tuple, weighted/canonical reference,
joint-center allocations, declarations, deallocations, and fingerprint mixing.

Run route checks, construction and Wannier90 MPI fixtures on 1/2/4/8 ranks,
production build, and `git diff --check`; then commit.

### Task 3: Re-run Si64

Run the unchanged Si64 input on 8 MPI ranks. Confirm Wannier90 convergence,
per-character alignment, inverse transform completion, post-gauge factored
point-cogroup proof, and final center-orbit validation. Record the next genuine
failure if one appears.

