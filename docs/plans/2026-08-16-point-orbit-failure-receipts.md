# Point-Orbit Failure Receipts Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Report the exact numerical stage that prevents Si64 point-orbit closure.

**Architecture:** Reuse internal scalar state and the existing message output. Preserve all success/failure decisions and add no workspace.

**Tech Stack:** Fortran 2008, MPI, Python source-contract tests, Si64 MPI-8 verification.

---

### Task 1: Add receipt contract RED

Modify `tests/dg/check_dg_overlapping_wannier_route.py` to require distinct
cover, Gram-rank, and cluster-leakage diagnostics. Run it and confirm failure,
then commit the RED.

### Task 2: Add stage-specific failure receipts

Modify `src/gs/dc/dg_overlapping_wannier_w90.f90` so each internal failure sets
the detailed message and the caller does not overwrite it. Run route checks,
Wannier90 MPI 1/2/4/8, production build, and `git diff --check`, then commit.

### Task 3: Re-run Si64

Run the unchanged Si64 MPI-8 input and record the first stage-specific receipt.
Use its numerical values to decide whether the next change belongs in the DC
source frame or in point-orbit clustering.

