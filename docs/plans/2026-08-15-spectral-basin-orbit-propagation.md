# Spectral Basin Orbit Propagation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build the complete row-owned trial frame from one diagonalized representative block per basin symmetry orbit.

**Architecture:** Derive deterministic generator words from basin maps, apply row-owned retained-space generator matrices to one block buffer, and write completed blocks directly into the final distributed frame. Verify complete rank with one distributed Gram reduction.

**Tech Stack:** Fortran 2008, MPI, existing spectral basin and row-owned representation fixtures.

---

### Task 1: Add representative propagation RED

Modify `tests/dg/test_dg_overlapping_wannier_eigenexa_mpi.f90` with a two-basin
orbit whose generator swaps two rank-two subspaces.  Supply only the first
basin's two representative vectors and require the complete identity trial
frame.  Run the EigenExa MPI fixture and observe the missing API failure.

### Task 2: Implement streamed orbit propagation

Modify `src/gs/dc/dg_overlapping_wannier_construction.f90`.  Validate collective
metadata and exact row ownership, derive parent generator words, propagate one
`Nstate x block_rank` buffer, fill row-owned output columns, and verify Gram-I.
Use checked extents, collective allocations, finite gates, fingerprinting, and
a conservative workspace receipt.

### Task 3: Verify and commit

Run EigenExa and construction fixtures on MPI 1/2/4/8, run `git diff --check`,
and commit only the source, focused fixture, and these plan documents.
