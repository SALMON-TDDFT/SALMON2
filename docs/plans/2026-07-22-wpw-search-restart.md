# WPW Search-history Restart Comparison Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Test whether accumulated LOBPCG search history causes the fixed-H residual plateau.

**Architecture:** Add a default-on input control, propagate an optional retain-history flag through the algebra step, and clear the search block after each update only when explicitly disabled.

**Tech Stack:** Fortran, MPI, Python source contracts, CMake.

---

### Task 1: Add and propagate the control

1. Add RED checks for declaration, namelist, default `y`, broadcast, log, and
   y/n validation.
2. Add RED matrix-free fixture coverage showing default history retention and
   explicit restart follow different, deterministic paths.
3. Implement an optional retain-history argument defaulting true.
4. In restart mode, clear `search` after the reduced update; change no other
   solver state or convergence rule.

### Task 2: Verify and compare

1. Run matrix-free, fixed-H, generalized-algebra, occupied-W tests and full build.
2. Request code review and resolve blockers.
3. Commit only intended changes.
4. Run fresh B=6 with preconditioner off and search history off.
5. Compare residual and metric-mode histories with Task 16 and record results.
