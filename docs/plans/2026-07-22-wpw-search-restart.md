# WPW Search-history Restart Comparison Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

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

## Result (2026-07-22)

Implemented in `54f334b`.  The default remains search-history retention; explicit
`yn_dg_wpw_search_history='n'` makes every post-update search block zero without
evaluating the discarded history update.  Focused matrix-free/fixed-H/generalized-
algebra/occupied-W tests, the MPI fixture, full MPI/EigenExa build, and code review
reported no blocking issue.

The fresh no-preconditioner complete-restart run is
`stage2d_wpw_runs/20260722_task17_search_restart_no_precondition_b6/run.log`.
It reproduced the same normalized B=6 seed and ended at the unchanged 160-inner
limit (`info=40`).  Against Task 16 (history retained), the state residuals were:

| inner | retained history occupied / extra | complete restart occupied / extra |
|---:|---:|---:|
| 32 | `1.9667E-03 / 4.0432E-03` | `4.1285E-03 / 7.0474E-03` |
| 96 | `3.6079E-04 / 1.1093E-03` | `1.1676E-03 / 1.8784E-03` |
| 160 | `2.5607E-04 / 2.0289E-04` | `8.0773E-04 / 1.2379E-03` |

Thus accumulated search history is not the plateau cause; removing it makes the
final occupied residual about 3.15 times and the extra residual about 6.10 times
worse.  Complete restart also keeps the reduced space less redundant (effective
rank 305 rather than 213 at inner 160), yet converges more slowly.  Search-space
rank loss alone therefore does not diagnose the plateau.  Retain history for now;
the next diagnosis should target the retained Rayleigh--Ritz/operator recurrence
or the quality of the correction directions, while keeping the physical basis,
cutoff, tolerance, and publication gates unchanged.
