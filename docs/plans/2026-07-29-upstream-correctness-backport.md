# Upstream Correctness Backport Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Backport applicable upstream correctness fixes without replacing the local DG/Wannier history.

**Architecture:** Treat upstream commits as patch references rather than merge
parents. Port and test each independent behavior, omit new dependencies, then
restore the in-progress Task 9 work.

**Tech Stack:** Fortran, MPI, CMake, Python focused source/runtime checks.

---

### Task 1: Preserve the in-progress Task 9 tree

Stash tracked and untracked changes, verify a clean HEAD, and record the stash
identity. Do not include Task 9 implementation in the upstream backport commit.

### Task 2: Backport DC input correctness fixes

Port upstream `ebfc4344` and `17bceedd`. Add focused input-contract tests, run
RED/GREEN verification, and confirm normal DC LCFO/EigenExa routing is unchanged.

### Task 3: Backport non-orthogonal-cell and NLCC fixes

Port upstream `901c3527`, `3d4a8781`, `a8951f9d`, `7299afdb`, and `0ccd23b6`
where their local counterparts exist. Test lower-triangular skew cells, full
lattice periodic replicas, Cartesian current output, and valid automatic
decompositions.

### Task 4: Backport TBmBJ correctness fixes

Port upstream `8daf081c` and `6342f489`. Test zero-density handling, full-cell
averaging, and the unbracketed BR root guard.

### Task 5: Backport RT restart fixes

Port upstream `891efc23` and `9a2e8904`. Test restart directory formatting and
preservation of the ground-state reference density.

### Task 6: Review, verify, and commit

Run focused tests, a clean-first overlay build, specification review, and code
quality review. Resolve Critical/Important findings and commit the backport with
upstream provenance.

### Task 7: Restore Task 9

Restore the exact saved Task 9 changes, verify status and hashes, then resume
the Si64 checkpoint-reuse diagnosis. Do not begin Task 10 until its gate passes.
