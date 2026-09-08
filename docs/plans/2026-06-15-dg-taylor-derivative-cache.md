# DG Taylor Derivative Cache Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Reduce Taylor4-PC overhead by reusing derivative row/map preparation across repeated Hamiltonian applications.

**Architecture:** Keep the numerical Taylor recurrence unchanged: each order still applies the current Hamiltonian-derived linear operator to the previous term. Cache only the derivative preparation data that depends on the fragment layout, local H blocks, spin, and required compact row set, not on the coefficient values.

**Tech Stack:** Fortran 90, existing DG fragment block-sparse runtime, MPI coefficient fetch path.

---

### Task 1: Cache Derivative Compact Row Metadata

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_derivative.f90`

**Steps:**
1. Add saved cache keys for spin, fragment count, local block counts, and nonlocal block counts.
2. Reuse `needed_row_ids`, `needed_row_pos`, and `output_row_needed` when the key matches.
3. Rebuild the cache only when the key changes.
4. Preserve all existing validity checks and fatal diagnostics.

### Task 2: Cache Fragment Compact Maps

**Files:**
- Modify: `src/rt/dg/rt_dg_integrator_derivative.f90`

**Steps:**
1. Reuse `frag_map_lid`, `frag_map_pos`, and `frag_map_count` after the row metadata cache is valid.
2. Invalidate the fragment map cache whenever the row metadata cache changes.
3. Keep fallback `build_compact_basis_map` behavior intact.

### Task 3: Verify Build

**Files:**
- Verify: existing build directory

**Steps:**
1. Run the narrowest available build command for the current tree.
2. Report compiler failures or warnings relevant to the touched file.
