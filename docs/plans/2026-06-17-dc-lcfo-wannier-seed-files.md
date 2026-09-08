# DC-LCFO Wannier Seed Files Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add the first DC-side Wannier90 export step by writing `.win` and `.eig` seed files from the Flux-LCFO diagonalization result.

**Architecture:** Keep Wannier setup in `&dc`, because the Wannier basis is generated during the GS/DC export phase. The first implementation only writes deterministic seed files; later steps will add SALMON-generated projection overlaps (`.amn`) and k-neighbor overlaps (`.mmn`) before calling Wannier90.

**Tech Stack:** Fortran namelist input, SALMON DC-LCFO-Flux output, Wannier90 text seed files.

---

### Task 1: Write Wannier Seed Files

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Steps:**
1. Add `yn_dc_lcfo_wannier`, `wannier_projection`, `wannier_num_wann`, `wannier_num_iter`, `wannier_dis_froz_max`, and `wannier_dis_win_max` to the `output` subroutine imports.
2. Add `write_wannier_seed_files` under `output`.
3. Write `<sysname>.win` in `data_dcdft/total/` with one Gamma point, unit cell, atom positions, `num_bands`, `num_wann`, windows, and projections.
4. Write `<sysname>.eig` from `esp_tot`.
5. Build with `USE_WANNIER90=ON`.

### Task 2: Keep The Interface Honest

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Steps:**
1. Print a concise `[DC-LCFO-WANNIER]` message when seed files are written.
2. Do not call Wannier90 until `.amn/.mmn` generation exists.
3. Keep DG-RT reading unchanged.

