# DG+PW Orthogonalized Mixed-Basis Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make `DG+PW` use only the orthogonalized mixed basis during runtime and compute fragment-plane-wave couplings from `core+buffer` support integrals.

**Architecture:** Keep raw fragment and PW objects only as ingredients for building the mixed generalized eigenproblem. After startup diagonalization, treat `mixed_transform` / `coef_mix` as the sole runtime representation and synchronize any raw backing arrays from it in one controlled place. Keep axis-wise support rules: fragmented axes use `core+buffer`, non-fragmented axes use fragment-box periodic wrapping.

**Tech Stack:** Fortran 90, MPI communicator helpers, DG fragment RT modules, plane-wave hybrid basis path, CMake build, SALMON regression runs.

---

### Task 1: Lock startup onto mixed-basis initialization

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-phi-box-cache-phase-a/src/rt/dg/rt_dg_fragment_hamiltonian.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-phi-box-cache-phase-a/src/rt/dg/rt_dg_plane_wave.f90`

**Step 1: Add a failing guard**

Require that `DG+PW` startup ends with `mixed_basis_ready=.true.` after mixed-basis setup.

**Step 2: Route startup through diagonalization**

Replace any raw `prepare_mixed_basis_startup` runtime fallback in the startup path with unconditional `diagonalize_mixed_basis(...)` for `use_plane_wave_basis` cases.

**Step 3: Run build**

Run: `cmake --build /Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-phi-box-cache-phase-a/build -j4`

**Step 4: Verify startup log**

Run a small `DG+PW` case and verify the log reports mixed-basis diagonalization and `mixed_basis_ready`.

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_fragment_hamiltonian.f90 src/rt/dg/rt_dg_plane_wave.f90
git commit -m "refactor: require mixed-basis startup for dg+pw"
```

### Task 2: Remove raw runtime current path for DG+PW

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-phi-box-cache-phase-a/src/rt/dg/rt_dg_observables.f90`

**Step 1: Write the failing assertion path**

Add a guarded branch so `DG+PW` current calculation errors if it would use raw fragment/PW current while `mixed_basis_ready` is expected.

**Step 2: Implement mixed-only runtime logic**

Update current evaluation to use only mixed-basis representations when plane waves are active and mixed basis is ready.

**Step 3: Include existing current components**

Preserve:
- paramagnetic contribution
- diamagnetic `A*rho`
- nonlocal current contribution

**Step 4: Build and run**

Run the same `DG+PW` smoke case and ensure no raw-path fallback is hit.

**Step 5: Commit**

```bash
git add src/rt/dg/rt_dg_observables.f90
git commit -m "refactor: use mixed-basis current path for dg+pw"
```

### Task 3: Remove raw runtime density path for DG+PW

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-phi-box-cache-phase-a/src/rt/dg/rt_dg_density_reconstruct.f90`

**Step 1: Add guard for mixed-basis requirement**

If plane waves are enabled and mixed basis is expected, fail instead of falling back to raw fragment/PW density handling.

**Step 2: Use mixed-basis reconstruction only**

Ensure density reconstruction uses `mixed_transform`, `coef_mix`, and cached PW coefficients consistently.

**Step 3: Run targeted smoke test**

Check that `Ne_raw` is finite and stable in a short `DG+PW` run.

**Step 4: Commit**

```bash
git add src/rt/dg/rt_dg_density_reconstruct.f90
git commit -m "refactor: require mixed-basis density reconstruction for dg+pw"
```

### Task 4: Make FP overlap integrals axis-wise support aware

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-phi-box-cache-phase-a/src/rt/dg/rt_dg_plane_wave.f90`
- Reference: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-phi-box-cache-phase-a/src/rt/dg/rt_dg_fragment.f90`

**Step 1: Define axis-wise support rule**

Use fragment support metadata so fragmented axes integrate over `core+buffer`, while non-fragmented axes wrap within the fragment box.

**Step 2: Update fragment-PW overlap**

Modify `compute_fragment_pw_overlap` to use the axis-wise support rule consistently.

**Step 3: Build and verify**

Check logs or probes to confirm support widths are correct by axis.

**Step 4: Commit**

```bash
git add src/rt/dg/rt_dg_plane_wave.f90 src/rt/dg/rt_dg_fragment.f90
git commit -m "refactor: use axis-wise support in dg+pw overlap integrals"
```

### Task 5: Make FP Hamiltonian couplings axis-wise support aware

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-phi-box-cache-phase-a/src/rt/dg/rt_dg_plane_wave.f90`

**Step 1: Update FP off-diagonal integral construction**

Evaluate fragment-plane-wave Hamiltonian couplings on `core+buffer` support for fragmented axes and fragment-box PBC for non-fragmented axes.

**Step 2: Keep coupling definitions aligned**

Ensure `S_mat_frag_pw` and `H_mat_frag_pw` are built on the same support rule.

**Step 3: Verify with short run**

Run a short `DG+PW` test and confirm no shape/normalization mismatch appears in mixed-basis assembly.

**Step 4: Commit**

```bash
git add src/rt/dg/rt_dg_plane_wave.f90
git commit -m "refactor: align dg+pw couplings with axis-wise support"
```

### Task 6: Align propagation with mixed-basis-only runtime

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-phi-box-cache-phase-a/src/rt/dg/rt_dg_integrator_derivative.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-phi-box-cache-phase-a/src/rt/dg/rt_dg_integrator_stage_update.f90`

**Step 1: Audit raw-path usage**

Identify where runtime propagation still uses raw fragment/PW coupling matrices directly.

**Step 2: Restrict runtime to mixed representation**

Ensure mixed-basis-ready `DG+PW` uses synchronized mixed-state propagation consistently.

**Step 3: Build and smoke test**

Run short `DG+PW` propagation to verify startup, step update, and derivative paths remain coherent.

**Step 4: Commit**

```bash
git add src/rt/dg/rt_dg_integrator_derivative.f90 src/rt/dg/rt_dg_integrator_stage_update.f90
git commit -m "refactor: align dg+pw propagation with mixed basis"
```

### Task 7: Add regression probes for DG+PW representation consistency

**Files:**
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-phi-box-cache-phase-a/src/rt/dg/rt_dg_observables.f90`
- Modify: `/Users/otobetoshihito/Library/CloudStorage/OneDrive-qst.go.jp/SALMON-v.2.2.2/.worktrees/codex/dg-phi-box-cache-phase-a/src/rt/dg/rt_dg_density_reconstruct.f90`

**Step 1: Add concise diagnostics**

Emit probe lines showing:
- mixed basis active
- mixed basis dimension
- whether raw runtime path was bypassed

**Step 2: Verify probe output**

Run a small `DG+PW` case and confirm the intended path is visible in the log.

**Step 3: Commit**

```bash
git add src/rt/dg/rt_dg_observables.f90 src/rt/dg/rt_dg_density_reconstruct.f90
git commit -m "chore: add dg+pw mixed-basis consistency probes"
```

### Task 8: Physics comparison against Full

**Files:**
- No source changes required

**Step 1: Run matched short comparisons**

Use a short `DG+PW` vs `Full` comparison at the same laser conditions and compare:
- `A(t)`
- current response
- `Ne_raw`

**Step 2: Run directional comparison**

Compare `x/y/z` polarization to ensure the representation cleanup did not introduce obvious directional artifacts.

**Step 3: Record results**

Summarize whether `DG+PW` is now internally consistent enough to become the new practical path for metallic/solid-like use.

**Step 4: Commit notes if needed**

```bash
git add docs/plans/2026-04-09-dg-pw-orthogonalized-mixed-basis-design.md docs/plans/2026-04-09-dg-pw-orthogonalized-mixed-basis-plan.md
git commit -m "docs: add dg+pw orthogonalized mixed-basis design and plan"
```
