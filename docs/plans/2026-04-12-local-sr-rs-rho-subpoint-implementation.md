# Local-SR RS Rho-Subpoint Probe Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a diagnostic-only rho-subpoint mode to the local-SR RS subdivision probe and run a convergence comparison for FHI, PSP, and VPS representative cases.

**Architecture:** Keep the existing production RS stress path unchanged. Extend the current diagnostic probe with a new rho sampling mode that can either reuse the cell-center density or interpolate density at each subpoint from an overlap-aware scalar box. Compare modes across separate runs.

**Tech Stack:** Fortran 90, SALMON stress diagnostics, MPI overlap exchange via `srg_scalar`, Wisteria batch verification.

---

### Task 1: Add the failing source checks

**Files:**
- Modify: `/Users/mizukitani/Documents/DFT/SALMON2-local-notes/feature-stress-tensor/tests/source/check_local_sr_rs_sampled_dump.sh`
- Test: same script

**Step 1: Write the failing test**

Add grep checks for:

- `mode_loc_sr_rs_subdiv_probe_rho`
- `mode_loc_sr_rs_subdiv_probe_rho = 'center'`
- `call fail_stress_input("mode_loc_sr_rs_subdiv_probe_rho must be`
- `probe_rho_mode`
- `trilinear`

**Step 2: Run test to verify it fails**

Run: `bash /Users/mizukitani/Documents/DFT/SALMON2-local-notes/feature-stress-tensor/tests/source/check_local_sr_rs_sampled_dump.sh`

Expected: FAIL because the new input and mode handling do not exist yet.

**Step 3: Write minimal implementation**

No implementation in this task.

**Step 4: Run test to verify it still fails for the expected reason**

Run the same command and confirm the missing grep target is the cause.

**Step 5: Commit**

Do not commit yet.

### Task 2: Add the rho probe input plumbing

**Files:**
- Modify: `/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/io/salmon_global.f90`
- Modify: `/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/io/inputoutput.f90`

**Step 1: Write the failing test**

Already done in Task 1.

**Step 2: Run test to verify it fails**

Use the same script from Task 1.

**Step 3: Write minimal implementation**

Add:

- `character(16) :: mode_loc_sr_rs_subdiv_probe_rho = 'center'`

Plumb it through:

- `/control/` namelist
- broadcast
- `variables.log`
- validation accepting only `center` and `trilinear`

**Step 4: Run test to verify partial progress**

Run: `bash /Users/mizukitani/Documents/DFT/SALMON2-local-notes/feature-stress-tensor/tests/source/check_local_sr_rs_sampled_dump.sh`

Expected: may still fail until stress-side usage is added.

**Step 5: Commit**

Do not commit yet.

### Task 3: Add diagnostic-only trilinear rho sampling

**Files:**
- Modify: `/Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor/src/common/stress.f90`

**Step 1: Write the failing test**

Reuse the source grep script, now expecting:

- stress-side mode variable use
- a trilinear helper
- dump header line for rho mode

**Step 2: Run test to verify it fails**

Run: `bash /Users/mizukitani/Documents/DFT/SALMON2-local-notes/feature-stress-tensor/tests/source/check_local_sr_rs_sampled_dump.sh`

Expected: FAIL before the helper and header are added.

**Step 3: Write minimal implementation**

In `calc_stress_loc_sr_rs`:

- accept optional `srg_scalar`
- build temporary total-density overlap box only when probe mode is `trilinear`
- use `update_overlap_real8(srg_scalar, mg, density_box)` when `info%if_divide_rspace`
- add helper to interpolate density at each subpoint
- keep `center` mode behavior unchanged
- write `probe_rho_mode` to the sampled dump header

**Step 4: Run test to verify it passes**

Run:

- `bash /Users/mizukitani/Documents/DFT/SALMON2-local-notes/feature-stress-tensor/tests/source/check_local_sr_rs_sampled_dump.sh`
- `bash /Users/mizukitani/Documents/DFT/SALMON2-local-notes/feature-stress-tensor/tests/source/check_local_sr_grad_sampled_dump.sh`
- `bash /Users/mizukitani/Documents/DFT/SALMON2-local-notes/feature-stress-tensor/tests/source/check_local_sr_stress_spline_backend.sh`
- `git -C /Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor diff --check`

Expected: PASS.

**Step 5: Commit**

```bash
git -C /Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor add \
  docs/plans/2026-04-12-local-sr-rs-rho-subpoint-design.md \
  docs/plans/2026-04-12-local-sr-rs-rho-subpoint-implementation.md \
  src/common/stress.f90 src/io/inputoutput.f90 src/io/salmon_global.f90
git -C /Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor commit -m "Add rho-subpoint mode to RS subdivision probe"
```

### Task 4: Push and rebuild on Wisteria

**Files:**
- Modify: none

**Step 1: Write the failing test**

Not applicable.

**Step 2: Run test to verify it fails**

Not applicable.

**Step 3: Write minimal implementation**

Push and rebuild:

```bash
git -C /Users/mizukitani/Documents/DFT/SALMON2/.worktrees/stress-tensor push origin feature/stress-tensor
```

On Wisteria:

```bash
git -C /work/go33/o33000/src/stress/SALMON2 pull --ff-only origin feature/stress-tensor
module load odyssey
cmake --build /work/go33/o33000/src/stress/SALMON2/build -j 48
cmake --install /work/go33/o33000/src/stress/SALMON2/build
```

**Step 4: Run verification**

Check installed binary mtime with:

```bash
stat -c '%y %n' /work/go33/o33000/src/stress/SALMON2/bin/salmon
```

**Step 5: Commit**

No new commit.

### Task 5: Run representative convergence campaigns

**Files:**
- Create under: `/work/go33/o33000/salmon/virial/pSi/test/`

**Step 1: Write the failing test**

Not applicable.

**Step 2: Run test to verify it fails**

Not applicable.

**Step 3: Write minimal implementation**

Prepare campaigns for:

- `si_pz_fhi_lda`
- `si_pz_psp_lda_str`
- `si_pz_vps_ca19`

For each case, run:

- `rho_mode=center`, `n=1,2,4`
- `rho_mode=trilinear`, `n=1,2,4`

Use `yn_out_loc_sr_rs_subdiv_probe='y'`.

**Step 4: Run verification**

Verify:

- `error.log` is empty
- `variables.log` contains the new mode
- sampled dump header contains `probe_rho_mode`
- `info.data` contains the probe pressure line

**Step 5: Commit**

No new commit.

### Task 6: Postprocess convergence tables

**Files:**
- Create under campaign root

**Step 1: Write the failing test**

Not applicable.

**Step 2: Run test to verify it fails**

Not applicable.

**Step 3: Write minimal implementation**

Generate tables containing:

- case
- rho mode
- `n`
- `P_loc_sr_rs`
- `P_loc_sr_grad`
- `P_loc_sr_rs_subdiv_probe`
- `probe-current`
- `n=2-n=1`
- `n=4-n=2`
- `r50/r90` for `abs(probe-current)`

Also generate a direct `center` vs `trilinear` compare for each `(case, n)`.

**Step 4: Run verification**

Check tables exist and values are internally consistent with `info.data`.

**Step 5: Commit**

No new commit.
