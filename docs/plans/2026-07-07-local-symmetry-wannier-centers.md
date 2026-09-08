# Local Symmetry Wannier Centers Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add diagnostics and then a gated generation path so local bond-center Wannier seeds preserve fragment-local structural symmetry.

**Architecture:** Start with read-only diagnostics in the DC local Wannier export path, because current evidence shows the asymmetry already exists in `local_wannier_basis.bin` and `buffer_periodic_wannier_basis.bin`. Once diagnostics locate the broken transformation, add an experimental symmetry-gauge branch behind a namelist option.

**Tech Stack:** Fortran 90 in `src/gs/dc/lcfo_flux.f90`, SALMON namelist plumbing in `src/io/inputoutput.f90` and `src/io/salmon_global.f90`, local MPI/root logging, existing binary export format.

---

### Task 1: Add a Static Diagnostic Test

**Files:**
- Test: one-off static command from repository root
- Inspect: `src/gs/dc/lcfo_flux.f90`

**Step 1: Write the failing test**

Run:

```bash
python3 - <<'PY'
from pathlib import Path
text = Path('src/gs/dc/lcfo_flux.f90').read_text()
assert 'diagnose_local_bond_center_orbit_closure' in text
assert 'diagnose_local_wannier_center_orbit_closure' in text
assert '[DC-LCFO-LOCAL-WANNIER-SYM]' in text
PY
```

Expected: FAIL because the diagnostic routines do not exist yet.

**Step 2: Implement minimal diagnostics**

Add two routines near the existing local Wannier export helpers:

- `diagnose_local_bond_center_orbit_closure(dc, ncenter, center_bohr)`
- `diagnose_local_wannier_center_orbit_closure(dc, ncenter, center_bohr, label)`

The first reports closure of the seed bond-center set under the fragment-local inversion candidates. The second reports closure of exported `wcenter_legacy` and BPW `wcenter`.

**Step 3: Run the static test**

Run the same Python command. Expected: PASS.

**Step 4: Commit**

```bash
git add src/gs/dc/lcfo_flux.f90
git commit -m "Diagnose local Wannier center symmetry"
```

### Task 2: Wire Diagnostics Into Local Wannier Export

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Step 1: Write the failing test**

Run:

```bash
python3 - <<'PY'
from pathlib import Path
text = Path('src/gs/dc/lcfo_flux.f90').read_text()
assert 'call diagnose_local_bond_center_orbit_closure' in text
assert 'call diagnose_local_wannier_center_orbit_closure(dc, nkeep_legacy, wcenter_legacy' in text
assert 'call diagnose_local_wannier_center_orbit_closure(dc, nkeep, wcenter' in text
PY
```

Expected: FAIL until calls are added.

**Step 2: Add calls**

In `write_local_wannier_seed`, call:

- seed closure diagnostic immediately after `build_local_bond_center_projection_map`
- local legacy center closure after `wcenter_legacy` is computed
- BPW center closure after `wcenter` is computed

Only run these calls on `dc%id_frag == 0`.

**Step 3: Build**

Run:

```bash
cmake --build build-mpi-eigenexa-wannier-lib -j 2
```

Expected: build succeeds.

**Step 4: Commit**

```bash
git add src/gs/dc/lcfo_flux.f90
git commit -m "Report local Wannier symmetry closure"
```

### Task 3: Regenerate a Short DC Export and Confirm Broken Layer

**Files:**
- Use: `/private/tmp/c64_bond_local_converged/inputfile_gs_bond_local_converged`
- Output: `/private/tmp/c64_bond_local_symdiag`

**Step 1: Run DC export with existing settings**

Use the current build and the existing C64 bond-center input. Keep `num_rgrid_buffer=6` and `wannier_projection='bond_centers'`.

**Step 2: Inspect logs**

Search:

```bash
rg -n "DC-LCFO-LOCAL-WANNIER-SYM|DC-LCFO-WANNIER|end SALMON|FATAL" run*.log
```

Expected:

- seed bond-center closure should be good for diamond fragments
- final `wcenter_legacy` and BPW `wcenter` closure should fail, matching the current binary-file diagnosis

### Task 4: Add Gated Symmetry-Gauge Option

**Files:**
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`

**Step 1: Write namelist static test**

Run:

```bash
python3 - <<'PY'
from pathlib import Path
all_text = '\\n'.join(Path(p).read_text() for p in [
  'src/io/salmon_global.f90',
  'src/io/inputoutput.f90',
  'src/gs/dc/lcfo_flux.f90',
])
assert 'yn_dc_lcfo_wannier_symmetry_gauge' in all_text
assert \"yn_dc_lcfo_wannier_symmetry_gauge = 'n'\" in all_text
PY
```

Expected: FAIL until namelist plumbing exists.

**Step 2: Add option**

Add `yn_dc_lcfo_wannier_symmetry_gauge` under `&dc`, default `'n'`.

**Step 3: Build**

Run:

```bash
cmake --build build-mpi-eigenexa-wannier-lib -j 2
```

Expected: build succeeds.

### Task 5: Implement Experimental Symmetry Gauge

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Step 1: Add failing static test**

Check for a helper named `apply_local_bond_center_symmetry_gauge`.

**Step 2: Implement minimal branch**

When `wannier_projection='bond_centers'` and `yn_dc_lcfo_wannier_symmetry_gauge='y'`:

- Work in the selected `wcoef_legacy` subspace.
- Build the projected position matrices.
- Diagonalize a small weighted position operator whose target centers are the seed bond centers.
- Rotate only within the selected local subspace.
- Recompute `wcenter_legacy`, `r_wann_legacy`, and BPW projected quantities after the rotation.

**Step 3: Verify**

Run the C64 DC export. Expected: final local center closure improves to below `1e-3` bohr for the local center set.

### Task 6: Compare RT HHG

**Files:**
- Use generated DC data from Task 5
- Use RT input based on `inputfile_rt_bond_local_laser_dt5_nt240_fullheigen_perwcenter`

**Step 1: Run init-only**

Expected:

- Full-H seed residual remains around `1e-12`
- `DG-FULL-H-SEED-XI` uses the per-Wannier center path

**Step 2: Run laser**

Use `dt=5, nt=240` first.

**Step 3: Analyze HHG**

Expected success criterion:

- H2 and H4 are dips relative to neighboring odd harmonics.
- Time-domain `Pz(+E)` and `Pz(-E)` satisfy odd response within numerical tolerance.

