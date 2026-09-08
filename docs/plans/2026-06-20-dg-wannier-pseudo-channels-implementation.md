# DG-Wannier Pseudo-Channel Projection Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Generalize DC-LCFO Wannier export from a diamond-only `C:sp3` path to a material-independent `pseudo_channels` projection path that supports `num_wann < num_proj` and enables Wannier-count convergence studies.

**Architecture:** Keep the current `C:sp3` export as a regression path.  Add a new `wannier_projection='pseudo_channels'` path that generates AO candidate projectors, compresses them to `num_wann` via projectability/SVD, writes standard Wannier90 seed files, and reuses the existing DG-Wannier+BPW diagnostics CSV workflow.

**Tech Stack:** Fortran DC-LCFO export (`src/gs/dc/lcfo_flux.f90`), SALMON namelist input, Wannier90 `.win/.eig/.amn/.mmn/.chk`, existing DG mixed-Z BPW diagnostics, MPI 8 / OMP 2 validation runs.

---

### Task 1: Add Projection Mode Detection

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/io/inputoutput.f90` only if validation messages need adjustment

**Step 1: Add a helper for pseudo-channel mode**

Add a helper near `is_c_sp3_projection`:

```fortran
logical function is_pseudo_channel_projection(text) result(enabled)
  character(*), intent(in) :: text
  character(len(text)) :: work
  work = adjustl(text)
  enabled = (trim(work) == 'pseudo_channels' .or. trim(work) == 'PSEUDO_CHANNELS')
end function is_pseudo_channel_projection
```

**Step 2: Keep existing C:sp3 validation**

Change the hard stop:

```fortran
if(.not. is_c_sp3_projection(trim(wannier_projection))) &
  stop "DC-LCFO Wannier export: only wannier_projection='C:sp3' is implemented."
```

to:

```fortran
if(.not. is_c_sp3_projection(trim(wannier_projection)) .and. &
   .not. is_pseudo_channel_projection(trim(wannier_projection))) &
  stop "DC-LCFO Wannier export: supported projections are C:sp3 and pseudo_channels."
```

**Step 3: Build**

Run:

```bash
cmake --build build-mpi-eigenexa-wannier-lib --target salmon -j 4
```

Expected: build succeeds.

**Step 4: Commit**

```bash
git add src/gs/dc/lcfo_flux.f90
git commit -m "feat: accept pseudo-channel Wannier projection mode"
```

### Task 2: Generate AO Candidate Metadata

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Step 1: Add AO candidate arrays**

Create a small metadata representation:

```fortran
integer, allocatable :: proj_atom(:)
integer, allocatable :: proj_l(:)
integer, allocatable :: proj_m(:)
integer :: nproj_raw
```

Use `l=0` for `s`, `l=1` for `p`, and `l=2` for `d`.

**Step 2: Add pseudo-channel candidate builder**

Initial implementation:

```fortran
subroutine build_pseudo_channel_projection_map(proj_atom, proj_l, proj_m, nproj)
```

For each atom:

```text
s: 1
p: 3
d: 5
```

For the first validation, generate `s+p+d` for every atom.  This avoids relying on incomplete channel metadata and is safe for C, Si, O, BN/hBN first tests.

**Step 3: Log projection count**

Print:

```text
[DC-LCFO-WANNIER] pseudo_channels candidates=<nproj_raw> target_wann=<wannier_num_wann>
```

**Step 4: Build**

Run:

```bash
cmake --build build-mpi-eigenexa-wannier-lib --target salmon -j 4
```

Expected: build succeeds.

**Step 5: Commit**

```bash
git add src/gs/dc/lcfo_flux.f90
git commit -m "feat: generate pseudo-channel AO candidates"
```

### Task 3: Implement Real AO Values

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Step 1: Add real AO evaluator**

Add:

```fortran
real(8) function pseudo_channel_projection_value(x, y, z, ia, l, m, sigma) result(val)
```

Use the existing minimum-image displacement convention from `c_sp3_projection_value`.

Basis:

```text
s        exp(-r2/2s2)
px,py,pz x,y,z times Gaussian
dxy      x*y times Gaussian
dyz      y*z times Gaussian
dzx      z*x times Gaussian
dx2-y2   (x*x-y*y) times Gaussian
dz2      (2*z*z-x*x-y*y) times Gaussian
```

**Step 2: Normalize by real-space norm**

Reuse the existing `norm_local/norm_sum` normalization pattern in `.amn` generation.

**Step 3: Build**

Run:

```bash
cmake --build build-mpi-eigenexa-wannier-lib --target salmon -j 4
```

Expected: build succeeds.

**Step 4: Commit**

```bash
git add src/gs/dc/lcfo_flux.f90
git commit -m "feat: evaluate pseudo-channel AO projections"
```

### Task 4: Allow `num_wann < num_proj`

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Step 1: Split raw projection count and exported projection count**

Keep:

```fortran
nproj_raw = number of AO candidates
nproj = wannier_num_wann
```

Only the compressed AMN columns are written to Wannier90.

**Step 2: Compute raw AMN**

Compute `amn_raw(nband_wann,nproj_raw)` in chunks when possible.

For a first implementation, if memory is acceptable for diamond:

```fortran
allocate(amn_raw(nband_wann,nproj_raw))
```

C64 with 320 bands and 576 raw `s+p+d` projectors is small enough for this dense diagnostic path.

**Step 3: Compress raw projectors**

Use projectability/Gram selection:

```fortran
q(ip) = sum(amn_raw(:,ip)**2)
```

Keep the top `wannier_num_wann` raw projectors for the first implementation.

This is the minimal version.  SVD compression can replace it in a later task.

**Step 4: Preserve C:sp3 exact behavior**

For `C:sp3`, keep the current cyclic-shell behavior and current stop when `nproj_csp3 > wannier_num_wann`.  Only `pseudo_channels` must allow `num_wann < num_proj`.

**Step 5: Build**

Run:

```bash
cmake --build build-mpi-eigenexa-wannier-lib --target salmon -j 4
```

Expected: build succeeds.

**Step 6: Commit**

```bash
git add src/gs/dc/lcfo_flux.f90
git commit -m "feat: compress pseudo-channel projections to num_wann"
```

### Task 5: Add Input Templates For Diamond Validation

**Files:**
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_gs_w90_pseudo_occ_only`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_gs_w90_pseudo_occ_lowc05`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_gs_w90_pseudo_occ_lowc10`

**Step 1: Copy the existing global input**

Use `inputfile_gs_w90_global` as the base.

**Step 2: Change projection and state counts**

For C64, `Nocc=128`.

```text
occ_only:
  nstate = 128
  wannier_num_bands = 128
  wannier_num_wann = 128

occ_lowc05:
  nstate = 192
  wannier_num_bands = 192
  wannier_num_wann = 192

occ_lowc10:
  nstate = 256
  wannier_num_bands = 256
  wannier_num_wann = 256
```

Set:

```text
wannier_projection = 'pseudo_channels'
num_wannier_cluster(1:3) = 1,1,1
```

**Step 3: Commit**

```bash
git add samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_gs_w90_pseudo_*
git commit -m "test: add diamond pseudo-channel Wannier inputs"
```

### Task 6: Run Diamond GS/Wannier Smoke Cases

**Files:**
- No source modifications

**Step 1: Run occ_only**

Run from a copied scratch directory or the sample directory:

```bash
OMP_NUM_THREADS=2 mpirun -np 8 /path/to/salmon < inputfile_gs_w90_pseudo_occ_only
```

Expected:

```text
[DC-LCFO-WANNIER] export bands=128 wann=128
[DC-LCFO-WANNIER] pseudo_channels candidates=...
[DC-LCFO-W90-GLOBAL] wrote 128 Wannier functions
```

**Step 2: Run occ_lowc05 and occ_lowc10**

Use the same command with the corresponding input files.

Expected: both finish and write `wannier90_global_basis.bin`, AA_R, and flux seed files.

**Step 3: Commit only if input fixes are needed**

Do not commit generated binary data.

### Task 7: Produce Diagnostics CSVs

**Files:**
- Use: `tools/dg_bpw_ecut_sweep.py`
- Use: `tools/plot_dg_bpw_sweep.py`

**Step 1: Run BPW Ecut sweep for each Wannier case**

Example:

```bash
python3 tools/dg_bpw_ecut_sweep.py \
  --allow-diagnostic-stop \
  --salmon build-mpi-eigenexa-wannier-lib/salmon \
  --input samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_rt_w90_mixed_fsum_diag_pw128_k2 \
  --workdir /path/to/case \
  --output /tmp/dg_wannier_pseudo/occ_only.csv \
  --log-dir /tmp/dg_wannier_pseudo/logs_occ_only \
  --case-label occ_only \
  --np 8 \
  --omp 2 \
  --sperp-tol 1e-1 \
  --ecut 0.5 0.75 1.0 1.5 2.0 3.0
```

**Step 2: Plot comparison**

Run:

```bash
python3 tools/plot_dg_bpw_sweep.py \
  /tmp/dg_wannier_pseudo/occ_only.csv \
  /tmp/dg_wannier_pseudo/occ_lowc05.csv \
  /tmp/dg_wannier_pseudo/occ_lowc10.csv \
  --out-dir /tmp/dg_wannier_pseudo/plots \
  --metrics np fsum_avg c_comm_sum herm selected_g2_max min_sperp max_sperp
```

Expected: SVG plots are generated.

### Task 8: Select The First Candidate

**Files:**
- No source modifications

**Step 1: Compare CSVs**

Selection criteria:

```text
occ_only -> lowc05 improves fsum_avg / C_comm_sum / herm
lowc05 -> lowc10 is small
min_sperp is not pathologically small
max_sperp is stable
Z Hermiticity remains <= 1e-12
```

**Step 2: Record the decision**

Append a short note to:

```text
docs/plans/dg_wannier_pseudo_channels.md
```

Include the selected case and reason.

**Step 3: Commit**

```bash
git add docs/plans/dg_wannier_pseudo_channels.md
git commit -m "docs: record diamond pseudo-channel Wannier selection"
```
