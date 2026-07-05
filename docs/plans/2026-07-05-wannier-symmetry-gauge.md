# Wannier Symmetry Gauge Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a general fragment-local symmetry-gauge preparation path for Wannier90-derived DG bases so each DG fragment keeps its local inversion symmetry with `WW` enabled.

**Architecture:** Implement a DC-side post-processing layer in `src/gs/dc/lcfo_flux.f90` that detects symmetry operations closing inside each owner fragment, diagnoses their action on the fragment-owned Wannier centers, and then rotates the imported Wannier data with a block-diagonal unitary in owner-fragment space. Start with fragment-local inversion as the enabled operation, but store operations through a general local affine `{R | tau}` representation. Do not use full-supercell space-group operations unless their restricted action closes inside the current DG fragment.

**Tech Stack:** Fortran, SALMON DC-LCFO Wannier import/export, Wannier90 checkpoint data, real-space grid wavefunction reconstruction, LAPACK/EigenExa eigensolver wrappers, existing DG binary seed files.

---

### Task 1: Add Fragment-Local Symmetry Operation Data Structures

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Step 1: Add local derived types near the Wannier constants**

Define:

```fortran
type :: t_wannier_symop
  integer :: owner_frag = 0
  integer :: rot(3,3) = 0
  real(8) :: origin_bohr(3) = 0.0d0
  real(8) :: tau_local(3) = 0.0d0
  real(8) :: atom_residual = 0.0d0
  character(32) :: label = ''
end type t_wannier_symop
```

**Step 2: Add helper routines**

Add pure helpers:

- `wrap_periodic_delta(delta, cell_length)`
- `local_distance2(delta)`
- `integer_det3(rot)`
- `matmul_int_real(rot, vec)`

**Step 3: Build**

Run:

```bash
cmake --build build-mpi-eigenexa-wannier-lib -j 2
```

Expected: build succeeds with no behavior change.

**Step 4: Commit**

```bash
git add src/gs/dc/lcfo_flux.f90
git commit -m "Add Wannier symmetry operation helpers"
```

### Task 2: Detect Fragment-Local Structure Symmetry Operations

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Step 1: Write the detector**

Add:

```fortran
subroutine detect_wannier_fragment_symops(dc, nsym, symops)
```

The first implementation must:

- collect atoms whose primary positions are inside each owner fragment core
- express positions relative to the fragment center or a candidate local origin
- always test identity for each fragment
- test inversion candidates around local origins generated from atom pairs inside the fragment
- accept an operation only when every atom owned by that fragment maps to an atom of the same species inside the same fragment
- reject operations that only close after mapping to another fragment, even if they are valid symmetries of the full crystal

Keep the internal representation general: `rot = -I` for inversion, `origin_bohr` and `tau_local` arbitrary.

**Step 2: Log diagnostics**

Print:

```text
[DC-LCFO-W90-SYM] fragment=... detected symop label=... atom_residual=...
```

**Step 3: Run a DC import-only smoke**

Use an existing Wannier90 output directory and run import-only. Expected:

- identity detected for each fragment
- local inversion detected for each C64 fragment that has a locally inversion-symmetric atom set
- no change yet to written basis

**Step 4: Commit**

```bash
git add src/gs/dc/lcfo_flux.f90
git commit -m "Detect fragment-local Wannier symmetries"
```

### Task 3: Diagnose Fragment-Owned Wannier Center Closure

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Step 1: Add center action routine**

Add:

```fortran
subroutine diagnose_fragment_wannier_center_symmetry(center_bohr, owner_frag, nsym, symops, &
    center_perm, center_shift, max_residual, rms_residual, closed)
```

For each operation and each center owned by `symop%owner_frag`, map `g r_i` in local coordinates, search the nearest center with the same owner, and record:

- permutation index
- residual distance
- whether the owner-fragment WF center set is closed under the operation

**Step 2: Log closure**

Expected log:

```text
[DC-LCFO-W90-SYM] fragment=... center closure label=inversion max=... rms=... closed=T/F
```

**Step 3: Verify on current C64 Wannier output**

Expected before gauge fixing:

- fragment-local inversion operation exists structurally
- owner-fragment center closure is imperfect for the current raw MLWF gauge

**Step 4: Commit**

```bash
git add src/gs/dc/lcfo_flux.f90
git commit -m "Diagnose fragment Wannier center symmetry closure"
```

### Task 4: Build Symmetry Representation In Each Fragment WF Subspace

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Step 1: Reconstruct real-space Wannier functions**

Use the existing DC wavefunction seed and `v_matrix` to reconstruct Wannier functions owned by one fragment on the fragment-local real-space grid. Only WFs in the same owner-fragment block enter this representation.

**Step 2: Apply a symmetry operation to each WF**

For each fragment-local grid point `r`, map to `g^{-1} r` inside the same fragment box and sample the transformed WF. Start with grid-point exact mappings for inversion; do not add interpolation until needed. If the mapped grid point leaves the fragment-local box, mark that operation as non-closing for the fragment.

**Step 3: Compute overlap matrix**

Compute:

```text
S_g(i,j) = < W_i^frag | g W_j^frag >
```

using the real-space volume element.

**Step 4: Unitarize**

Use polar decomposition through Hermitian diagonalization:

```text
D_g = S_g (S_g^dagger S_g)^(-1/2)
```

If the smallest eigenvalue is below tolerance, print a warning and do not rotate that fragment.

**Step 5: Commit**

```bash
git add src/gs/dc/lcfo_flux.f90
git commit -m "Build fragment Wannier symmetry representation"
```

### Task 5: Apply Fragment-Local Inversion Gauge Symmetrization

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`

**Step 1: Compute the target inversion representation**

From center closure inside one owner-fragment block, build the permutation representation `P_inv`.

**Step 2: Find a unitary correction**

Construct a unitary rotation `Q` that makes the measured representation closer to `P_inv`. Start with inversion-only projection:

```text
Q_f,new = polar( Q_f + D_f,inv Q_f P_f,inv^dagger )
```

Iterate until `||D_inv Q - Q P_inv||` is below tolerance or max iterations.

**Step 3: Rotate all Wannier data**

Assemble `Q = direct_sum_f Q_f` over owner fragments, then apply:

```text
V' = V Q
AA'_a = Q^dagger AA_a Q
center'_a(i) = real(AA'_a(i,i))
seed' = Q^dagger seed
```

Then recompute owner fragments from `center'`.

**Step 4: Log final residuals**

Print before/after:

- fragment-local center closure max/rms
- fragment-local representation residual
- `AA_R` Hermiticity residual

**Step 5: Commit**

```bash
git add src/gs/dc/lcfo_flux.f90
git commit -m "Symmetrize fragment Wannier gauge for inversion"
```

### Task 6: Wire A Namelist Switch

**Files:**
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`

**Step 1: Add setting**

Add:

```fortran
character(32) :: dg_wannier_fragment_symmetry_gauge
```

Accepted values:

- `none`
- `diagnose`
- `local_inversion`

Default: `diagnose` while developing, or `none` if preserving strict historical behavior is required.

**Step 2: Add variable log output**

Write the chosen mode to `variables.log`.

**Step 3: Validate**

Reject unknown strings in `inputoutput.f90`.

**Step 4: Commit**

```bash
git add src/io/salmon_global.f90 src/io/inputoutput.f90 src/gs/dc/lcfo_flux.f90
git commit -m "Add Wannier symmetry gauge control"
```

### Task 7: Regression Runs

**Files:**
- Use: current C64 Wannier/DC output
- Use: current C64 RT laser input

**Step 1: Rebuild**

```bash
cmake --build build-mpi-eigenexa-wannier-lib -j 2
```

Expected: build succeeds.

**Step 2: Regenerate symmetric DC Wannier binary files**

Run DC import-only with:

```fortran
dg_wannier_fragment_symmetry_gauge = 'local_inversion'
```

Expected:

- symmetry logs show fragment-local inversion detected
- per-fragment center and representation residuals improve
- `wannier90_global_basis.bin` and `wannier_flux_eigen_seed.bin` are rewritten

**Step 3: Run RT with WW enabled**

Use:

```fortran
yn_dg_mixed_z_include_ww = 'y'
```

Expected:

- `||Z_WW||_F` is nonzero
- Pz remains inversion-symmetric under +/-E
- H2/H1 and H4/H1 are near the previous `include_ww=n` control

**Step 4: Commit final verification notes if scripts or docs are added**

```bash
git add docs/plans/2026-07-05-wannier-symmetry-gauge*.md
git commit -m "Document Wannier symmetry gauge plan"
```
