# Occupied-W Localization Diagnostic Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Measure the runtime projected occupied-W widths with the same discrete Gamma-point spread convention used by the existing Wannier90 path, without changing the occupied gauge or production gates.

**Architecture:** Each fragment computes link matrices for its owned 16 W rows from one unique canonical-cell preimage contained in P. Existing reciprocal-neighbor weights feed a standalone spread kernel. MPI communicates only small matrices and diagnostics; localization is deferred until the measured distribution is reviewed.

**Tech Stack:** Fortran 2008, MPI, existing Wannier90 reciprocal-shell helpers, DG Python/Fortran test drivers, CMake.

---

### Task 1: Exact discrete spread kernel

**Files:**
- Modify: `src/gs/dc/lcfo_wannier_sawf_seed.f90`
- Modify: `tests/dg/test_dg_wpw_core_wannier_seed_mpi.f90`
- Modify: `tests/dg/run_dg_wpw_core_wannier_seed_mpi.py`

**Step 1: Write failing tests**

Add analytic diagonal-link fixtures for localized W functions, including
positive/negative centers and a center across the phase seam. Test
`q=Im log(M/N)`, `rbar=-sum_b w*b*q`, and
`Omega=sum_b w*(1-|M/N|^2+(q+b.rbar)^2)`, plus units and invariance under a
common input scaling after normalization. Require `abs(N-1)<=1e-8` for
production-valid W while allowing the pure formula fixture to demonstrate
scaling invariance. Add
inconsistent +/-b conjugate pairs, invalid shapes or reciprocal weights,
non-finite input, and undefined-center cases. A fixed-b link is not required
to be Hermitian.

**Step 2: Verify RED**

Run: `python3 tests/dg/run_dg_wpw_core_wannier_seed_mpi.py`

Expected: compile failure because the spread API is absent.

**Step 3: Implement the minimal kernel**

Expose a routine accepting diagonal `M(n,b)`, reciprocal vectors, and weights.
For the orthorhombic Gamma shell construct `b=+/- (2*pi/L_a)e_a` and
`w=1/(2*|b_a|^2)`, then verify
`sum_b w*b_alpha*b_beta=delta_alpha_beta`. When an `nnkp` file is available,
verify its G vectors are the same +/- Cartesian shell. Implement total per-W
spread only and return atomic-unit quantities plus an explicit center-validity
mask. Do not implement an incomplete 16-W `Omega_I/Omega_OD` decomposition.

**Step 4: Verify GREEN**

Run the focused core-W seed driver. Expected: PASS with analytic values and rejection
cases.

**Step 5: Commit**

```bash
git add src/gs/dc/lcfo_wannier_sawf_seed.f90 tests/dg/test_dg_wpw_core_wannier_seed_mpi.f90 tests/dg/run_dg_wpw_core_wannier_seed_mpi.py
git commit -m "feat: compute discrete occupied Wannier spreads"
```

### Task 2: Unique canonical-cell extraction from P

**Files:**
- Modify: `src/gs/dc/dg_wpw_occupied_w_basis.f90`
- Modify: `tests/dg/test_dg_wpw_occupied_w_basis.f90`
- Modify: `tests/dg/run_dg_wpw_occupied_w_basis.py`

**Step 1: Write failing tests**

For core 12 and total cell 24, test B=6 and B=10 cyclic/unwrapped P arrays.
Require one deterministic P index for every canonical grid point, identical
canonical values for both B, and no alias sums.  Add B<6 incomplete coverage,
inconsistent duplicate images, invalid origin/shape, and 3-D seam cases.

**Step 2: Verify RED**

Run: `python3 tests/dg/run_dg_wpw_occupied_w_basis.py`

Expected: compile failure because canonical-cell extraction is absent.

**Step 3: Implement minimal extraction**

Add a helper that selects a half-open, total-cell-sized window from P, maps it
bijectively to canonical coordinates, and optionally checks all additional P
images against the selected values.  Return failure on missing points,
duplicates in the selected window, or periodic inconsistency. Extra P
preimages outside the selected window are allowed, checked for equality, and
never summed. Use `100*epsilon(1d0)*max(1,maxabs)` as the duplicate-image
consistency tolerance.

**Step 4: Verify GREEN**

Run the occupied-W basis driver. Expected: PASS for B=6/B=10 and all seam
fixtures.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_wpw_occupied_w_basis.f90 tests/dg/test_dg_wpw_occupied_w_basis.f90 tests/dg/run_dg_wpw_occupied_w_basis.py
git commit -m "feat: extract one physical cell from occupied W buffer"
```

### Task 3: Fragment-local link-matrix assembly

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90:4432-4640`
- Modify: `tests/dg/test_dg_wpw_core_wannier_seed_mpi.f90`
- Modify: `tests/dg/run_dg_wpw_core_wannier_seed_mpi.py`

**Step 1: Write failing MPI tests**

Build a synthetic 16-W fragment block split over ranks. Require the 16
diagonal link sets
`M_n(b)=sum_r |W_n(r)|^2*exp(-i*b.r)*hvol`, identical to a serial reference,
invariant between B=6 and B=10, and collective failure for missing
canonical coverage or inconsistent images.

**Step 2: Verify RED**

Run: `python3 tests/dg/run_dg_wpw_core_wannier_seed_mpi.py`

Expected: failure because production link assembly is absent.

**Step 3: Implement rank-local assembly**

Only the fragment representative (`dc%id_frag==0`) extracts one canonical cell
from its replicated `buffer_values` and assembles the 16 diagonal link sets;
other ranks contribute no grid sum. Collect the resulting 16 per-W records per
fragment once on `dc%icomm_tot`. Keep stable global W IDs and never mix rows
across fragments. Synchronize validation before every collective. The MPI
fixture must explicitly provide replicated P on multiple ranks and prove the
result is not multiplied by the rank count. Serial/MPI agreement uses
`500*epsilon(1d0)*max(1,maxabs)`.

**Step 4: Verify GREEN**

Run the MPI driver with its normal rank counts. Expected: PASS and exact serial
agreement within the named floating-point tolerance.

**Step 5: Commit**

```bash
git add src/gs/dc/lcfo_flux.f90 tests/dg/test_dg_wpw_core_wannier_seed_mpi.f90 tests/dg/run_dg_wpw_core_wannier_seed_mpi.py
git commit -m "feat: assemble fragment occupied Wannier links"
```

### Task 4: Production diagnostics before the existing buffer gate

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90:4432-4640`
- Modify: `tests/dg/check_dg_wpw_fixed_h_relaxation.py`
- Modify: `tests/dg/check_dg_wpw_production_consumer.py`

**Step 1: Write failing source-contract tests**

Require per-W stable ID, fragment, center validity, `Omega_A2`, and `width_A`,
followed by global min/mean/median/p90/max and count
above 1.2 Angstrom.  Require the diagnostic before `buffer_sufficiency`, forbid
changing that gate, and forbid dense 128-by-128 W rotation/gather.

**Step 2: Verify RED**

Run both Python contract tests. Expected: FAIL on the missing diagnostic.

**Step 3: Add output and collective validation**

Convert Bohr squared to Angstrom squared only for output.  Preserve undefined
center records and stop collectively only after diagnostics are emitted.
Leave W values, gradients, density, gate tolerance, and publication state
unchanged.

**Step 4: Verify GREEN and regression tests**

Run:

```bash
python3 tests/dg/check_dg_wpw_fixed_h_relaxation.py
python3 tests/dg/check_dg_wpw_production_consumer.py
python3 tests/dg/run_dg_wpw_occupied_w_basis.py
python3 tests/dg/run_dg_wpw_occupied_w_basis_mpi.py
python3 tests/dg/run_dg_wpw_core_wannier_seed_mpi.py
```

Expected: PASS, with no changed bootstrap acceptance behavior.

**Step 5: Commit**

```bash
git add src/gs/dc/lcfo_flux.f90 tests/dg/check_dg_wpw_fixed_h_relaxation.py tests/dg/check_dg_wpw_production_consumer.py
git commit -m "feat: report occupied Wannier localization"
```

### Task 5: Fresh Si64 measurement and decision checkpoint

**Files:**
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`
- Modify: `docs/plans/2026-07-22-wpw-occupied-w-localization.md`

**Step 1: Verify before physical execution**

Run the focused suite, `cmake --build build-mpi-eigenexa -j4`, and
`git diff --check`. Expected: all pass.

**Step 2: Run fresh B=6 and B=10 diagnostics**

Use new directories and the validated state counts. Preserve all prior runs.
Both cases must report the same physical-cell spread distribution within
`1e-10 Angstrom` for centers modulo the cell and `1e-10 Angstrom^2` for spread;
otherwise stop at periodic extraction/link assembly investigation.

**Step 3: Review the measurement**

Compare the full exact distribution with the historical approximately 1.2
Angstrom value.  Do not use an automatic post-hoc threshold.  Record min,
median, p90, mean, max, widest W identity, and B=6/B=10 discrepancy.  If the W
are excessive, write and review the separate 16-by-16 equal-occupation
fragment-local MV localization plan specified by the design boundary.

**Step 4: Request findings-first code review**

Review formula equivalence with the existing Wannier90 path, reciprocal
weights, canonical P coverage, duplicate-image checks, MPI collective order,
units, and unchanged production behavior. Fix all Critical/Important findings
and repeat review.

**Step 5: Commit the checkpoint**

```bash
git add docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md docs/plans/2026-07-22-wpw-occupied-w-localization.md
git commit -m "docs: record occupied Wannier spread measurement"
```
