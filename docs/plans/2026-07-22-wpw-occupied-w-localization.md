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

## 2026-07-22 execution checkpoint

Tasks 1-4 passed focused serial/MPI tests, the full MPI/EigenExa build, and
findings-first review with no remaining Critical or Important findings.

Fresh B=6 (`20260722_task7_spread_b6/run.log`) converged DC-SCF at iteration 87
with `diff=9.8498E-10` and charge 256. All 128 centers were valid. The exact
per-W discrete widths were:

- minimum `1.26877 Angstrom`;
- mean `1.28540 Angstrom`;
- median/p90/maximum `1.29094 Angstrom`;
- 128 rows above the historical approximate `1.2 Angstrom` reference.

This is close to, rather than qualitatively broader than, the historical
Wannier90 result. The unchanged outer-shell gate still failed at ratio
`7.2686E-03`.

Fresh B=10 (`20260722_task7_spread_b10/run.log`) converged DC-SCF at iteration
85 with `diff=9.4464E-10` and charge 256, then failed collectively at
`occupied_w_link_assembly`. P contains multiple preimages of canonical points
when its extent is 32 and the physical cell is 24, so this failure is
consistent with a duplicate-image inconsistency, but the current log does not
yet isolate extraction from numerical link assembly or report the failing
coordinate and mismatch. No image was summed or silently selected, and no
B=10 spread was reported.

Task 5 therefore stops at its specified periodic extraction/link-assembly
investigation boundary. The evidence does not justify implementing
localization: B=6 W widths are already close to 1.2 Angstrom. The next root
cause investigation must instrument the extraction/link boundary and determine
whether the fragment projected-W P field violates the assumed physical-cell
image identity for B>6.

### Instrumented B=10 result

The extraction and link stages were separated and the extraction failure now
reports its largest duplicate-image mismatch. Fresh run
`20260722_task8_image_diag_b10/run.log` again converged at iteration 85
(`diff=9.4464E-10`, charge 256). Every fragment reported a mismatch between P
coordinates `[1,16,19]` and `[25,16,19]`, which differ by exactly the physical
24-point period in x and map to the same canonical point. The maximum absolute
W-value difference was `3.92919E-02` on every fragment. The affected local W
row and canonical coordinate transform with the fragment as expected.

The run then stopped collectively at the newly isolated
`occupied_w_canonical_cell_extraction` stage. This directly proves that the raw
B=10 fragment projected-W P values are not periodic under the 24-point
physical-cell translation. The discrepancy is many orders of magnitude above
the `100*epsilon*scale` image-consistency tolerance and is not roundoff.

Static root-cause tracing reaches `src/gs/dc/dcdft.f90:init_fragment`. The
fragment calculation replaces its cell/grid by `domain+2B`; hence B=10
fragment orbitals obey a 32-point fragment-cell periodicity. In contrast,
`dc%jxyz_tot` folds those 32 storage points onto the 24-point physical cell.
Storage indices 15 and 23 can therefore map to the same physical grid point
without being equal fragment-cell points. The mismatch exists at the raw
fragment-wavefunction/physical-grid mapping boundary, before projected-W
transformation or P reordering.

The corrective design must construct the WPW physical-periodic P from a unique
representative for each `dc%jxyz_tot` value and populate every P image from
that representative before differentiation. It must not require the raw
fragment SCF orbital, whose boundary condition belongs to the larger fragment
box, to be 24-point periodic.
