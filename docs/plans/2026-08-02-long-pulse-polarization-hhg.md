# Long-Pulse Polarization HHG Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Produce and validate reproducible long-pulse, large-time-step HHG spectra derived strictly from overlapping-Wannier polarization.

**Architecture:** Keep the accepted V3 generalized-eigenvalue exponential coefficient RT unchanged and drive it in length gauge with a ten-cycle pulse.  Separate qualitative finite-cell morphology gates from quantitative diagnostics, verify `dt=2.0 a.u.` against `dt=1.0 a.u.`, and generate publication-ready plots from the polarization spectra.

**Tech Stack:** Fortran 2008, MPI, Python 3, NumPy, Matplotlib, existing SALMON spectral analyzer.

---

### Task 1: Freeze the long-pulse and finite-cell acceptance RED

**Files:**
- Modify: `tests/dg/test_si64_harmonic_morphology.py`
- Modify: `tests/dg/check_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/run_si64_overlapping_wannier_response_hhg.py`

**Step 1: Write RED assertions**

Require manifests to report a ten-cycle pulse, at least two post-pulse cycles,
`dt=2.0 a.u.` for production, and `dt=1.0 a.u.` for the convergence case.
Require spectra to use polarization.  Keep H2/H3/H4 local morphology and axis
covariance as gates; report absolute powers and strong/weak ratios without
thresholding them.

**Step 2: Run RED**

Run:

```bash
python3 tests/dg/test_si64_harmonic_morphology.py
python3 tests/dg/check_si64_overlapping_wannier_response_hhg.py \
  /tmp/si64-exact-response-projector-fix \
  --displaced-root /tmp/si64-exact-response-displaced
```

Expected: FAIL because the existing matrices use the short pulse and old time steps.

**Step 3: Implement the minimal runner/checker policy**

Set the physical duration once in the runner and derive `nt` from each `dt`.
Store pulse cycles, post-pulse cycles, duration, and `dt` in every laser
manifest.  Remove only finite-cell quantitative intensity gates; do not weaken
provenance, symmetry, morphology, or time-step gates.

**Step 4: Focused verification and commit**

Run the morphology test, Python compilation, route/removal contracts, and
`git diff --check`.  Perform specification and code-quality reviews, resolve
all Critical and Important findings, then commit.

### Task 2: Configure the ten-cycle pulse and large time step

**Files:**
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_laser_weak.in`
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_laser_hhg.in`
- Modify: `tests/dg/run_si64_overlapping_wannier_response_hhg.py`

**Step 1: Verify the parameter RED from Task 1**

Confirm it fails specifically on pulse duration or time-step metadata.

**Step 2: Implement parameters**

For `omega=0.05696 a.u.`, use ten carrier periods for `tw1` and a total record
of twelve periods.  Run the main laser cases at `dt=2.0 a.u.` and the x
convergence case at `dt=1.0 a.u.` over the identical physical duration.

**Step 3: Focused verification and commit**

Run Python fixture checks and verify that rendered inputs have identical end
times within one fine step.  Review, resolve Critical/Important findings, and commit.

### Task 3: Generate reproducible polarization-HHG figures

**Files:**
- Add: `tools/plot_overlapping_wannier_hhg.py`
- Add: `tests/dg/test_plot_overlapping_wannier_hhg.py`

**Step 1: Write plotting RED**

Create synthetic spectra with H2 dip, H3 peak, and H4 slope.  Require PNG and
PDF output, log scaling, integer-harmonic guides, source labels, and a JSON
sidecar containing input hashes and plotted limits.

**Step 2: Run RED**

```bash
python3 tests/dg/test_plot_overlapping_wannier_hhg.py
```

Expected: FAIL because the plotting tool is absent.

**Step 3: Implement the plotter**

Plot ideal x/y/z in the first panel.  If a displaced root is supplied, compare
ideal and fixed-displaced x in a second panel.  Read only `hhg-spectrum.tsv`
files whose summaries declare `spectrum_source=polarization`.

**Step 4: Focused verification and commit**

Run the plotting test, render a fixture figure, inspect it visually, run
`git diff --check`, complete both reviews, resolve findings, and commit.

### Task 4: Run ideal and fixed-displaced production evidence

**Files:**
- Modify: `docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md`
- Modify: `docs/plans/2026-07-31-obsolete-dg-route-removal-results.md`

**Step 1: Run a clean-first parent-prerequisite overlay build**

Create a fresh `git archive HEAD` tree, apply the parent prerequisite and
current diff overlays, configure MPI, ScaLAPACK, EigenExa, and spglib, and run
`cmake --build <build> --clean-first -j1`.  Record options and binary SHA-256.

**Step 2: Run focused tests**

Run overlapping-Wannier nonlocal, observable, symmetry, construction, solver,
SCF, checkpoint, and RT runners on 1/2/4/8 ranks plus route/removal contracts,
normal DC LCFO+EigenExa, Python tests, and `git diff --check`.

**Step 3: Run genuine production matrices**

Generate fresh ideal and fixed-displaced GS checkpoints with the accepted
binary.  Run field-off, impulse, long weak laser, long HHG x/y/z, amplitude
comparison, and `dt=2.0/1.0` convergence cases.  Reject forbidden route markers.

**Step 4: Plot and inspect**

Generate PNG/PDF/JSON artifacts, inspect the PNG, and record H2/H3/H4
morphology, axis covariance, time-step difference, hashes, runtime, and memory.

**Step 5: Reviews and commit**

Perform full specification and code-quality reviews.  Resolve every Critical
and Important finding and rerun affected verification before committing results.

### Task 5: Final branch verification and publication

**Files:**
- Modify: `docs/plans/2026-07-31-obsolete-dg-route-inventory.md`
- Modify: `docs/plans/2026-07-31-obsolete-dg-route-removal-results.md`

**Step 1: Repeat clean-first overlay acceptance**

Repeat Task 4 build and acceptance commands from a new temporary directory.

**Step 2: Audit routes and physics scope**

Confirm the only DG route is buffered exact-symmetry overlapping-Wannier GS to
V3 to generalized-eigenvalue Exp coefficient RT.  Confirm normal SALMON and
normal DC LCFO+EigenExa remain available.  Confirm no lattice dynamics or
dephasing claim is made for the fixed-displaced snapshot.

**Step 3: Final reviews, commit, and push**

After verification-before-completion, commit the audit, push the branch to
both `origin` and `upstream`, and verify both remote refs equal local HEAD.

