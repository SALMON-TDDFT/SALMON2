# Wannier+PW Exp Full TDDFT Validation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make one namelist-driven global Wannier+PW exponential propagation path reproduce Full TDDFT `Pz(t)` within 5 percent relative RMS at `dt=2 a.u.` under a `1e11 W/cm^2`, `1.55 eV`, `10 fs`, `sin^2`, z-polarized pulse.

**Architecture:** Promote the existing reduced Wannier+BPW exponential path from environment-controlled diagnostics to an explicit production mode selected through `&t_dg`. Use one fixed global Wannier+PW basis with the complete DG interface operator, solve its generalized eigenproblem for the initial occupied states, and propagate with an `S`-metric midpoint exponential while updating the TDDFT Hartree and exchange-correlation potential. Keep MPI coefficient ownership and communication outside the first physics comparison, then use one manifest-driven analysis script to compare polarization and its identically differentiated current.

**Tech Stack:** Fortran 2008, MPI, LAPACK/ScaLAPACK or EigenExa through the existing build, Python 3 standard library plus NumPy/Matplotlib for analysis, CMake/CTest-style source checks.

---

### Task 1: Freeze the production-mode input contract

**Files:**
- Create: `tests/dg/check_wpw_exp_production_input.py`
- Modify: `src/io/salmon_global.f90:157-201`
- Modify: `src/io/inputoutput.f90:317-355,849-887,1497-1528,2490-2545,3179-3230,3615-3635`

**Step 1: Write the failing source-contract test**

Add a Python test that reads `salmon_global.f90` and `inputoutput.f90` and asserts that:

```python
required = {
    "yn_dg_wpw_exp_production": "character(1)",
    "yn_dg_wpw_exp_trace": "character(1)",
}
```

Each variable must be declared, included in the `&t_dg` namelist, initialized,
broadcast, written to `variables.log`, and checked with `yn_argument_check`.
The test must also assert that production mode requires `yn_dg_fragment_rt='y'`,
`yn_dg_length_gauge='y'`, `yn_plane_wave_basis='y'`,
`time_integrator_dg_fragment='expdiag'`, and global coefficient ownership.

**Step 2: Run the test and verify failure**

Run:

```bash
python3 tests/dg/check_wpw_exp_production_input.py
```

Expected: FAIL because the two production controls do not exist.

**Step 3: Add the minimal namelist controls**

Declare and plumb the two variables through the normal SALMON input lifecycle.
Use defaults:

```fortran
yn_dg_wpw_exp_production = 'n'
yn_dg_wpw_exp_trace = 'n'
```

Add explicit fatal guards for unsupported combinations. Do not silently fall
back to local BPW propagation.

**Step 4: Run the focused test**

Run `python3 tests/dg/check_wpw_exp_production_input.py`.
Expected: PASS.

**Step 5: Commit**

```bash
git add tests/dg/check_wpw_exp_production_input.py src/io/salmon_global.f90 src/io/inputoutput.f90
git commit -m "feat: add Wannier PW Exp production controls"
```

### Task 2: Add numerical tests for generalized eigenstates and S-metric Exp

**Files:**
- Create: `tests/dg/test_wpw_generalized_exp.py`
- Create: `tests/dg/fixtures/wpw_generalized_exp_driver.F90`
- Create: `src/rt/dg/rt_dg_generalized_exp.f90`
- Modify: `src/rt/dg/CMakeLists.txt`
- Modify: `tests/dg/CMakeLists.txt`

**Step 1: Write the failing numerical tests**

Construct deterministic small complex Hermitian matrices `H`, positive-definite
overlaps `S`, and position matrices `Z`. The Fortran driver and Python oracle
must test:

```text
H C = S C epsilon
C^H S C = I
U^H S U = S
C(dt) = X exp[-i X^H H X dt] X^H S C(0)
```

Add a negative test that applies ordinary Euclidean projection/exponentiation
when `S /= I` and must exceed the accepted error. Add another negative test that
removes a nonzero DG face block and must fail the eigen-residual or reference
eigenvalue comparison.

**Step 2: Run the tests and verify failure**

Run:

```bash
python3 -m unittest tests.dg.test_wpw_generalized_exp -v
```

Expected: FAIL because the numerical driver and production generalized-Exp
helper do not exist.

**Step 3: Implement a reusable S-metric helper**

Create or extract a production helper that diagonalizes `S`, rejects null or
ill-conditioned metric directions according to a documented threshold, builds
`X=S^{-1/2}`, diagonalizes `X^H H X`, and applies the exponential with the
correct forward and backward metric transforms. Return generalized residual,
S-orthonormality, and S-unitarity diagnostics.

**Step 4: Run the numerical tests**

Run `python3 -m unittest tests.dg.test_wpw_generalized_exp -v`.
Expected: PASS, including both negative controls.

**Step 5: Commit**

```bash
git add tests/dg/test_wpw_generalized_exp.py tests/dg/fixtures/wpw_generalized_exp_driver.F90 tests/dg/CMakeLists.txt src/rt/dg/rt_dg_generalized_exp.f90 src/rt/dg/CMakeLists.txt
git commit -m "test: validate generalized Wannier PW exponential"
```

### Task 3: Route production mode through the reduced mixed-space Exp update

**Files:**
- Create: `tests/dg/check_wpw_exp_production_route.py`
- Modify: `src/rt/dg/rt_dg_integrator_expdiag.f90:1-160,802-986,2277-2593`
- Modify: `src/rt/dg/rt_dg_fragment_types.f90`

**Step 1: Write the failing route test**

The source test must assert that production enablement comes only from
`yn_dg_wpw_exp_production`, initialization requires
`yn_dg_full_h_eigen_seed='y'`, and tracing comes from
`yn_dg_wpw_exp_trace`. It must reject calls to `get_environment_variable`
for `SALMON_DG_WPW_REDUCED_EXPDIAG`,
`SALMON_DG_WPW_REDUCED_INIT_PROJECT`, and the reduced Exp trace controls. It
must also reject `initialize_wpw_reduced_self_projection` from the production
branch.

It must also assert a production call to a clearly named routine such as:

```fortran
call propagate_wpw_exp_production(state_first, state_last, dt, itt, E_mid, Ac_ham)
```

and require a fatal stop if the mixed Hamiltonian, mixed position operator, or
initial mixed coefficients are unavailable.

**Step 2: Verify the test fails**

Run `python3 tests/dg/check_wpw_exp_production_route.py`.
Expected: FAIL on the environment-variable reads and missing production routine.

**Step 3: Extract the production routine**

Rename/refactor the scientifically required part of
`dryrun_wpw_reduced_expdiag` into `propagate_wpw_exp_production`. Keep optional
diagnostics separate. The production routine must:

1. require occupied initial coefficients from the generalized eigensystem of
   the complete zero-field DG Hamiltonian, including interface/flux terms;
2. construct or reuse the Hermitian mixed Hamiltonian;
3. add the midpoint length-gauge field term using the complete WW, WP/PW, PP,
   and required DG-interface position blocks;
4. transform the generalized problem with the validated `S^{-1/2}` helper;
5. apply the `S`-metric exponential to every occupied mixed coefficient vector;
6. write the propagated Wannier and BPW-perpendicular coefficients back;
7. stop rank-locally on non-finite coefficients before any collective reduction.

Delete the three required environment-variable switches after their namelist
replacements are active. Preserve purely diagnostic environment variables only
when they do not alter the scientific result.

**Step 4: Run focused tests and compile**

Run:

```bash
python3 tests/dg/check_wpw_exp_production_input.py
python3 tests/dg/check_wpw_exp_production_route.py
cmake --build build-mpi-eigenexa -j 2
```

Expected: both tests PASS and the build completes successfully.

**Step 5: Commit**

```bash
git add tests/dg/check_wpw_exp_production_route.py src/rt/dg/rt_dg_integrator_expdiag.f90 src/rt/dg/rt_dg_fragment_types.f90
git commit -m "feat: promote mixed Wannier PW exponential propagation"
```

### Task 4: Make the fixed-basis full-DG eigensystem the mandatory initial state

**Files:**
- Create: `tests/dg/check_wpw_exp_invariants.py`
- Modify: `src/rt/dg/rt_dg_fragment.f90`
- Modify: `src/rt/dg/rt_dg_integrator_expdiag.f90`
- Modify: `src/rt/dg/rt_dg_observables.f90`

**Step 1: Write the failing invariant test**

Require production-mode checks for:

```text
max occupied residual ||H_DG C - S C epsilon||
S-orthonormality ||C^H S C - I||
max_abs(H-H^H)
occupied mixed-space norm before/after Exp
maximum coefficient magnitude
Pz at the initial and first propagated steps
```

The test must require the full Wannier+PW DG Hamiltonian, including DG
interface/flux blocks, and the explicit overlap matrix in the generalized
eigenproblem. It must reject initialization by projection alone and require the
lowest occupied physical full-DG eigenvectors to populate the RT coefficients.
It must require that the complete DG face/flux blocks are assembled after the
fixed Wannier+PW basis is finalized. A test configuration with a required face
block removed must fail.

The test must assert that non-Hermiticity, non-finite values, or a configurable
hard norm failure stops before MPI reduction. Normal runs should print only a
compact first-step summary; per-step output is enabled by the trace namelist.

**Step 2: Verify failure**

Run `python3 tests/dg/check_wpw_exp_invariants.py`.
Expected: FAIL because the production-specific invariant block is absent.

**Step 3: Implement the eigenseed and invariant checks**

Reuse and complete the existing `yn_dg_full_h_eigen_seed` generalized-overlap
path. Implement the following fixed-potential initialization:

```text
read and freeze V_eff from the converged reference GS
construct and freeze global Wannier+PW basis Phi
rebuild S[Phi]
rebuild every volume, nonlocal, W-PW, PW-PW, and DG face/flux block
solve H_DG[Phi] C = S[Phi] C epsilon
copy the occupied full-DG eigenvectors into the RT coefficients
```

The effective potential stays fixed throughout eigenseed construction. The
Wannier+PW basis is not regenerated from the DG eigenvectors in this milestone.
A coefficient or eigenvector update in the unchanged basis must not rebuild the
DG surface matrix. Sort the physical eigenpairs and use
tolerances based on machine precision and matrix scale for Hermiticity,
eigen-residual, and S-orthonormality. Record norm drift without renormalizing the
propagated state. A violation must reveal the actual error rather than conceal
it.

Acceptance must include all of:

```text
max_i ||H_DG C_i-S C_i epsilon_i||
||C^H S C-I||_F
field-off ||Q(t)-Q(0)||_F / ||Q(0)||_F
```

Add a field-off stationarity smoke test: after removing the analytically
expected eigenphases, `Pz`, density, and occupied-subspace projector must remain
stationary within numerical tolerance.

**Step 4: Test and build**

Run the three `check_wpw_exp_*.py` tests and
`cmake --build build-mpi-eigenexa -j 2`.
Expected: PASS.

**Step 5: Commit**

```bash
git add tests/dg/check_wpw_exp_invariants.py src/rt/dg/rt_dg_fragment.f90 src/rt/dg/rt_dg_integrator_expdiag.f90 src/rt/dg/rt_dg_observables.f90
git commit -m "feat: seed Wannier PW Exp from full DG eigenstates"
```

### Task 5: Validate the complete length-gauge position operator

**Files:**
- Create: `tests/dg/check_wpw_position_operator.py`
- Modify: `src/rt/dg/rt_dg_fragment.f90`
- Modify: `src/rt/dg/rt_dg_mixed_fsum_diagnose.f90`

**Step 1: Write failing position-operator tests**

Require the production position operator to record WW, WP/PW, PP, and DG
interface contributions separately. Test Hermiticity after transformation to
the orthonormal metric, the existing momentum-position consistency diagnostic,
and odd response under `E_z -> -E_z`. A negative fixture that suppresses a
nonzero interface contribution must fail.

**Step 2: Run and verify failure**

Run `python3 tests/dg/check_wpw_position_operator.py`.
Expected: FAIL while neighbor/interface projection is incomplete.

**Step 3: Complete and validate the operator**

Implement the missing BLAS-style projection of required neighbor position/flux
blocks. Do not repair non-Hermiticity by unconditional symmetrization before
reporting the unsymmetrized residual.

**Step 4: Test and commit**

Run the focused test and build, then commit the operator and test together.

### Task 6: Add midpoint self-consistent TDDFT potential updates

**Files:**
- Create: `tests/dg/check_wpw_exp_midpoint_scf.py`
- Modify: `src/rt/dg/rt_dg_integrator_expdiag.f90`
- Modify: `src/rt/dg/rt_dg_density_hamiltonian_update.f90`

**Step 1: Write failing numerical and route tests**

Require production field-on propagation to reconstruct density, Hartree, and
exchange-correlation potentials inside a midpoint predictor-corrector. Assert
that DG kinetic/surface blocks are not rebuilt when the basis is unchanged, but
density-dependent potential matrix blocks are rebuilt. A negative test with the
potential update disabled must be labeled independent-particle and must not be
accepted as Full TDDFT mode.

**Step 2: Verify failure**

Run `python3 tests/dg/check_wpw_exp_midpoint_scf.py`.
Expected: FAIL because the production Exp path does not yet close the midpoint
density/potential loop.

**Step 3: Implement the midpoint predictor-corrector**

At each step, predict with the input Hamiltonian, form the midpoint density,
update `V_H+V_xc`, reconstruct only the potential-dependent matrix blocks, and
repeat until a documented density or Hamiltonian tolerance is reached. Stop on
non-convergence; do not silently accept the last iterate.

**Step 4: Test, build, and commit**

Run the focused numerical tests and `cmake --build build-mpi-eigenexa -j 2`, then
commit the midpoint implementation and tests together.

### Task 7: Create reproducible Full TDDFT and Wannier+PW input families

**Files:**
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_full_tddft_i1e11_dt2`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_global_wpw_exp_pw16_i1e11_dt2`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_global_wpw_exp_pw64_i1e11_dt2`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_global_wpw_exp_pw128_i1e11_dt2`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/stage2d_wpw_exp_manifest.tsv`
- Create: `tests/dg/check_stage2d_wpw_exp_inputs.py`

**Step 1: Write the failing input-consistency test**

Parse the four inputs and assert identical system geometry, GS provenance,
`dt=2`, pulse energy, duration, envelope, intensity, polarization direction,
total simulated time, and volume convention. Assert that only the declared
Wannier/PW basis fields differ across the three reduced-space inputs.

The manifest columns are:

```text
path basis_id WF_count PW_count PW_cutoff_or_shell dt propagator_kind observable_source gauge volume_normalization
```

**Step 2: Verify failure**

Run `python3 tests/dg/check_stage2d_wpw_exp_inputs.py`.
Expected: FAIL because the Stage 2D inputs and manifest are absent.

**Step 3: Add the input family and manifest**

Derive all four inputs from the same verified converged GS, not from a smoke or
unconverged GS. Set `yn_dg_wpw_exp_production='y'` and
`yn_dg_full_h_eigen_seed='y'` for the three reduced-space inputs. The imported
GS supplies the density/potential and basis provenance; it does not directly
supply the occupied RT coefficients. Those coefficients must be replaced by
the converged occupied eigenstates of the complete DG Wannier+PW Hamiltonian.
Record that the basis is fixed, the full-DG eigenseed residual, overlap metric
threshold, midpoint SCF tolerance, and exact restart provenance in comments and
the manifest.

**Step 4: Run the parser test and 20-step smoke runs**

Run the parser test, then run each input with `nt=20` copies in a scratch
directory. Expected: no guard failure, finite `Pz`, and negligible zero-field
drift in a matching field-off smoke case.

**Step 5: Commit**

```bash
git add tests/dg/check_stage2d_wpw_exp_inputs.py samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_* samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/stage2d_wpw_exp_manifest.tsv
git commit -m "test: add Stage 2D Wannier PW Exp inputs"
```

### Task 8: Implement the Pz-primary comparison and derived Jz

**Files:**
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/compare_stage2d_wpw_exp.py`
- Create: `tests/dg/test_compare_stage2d_wpw_exp.py`
- Create: `tests/dg/fixtures/stage2d_full_pz.data`
- Create: `tests/dg/fixtures/stage2d_wpw_pz.data`

**Step 1: Write failing unit tests**

Test that the analysis tool:

- aligns identical time samples and rejects extrapolation;
- removes only a documented common initial polarization offset;
- computes `rms(P_wpw-P_full)/rms(P_full)`;
- uses the same second-order centered difference for both currents, with
  second-order one-sided endpoint differences;
- reports PASS only when `rel_rms <= 0.05`;
- does not use the current as a substitute for the polarization gate.

**Step 2: Verify failure**

Run `python3 -m unittest tests.dg.test_compare_stage2d_wpw_exp -v`.
Expected: FAIL because the comparison module does not exist.

**Step 3: Implement the analysis tool**

Accept the manifest plus Full TDDFT polarization data. Write:

```text
stage2d_wpw_exp_summary.tsv
stage2d_wpw_exp_waveforms.tsv
stage2d_wpw_exp_comparison.png
```

The summary must include basis metadata, `Pz_rel_rms`, `Jz_rel_rms`, maximum
absolute errors, zero-field drift when available, and PASS/FAIL.

**Step 4: Run unit tests**

Run `python3 -m unittest tests.dg.test_compare_stage2d_wpw_exp -v`.
Expected: PASS.

**Step 5: Commit**

```bash
git add tests/dg/test_compare_stage2d_wpw_exp.py tests/dg/fixtures/stage2d_* samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/compare_stage2d_wpw_exp.py
git commit -m "feat: compare Wannier PW polarization with Full TDDFT"
```

### Task 9: Run the convergence gate and document the accepted basis

**Files:**
- Modify: `doc/NOTE_DG.md:295-320,420-535`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/stage2d_wpw_exp_summary.tsv`

**Step 1: Build and run the reference matrix**

Build with:

```bash
cmake --build build-mpi-eigenexa -j 2
```

Run the Full TDDFT reference and PW16/PW64/PW128 Wannier+PW Exp cases through
the full pulse and agreed post-pulse interval. Do not interpret a run until the
actual input, GS provenance, pulse mode, and output normalization are recorded.

**Step 2: Run the comparison gate**

Run:

```bash
python3 samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/compare_stage2d_wpw_exp.py \
  --manifest samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/stage2d_wpw_exp_manifest.tsv \
  --full-pz <full-tddft-polarization-file> \
  --threshold 0.05
```

Expected: at least one converged Wannier+PW basis reports `PASS` with
`Pz_rel_rms <= 0.05`. If none passes, stop and classify the discrepancy using
the design document before changing fragment/distributed code.

**Step 3: Run all focused regression tests**

Run:

```bash
python3 tests/dg/check_wpw_exp_production_input.py
python3 tests/dg/check_wpw_exp_production_route.py
python3 tests/dg/check_wpw_exp_invariants.py
python3 -m unittest tests.dg.test_wpw_generalized_exp -v
python3 tests/dg/check_wpw_position_operator.py
python3 tests/dg/check_wpw_exp_midpoint_scf.py
python3 tests/dg/check_stage2d_wpw_exp_inputs.py
python3 -m unittest tests.dg.test_compare_stage2d_wpw_exp -v
cmake --build build-mpi-eigenexa -j 2
```

Expected: all tests PASS and the build succeeds.

**Step 4: Document the accepted configuration**

Update `NOTE_DG.md` with the production namelist, GS provenance, accepted
Wannier/PW counts, comparison window, `Pz_rel_rms`, derived `Jz_rel_rms`, and
known limitations. State explicitly that local symmetry is not constrained and
that distributed fragment validation is the next milestone.

**Step 5: Commit**

```bash
git add doc/NOTE_DG.md samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/stage2d_wpw_exp_summary.tsv
git commit -m "docs: record Wannier PW Exp Full TDDFT validation"
```
