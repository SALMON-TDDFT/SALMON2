# Wannier+PW Exp Full TDDFT Validation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make one namelist-driven global Wannier+PW exponential propagation path reproduce Full TDDFT `Pz(t)` within 5 percent relative RMS at `dt=2 a.u.` under a `1e11 W/cm^2`, `1.55 eV`, `10 fs`, `sin^2`, z-polarized pulse.

**Architecture:** Promote the existing reduced Wannier+BPW exponential path from environment-controlled diagnostics to an explicit production mode selected through `&t_dg`. Keep the first validation global so fragment ownership and communication are outside the physics comparison, then use one manifest-driven analysis script to compare polarization and its identically differentiated current.

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

Expected: FAIL because the three production controls do not exist.

**Step 3: Add the minimal namelist controls**

Declare and plumb the three variables through the normal SALMON input lifecycle.
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

### Task 2: Route production mode through the reduced mixed-space Exp update

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
3. add the midpoint length-gauge field term using the complete WW, WP/PW, and PP position blocks;
4. diagonalize the Hermitian matrix;
5. apply `exp(-i H_mid dt)` to every occupied mixed coefficient vector;
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

### Task 3: Make the self-consistent full-DG eigensystem the mandatory initial state

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
It must also require a basis-generation counter and DG-operator-generation
counter with equality at diagonalization: a basis refresh followed by reuse of
old surface blocks is a fatal error.

The test must assert that non-Hermiticity, non-finite values, or a configurable
hard norm failure stops before MPI reduction. Normal runs should print only a
compact first-step summary; per-step output is enabled by the trace namelist.

**Step 2: Verify failure**

Run `python3 tests/dg/check_wpw_exp_invariants.py`.
Expected: FAIL because the production-specific invariant block is absent.

**Step 3: Implement the eigenseed and invariant checks**

Reuse and complete the existing `yn_dg_full_h_eigen_seed` generalized-overlap
path. Implement the following fixed-potential initialization loop:

```text
read and freeze V_eff from the converged reference GS
construct trial global Wannier+PW basis Phi(0)
for k = 0, 1, ...:
    rebuild S[Phi(k)]
    rebuild every volume, nonlocal, W-PW, PW-PW, and DG face/flux block
    solve H_DG[Phi(k)] C(k) = S[Phi(k)] C(k) epsilon(k)
    form the occupied projector Q(k) in the S metric
    if the basis definition is fixed: accept after residual checks
    otherwise construct Phi(k+1) from Q(k), align its gauge to Phi(k), and test convergence
end
copy the converged occupied full-DG eigenvectors into the RT coefficients
```

The effective potential stays fixed throughout this loop. Basis refresh is the
only event that invalidates and rebuilds the zero-field DG matrix. A coefficient
or eigenvector update in an unchanged basis must not rebuild the DG surface
matrix. Sort the physical eigenpairs and use
tolerances based on machine precision and matrix scale for Hermiticity,
eigen-residual, and S-orthonormality. Record norm drift without renormalizing the
propagated state. A violation must reveal the actual error rather than conceal
it.

Convergence must include all of:

```text
||Q(k)-Q(k-1)||_F / ||Q(k)||_F
max_i |epsilon_i(k)-epsilon_i(k-1)|
||H_DG(k)-T^H H_DG(k-1) T||_F / ||H_DG(k)||_F
max_i ||H_DG C_i-S C_i epsilon_i||
```

Here `T` is the gauge/alignment transformation between successive bases. For a
fixed global Wannier+PW basis, report `basis_iterations=1` and do not pretend
that repeated diagonalization is a self-consistency cycle.

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

### Task 4: Create reproducible Full TDDFT and Wannier+PW input families

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
Record whether the basis was fixed or iterated, the number of basis iterations,
the final operator generation, and exact restart provenance in comments and the
manifest.

**Step 4: Run the parser test and 20-step smoke runs**

Run the parser test, then run each input with `nt=20` copies in a scratch
directory. Expected: no guard failure, finite `Pz`, and negligible zero-field
drift in a matching field-off smoke case.

**Step 5: Commit**

```bash
git add tests/dg/check_stage2d_wpw_exp_inputs.py samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_stage2d_* samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/stage2d_wpw_exp_manifest.tsv
git commit -m "test: add Stage 2D Wannier PW Exp inputs"
```

### Task 5: Implement the Pz-primary comparison and derived Jz

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

### Task 6: Run the convergence gate and document the accepted basis

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
