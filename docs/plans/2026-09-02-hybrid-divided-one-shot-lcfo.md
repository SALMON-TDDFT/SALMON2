# Hybrid Divided-SCF and One-Shot LCFO Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Make fragment-local WF+PW density SCF followed by one complete LCFO solve the production Hybrid route, with explicit optional LCFO refinement and complete-v3 RT publication.

**Architecture:** Reuse the already implemented divided route, but replace its non-authoritative convergence, stale occupations, unchecked electron count, and ordinary-`hpsi` operator with shared conventional-DC semantics and the fixed broken-volume/SIPG payload.  Keep the repeated complete-Hybrid continuation as an explicit oracle.  Extract its final complete solve, dynamic energy-window certification, localized RT gauge, and version-3 checkpoint publication so the divided route can use them once by default or exactly `N+1` times when the user requests `N` refinements.

**Tech Stack:** Fortran 2008, SALMON DC/LCFO, MPI, ScaLAPACK, Wannier90, CMake, standalone Python source-contract tests, linked Fortran MPI fixtures.

---

The worktree is intentionally dirty.  Preserve every existing modification and
all verification directories.  Never use `git add -A`, `git commit -a`,
`git reset`, `git checkout --`, or `git clean`.  Use `git add -p` for every
already modified file, inspect `git diff --cached`, and commit only the current
task's hunks.  Do not create another worktree.

The currently running eight-rank continuation job is a reference calculation.
Do not start another heavy Si64 calculation concurrently and do not remove or
overwrite any of its output.  Source-level and small MPI tests may proceed;
wait for the reference process before the next full Si64 run.

The production scope remains the accepted Gamma, non-SOI, PZ-LDA, gapped
Hybrid route.  Do not broaden theory scope in this plan.

### Task 1: Share the authoritative conventional-DC density convergence metric

**Files:**

- Create: `src/gs/dc/dc_scf_convergence.f90`
- Create: `tests/dg/test_dc_scf_convergence_mpi.f90`
- Create: `tests/dg/run_dc_scf_convergence_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/scf_iteration_dft.f90:497-537`
- Modify: `src/gs/dc/dg_hybrid_divided_scf.f90:82-101`

**Step 1: Write the failing MPI fixture**

For a distributed density difference with known absolute sum and square sum,
require the conventional definitions:

```fortran
call reduce_dc_density_convergence(comm,'rho_dne',local_abs,local_square,&
  hvol,nelec,ngrid,value,ok,message)
call require(ok .and. abs(value-hvol*global_abs/nelec)<1d-14,&
  'rho_dne normalization changed')

call reduce_dc_density_convergence(comm,'norm_rho',local_abs,local_square,&
  hvol,nelec,ngrid,value,ok,message)
call require(ok .and. abs(value-global_square)<1d-14,&
  'norm_rho normalization changed')

call reduce_dc_density_convergence(comm,'norm_rho_dng',local_abs,local_square,&
  hvol,nelec,ngrid,value,ok,message)
call require(ok .and. abs(value-global_square/real(ngrid,real64))<1d-14,&
  'norm_rho_dng normalization changed')
```

Also require collective rejection of unsupported modes, non-finite
accumulators, non-positive volume/electron/grid counts, and rank-disagreeing
scalar controls.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dc_scf_convergence_mpi.py`

Expected: compile failure because `dc_scf_convergence` does not exist.

**Step 3: Implement the shared reducer**

Expose only this narrow API:

```fortran
module dc_scf_convergence
  use,intrinsic::iso_fortran_env,only:real64
  implicit none
  private
  public::reduce_dc_density_convergence
contains
  subroutine reduce_dc_density_convergence(comm,mode,local_abs,local_square,&
      hvol,electron_count,global_point_count,value,ok,message)
    integer,intent(in)::comm,global_point_count
    character(*),intent(in)::mode
    real(real64),intent(in)::local_abs,local_square,hvol,electron_count
    real(real64),intent(out)::value
    logical,intent(out)::ok
    character(*),intent(out)::message
```

Use one collective reduction of `[local_abs,local_square]`, then apply exactly
the formulas in `scf_iteration_dft.f90`.  Do not take square roots and do not
replace `rho_dne` by a maximum norm.

**Step 4: Route both implementations through the helper**

Keep each caller's existing local grid loop, but send its local absolute and
square sums to the shared reducer.  `norm_pot` and `norm_pot_dng` remain in the
ordinary SCF implementation and are outside the divided-route scope.

**Step 5: Run GREEN and protected checks**

Run:

```text
python3 tests/dg/run_dc_scf_convergence_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
cmake --build build-hybrid-commit -j2
```

Expected: MPI 1/2/4 ranks PASS, divided-driver tests PASS, and the build
completes.

**Step 6: Commit only Task 1**

```text
git add src/gs/dc/dc_scf_convergence.f90 tests/dg/test_dc_scf_convergence_mpi.f90 tests/dg/run_dc_scf_convergence_mpi.py
git add -p src/gs/dc/CMakeLists.txt src/gs/scf_iteration_dft.f90 src/gs/dc/dg_hybrid_divided_scf.f90
git diff --cached --check
git diff --cached
git commit -m "fix(dc): share divided density convergence semantics"
```

Checkpoint after this task: report the exact RED and GREEN evidence before
continuing.

### Task 2: Enforce divided-SCF electron count and terminal potential refresh

**Files:**

- Modify: `src/gs/dc/dg_hybrid_divided_scf.f90`
- Modify: `tests/dg/test_dg_hybrid_divided_scf_mpi.f90`
- Modify: `src/gs/main_dft.f90:3613-3622`

**Step 1: Extend the fixture RED**

Add `core_weights`, `expected_electron_count`, and `electron_tolerance` to the driver fixture.
Require collective rejection when any density callback returns NaN or when

```fortran
abs(electron_count-expected_electron_count) > electron_tolerance
```

Require the rejected path not to call the mixer or publish a converged density.
Independently integrate `sum(core_weights*new_density)` across the total
communicator and reject disagreement with either the callback count or the
target count.  Apply the same gate to the mixed density.
On convergence, require one terminal `update_total_potential(new_density)`
after the accepted density callback, so the returned density and the potential
used by final LCFO have the same potential epoch.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py`

Expected: fixture compile failure on the new driver arguments or an assertion
showing that the bad electron count was accepted.

**Step 3: Extend the driver contract**

Use this signature tail:

```fortran
...,mix_dc_density,maximum_iterations,core_weights,expected_electron_count,electron_tolerance,&
converged_density,iterations,convergence_value,electron_defect,ok,message)
```

Validate all scalars collectively before iteration.  After every core-density
callback, check finiteness and electron count before calculating convergence.
When convergence passes, refresh the potential from `new_density`; publish
nothing if that refresh fails.

**Step 4: Update the production call and run GREEN**

Pass `dc%elec_num_tot` and `dg_dc_gs_electron_count_tolerance` from
`main_dft.f90`.  Emit the final electron defect in the divided-SCF receipt.

Run:

```text
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/run_dg_dc_seed_state_mpi.py
```

Expected: all PASS for their configured rank sets.

**Step 5: Commit**

```text
git add -p src/gs/dc/dg_hybrid_divided_scf.f90 tests/dg/test_dg_hybrid_divided_scf_mpi.f90 src/gs/main_dft.f90
git diff --cached --check
git diff --cached
git commit -m "fix(dg): gate divided density electron count"
```

### Task 3: Separate fragment eigensolving from occupation-owned density

**Files:**

- Modify: `src/gs/dc/dg_hybrid_fragment_solver.f90`
- Modify: `tests/dg/test_dg_hybrid_fragment_solver_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_fragment_solver_mpi.py`

**Step 1: Write RED split-phase tests**

Require a spectrum phase that returns coefficients/eigenvalues/core norms
without accepting occupations, followed by a density phase that accepts the
current iteration's occupations.  Verify that the compatibility wrapper gives
the same density as the two explicit phases.

```fortran
call solve_dg_hybrid_fragment_spectrum(...,coefficients,eigenvalues,&
  core_norms,residual,orthogonality,ok,message)
call reconstruct_dg_hybrid_fragment_density(basis,coefficients,occupations,&
  core_mask,point_weights,density,electron_count,ok,message)
```

Reject negative/greater-than-spin-degeneracy occupations, extent mismatches,
non-finite coefficients, and duplicated/non-unique core ownership.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_fragment_solver_mpi.py`

Expected: compile failure because the two new public routines are absent.

**Step 3: Extract without changing algebra**

Move the existing generalized solve into
`solve_dg_hybrid_fragment_spectrum`; move density reconstruction and electron
integration into `reconstruct_dg_hybrid_fragment_density`.  Keep
`solve_dg_hybrid_fragment_basis` as a thin compatibility wrapper until all
callers are migrated.

**Step 4: Run GREEN**

Run:

```text
python3 tests/dg/run_dg_hybrid_fragment_solver_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
```

Expected: PASS for 1/2/4 ranks.

**Step 5: Commit**

```text
git add -p src/gs/dc/dg_hybrid_fragment_solver.f90 tests/dg/test_dg_hybrid_fragment_solver_mpi.f90 tests/dg/run_dg_hybrid_fragment_solver_mpi.py
git diff --cached --check
git diff --cached
git commit -m "refactor(dg): separate fragment spectrum and density"
```

### Task 4: Reuse one common DC chemical potential and current occupations

**Files:**

- Create: `src/gs/dc/dc_fragment_occupation.f90`
- Create: `tests/dg/test_dc_fragment_occupation_mpi.f90`
- Create: `tests/dg/run_dc_fragment_occupation_mpi.py`
- Modify: `src/gs/occupation_kernel.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/dc/dcdft.f90:743-933`
- Modify: `src/gs/main_dft.f90:4911-4939`

**Step 1: Write the common-kernel RED test**

Build two gapped fragment spectra with unequal core norms.  Require one common
chemical potential, occupations in `[0,wspin]`, and

```fortran
sum(occupation(fragment,state)*core_norm(fragment,state)) == Ne
```

within the existing electron tolerance.  Test zero temperature, finite
temperature with an explicitly supplied guard-state tail, degeneracy at the
Fermi edge, decomposition invariance, and collective rejection when the
available spectrum cannot carry `Ne`.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dc_fragment_occupation_mpi.py`

Expected: compile failure because `dc_fragment_occupation` is absent.

**Step 3: Generalize the existing authoritative numerical kernel**

Add a weighted-state entry point to `occupation_kernel.f90` and make the
existing `solve_spectrum_occupations` a shape adapter around it:

```fortran
subroutine solve_weighted_state_occupations(eigenvalues,state_weights,&
    electron_target,electronic_temperature,maximum_occupation,&
    occupations,chemical_potential,electron_count,ok,message)
```

Preserve the existing Fermi function, bracketing, exact zero-temperature
fallback, and stopping semantics.  A state weight is the unique-core norm of
that fragment state; it is not a fragment-wide k-point weight.

**Step 4: Add the distributed fragment adapter**

Expose a routine independent of SALMON structure types:

```fortran
subroutine determine_dc_fragment_occupations(comm,energies,core_norms,&
    representative_mask,temperature,wspin,expected_electrons,tolerance,&
    chemical_potential,occupations,electron_count,ok,message)
```

Only one representative per fragment contributes to the total communicator;
broadcast the result within each fragment communicator.  Do not introduce a
per-fragment chemical potential.

**Step 5: Route conventional and Hybrid callers through it**

Make `ne2mu_dcdft` use the extracted kernel after its existing spectrum/core
norm gathering.  Change `solve_dg_hybrid_divided_fragments` to:

1. solve all fragment spectra;
2. gather representative eigenvalues and core norms;
3. determine the common chemical potential and occupations;
4. reconstruct the local core density with those occupations.

Delete the stale `fragment_occupations=system%rocc(...)` assignment from the
production path.

**Step 6: Run GREEN and conventional protection**

Run:

```text
python3 tests/dg/run_dc_fragment_occupation_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_solver_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
cmake --build build-hybrid-commit -j2
```

Expected: all PASS and the build completes.

**Step 7: Commit**

```text
git add src/gs/dc/dc_fragment_occupation.f90 tests/dg/test_dc_fragment_occupation_mpi.f90 tests/dg/run_dc_fragment_occupation_mpi.py
git add -p src/gs/occupation_kernel.f90 src/gs/dc/CMakeLists.txt src/gs/dc/dcdft.f90 src/gs/main_dft.f90
git diff --cached --check
git diff --cached
git commit -m "fix(dc): share fragment occupation policy"
```

### Task 5: Build fragment self blocks from the fixed broken-volume/SIPG payload

**Files:**

- Create: `src/gs/dc/dg_hybrid_divided_operator.f90`
- Create: `tests/dg/test_dg_hybrid_divided_operator_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_divided_operator_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/main_dft.f90:3440-3612`

**Step 1: Write a RED two-fragment operator fixture**

Create a row-distributed synthetic payload with:

- broken-volume kinetic rows;
- a nonlocal projector whose support reaches a fragment not identified only by
  a Cartesian face;
- local-potential rows;
- SIPG self and cross-fragment rows; and
- a nonidentity metric.

Require `extract_dg_hybrid_fragment_self_block` to return the exact `H_ff` and
`S_ff`, including all diagonal/self pieces of the SIPG and nonlocal operators.
Require `compose_dg_hybrid_complete_rows` to reconstruct the direct full
reference with every cross-fragment contribution exactly once.  Check
Hermiticity and rank-decomposition invariance.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_divided_operator_mpi.py`

Expected: compile failure because `dg_hybrid_divided_operator` is absent.

**Step 3: Implement narrow matrix APIs**

```fortran
subroutine extract_dg_hybrid_fragment_self_block(comm_fragment,fragment_id,&
    row_ids,basis_fragment,fixed_payload,local_rows,hff,sff,ok,message)

subroutine compose_dg_hybrid_complete_rows(comm,row_ids,fixed_payload,&
    local_rows,hamiltonian_rows,operator_fingerprint,ok,message)
```

Use the frozen basis directory and row IDs; never infer a fragment from a
rank-local offset.  The fixed payload contributes kinetic, nonlocal, metric,
and SIPG rows.  Only `local_rows` changes with density.

**Step 4: Share payload construction between divided and reference routes**

Move the already accepted basis-directory, interior, nonlocal, face, and
`freeze_dg_hybrid_variational_payload` construction before the
`yn_dg_hybrid_continuation_scf`/`yn_dg_hybrid_divided_scf` branch.  Keep one
immutable payload and one fingerprint contract for both routes.

**Step 5: Run GREEN**

Run:

```text
python3 tests/dg/run_dg_hybrid_divided_operator_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_support_redistribution_mpi.py
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/check_dg_hybrid_continuation_route.py --interface-only
```

Expected: all PASS.

**Step 6: Commit**

```text
git add src/gs/dc/dg_hybrid_divided_operator.f90 tests/dg/test_dg_hybrid_divided_operator_mpi.f90 tests/dg/run_dg_hybrid_divided_operator_mpi.py
git add -p src/gs/dc/CMakeLists.txt src/gs/main_dft.f90
git diff --cached --check
git diff --cached
git commit -m "feat(dg): project fixed DG fragment self blocks"
```

### Task 6: Connect the production divided loop to the fixed DG operator

**Files:**

- Modify: `src/gs/main_dft.f90:3613-3668,4867-4995`
- Modify: `tests/dg/check_dg_hybrid_divided_dc_controls.py`
- Modify: `tests/dg/check_dg_hybrid_divided_lcfo_route.py`
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in`
- Modify: `src/io/inputoutput.f90:3130-3170`

**Step 1: Make the route checks RED**

Require the divided branch to use
`extract_dg_hybrid_fragment_self_block` and the split spectrum/occupation/
density phases.  Forbid production calls to
`apply_dg_hybrid_divided_fragment_hpsi` and the identity metric callback.
Require `yn_dg_hybrid_divided_scf`, `yn_dg_hybrid_continuation_scf`, and
`yn_dg_hybrid_scf` to be mutually exclusive.

**Step 2: Run RED**

Run:

```text
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
```

Expected: FAIL because the old divided callback still calls ordinary `hpsi`
with an identity metric and flag exclusivity is absent.

**Step 3: Wire the divided callbacks**

At each divided iteration:

1. update the total potential from the current unique-core density;
2. assemble the density-dependent local rows;
3. extract and solve each fragment `H_ff C_f=S_ff C_f epsilon_f`;
4. determine current occupations with the common DC chemical potential;
5. reconstruct unique-core density; and
6. use the shared DC convergence/mixing routines.

The local loop must contain no complete generalized eigensolve.  The old
ordinary-`hpsi` callbacks may remain only for isolated compatibility fixtures;
they are not production entry points.

**Step 4: Run GREEN and build**

Run:

```text
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/run_dg_hybrid_divided_operator_mpi.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
cmake --build build-hybrid-commit -j2
```

Expected: all PASS and no complete eigensolve appears inside the divided SCF
loop.

**Step 5: Commit**

```text
git add tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in
git add -p src/gs/main_dft.f90 src/io/inputoutput.f90 tests/dg/check_dg_hybrid_divided_dc_controls.py tests/dg/check_dg_hybrid_divided_lcfo_route.py
git diff --cached --check
git diff --cached
git commit -m "feat(dg): run production fragment DG density SCF"
```

### Task 7: Extract one terminal LCFO certification and complete-v3 publication

**Files:**

- Create: `src/gs/dc/dg_hybrid_lcfo_finalization.f90`
- Create: `tests/dg/test_dg_hybrid_lcfo_finalization_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_lcfo_finalization_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/main_dft.f90:5060-5730`
- Modify: `tests/dg/check_dg_hybrid_divided_lcfo_route.py`
- Modify: `tests/dg/check_dg_hybrid_continuation_route.py`

**Step 1: Write RED finalization tests**

Use a small complete eigensystem with an occupied block and a degenerate
unoccupied boundary cluster.  Require one call to the complete generalized
solver, current-spectrum occupations, extension from the requested cutoff to
the first symmetry-closing cluster, physical occupied/window symmetry checks,
unconstrained construction WFs, and a complete version-3 payload.  Require
construction-basis nonclosure to remain diagnostic.  Reject missing operator,
density, selection, pseudopotential, or ownership fingerprints.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_lcfo_finalization_mpi.py`

Expected: compile failure because the finalization module is absent.

**Step 3: Extract the terminal logic**

Expose a routine that performs no density iteration:

```fortran
subroutine finalize_dg_hybrid_lcfo_once(comm,operator_epoch,row_ids,&
    hamiltonian_rows,fixed_payload,basis_values,grid_ids,grid_weights,&
    physical_symmetry,selection,provenance,solve_complete,&
    result,lcfo_density,checkpoint_payload,receipt,ok,message)
```

Move, without weakening, the current continuation logic for:

- `solve_dg_hybrid_generalized_complete_once`;
- `derive_dg_hybrid_occupation_policy`;
- occupied density/projector reconstruction;
- requested/certified/proof energy-window selection;
- physical occupied/window symmetry certification;
- certified RT basis localization;
- component and energy receipts; and
- complete-v3 payload authentication/publication preparation.

Do not move the continuation stage scheduler, lambda trials, density mixing, or
repeated candidate loop.  Publication remains a separate final call so an
optional refinement does not write intermediate checkpoints.

**Step 4: Route both algorithms through the terminal routine**

The divided route composes `H[rho_divided]` and calls the routine once.  The
reference continuation calls it only after its accepted final refresh.  The
default divided path must not call the old version-2 occupied-checkpoint
publisher.

**Step 5: Run GREEN**

Run:

```text
python3 tests/dg/run_dg_hybrid_lcfo_finalization_mpi.py
python3 tests/dg/run_dg_hybrid_occupation_policy_mpi.py
python3 tests/dg/run_dg_hybrid_low_energy_symmetry_mpi.py
python3 tests/dg/run_dg_hybrid_certified_rt_basis_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
```

Expected: all PASS; zero-refinement divided fixtures report exactly one
complete eigensolve and version 3.

**Step 6: Commit**

```text
git add src/gs/dc/dg_hybrid_lcfo_finalization.f90 tests/dg/test_dg_hybrid_lcfo_finalization_mpi.f90 tests/dg/run_dg_hybrid_lcfo_finalization_mpi.py
git add -p src/gs/dc/CMakeLists.txt src/gs/main_dft.f90 tests/dg/check_dg_hybrid_divided_lcfo_route.py tests/dg/check_dg_hybrid_continuation_route.py
git diff --cached --check
git diff --cached
git commit -m "feat(dg): finalize one-shot LCFO into certified RT space"
```

### Task 8: Add explicit user-controlled LCFO refinement epochs

**Files:**

- Create: `src/gs/dc/dg_hybrid_lcfo_refinement.f90`
- Create: `tests/dg/test_dg_hybrid_lcfo_refinement_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_lcfo_refinement_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/io/salmon_global.f90:480-490`
- Modify: `src/io/inputoutput.f90:630-645,1155-1170,1885-1905,2955-2975,3130-3170`
- Modify: `src/gs/main_dft.f90:3613-3670`
- Modify: `src/rt/dg/rt_dg_hybrid_checkpoint.f90`
- Modify: `src/rt/dg/rt_dg_hybrid_initialization.f90:640-675`
- Modify: `src/rt/dg/rt_dg_hybrid_stationarity.f90`
- Modify: `src/rt/main_tddft.f90:330-405`
- Modify: `tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_initialization_mpi.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_stationarity_mpi.f90`
- Modify: `tests/dg/check_dg_hybrid_localization_first_inputs.py`
- Modify: `tests/dg/check_dg_hybrid_divided_lcfo_route.py`

**Step 1: Write RED controller and input tests**

Require `dg_hybrid_lcfo_refine_steps=0` by default and broadcast it on all
ranks.  Reject negative values, positive values unless the divided route is
selected, and simultaneous divided/continuation/full-Hybrid route flags.  A
callback-counting fixture must require:

- `0` refinements: one complete solve;
- `2` refinements: three complete solves, two density mixes, and two potential
  updates;
- fixed metric/kinetic/nonlocal/SIPG fingerprints at every epoch;
- a new local/total operator fingerprint after each density update;
- no intermediate publication; and
- final coefficients and eigenpair receipt from the final operator epoch.

**Step 2: Run RED**

Run:

```text
python3 tests/dg/run_dg_hybrid_lcfo_refinement_mpi.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
```

Expected: compile/source-contract failure because the input and controller are
absent.

**Step 3: Implement the explicit state machine**

```fortran
call solve_once(epoch=0)
do refinement=1,refine_steps
  call measure_shared_dc_density_difference(rho_potential,rho_lcfo,delta_rho)
  call mix_density(rho_potential,rho_lcfo,rho_next)
  call update_potential(rho_next)
  call assemble_local_and_complete_rows(rho_next,epoch=refinement)
  call solve_once(epoch=refinement)
enddo
call publish_final_epoch_only()
```

Record route code, requested/completed refinement count, final `Delta rho`,
potential-density fingerprint, and final operator epoch in the existing
version-3 continuation receipt/fingerprint fields.  Use the two existing v3
density arrays without changing the serialized layout:

- `payload%density`: the density that generated the final stored operator;
- `payload%rt_space%density`: the density reconstructed from the final LCFO
  state and used to initialize RT.

Update authentication so the reconstructed certified density must equal only
`rt_space%density`; authenticate the potential density through the overall
payload fingerprint and epoch receipt.  Existing continuation v3 files, where
the two densities are equal, remain valid.  Do not silently add a refinement
based on `Delta rho`.

At RT startup, first validate the stored eigenpair against the stored operator
and its potential-density epoch.  Then reconstruct the orbital density, verify
`rt_space%density`, update the TD Hamiltonian from that physical density, and
emit the initial density/Hamiltonian change.  Split zero-field stationarity
into measurement and enforcement:

- continuation-reference checkpoints retain the existing strict tolerance
  gate;
- divided one-shot/refined checkpoints require finite authenticated receipts
  and report drift without a magnitude gate; and
- non-finite drift or any broken fingerprint/epoch remains fatal in every
  mode.

Do not implement measurement mode by passing artificially large tolerances.

**Step 4: Run GREEN**

Run:

```text
python3 tests/dg/run_dg_hybrid_lcfo_refinement_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_finalization_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/run_rt_dg_hybrid_initialization_mpi.py
python3 tests/dg/run_rt_dg_hybrid_stationarity_mpi.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
cmake --build build-hybrid-commit -j2
```

Expected: all PASS and the default route still counts one solve.

**Step 5: Commit**

```text
git add src/gs/dc/dg_hybrid_lcfo_refinement.f90 tests/dg/test_dg_hybrid_lcfo_refinement_mpi.f90 tests/dg/run_dg_hybrid_lcfo_refinement_mpi.py
git add -p src/gs/dc/CMakeLists.txt src/io/salmon_global.f90 src/io/inputoutput.f90 src/gs/main_dft.f90 src/rt/dg/rt_dg_hybrid_checkpoint.f90 src/rt/dg/rt_dg_hybrid_initialization.f90 src/rt/dg/rt_dg_hybrid_stationarity.f90 src/rt/main_tddft.f90 tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90 tests/dg/test_rt_dg_hybrid_initialization_mpi.f90 tests/dg/test_rt_dg_hybrid_stationarity_mpi.f90 tests/dg/check_dg_hybrid_localization_first_inputs.py tests/dg/check_dg_hybrid_divided_lcfo_route.py
git diff --cached --check
git diff --cached
git commit -m "feat(dg): add explicit LCFO refinement epochs"
```

### Task 9: Validate Si64 one-shot GS and separate-directory zero-field RT

**Files:**

- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in`
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_dg_zero_field_rt.in`
- Modify: `tests/dg/run_dg_hybrid_si64_divided_lcfo.py`
- Modify: `tests/dg/run_dg_hybrid_si64_continuation_rt.py` only to share parser helpers when duplication would otherwise occur

**Step 1: Extend parser-only acceptance RED**

The runner must require:

- a reused conventional DC seed and no `DC #SCF =` lines in the reused GS;
- the same eight ranks and exact rank--fragment mapping;
- unconstrained WF localization and a dynamic, non-hard-coded retained rank;
- divided-SCF convergence, common chemical potential, and per-iteration
  electron receipts;
- complete SIPG/nonlocal final operator receipts;
- exactly one complete eigensolve when refinement count is zero;
- requested/certified/proof energy-window receipts;
- physical occupied/window/density symmetry receipts;
- a complete-v3 checkpoint and matching GS/RT fingerprints;
- `Delta rho` recorded as a diagnostic; and
- finite zero-field stationarity measurements for every requested RT step,
  without applying the continuation-reference magnitude gate to divided mode.

Add negative parser fixtures for a hidden second solve, stale density/operator
epoch, basis-closure used as an acceptance gate, rank/mapping mismatch, and an
old version-2 occupied checkpoint.

**Step 2: Run parser RED**

Run: `python3 tests/dg/run_dg_hybrid_si64_divided_lcfo.py --parser-only`

Expected: FAIL until the new receipts are parsed and required.

**Step 3: Implement the runner and fresh-directory policy**

Reuse the already generated final canonical DC seed; do not repeat
conventional DC.  Create separate fresh directories for one-shot GS and RT.
Hash the seed before and after each run and preserve all logs.  Never overwrite
the current continuation-oracle directory.

Record the one-shot zero-field drift even when it exceeds the strict
continuation-reference tolerance.  Run a second fresh GS with a small explicit
refinement count and quantify the improvement.  Non-finite drift remains a
failure.  Do not turn the observation into automatic fallback logic.

**Step 4: Run focused GREEN before Si64**

Run:

```text
python3 tests/dg/run_dc_scf_convergence_mpi.py
python3 tests/dg/run_dc_fragment_occupation_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/run_dg_hybrid_divided_operator_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_finalization_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_refinement_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/run_rt_dg_hybrid_initialization_mpi.py
python3 tests/dg/run_rt_dg_hybrid_stationarity_mpi.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
cmake --build build-hybrid-commit -j2
```

Expected: all PASS.

**Step 5: Run fresh eight-rank one-shot and RT**

Run only after the current heavy reference job has ended:

```text
OMP_NUM_THREADS=1 python3 tests/dg/run_dg_hybrid_si64_divided_lcfo.py --mpi-ranks 8 --binary build-hybrid-commit/salmon --result-dir verification-si64-divided-one-shot-lcfo-20260902
```

Expected: reused DC seed, converged divided SCF, one final LCFO solve, complete
version-3 checkpoint, and a separate zero-field RT result.  Preserve failure
evidence if any gate fails.

**Step 6: Compare with the reference oracle**

Record energy, band gap, occupied-projector difference, density difference,
certified rank, symmetry defects, wall time, and peak RSS.  The comparison is
validation evidence, not a post-LCFO density convergence gate.

**Step 7: Commit the runner contract**

```text
git add tests/dg/run_dg_hybrid_si64_divided_lcfo.py
git add -p tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_dg_zero_field_rt.in tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in tests/dg/run_dg_hybrid_si64_continuation_rt.py
git diff --cached --check
git diff --cached
git commit -m "test(dg): certify divided one-shot LCFO and RT"
```

### Task 10: Run protected regressions, review, and finish the branch

**Files:**

- Review all files changed in Tasks 1-9
- Preserve all untracked verification output

**Step 1: Run all focused and protected verification fresh**

Run:

```text
git diff --check
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
python3 tests/dg/run_dg_dc_seed_checkpoint_mpi.py
python3 tests/dg/run_dg_dc_seed_state_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_solver_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_mpi.py
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/run_dg_hybrid_occupation_policy_mpi.py
python3 tests/dg/run_dg_hybrid_low_energy_symmetry_mpi.py
python3 tests/dg/run_dg_hybrid_certified_rt_basis_mpi.py
python3 tests/dg/run_dg_hybrid_continuation_controller_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/run_rt_dg_hybrid_initialization_mpi.py
python3 tests/dg/run_rt_dg_hybrid_stationarity_mpi.py
python3 tests/dg/run_rt_dg_hybrid_length_gauge_mpi.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
cmake --build build-hybrid-commit -j2
ctest --test-dir build-hybrid-commit --output-on-failure
```

Expected: all standalone checks PASS, the build completes, and CTest has no
failures.

**Step 2: Request code review**

Invoke `@superpowers:requesting-code-review`.  Review especially:

- exact conventional-DC convergence equivalence;
- common chemical potential and electron count;
- SIPG/nonlocal terms exactly once;
- no complete solve inside fragment SCF;
- zero-refinement solve count exactly one;
- physical-space rather than individual-WF symmetry acceptance;
- version-3 density/operator epoch provenance; and
- exact rank-count/rank--fragment DC-seed reuse.

**Step 3: Resolve findings rigorously**

Use `@superpowers:receiving-code-review`.  For each Critical or Important
finding, reproduce it with a RED test, apply the minimal fix, and rerun every
affected check.  Commit review fixes separately.

**Step 4: Verify before any completion claim**

Invoke `@superpowers:verification-before-completion` and rerun affected commands
fresh.  A preserved prior log is evidence for comparison, not a substitute for
fresh verification of changed code.

**Step 5: Finish without automatic merge**

Invoke `@superpowers:finishing-a-development-branch`.  Present the integration
options and do not merge, delete the worktree, clean output, or push unless the
user explicitly requests it.
