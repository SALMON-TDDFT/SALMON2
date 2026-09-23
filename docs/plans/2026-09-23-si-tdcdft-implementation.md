# Si TDCDFT Implementation Plan

> Execute task by task using the executing-plans skill after design approval.

**Goal:** Add an opt-in macroscopic LRC xc vector potential and validate Si impulse and pump–probe response.

**Architecture:** Store xc history separately from classical electromagnetic fields in s_rt. Couple the effective potential consistently to propagation, physical current and nonlocal projectors; reuse existing Si ground-state and delayed impulse facilities.

**Tech Stack:** Fortran, CMake, existing SALMON test suite; Python for response analysis.

Status: user-approved implementation on TDCDFT. Tasks 1–5 implemented and tested; tasks 6–7 have inputs and analysis tests, but physical convergence and pump–probe validation remain pending. See ../tdcdft-validation.md for actual evidence.

## Task 1: Establish conventions and reference tests

Read src/common/density_matrix.f90, src/rt/em_field.f90, src/rt/time_evolution_step.f90, src/rt/initialization_rt.f90 and the optical Fourier-analysis implementation. Trace the charge-current sign, A/c normalization, timestep indexing and nonlocal-current phase. Record the derivation and exact discrete equation in the design before implementation.

Build an unmodified serial reference out of source using configure.py/CMake supported options. Run testsuites/111_bulk_Si_gs_dp and 112_bulk_Si_rt_response_dp through their existing CMake harness. Preserve their small reference outputs outside version control. Inspect test registration and verification scripts before adding tests.

## Task 2: Implement the xc history update with tests

Create src/rt/tdcdft_lrc.f90 and register it in src/rt/CMakeLists.txt. First add a failing standalone Fortran test under testsuites/ for zero drive, constant drive with analytic quadratic Axc, and smooth drive with second-order timestep convergence. Use explicit arguments, no mutable global state. Verify failure before implementation and success afterward. The discrete sign and initial half step must follow Task 1, not an assumed current convention.

## Task 3: Input and state lifecycle

Modify src/io/salmon_global.f90, src/io/inputoutput.f90, src/common/structures.f90 and src/rt/initialization_rt.f90. Add opt-in mode and nonnegative alpha, defaults, normalization, broadcasting, logging and compatibility checks. Test omitted input, disabled mode, alpha=0, finite positive alpha, invalid values and unsupported execution modes. Check all allocation/deallocation and restart initialization paths. Update the corresponding SALMON-DOCS manual source in a separate checkout and report its changed files.

## Task 4: Couple effective potential consistently

Modify src/rt/time_evolution_step.f90 and src/rt/em_field.f90 only where needed. Use the same effective field in midpoint propagation and endpoint physical-current evaluation, refreshing nonlocal projectors at the appropriate time. Preserve classical Ac_tot output and response normalization. Add a small Si test proving disabled/alpha=0 reference equivalence and nonzero response for positive alpha, plus impulse-amplitude and timestep checks. Reject unimplemented propagators instead of silently omitting Axc.

## Task 5: Output and restart

Modify src/io/write.f90 and src/io/checkpoint_restart.f90, following existing root-only and filesystem helpers. Add separately identified xc output with units and symmetric optional restart data. Compare uninterrupted and split runs for current and xc history; test incompatible restart settings. Preserve legacy behavior for disabled mode. Run representative serial and MPI cases.

## Task 6: Si linear-response examples and analysis

Create samples/exercise_si_tdcdft/ with documented GS/RT inputs based on the existing Si samples, alpha variants, and an analysis script. Add synthetic-spectrum tests for peak extraction and numerical integration; test that incompatible grids/windows are rejected. Produce alpha=0 versus finite-alpha smoke results, then carry out k-grid, real-grid, timestep, duration and impulse-amplitude convergence before making quantitative claims. Report E1/E2 peak positions and intensities with consistent broadening and explicitly separate any energy-axis alignment.

## Task 7: Laser-excited response

Add pump-only and delayed weak-probe examples using ae_shape2='impulse', documenting the verified delay convention. Add a difference-current analysis tool and tests using synthetic pump cancellation and a known linear impulse response. Normalize to the actual probe field and shift time origin to the probe. Verify no-pump agreement with Task 6 and weak-probe convergence before an intensity/delay scan. Save windows, broadenings and pump parameters with every result. Do not interpret fixed-alpha results as density-dependent screening.

## Task 8: Review and completion

Check formatting and diff, build with gfortran, run existing and new relevant tests, and review effective-potential consistency and restart compatibility. Document tested options and untested accelerator/library paths. Commit approved changes on TDCDFT; do not push unless requested. Distinguish working implementation/smoke tests from converged scientific results in the final report.
