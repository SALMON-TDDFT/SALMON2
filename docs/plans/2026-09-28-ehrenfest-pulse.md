# Real-space Ehrenfest finite pulse implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a bounded finite-duration laser to the existing DC-initialized real-space PBEh Ehrenfest route and verify energy against external work.

**Architecture:** Keep real-space Taylor4/MLWF/ACE and midpoint ionic geometry. At the endpoint, evaluate kinetic/nonlocal energy, electronic current, force and completed ionic velocity consistently. For the imposed transverse Acos2 field, use centered endpoint E=-dA/dt for ionic forcing and diagnostic current/work. Occupations are not thermalized; LCFO only prepares the initial mesh state.

**Tech Stack:** Native Fortran mesh TDDFT, MPI, Python executable regression tests.

Spec: user-approved real-space Ehrenfest continuation from 426259e7. Full EXX support, occupied-only, NVE and k-only native exchange limitations remain.

### Task 1: Pulse acceptance and work oracle

Extend testsuites/unit_pbeh_rvv10/test_ehrenfest.py with Acos2 inputs using theory=tddft_pulse, a finite pulse followed by field-free propagation, dt/dt2/dt4. Independently integrate volume*(Jion-Jmatter).E with endpoint trapezoidal quadrature; compare to Eall+Tion and refine trajectory/current. Observe current input rejection first. Preserve an unsupported pulse-shape rejection case.

### Task 2: Correct time locations

Modify src/io/inputoutput.f90 and src/xc/hse_native.f90 to permit only the bounded pulse route with positive finite dt, omega and width, nonnegative start, transverse linear polarization and no second pulse/checkpoints. Modify src/rt/time_evolution_step.f90 to refresh endpoint A before energy evaluation and output ion current after completing the velocity step. Modify src/rt/initialization_rt.f90 and src/rt/em_field.f90 for centered endpoint fields on this new pulse route and initial direct ionic field force. Other functionals/projected LCFO behavior remains unchanged. Run pulse tests and inspect convergence; no tolerance relaxation to hide timing errors.

### Task 3: Verification and evidence

Add targeted endpoint checks (A/E timing, ionic current against output velocities, zero-amplitude behavior, invalid pulse controls). Run impulse/water and input regressions, HSE ON/OFF builds, one independent review. Save measured work/energy errors and valid sample input. Do not claim finite-pulse correctness from an acceptance-only test, or giant-system scalability from tiny native runs. Commit locally after verification.
