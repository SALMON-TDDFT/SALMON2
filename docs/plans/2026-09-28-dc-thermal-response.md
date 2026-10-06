# DC finite-temperature response implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Establish charge-conserving finite-temperature occupation and entropy derivatives for the approved DC free-energy force response.

**Architecture:** A standalone real-array Fortran module solves weighted Fermi occupations at fixed electron count and temperature. It exposes directional derivatives with respect to eigenvalues and core weights. PBEh DC uses the solver; the static diagnostic uses the same entropy definition. This is the occupation block of the coupled response, not a complete nuclear force.

**Tech Stack:** Fortran, Python unittest, MPI SALMON builds with HSE enabled and disabled.

**Spec:** `docs/plans/2026-09-27-dc-force-response.md`; user selected fixed finite temperature on 2026-09-28.

## Global constraints

Temperature is kBT in Hartree. Weights contain core norm and k-point weight; spin degeneracy is explicit. Preserve other functionals and zero-temperature paths. Never enable DC-MD based only on this block. Reject invalid input, insufficient capacity and unconverged chemical potential explicitly. Saturated occupations can satisfy charge numerically but must not yield a falsely certified response.

### Task 1: Weighted occupation and response kernel

Files: create `src/gs/dc/dc_thermal.f90`, `testsuites/653_functional/dc_thermal_probe.f90`, `testsuites/653_functional/test_dc_thermal.py`.

1. Write a standalone probe of weighted charge, entropy, analytic fixed-N directional derivatives versus independent central differences, common energy shifts, zero weights, degenerate states, extreme tails, insufficient capacity, invalid inputs and saturated response rejection.
2. Run `python3 -m unittest discover -s testsuites/653_functional -p test_dc_thermal.py -v`. Expected: missing production module failure.
3. Implement `solve_dc_thermal(e,w,T,g,N,mu,f,ts,status)` and `response_dc_thermal(e,w,T,g,mu,de,dw,dmu,df,dts,status)`. Solve a monotone bracket by bisection, compute stable Fermi tails, and enforce sum(g*(w*df+f*dw))=0. Use dTS=g*sum(T*s(f)*dw+(e-mu)*w*df).
4. Run the same test. Expected: PASS.

### Task 2: Native finite-temperature PBEh DC integration

Files: modify `src/gs/dc/CMakeLists.txt`, `src/gs/dc/dcdft.f90`, `src/gs/dc/dc_force.f90`, `testsuites/653_functional/test_dc_force.py`.

1. Add a native failure test for insufficient fragment-state capacity. Run against old binary. Expected: no explicit occupation failure.
2. Register module, call weighted solver for positive-T PBEh DC, preserve legacy paths, use shared entropy in force diagnostic. Status failures print a useful error on representative rank then stop collectively.
3. Build HSE ON MPI and HSE OFF serial. Expected: both succeed.
4. Run DC force tests including native water limit and MPI parity; unit thermal/projector and input guard regressions. Expected: PASS.

### Task 3: Review and evidence

Update accepted choice, equations, limits and verification evidence. Run one independent review of changes against `51d9b70c` and address important findings. Save local commit; no push or merge.

## Review focus

Capacity near full occupation, division by small charge susceptibility, spin/k/core weight factors, entropy weight derivative, replicated MPI failures, and distinction between occupation response and the still-missing coupled orbital/density response.
