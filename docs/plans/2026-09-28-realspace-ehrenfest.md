# Real-space PBEh Ehrenfest implementation plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Connect DC-prepared PBEh occupied mesh wavefunctions to SALMON real-time Ehrenfest propagation.

**Architecture:** Existing conventional-from-DC reconstruction populates spsi%zwf without configuring lcfo_rt_basis. Native Taylor4/ACE propagates it; midpoint nuclear pseudopotentials and endpoint forces connect to existing Verlet. No additional basis or thermal update.

**Tech Stack:** Native Fortran real-space kernels, MPI, Python native tests.

Spec: 2026-09-28-dc-ehrenfest-design.md, user explicit mesh representation correction.

### Task 1: Native rejection test then connection

Create testsuites/653_functional/test_ehrenfest.py: generate a 2-fragment DC initial state, reconstruct into one native full-grid rank without LCFO projection, impulse and moving NVE ions. Observe existing PBEh input rejection first. Modify inputoutput.f90 to permit precisely this metadata-checked route, retaining checkpoint, finite-support, projected-LCFO MD, fractional RT occupation and finite-pulse guards. Modify hse_native.f90 to admit the corresponding ionic extension. Require per-step pseudopotential updates and fresh per-step energy. Build and rerun.

### Task 2: Geometry and force consistency

Modify time_evolution_step.f90 to evaluate ionic potentials/nonlocal projectors at midpoint positions for electron propagation, then rebuild at endpoint positions before force/energy. Reuse existing real-space Taylor4 and MD steps. Add trajectory/energy timestep refinement, occupation/charge and native nonlocal-projector tests. Preserve failure gates and old functionals' behavior. No projected LCFO MD changes.

### Task 3: Verification and review

Build HSE ON/OFF, run native and input regressions. Independent review checks state representation, midpoint potential/endpoint force geometry, initial forces, exchange cache validity, metadata and misleading giant-system claims. Record actual metrics and restrictions; commit only if tests validate the newly enabled path.
