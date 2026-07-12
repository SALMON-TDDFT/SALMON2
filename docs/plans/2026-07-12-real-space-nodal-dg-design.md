# Real-Space Nodal DG Time Propagation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a matrix-free real-space DG prototype that propagates core grid values with Taylor expansion while treating fragment buffers as neighbor halos.

**Architecture:** Each fragment owns only its core real-space degrees of freedom. Six face descriptors retain `(axis, side)` identity even when two periodic faces map to the same neighboring fragment. Every Taylor Hamiltonian action refreshes the required halo values, applies the velocity-gauge finite-difference Hamiltonian on the core, and later adds DG numerical-flux and nonlocal-pseudopotential contributions.

**Tech Stack:** Fortran 2008, SALMON real-space stencils, MPI fragment communicators, Python structural regression tests, CMake.

---

### Task 1: Nodal storage and face topology

**Files:**
- Create: `src/rt/dg/rt_dg_nodal_types.f90`
- Modify: `src/rt/CMakeLists.txt`
- Test: `tests/dg/check_nodal_real_space_types.py`

1. Write a failing structural test requiring core-only wavefunction storage and six face descriptors with explicit axis and side.
2. Run the test and confirm it fails because the module is absent.
3. Add minimal nodal types without wiring them into the legacy basis route.
4. Run the structural test and compile SALMON.

### Task 2: Core/halo indexing and periodic face identity

**Files:**
- Create: `src/rt/dg/rt_dg_nodal_halo.f90`
- Test: `tests/dg/check_nodal_face_halo.py`

1. Write a failing test for six distinct faces and periodic `+/-` faces mapping to the same neighbor without being collapsed.
2. Implement face initialization and local pack/unpack bounds.
3. Add rank-local bounds and duplicate-face diagnostics before MPI communication.
4. Run the test and compile.

### Task 3: Matrix-free field-off Hamiltonian action

**Files:**
- Create: `src/rt/dg/rt_dg_nodal_hamiltonian.f90`
- Test: `tests/dg/check_nodal_hamiltonian_action.py`

1. Write a failing numerical test comparing a decomposed periodic finite-difference Laplacian with the undecomposed reference.
2. Implement core stencil application using refreshed face halos.
3. Add local potential multiplication; leave nonlocal pseudopotential as an explicit hard-stop until implemented.
4. Verify Hermiticity and decomposition independence.

### Task 4: Taylor propagation and input guards

**Files:**
- Create: `src/rt/dg/rt_dg_nodal_taylor.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `src/rt/dg/rt_dg_iteration.f90`
- Test: `tests/dg/check_nodal_taylor_route.py`

1. Add a default-off namelist switch for the nodal real-space route.
2. Require Taylor propagation and `dt <= 0.02 au` when enabled.
3. Refresh halos for every Taylor Hamiltonian action, not only once per RT step.
4. Test norm and field-off stationarity on a small periodic system.

### Task 5: Velocity gauge and physical validation

**Files:**
- Modify: `src/rt/dg/rt_dg_nodal_hamiltonian.f90`
- Test: `tests/dg/check_nodal_velocity_gauge.py`

1. Add the volume operator `-i A.grad + A^2/2` using the existing SALMON stencil.
2. Add the covariant normal derivative to the DG face flux.
3. Add the velocity-gauge nonlocal pseudopotential action.
4. Compare field-off drift and laser current with conventional Full TDDFT at `dt=0.02 au`.

### Task 6: Controlled O(N) localization

1. Establish the all-occupied-orbital nodal route as the correctness reference.
2. Add fragment-local orbital support radii and truncation diagnostics behind a separate switch.
3. Verify error versus buffer/support radius before claiming weak scaling or O(N) behavior.
