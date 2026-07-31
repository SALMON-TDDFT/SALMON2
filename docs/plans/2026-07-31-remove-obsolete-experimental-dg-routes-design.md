# Obsolete Experimental DG Route Removal Design

## Goal

Reduce maintenance cost by retaining only:

1. the established non-experimental SALMON paths, including normal DC
   LCFO plus EigenExa; and
2. the accepted buffer-supported, symmetry-preserving
   overlapping-Wannier ground-state path followed by V3-checkpoint-backed
   generalized-eigenvalue exponential coefficient RT.

Remove the other experimental DG execution paths, their input surface,
dispatch code, implementation modules, samples, and tests. Do not retain
fallbacks, automatic promotion, or compatibility shims for removed DG
options.

## Retained Architecture

### Established SALMON paths

Normal ground state, conventional RT, normal DC, LCFO, EigenExa, and
Wannier90 integration remain unchanged. Their input compatibility,
dispatch order, outputs, and numerical behavior are outside the removal
scope.

### Overlapping-Wannier DG path

The retained DG path consists of:

1. conventional DC fragment candidate generation;
2. periodic buffer-box, symmetry-preserving complete-s+p Wannier
   construction;
3. row-owned overlap, Hamiltonian, position, velocity, density, and
   coefficient state;
4. atomic V3 route-checkpoint publication and reuse; and
5. fixed-basis coefficient RT solving
   \(iS\dot C=H(t)C\) by the generalized-Hermitian exponential.

This path must not enter direct real-space DG, LCFO diagonalization,
EigenExa, WPW, conventional checkpoint publication, DG-Fragment/Nodal RT,
or conventional orbital RT.

Code with a legacy DG name is not retained merely for compatibility. If
the accepted path requires it, move or rename it into the
overlapping-Wannier implementation. Otherwise remove it.

## Removal Scope

### Direct real-space DG ground state

Remove:

- `yn_dg_dc_local_periodic`;
- direct SIPG ground-state continuation;
- direct-DG checkpoint and handoff state;
- direct-DG ground-state solver and SALMON adapter; and
- route-specific diagnostics, samples, and tests.

### DG-Fragment and Nodal RT

Remove:

- `yn_dg_fragment_rt`, `yn_dg_nodal_rt`, and their supporting input
  options;
- LCFO-to-DG seed projection and validation;
- DG-Fragment coefficient, density, Hamiltonian, lifecycle, halo,
  observable, and integrator implementations;
- Nodal state, preparation, ground-state relaxation, Hamiltonian, halo,
  nonlocal, current, density, and propagation implementations; and
- their samples, runners, and focused tests.

### WPW ground state and RT

Remove:

- `yn_dg_wpw_*` production and checkpoint-RT options;
- WPW fixed-operator, complement, SCF, density, projector, checkpoint,
  halo, face, quadrature, sparse, and matrix-free implementations;
- WPW-specific LCFO adapters and export consumers; and
- their samples, runners, diagnostics, and focused tests.

### Obsolete experimental option families

Remove route-specific mixed-z, full-H seed, adaptive-DG-basis,
experimental length-gauge, projected-H, delta-H, and automatic
LCFO/WPW/DG promotion controls when they have no retained non-DG consumer.

Normal LCFO, EigenExa, and Wannier90 controls are retained. A Flux or
export facility is removed only when dependency analysis proves it serves
only a removed DG consumer.

## Dependency-Driven Removal

Removal proceeds from entry points inward:

1. Add source-contract RED tests that fail while removed flags, dispatches,
   and CMake entries remain.
2. Remove obsolete inputs and dispatch branches.
3. Remove route-specific samples, runners, and tests.
4. Remove CMake registrations and implementation modules that become
   unreferenced.
5. Move genuinely shared code required by overlapping-Wannier into its
   namespace.
6. Use symbol search, build dependency data, and link verification to
   prove that no obsolete route remains.

Removed inputs are not silently reinterpreted. They fail as unknown
namelist entries.

## Verification

Every removal task requires:

- a RED source-contract or build test;
- focused verification for the affected boundary;
- overlapping-Wannier MPI fixtures on 1, 2, 4, and 8 ranks;
- specification and code-quality review;
- resolution of all Critical and Important findings; and
- `git diff --check`.

The final gate requires:

- no removed DG flags, dispatches, sources, CMake entries, samples, tests,
  or stale documentation;
- a clean-first parent-prerequisite overlay build with MPI, ScaLAPACK, and
  EigenExa enabled;
- unchanged normal-DC Si64 LCFO plus EigenExa behavior;
- a fresh Si64 overlapping-Wannier ground-state checkpoint publication
  and accepted reuse;
- short coefficient RT and restart-split byte identity; and
- no marker indicating direct DG, WPW, DG-Fragment/Nodal RT, normal
  checkpoint publication, or conventional RT from the retained route.

## Failure Handling

No removed route may fall back to another path. The retained
overlapping-Wannier path continues to reject stale checkpoint
fingerprints, invalid metrics or operators, incomplete periodic tails,
and incompatible route settings collectively. Normal SALMON paths retain
their existing error behavior.

