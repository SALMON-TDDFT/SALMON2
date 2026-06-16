# DG Flux Real-Space RT Design

## Goal

Use DC-LCFO-Flux eigenstates as the initial occupied real-space wavefunctions, then propagate those wavefunctions with the Flux-inclusive real-space Hamiltonian.

This changes the main RT path away from DG coefficient propagation. The Taylor expansion remains a Taylor approximation to `exp(-i H dt)` with a fixed Hamiltonian during each time step. Flux is part of the Hamiltonian application, not an update performed inside the Taylor order loop.

## Motivation

The current coefficient-space DG-RT path is efficient, but it exposes two physical and numerical risks:

- The RT basis truncation can underrepresent the paramagnetic response and f-sum strength, especially for diamond where absorption extends to high energy.
- The Hamiltonian used for propagation must match the DC-LCFO-Flux eigenstates. Recomputing or changing Flux during Taylor orders breaks the intended exponential propagation.

The new path keeps the physically important Flux Hamiltonian, but removes the coefficient-space truncation from the RT dynamics by reconstructing occupied wavefunctions in real space.

## Recommended Architecture

Add a new DG-Flux real-space RT mode:

1. Run DC with `yn_dc_lcfo='y'`, `yn_dc_lcfo_flux='y'`, and `yn_dc_lcfo_diag='y'`.
2. Read the Flux-inclusive LCFO eigenvectors from `data_dcdft/fragments/*/wavefunctions.bin`.
3. Reconstruct only occupied real-space orbitals from the fragment basis data.
4. Store the result in the usual RT orbital container, preferably `spsi%zwf`, because the impulse and time evolution are complex.
5. Use the standard real-space RT time evolution machinery.
6. Add the DG Flux boundary correction to every real-space Hamiltonian application used by RT.

This preserves the existing RT density, potential, observable, current, and Taylor infrastructure as much as possible.

## Existing Code Anchors

- `src/rt/initialization_rt.f90`
  - The DG-from-DC path currently skips conventional real-space wavefunction reconstruction.
  - The conventional path already calls `init_conventional_from_dcdft`.

- `src/gs/dc/lcfo.f90`
  - `init_conventional_from_dcdft` already reconstructs real-space wavefunctions from DC-LCFO files.
  - This routine should be extended or split so the Flux RT path can reconstruct occupied states from the Flux seed without enabling coefficient-space DG-RT propagation.

- `src/gs/conjugate_gradient.f90`
  - `apply_dc_lcfo_flux_hpsi_rwf` already applies the Flux correction to a real-valued `hpsi`.
  - RT needs the same operator semantics for complex `zwf`.

- `src/rt/taylor.f90` and `src/rt/time_evolution_step.f90`
  - These should remain the high-level propagation path.
  - The Flux correction should enter through the Hamiltonian application, so every Taylor order sees the same operator for that time step.

## Data Flow

### Initialization

The new mode should enter after RT allocation and before density/potential initialization.

For each fragment owned by the current process group:

1. Read fragment metadata and coefficient matrix from `wavefunctions.bin`.
2. Read grid mapping from `rgrid_index.bin`.
3. Read the same basis representation used by Flux DC, including buffer basis data when needed.
4. Select occupied columns only.
5. Build real-space orbitals:

   ```text
   psi_occ(r, io) = sum_lambda basis_lambda(r) * coef(lambda, io)
   ```

6. Write into `spsi_in%zwf`.

Only occupied orbitals should be propagated. The virtual manifold remains represented implicitly by the real-space grid Hamiltonian, which is the key difference from coefficient-space DG-RT.

### Hamiltonian Application

The ordinary real-space `hpsi` call should compute the local kinetic, local potential, nonlocal pseudopotential, and field terms. Then the DG Flux correction should add the boundary coupling contribution:

```text
H_flux psi = H_realspace psi + boundary_flux_terms(psi_neighbor_faces)
```

The correction must be applied once per Hamiltonian application. It must not update the Flux Hamiltonian during individual Taylor orders.

The existing `apply_dc_lcfo_flux_hpsi_rwf` should be refactored into a reusable DG Flux Hamiltonian helper and extended to complex `zwf`.

## Mode Selection

Introduce an explicit RT mode rather than overloading the old coefficient-space DG flag behavior.

Recommended input switch:

```text
yn_dg_flux_realspace_rt = 'y'
```

When enabled:

- Require `yn_conventional_from_dcdft='y'` or an equivalent DC-from-file initialization path.
- Require DC Flux seed files from a GS run with `yn_dc_lcfo_flux='y'` and `yn_dc_lcfo_diag='y'`.
- Disable coefficient-space DG-RT propagation for that run.
- Stop with a clear error if spin-orbit, nonorthogonal cell, optimized fragment geometry, or unsupported k-point settings are requested before those cases are implemented.

## Expected Physics

For the diamond impulse test:

- `J_tot,z` starts from the diamagnetic kick response.
- The paramagnetic current cancels the initial diamagnetic current.
- `J_tot,z` then oscillates around zero with damping.
- A persistent oscillation around `J_dia` indicates insufficient paramagnetic response or an inconsistent Hamiltonian/initial-state pairing.

The first target is qualitative correctness of this impulse response. Drift and long-time conservation diagnostics can be added after the zero-centered response is recovered.

## Validation Plan

Start with the smallest reproducible diamond case already used in the DG-RT tests.

1. Confirm DC-LCFO-Flux converges and writes Flux-inclusive eigenvectors.
2. Reconstruct occupied real-space wavefunctions.
3. Verify norm and orthogonality after reconstruction.
4. Apply one Hamiltonian operation and compare real-space Flux expectation values against the DC eigenvalues for occupied states.
5. Run a short impulse RT.
6. Check that `J_tot,z` drops from the initial kick toward zero instead of oscillating around `J_dia`.
7. Run a longer impulse RT only after the short test is qualitatively correct.

## Risks

- The current Flux helper only supports `rwf`; RT needs a complex `zwf` implementation.
- The current helper lives in GS conjugate-gradient code. It should be moved or wrapped so RT can use it without coupling RT to CG internals.
- The conventional real-space Hamiltonian may already apply periodic stencil coupling. The Flux mode must avoid double-counting or mixing full periodic boundary coupling with DG Flux boundary coupling.
- The initial reconstruction must use the same Flux seed basis and coefficients as the DC diagonalization; otherwise the initial state and propagation Hamiltonian will not match.
- Full support for spin-orbit and nonorthogonal cells should be deferred until the scalar orthogonal path is working.
