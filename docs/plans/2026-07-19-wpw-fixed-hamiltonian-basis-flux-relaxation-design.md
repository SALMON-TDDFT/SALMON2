# WPW Fixed-Hamiltonian Basis/Flux Relaxation Design

## Purpose

Construct the first Si64 WPW checkpoint without conflating three effects:

1. changing from the converged Full-TDDFT/DC representation to the W+P basis;
2. adding the DG interface bilinear form across neighboring fragments;
3. changing the physical density, Hartree potential, or XC potential.

The first checkpoint is a representation and interface-relaxation checkpoint, not a new self-consistent electronic structure calculation.

## Authoritative reference state

The converged DC+LCFO calculation supplies the frozen physical state:

- total density and electron count;
- ionic/local potential;
- Hartree and XC potentials;
- nonlocal pseudopotential data;
- all non-flux Hamiltonian components;
- occupations and spin convention;
- fragment-local LCFO basis and its provenance.

These quantities remain immutable throughout fixed-Hamiltonian relaxation. The implementation must snapshot their checksums/norms before relaxation and verify them after every accepted orbital update.

## Operator decomposition

The WPW generalized eigenproblem is written as

\[
  (H_0^{\mathrm{DC}} + \lambda_{\mathrm{interface}} H_{\mathrm{interface}}) C = S C E.
\]

`H0_DC` contains every volume and pseudopotential contribution evaluated from the converged reference density. `H_interface` contains the complete DG face bilinear form on self and cross-fragment traces.

The current stored representation does not support a safe independent `H_neighbor`/`H_flux` split: neighbor transport is an ownership/halo mechanism, while `ww_face_self`, `ww_cross_value`, and `wp_h_face` already contain assembled DG interface terms. Therefore the first implementation uses one continuation parameter, `lambda_interface`. A later split into consistency, symmetry, and penalty lambdas is permitted only after those raw terms are published separately and a dense oracle proves that their sum reproduces the existing face operator exactly.

The initial component map is:

| Frozen `H0_DC` | Continued `H_interface` | Transport only, never scaled |
|---|---|---|
| `ww_kinetic`, `ww_potential`, projector `ww_nonlocal` | `ww_face_self`, `ww_cross_value` | support IDs, halo records, owner schedules |
| `wp_h_volume`, `wp_h_nonlocal` | `wp_h_face` | candidate routing and canonical trace delivery |
| `pp_h_volume`, `pp_h_nonlocal` | none | P-owner redistribution |

The decomposition must be explicit in provenance and diagnostics. A potential rebuild from a reconstructed WPW density is forbidden in this stage.

## Occupied-subspace handoff

The relaxation must not begin from a random deterministic seed. Export the converged DC+LCFO occupied subspace, project it into the W+P metric, and solve the rank-revealing S-orthonormalization problem. The projected occupied projector is the authoritative initial state. Extra retained states may use deterministic complement vectors only after projection against the occupied subspace.

Record separately:

- projection loss before relaxation;
- occupied-projector defect after S-orthonormalization;
- reconstructed-density error of the projected, unrelaxed state;
- the additional change caused by fixed-H variational relaxation.

Coefficient-by-coefficient comparison is forbidden because gauge rotations and degenerate-state mixing are physically irrelevant. Compare the S-metric occupied projector, density, charge, and invariant subspace traces.

## Relaxation state machine

### Stage 0: Freeze and validate

Require converged, finite DC energy and charge. Snapshot the reference total-grid density/potential fields and the projected WW/WP/PP volume/nonlocal blocks. Record fingerprints for the input, atom coordinates, pseudopotential, grid, fragmentation, basis/window selection, metric, occupations, and operator decomposition.

### Stage 1: Zero-interface representation check

Build the W+P overlap and frozen non-interface Hamiltonian with `lambda_interface=0`. Project the converged occupied subspace into W+P before any minimization. First qualify the raw projection, then solve the fixed generalized eigenproblem from that projected state. Reconstruct density only for comparison; do not feed it back into H.

This stage measures basis/window truncation error independently of DG interfaces.

### Stage 2: DG-interface continuation

Introduce the complete bounded DG interface operator using `lambda_interface` from 0 to 1. At each accepted continuation point, solve the fixed generalized eigenproblem using the previous occupied subspace as the initial guess. Reject a step and reduce its size if the subspace merit, S-orthogonality, charge, Hermitian defect, or frozen-field invariants exceed their acceptance bounds. Only the mapped interface blocks may change.

### Stage 3: Checkpoint qualification

At `lambda_interface=1`, require a converged occupied invariant subspace and compare the reconstructed WPW density and frozen-H0 energy components against the Full-TDDFT/DC reference. Publish the checkpoint only after all invariants and provenance checks pass.

## Orbital solver and cooling

Use nonlinear conjugate-gradient/subspace iteration on the fixed generalized Rayleigh quotient. Density mixing is disabled because the density does not update H in this stage.

The orbital step has an explicit `orbital_mix`. Acceptance is defined for the occupied invariant subspace, not for individual eigenvectors. Its merit combines the occupied trace Rayleigh quotient, maximum occupied generalized residual, S-orthogonality, and occupied-projector change. Degenerate clusters are locked and rotated as subspaces rather than state by state.

The cooling rules are:

- begin conservatively;
- reduce it after rejected steps or residual growth;
- allow limited growth only after consecutive accepted steps;
- restart the CG direction when it ceases to be a descent direction or the S-metric rank changes;
- reorthonormalize in S after every trial step.

Continuation cooling and orbital cooling are separate controls. The previous accepted eigenspace is always the restart point after rejection.

## Invariants and failure policy

Every accepted step must preserve:

- frozen total-grid density and potential fingerprints;
- finite H, S, eigenvalues, and coefficients;
- exact electron count within tolerance;
- bounded S-orthogonality error;
- consistent operator epoch and ownership fingerprint;
- Hermitian WW/WP/PP block conventions;
- no mutation of volume/nonlocal blocks when only `lambda_interface` changes;
- unchanged occupations and occupied-state count;
- unitary/gauge-invariant subspace comparison across degenerate clusters.

Any invariant failure rejects the trial transactionally. Exhausted cooling or continuation limits fail closed and do not publish a checkpoint.

## Diagnostics and acceptance criteria

Report separately for raw projection, zero-interface relaxation, and interface continuation:

- interface lambda and accepted/rejected step counts;
- orbital mix and CG restart count;
- generalized residual, S-orthogonality, and subspace/projector change;
- eigenvalue min/max, occupied gap, locked clusters, and occupied trace Rayleigh quotient;
- reconstructed charge and density difference from the frozen reference;
- decomposed Rayleigh contributions from H0 and the complete DG interface operator;
- Hermitian defects and block norms;
- proof that frozen density/potential fingerprints did not change.

Two tolerance profiles are mandatory: a diagnostic preflight profile and a stricter publication profile. Both must be recorded in the manifest before a run. Density/projector publication limits are referenced to the measured zero-interface projection error: interface relaxation may not exceed a configured multiple of that baseline plus an explicit absolute floor. Charge, residual, S-orthogonality, Hermitian defect, and frozen-fingerprint tolerances remain independent absolute controls.

The first production checkpoint is accepted only when the fixed full-interface operator is converged and the representation errors satisfy the publication profile. It does not claim a new WPW self-consistent ground state.

## Deferred work

Density/potential self-consistency is a later, separately enabled stage. If required, it begins from the qualified fixed-H checkpoint and uses slow density mixing. Field-off RT and laser smoke remain blocked until the fixed-H checkpoint passes.
