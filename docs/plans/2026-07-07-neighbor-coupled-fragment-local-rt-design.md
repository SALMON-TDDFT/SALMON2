# Neighbor-Coupled Fragment-Local RT Design

## Goal

Make the fragment-local DG real-time propagation physically comparable to the global DG reference by including the neighbor and cross-fragment couplings that restore the assembled Hamiltonian symmetry.

## Motivation

Recent C64 local-Wannier laser tests showed that forcing fragment-local Wannier centers into an inversion-closed set strongly reduced the polarization offset and even harmonics, but did not remove the nonphysical even-order response. Direct diagnostics then showed that, in the projected local seed basis, the fragment-local overlap and Hamiltonian blocks are not individually inversion symmetric.

That result does not by itself mean the DG/DC Hamiltonian is broken. In a DG decomposition, a local block can be asymmetric while the full assembled operator is symmetric after neighbor and cross-fragment terms are included. The current fragment-local time evolution is therefore likely cutting the compensation terms that the global route keeps.

## Chosen Approach

Do not enforce artificial symmetry on each isolated fragment block. Instead, propagate each owner fragment together with its face-neighbor environment:

```text
C_env = [C_owner, C_neighbor]

H_env =
  [ H_owner,owner    H_owner,neighbor ]
  [ H_neighbor,owner H_neighbor,neighbor ]

R_env =
  [ R_owner,owner    R_owner,neighbor ]
  [ R_neighbor,owner R_neighbor,neighbor ]

H_eff(t) = H_env - E(t) R_env
```

For the first correctness route, include the six face neighbors only. Use a dense diagonalization/exponential backend for this environment block, because the current priority is physics validation rather than weak-scaling performance.

After propagating the environment block, write back only the owner-fragment component. Neighbor coefficients are input data for the local environment propagation, not owner data to be double-counted.

## Data Flow

At every RT step:

1. Synchronize or gather the current coefficients for the owner and face-neighbor fragments.
2. Build an environment basis with stable global IDs for owner and neighbor slots.
3. Assemble `S_env`, `H_env`, and the field-coupling position block `R_env`.
4. Orthonormalize or generalized-diagonalize the environment block consistently with the existing DG basis metric.
5. Apply the dense exponential propagation for the current time step.
6. Write back the propagated owner component only.
7. Recompute observables from the globally owned coefficients.

The global mixed-Z backend remains the reference. The new route is an accuracy/debug backend until the physical diagnostics pass.

## Required Audits

Before implementing the field-on route, audit whether the cross position blocks `R_owner,neighbor` and `R_neighbor,owner` are already exported or reconstructable from existing BPW/Wannier data. If they are not available, add an explicit export/reconstruction step rather than silently using only owner-local position blocks.

Also audit whether the existing face-neighbor enumeration is consistent between the Hamiltonian, BPW window, and coefficient-storage paths. The first implementation should reuse the existing face-neighbor helper rather than introduce a second neighbor definition.

## Diagnostics

The new backend should print or expose:

- owner and environment dimensions
- face-neighbor fragment IDs
- maximum Hermiticity error for `S_env`, `H_env`, and `R_env`
- condition number or smallest eigenvalue of the environment metric
- coefficient difference against the global reference in diagnostic mode
- field-off polarization/norm drift
- laser `Pz` odd-response check from `+E` and `-E`
- HHG ratios including `H2/H1`, `H3/H1`, and `H4/H1`

## Acceptance Criteria

The neighbor-coupled local backend is usable when:

- field-off propagation does not increase drift relative to the current fragment-local route
- short laser `Pz` follows the global reference phase and amplitude better than isolated fragment-local propagation
- even harmonics become dips for centrosymmetric C64/Si64 tests, within the same analysis used for the global reference
- `+E` and `-E` runs satisfy odd-response symmetry in the time-domain polarization
- the route remains controlled by namelist options, not environment-variable-only switches

## Out Of Scope

- Removing the existing global reference backend
- Optimizing communication or memory use
- Making the dense environment propagation the final weak-scaling implementation
- Requiring each isolated fragment Hamiltonian block to be inversion symmetric
