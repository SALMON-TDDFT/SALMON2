# Wannier90 Seed–DMN Gauge Consistency Design

## Problem

SALMON constructs a localized spatial trial frame and passes its overlap with the retained orthonormal frame as the Wannier90 `A` matrix.  The current restored DMN path instead publishes the retained fixed-center representation with an identity AMN.  The two inputs therefore describe different initial gauges.  With the corrected Bohr-to-Angstrom geometry this mismatch leaves the Si64 localization trajectory oscillating rather than converging.

## Design

Compute the Gamma-point raw overlap `A=<retained|trial>` once, collectively, before Wannier90 setup.  Because a localized trial frame is not generally orthonormal, replace the raw overlap by its closest unitary polar factor `Q=U V^H` from `A=U Sigma V^H`.  Use that exact `Q` in all three places:

1. DMN `u_matrix`/AMN input.
2. DMN target representation `D_trial=Q^H D_retained Q`.
3. Wannier90 library `A` input.

The retained representation remains the already validated fixed-center representation.  The physical lattice conversion, symmetry catalog, Wannier90 optimizer, complex-gauge canonicalization, and post-Wannier validation remain unchanged.

To avoid repeating the expensive spatial overlap, the existing M/A assembler accepts the validated precomputed A matrix and skips its A accumulation.  Its M assembly and memory bounds remain unchanged.

## Contracts

- Dimensions, weights, finite payloads, and coordinator allocation are collectively validated.
- The raw overlap must have full numerical rank; otherwise the seed cannot define the complete retained gauge and is rejected.
- Q is explicitly checked unitary and is present only on the coordinator, matching the existing Wannier90 matrix ownership.
- The precomputed A fingerprint is bound to the existing M/A input fingerprint.
- The production route checks that DMN and Wannier90 consume the same A payload.
- Failures reject collectively before shape-dependent communication.

## Verification

- Focused MPI 1/2/4/8 test of distributed overlap assembly followed by the unitary-polar characterization `Q^H Q=I`, with `Q^H A` Hermitian positive.
- Test that precomputed-A reuse produces the same full M/A result as direct assembly.
- Route test requiring `Q^H D Q` for DMN and post-Wannier covariance.
- Existing complex-gauge, symmetry, build, and diff checks.
- Si64 run with correct geometry; acceptance requires Wannier90 convergence before downstream symmetry and Hybrid-SCF gates.
