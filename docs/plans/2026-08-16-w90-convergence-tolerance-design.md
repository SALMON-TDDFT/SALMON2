# Wannier90 Convergence Tolerance Design

## Context

The Si64 direct-retained-frame run reaches a stationary Wannier90 solution but
does not emit the explicit convergence message.  Its total spread is about
6947 Angstrom squared and its final five iteration changes are between
approximately 0.9e-12 and 4.5e-12.  The generated input currently requests an
absolute `conv_tol` of `1e-12` with a five-iteration window, which is at the
double-precision resolution of a spread of this magnitude.

## Decision

Generate `conv_tol = 1.d-10` and retain `conv_window = 5` and `num_iter = 200`.
SALMON will continue to require Wannier90's explicit convergence message and
will continue to reject iteration-limit exhaustion.  It will not infer
convergence from a final state or parse its own plateau criterion.

## Verification

Add a focused source-route assertion for the generated tolerance so the input
contract cannot silently regress.  Run the Wannier90 MPI fixture on 1, 2, 4,
and 8 ranks, the route checks, and a production build.  Then rerun Si64 and
confirm that Wannier90 emits its convergence message before iteration 200 and
that SALMON proceeds beyond the previous validation failure.

