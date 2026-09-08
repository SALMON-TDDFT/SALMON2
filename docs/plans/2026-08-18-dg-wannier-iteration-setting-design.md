# DG Wannier Iteration Setting Design

## Problem

The DG overlapping-Wannier production route writes `num_iter = 400` into the
Wannier90 input and validates the resulting log against the same literal.  This
bypasses SALMON's existing `wannier_num_iter` input.  The Si64 calculation was
still reducing its spread at iteration 400, so Wannier90 completed normally but
SALMON correctly rejected the unconverged result.

## Design

Use the existing `wannier_num_iter` setting as the single source of truth.  Pass
it from `main_dft` into the Gamma-library setup routine and the Gamma-library run
routine.  The setup routine writes that value to the `.win` file; the run routine
uses the identical value when validating convergence.  Reject non-positive
iteration limits through the existing collective contracts.

Do not add a second DG-specific setting, change Wannier90's convergence
algorithm, weaken `conv_tol`, or replace the limit with another literal.  The
Si64 verification input will explicitly request 1000 iterations.

## Verification

The route checker must fail if either production call omits
`wannier_num_iter`, if a fixed `num_iter = 400` remains, or if convergence
validation uses a literal limit.  Then rebuild SALMON and rerun the focused route
checks before starting a new Si64 run.
