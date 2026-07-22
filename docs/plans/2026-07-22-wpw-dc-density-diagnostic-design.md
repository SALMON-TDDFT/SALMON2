# WPW DC-density Diagnostic Design

## Context

The normalized B=6 occupied-W bootstrap passes the structural representation
checks: the source density is reproduced by the W+P projection to
`7.7960E-11`, the normalization residual is `1.3946E-15`, and the projected
charge error is `5.9380E-11`.  The route nevertheless stops before fixed-H
because the fragment-source density differs from the converged DC density by
`6.0445E-02`, and that approximation diagnostic is compared with the numerical
SCF residual tolerance `1E-8`.

The projected Wannier construction is an approximate representation.  The DC
density comparison measures its physical quality, but a finite difference
does not prove corrupt data or an invalid W+P projection.  It must therefore
not share the fatal path used for nonfinite data, invalid norm, rank loss,
projection failure, normalization failure, or charge failure.

## Decision

Keep all structural source-to-W+P qualifications fatal and unchanged.  Treat
the relative source-to-DC density residual as a nonblocking diagnostic:

- a finite residual with a positive finite DC-density denominator is valid;
- a valid residual above `dg_wpw_scf_residual_tolerance` emits
  `[DG-WPW-DC-DENSITY-WARNING]` and continues;
- an invalid/nonfinite numerator or denominator, or a nonpositive denominator,
  remains a collective failure;
- retain the configured tolerance and reported residual so buffer and basis
  convergence can be assessed from physical runs;
- do not weaken fixed-H convergence, frozen-state, checkpoint-publication, or
  RT-entry gates.

This separates approximation quality from representation integrity and lets
the actual fixed-H calculation determine whether the normalized Wannier basis
is adequate.

## Verification

Add a pure classifier test covering below/equal/above tolerance and invalid
inputs.  Add a source-contract assertion that the valid above-tolerance case
prints a warning rather than setting the seed failure code.  Run the focused
fixed-H and occupied-W MPI tests, build the MPI target, and repeat the B=6
physical run in a fresh directory.
