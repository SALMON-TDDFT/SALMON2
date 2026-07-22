# WPW Search-space Metric-mode Diagnostic Design

## Decision

Diagnose the LOBPCG search space `Z=[Q,Z_R,P]` at the existing selected
iterations.  Using the already assembled reduced metric `S_Z=Z^dagger S Z`
and a new overlap `B_R=Z^dagger R`, diagonalize `S_Z` with the same cutoff
convention as the reduced eigensolver and report:

- full minimum/maximum metric eigenvalue and effective retained rank;
- smallest retained eigenvalue divided by the maximum;
- fraction of `|V^dagger B_R|^2` lying in discarded metric modes, separately
  for occupied and extra residual columns;
- the corresponding fractions in the lowest retained decade above cutoff.

This diagnostic distinguishes search-space redundancy from a residual that
stalls in well-conditioned retained directions.  It does not change the
metric cutoff, retained subspace, reduced eigenproblem, convergence test,
preconditioner, or publication boundary.  Nonfinite or inconsistent
diagnostic inputs remain fatal.

The result concerns LOBPCG search-space conditioning.  It must not by itself
be described as proof that the physical occupied-W span is deficient.

## Interpretation

- Large discarded/near-cutoff residual fractions: search construction is
  feeding redundant metric directions; revise orthogonalization/restart.
- Small fractions with stalled residual: the retained Rayleigh-Ritz update or
  the underlying operator/basis conditioning is responsible.
- A well-conditioned no-preconditioned search space with continued stagnation
  motivates a direct/dense diagnostic solve before changing the Wannier basis.
