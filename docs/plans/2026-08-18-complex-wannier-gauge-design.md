# Complex Wannier Gauge Acceptance Design

## Decision

SALMON will accept a finite complex Wannier90 transformation when it is unitary,
does not increase the gauge-dependent spread beyond tolerance, and passes the
existing downstream symmetry and operator checks.  The transformation need not
be element-wise real at the Gamma point.

Realness is a gauge choice, not an observable.  Initial-state equivalence is
established by preserving the represented subspace and consistently transforming
the Hamiltonian, transition operators, values, and gradients.  Character-sector
pairing and point-cogroup validation remain responsible for symmetry.

The change removes only the early `max(abs(aimag(transform)))` rejection.  It
does not weaken finite-value, unitarity, spread, character, cocycle, or
post-gauge checks.
