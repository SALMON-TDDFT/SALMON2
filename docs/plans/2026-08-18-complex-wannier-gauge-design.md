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

Both the immediate result validator and the production transform applicator
must accept complex unitary matrices.  The applicator fixes only each output
column's arbitrary scalar phase: if the deterministic maximum-amplitude pivot
is `z`, it multiplies the transform column, spatial values, and gradients by
`conjg(z)/abs(z)`.  This makes the pivot positive real without changing the
represented orbital, its observables, or any internal unitary mixing.

The change does not weaken finite-value, unitarity, spread, character, cocycle,
or post-gauge checks.  Character-sector sewing may later construct a real final
orbit basis where required by the real generalized eigensolver; that is a
separate, symmetry-aware operation rather than an entry condition on the raw
Wannier90 matrix.
