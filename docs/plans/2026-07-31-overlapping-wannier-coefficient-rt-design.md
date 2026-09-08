# Overlapping-Wannier Coefficient RT Design

## Scope

Propagate occupied coefficients in the accepted, fixed overlapping-
Wannier basis. The route loads only its own accepted checkpoint and never
falls through to conventional orbital RT, direct real-space DG, LCFO,
EigenExa, or WPW.

## Checkpoint V3

The ground-state publication is extended to store row-owned matrices:

- metric `S`;
- field-free Hamiltonian `H0`;
- three wrapped-position matrices `X`;
- three velocity matrices `V`;
- nonlocal velocity components where required by the validated observable
  convention.

The manifest retains coefficients, occupations, basis generation,
geometry generation, basis fingerprint, and immutable operator
fingerprint. It adds matrix fingerprints and an explicit gauge convention.
RT rejects V2, incomplete row ownership, mismatched generations, stale
operator fingerprints, or failed GS acceptance gates.

An RT restart adds time, step, coefficients, field convention, and the
immutable V3 basis/operator provenance. It never republishes a GS
checkpoint.

## Propagation

For the midpoint Hamiltonian,

```text
Hmid = H0 + sum_axis E_mid(axis) X(axis)
```

in length gauge, or the separately validated velocity-gauge expression,
solve

```text
Hmid U = S U epsilon,  U^dagger S U = I
```

and propagate

```text
C_next = U exp(-i epsilon dt) U^dagger S C_current.
```

This is exact for a time-independent Hamiltonian and preserves the
S metric. A field-free eigensystem may be cached and reused. A
time-dependent field uses the midpoint Hamiltonian and is diagonalized
for that step. Crank--Nicolson is not an automatic fallback.

The initial implementation uses the existing dense LAPACK generalized
Hermitian eigensolver. Production assembly respects the established row
ownership and must not persist an additional full matrix on every rank.

## Gates

Focused tests require:

- `C†SC` conservation;
- field-free energy conservation;
- covariance under coefficient phase and retained-basis unitary rotation;
- length-gauge position and velocity-matrix coupling;
- one-shot versus restart-split identity;
- collective rejection of stale basis, operator mismatch, incomplete
  rows, non-Hermitian operators, and tail-generation mismatch;
- identical fingerprints and numerical results on 1, 2, 4, and 8 ranks.

The basis, metric, and observable matrices remain fixed. V3 validates that
every retained tail contains every physical periodic-grid ID at least once;
overlapping periodic buffers may contain duplicate IDs.
Consequently a fixed-basis coefficient combination cannot escape the
periodic support: a representational tail escape is structurally zero,
while any request to update/adapt the basis terminates this route
explicitly.
