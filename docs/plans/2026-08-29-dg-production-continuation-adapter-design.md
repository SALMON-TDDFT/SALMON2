# DG Production Continuation Adapter Design

## Purpose

The continuation controller already enforces the lambda transaction and
acceptance protocol, but its callback bundle has no explicit production
context.  The existing callbacks used by divided SCF solve independent
fragments and cannot represent the complete DG generalized eigenproblem.
Task 7b therefore needs a tested adapter boundary before it can be wired into
`main_dft`.

## State ownership

The callback bundle shall own one explicit polymorphic backend object.  The
controller invokes deferred backend procedures through that object instead of
calling context-free `nopass` procedure pointers.  No module-global singleton
is permitted: a singleton would mix independent calculations, make rollback
state implicit, and prevent deterministic unit tests.

The production backend owns references or immutable copies of:

- the frozen effective WF-block and PW-packet selections and their actions;
- the distributed WF+PW fragment basis and ownership;
- the fixed DG metric and complete coefficient-independent SIPG face blocks;
- the density-dependent volume-operator workspace;
- the current occupied coefficients, occupations, eigenvalues, projector,
  density, and interface trace;
- operator, basis, density, and state fingerprints and epochs.

## Callback sequence

For every inner iteration the existing controller calls the backend in this
order.

1. Rebuild the volume operator from the supplied density and combine it with
   the frozen complete interface operator multiplied by the single global
   lambda.
2. Solve the complete distributed generalized eigenproblem.
3. Rebuild the occupied S-metric projector from the solved coefficients.
4. Reconstruct density and all interface traces from the same occupied
   subspace.
5. Evaluate coefficient-space, real-space, density, interface, metric,
   electron-number, and optional occupied-subspace symmetry residuals.
6. Mix density only.  Interface traces are refreshed observables and are not
   stale Hamiltonian inputs for coefficient-independent SIPG faces.
7. Produce an acceptance receipt only from the completely refreshed state.

Rollback state remains in the continuation controller.  The backend must be a
pure function of the accepted density, frozen basis/operator payload, lambda,
and its current solve result; it must not retain an uncommitted density as an
authoritative state.

## Symmetry

WF closure is performed on fragment/WF blocks, not on individual eigenvectors.
The production fragment action supplies the block permutation.  A general
complex Wannier representation remains an internal gauge action within each
block and is not converted to an orbital permutation.  Identity-only systems
use one explicit identity block action.  The optional symmetry residual is an
occupied-projector residual; individual states are not required to transform
as one-dimensional irreducible representations.

## Supported first production scope

The first backend accepts the scope already enforced by
`build_dg_hybrid_scope_receipt`: periodic, spin-unpolarized, no spin orbit,
no DFT+U, no exact exchange/HSE, no fixed/history-dependent functional, and
supported density-only XC types.  Unsupported selectors fail before backend
construction and leave protected routes unchanged.

## Testing

The callback API refactor is first tested against the existing deterministic
continuation fixture.  A separate MPI production-adapter fixture then uses a
small distributed Hermitian generalized problem with nonzero cross-fragment
SIPG blocks.  It verifies operation order, one uniform lambda on all faces,
coefficient-independent face reuse, density-dependent volume refresh,
projector gauge invariance, fresh trace epochs, rollback reproducibility, and
fully refreshed lambda-one acceptance.  Existing continuation controller,
acceptance, SIPG, face-trace, divided-SCF, and protected-route tests remain
green.
