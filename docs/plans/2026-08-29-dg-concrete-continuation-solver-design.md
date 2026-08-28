# Concrete WF+PW DG Continuation Solver Design

## Decision

Implement the production continuation as one concrete solver, not as a generic
callback controller wrapped by a production adapter.  The solver owns one
explicit state transition and uses the existing SALMON numerical routines
directly.  This keeps every MPI synchronization point, physical state update,
rollback, and acceptance decision in one auditable path.

The earlier production-adapter design is superseded.  It remains in history as
design evidence but is not an implementation target.

## Three data objects

The production route has exactly three conceptual objects.

### Immutable catalog

The catalog is frozen before lambda-zero iteration and contains:

- effective, symmetry-closed WF blocks and PW packets, including the actions
  and the requested-to-effective selection record;
- distributed basis row IDs, exactly-one ownership, basis values, and grid
  distribution;
- the fixed DG metric and its independent sparse graph;
- the canonical physical-face list, face ownership, quadrature, normals,
  penalty data, and the complete coefficient-independent SIPG interface
  blocks;
- supported-scope, symmetry-analysis, cutoff, selection, and basis provenance;
- decomposition-independent fingerprints computed from canonical global IDs
  and values after duplicate and missing-ID checks.

There is one scalar lambda for the complete catalog.  A face-local or
fragment-local lambda is not representable by the solver interface.

### Mutable fixed-point state

One state value contains all quantities belonging to the same iterate:

- input and output density;
- density-dependent volume operator and complete
  `H_volume + lambda H_interface`;
- occupied coefficients, occupations, and eigenvalues;
- the basis-space occupied map
  `Q_occ = C_occ C_occ^dagger S` and the occupation density matrix
  `Gamma_occ = C f C^dagger`;
- gauge-invariant interface traces reconstructed from `Gamma_occ`;
- coefficient-space and real-space residuals, electron number, Hermiticity,
  metric, projector, and optional physical-symmetry residuals;
- lambda, iteration counters, and consistent epochs/fingerprints.

The state never uses `C^dagger S C` as the occupied projector.  Raw
eigenvectors are not mixed.  Interface traces are refreshed observables, not
independently mixed Hamiltonian boundary data, because the complete SIPG
operator is fixed for the catalog.

### Accepted checkpoint in memory

The solver keeps one deep copy of the last accepted mutable state.  A rejected
trial restores this complete copy transactionally.  No rejected density,
trace, operator, occupation, epoch, or mixing history may survive.  This
in-memory checkpoint is distinct from the final GS-to-RT file.

## Concrete solver loop

The initial input density is exactly the converged DC total density and is
fingerprint-checked before any WF+PW density reconstruction.  Lambda zero is
then converged as a finite-basis projected volume fixed point; equality to the
DC seed is provenance, not proof of lambda-zero convergence.

For each trial lambda and inner iteration, all ranks execute the same ordered
collective path:

1. establish collective validity of the input state;
2. build the volume operator from the current input density;
3. establish collective success before entering the next collective kernel;
4. add the same lambda times every frozen SIPG face block;
5. solve the distributed generalized eigenproblem;
6. build `Q_occ` and `Gamma_occ` from that solve;
7. reconstruct output density and every interface trace from `Gamma_occ`;
8. evaluate the independent residual channels from this one state;
9. run expensive real-space and symmetry gates only for a candidate that has
   passed the inexpensive gates;
10. accept the stage, or mix density only and continue.

Every rank-local numerical failure is converted to a communicator-wide result
before any later collective is entered.  The concrete solver therefore owns
the collective schedule; numerical helpers must not hide unmatched
collectives behind a rank-local early return.

Acceptance is a pure evaluation of the current, fully refreshed state.  It
does not accept cached booleans and cannot be called after restoring older
density or trace fields.  Electron number, Hermiticity, symmetry, and every
numeric residual are recomputed or read from the same state epoch.

The adaptive lambda state machine is deliberately small.  Accepted stages may
increase the bounded step; failed convergence, residual growth, or occupied
subspace discontinuity restores the accepted state and shrinks it.  Gap size
informs the proposal and degeneracy clustering but is not by itself a
material-specific rejection threshold.

## Final lambda-one refresh

After apparent convergence at lambda one, the solver performs one explicit
unmixed refresh from the converged occupied state:

1. reconstruct density and traces;
2. rebuild the volume operator;
3. combine it with the full interface operator;
4. solve the generalized problem;
5. reconstruct projector, density, and traces again;
6. evaluate every final residual and provenance gate once on that state.

The solver publishes this state directly.  It must not restore pre-refresh
density or trace values and then repeat acceptance.

## Symmetry and general systems

The same path supports crystals, liquids, defects, interfaces, and surfaces.
An authoritative analysis that finds no nonidentity operation produces the
explicit identity group.  A failed or missing analysis is not converted to
identity-only.  For a nontrivial actual group, retained-basis closure,
operator covariance, and covariance of the complete occupied subspace are
required.  Individual eigenvectors are never required to be symmetric.

RT driven-state symmetry is outside this ground-state solver.  At GS-to-RT
handoff only the accepted zero-field operator, metric, basis, and complete
occupied subspace are checked.

## Production integration boundary

`main_dft` performs only route selection, scope validation, construction of
the frozen catalog, invocation of the concrete solver, and final checkpoint
publication.  It does not implement the inner SCF loop.

The solver may call a small number of existing concrete SALMON routines for
volume assembly, distributed diagonalization, density reconstruction, and
real-space action.  It does not introduce an abstract backend, `class(*)`
context, procedure-pointer table, or public callback protocol.  If a test
needs a small deterministic problem, it supplies concrete arrays to the same
solver path rather than manually invoking internal phases.

Protected DC+LCFO/Wannier90, overlapping-Wannier, ordinary GS, and ordinary RT
branches are unchanged.

The concrete solver and its production catalog connection are implemented as
one task.  They must not be separated by a temporary generic backend.  The
catalog type is defined before the solver test and contains the actual matrix,
basis, face, selection, and ownership payload needed by the solver.  The
`main_dft` branch constructs that type through existing SALMON routines; the
solver never guesses how to rebuild missing production data.

## Tests

The principal MPI fixture runs the complete solver loop on a small
nonorthogonal two-fragment problem with nonzero cross-fragment SIPG blocks.  It
must cover:

- 1, 2, 4, and 8 rank decomposition independence;
- exact DC seed-density provenance;
- every physical face receiving the same lambda exactly once;
- a nontrivial metric and the basis-space occupied projector;
- phase and degenerate occupied-space rotation invariance;
- independent rejection by each residual and physical gate;
- rank-local kernel failure without deadlock;
- complete transactional rollback;
- final lambda-one refresh with no stale state;
- identity-only and nontrivial-group cases.

Focused algebra, SIPG, selection, and face-trace tests remain useful, but they
do not replace this full-path fixture.  Si64 remains the first production
acceptance calculation with eight MPI ranks, `OMP_NUM_THREADS=1`, and no time
cutoff; it is not embedded as a solver assumption.

## Deliberately omitted mechanisms

The first production solver does not add independent trace mixing, raw
eigenvector mixing, a polymorphic backend, plugin callbacks, a second
acceptance receipt layer, face-local continuation, basis changes inside a
lambda stage, or RT driven-state symmetry acceptance.  None is required to
establish the complete self-consistent DG ground state.
