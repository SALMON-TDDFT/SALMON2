# WPW Production API Hardening Design

## Scope

This checkpoint hardens the distributed `windowed_kg` operator infrastructure
before adding the quadrature assembler or production SCF consumer.  It does not
change Tasks 0--4 and does not connect the dense fixed-basis oracle to
production.

## Decisions

1. Callback contexts are created only through a binder that associates
   explicitly persistent `TARGET` storage.  The binder validates every target
   and records a bound state; callbacks reject an unbound context.
2. WW Hamiltonian inputs have an explicit component contract: the real WW
   blocks contain local, kinetic, and SIPG terms only, while the complex WW
   blocks contain the nonlocal pseudopotential contribution only.  The binder
   requires these distinct component tags before H application.
3. Sparse candidates carry provenance.  WP candidates may be volume or face
   contributions; PP candidates must be volume contributions.  The builder
   fails closed on a PP face or penalty candidate.
4. Candidate generation is rank-local.  The builder accepts only candidates
   whose output row/column is owned by the calling rank; it no longer scans a
   globally replicated candidate list and filters it afterward.
5. MPI communication failures are fatal for the operator communicator.  MPI
   collective errors cannot be made safely recoverable by a later collective,
   so the exchange layer calls `MPI_Abort` after reporting the failed phase
   instead of returning one rank into a mismatched collective sequence.

## Validation boundary

Each contract receives a regression test that is observed RED before the
minimal implementation.  The existing focused tests, Fortran numerical
fixtures, 3-rank owner exchange fixture, 2-rank dense-equivalence fixture, and
MPI/EigenExa SALMON build must pass before the checkpoint commit.

The next checkpoint remains the real-grid/canonical-face quadrature assembler.
Production occupied-subspace SCF integration follows only after that assembler
can supply support-local WP/PP COO entries.
