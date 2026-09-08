# WF+PW DG Production Interface Traces Design

## Purpose

The production WF+PW continuation route must assemble the complete SIPG
interface operator from the same frozen basis used by the generalized
eigenproblem.  A fragment basis containing only buffer-grid values is not a
sufficient production contract: the assembler also needs canonical faces,
the two adjacent fragment traces, and normal derivatives.

## Frozen interface payload

For each canonical face, store one immutable record containing the global face
ID, minus and plus fragment IDs, periodic image shift, canonical owner, normal,
quadrature weights, face point IDs, and the global basis IDs participating on
both sides.  Store basis values and derivatives with respect to the one
canonical normal.  The plus-side derivative is converted from its outward
normal when the record is built, not later in the SIPG algebra.

The payload fingerprint covers topology, ownership, geometry, basis IDs,
values, derivatives, and quadrature weights.  It is collective and must be
identical on all ranks.  Every physical interface has exactly one canonical
owner.  Periodic faces use the same representation as internal faces.

## Construction and data flow

The existing production fragment-basis construction remains responsible for
WF+PW materialization.  After symmetry closure has produced the effective WF
and PW selections, a separate builder reconstructs gradients with SALMON's
existing production stencil, identifies canonical neighboring fragment
faces, and samples both traces.  It rejects missing neighbors, duplicate
owners, incomplete quadrature correspondence, nonfinite traces, and a basis
selection that is not closed under the accepted group action.

The resulting immutable face payload is passed to
`assemble_dg_hybrid_sipg_face`.  Complete face blocks are assembled once and
remain coefficient-independent during a continuation attempt.  The
continuation controller applies one uniform scalar lambda to all blocks.
Occupied interface density matrices are derived separately from the current
occupation kernel and are acceptance diagnostics; they do not rebuild or mix
the fixed kinetic face blocks.

## Scope and protected paths

Only the explicit `yn_dg_hybrid_continuation_scf='y'` route constructs this
payload.  The existing divided one-shot, overlapping-Wannier, DC+LCFO, and
ordinary GS paths retain their present calls and data structures.  The common
fragment basis type may gain an optional interface component, but legacy
callers are valid when it is absent and never execute continuation validation.

## Tests

An MPI fixture uses two irregular fragments and a periodic image face.  It
checks canonical ownership, normal orientation, trace and derivative values,
collective fingerprints, nonzero cross-fragment SIPG blocks, Hermiticity, and
rejection of duplicate or incomplete faces.  A source contract then requires
the production continuation branch to build the frozen payload before
initializing and running the coupled continuation.
