# DG Continuation Review Remediation Design

## Purpose

Correct the production WF+PW DG continuation before extending its checkpoint.
The remediation must preserve the explicit hybrid route and every legacy
DC+LCFO/Wannier90 route while making the procedure valid for insulating,
metallic, liquid, defect, surface, and interface systems.

## Fixed physical decomposition

The continuation Hamiltonian remains

\[
H_\lambda[\rho]=T_{\rm broken}+V_{\rm NL}+V_{\rm local}[\rho]
  +\lambda H_{\rm SIPG}.
\]

`T_broken`, `V_NL`, the metric, the retained basis, and `H_SIPG` are frozen
for one continuation attempt.  Only `V_local[rho]` changes in the inner SCF.
The exact converged DC density is the first Hamiltonian input.  Density is the
only mixed nonlinear variable.

The nonlocal projector coefficient has two deliberately distinct forms:

- projected matrix assembly uses `hvol * rinv_uvu`;
- pointwise strong action uses `rinv_uvu`.

They must not share one ambiguously named strength array.

## Setup communication

Interior basis data are materialized once before continuation.  Requests are
grouped by basis owner, fragment, and destination rank.  A destination sends
its required physical point IDs once; the owner returns all requested basis
columns with value, three Cartesian derivatives, and kinetic stencil action in
one packed response.  This remains point-to-point communication and removes
the previous basis-by-basis, all-peer request loop.

Projectors use the physical key `(atom_id, projector_ordinal)`.  Keys are
sorted deterministically, assigned one canonical owner, and partial overlaps
are reduced to that owner.  Only the complete overlaps needed for distributed
matrix rows and locally supported strong action are transferred.  No rank
stores every fragment copy and no quadratic duplicate search is used.

Every rank-local validation sets a local failure flag.  Before any subsequent
collective, the flag is reduced on the same communicator.  All ranks either
continue or return together.

## Density and potential update

The continuation reuses the established DC density redistribution and Hartree
FFT layout.  It does not allocate and all-reduce a complete real-space density
on every rank.  Hartree is computed from the distributed total density.  XC is
evaluated fragment-locally, including the required halo, and the resulting
local potential is projected without recomputing the frozen kinetic matrix.
SCF work arrays are allocated outside the inner loop.

## Inner and candidate gates

Every inner iteration performs only the work needed to update the fixed point:

1. distributed density to Hartree/XC potential;
2. local-potential matrix update;
3. generalized eigensolve;
4. occupation kernel, occupied projector, density, and interface traces;
5. coefficient, density, trace, orthogonality, electron-number, and projector-
   change residuals;
6. density mixing or candidate evaluation.

Hamiltonian covariance, occupied-projector covariance, Hermiticity, and the
reconstructed real-space DG residual run only after the inexpensive gates
pass.  Identity-only symmetry still executes the common validation boundary,
but its prepared identity action returns an exact zero covariance defect
without dense matrix communication.

The final lambda-one refresh is a separate, unmixed solve-and-check operation
with its own evaluation opportunity; it does not consume the last ordinary
SCF iteration through an in-loop `cycle`.

## Complete reconstructed DG residual

The volume channel applies

\[
T_{\rm broken}+V_{\rm NL}+V_{\rm local}[\rho]
\]

pointwise to each reconstructed occupied orbital using the SALMON Cartesian
stencil and unweighted nonlocal strong coefficient.  SIPG is not represented
as a fictitious volume field.  Its consistency, adjoint-consistency, and
penalty functionals are applied to the reconstructed face traces and compared
with the corresponding generalized eigen-equation boundary functional.  The
three normalized face-action residuals are reported independently.  `R_T`,
which measures change of interface observables between iterations, remains a
fixed-point residual but cannot substitute for any face-action residual.

A candidate passes only when the volume residual and all three face-action
residuals meet their stage tolerances.  The lambda-one refresh repeats them at
the final tolerances.

## Occupations and crossings

There is no positive HOMO--LUMO-gap acceptance gate.  The solver obtains all
states required by the configured occupation/smearing model.  The
occupation-weighted density matrix and the occupied or thermally active
`S`-metric projector define state continuity.  A meaningful gap may reduce the
next lambda step when it shrinks, but a zero gap does not itself reject a
stage.  Crossings are accepted when the occupation kernel and tracked
subspace remain continuous and all physical residuals converge.

## Mixing and state ownership

The controller owns the single density damping parameter.  The production
driver does not duplicate it as a local constant.  Rejection restores the
accepted density and controller histories; coefficients and derived caches
are invalidated and rebuilt by the next solve.  Obsolete continuation-specific
mixing or callback paths are removed once no caller remains.

## Checkpoint simplification

Task 8 resumes only after this remediation passes.  The checkpoint stores one
authoritative metric graph and does not serialize a second operator-side copy
of the metric.  The operator-union graph stores Hamiltonian components and
position operators.  Its graph may differ from the metric graph.  The complete
payload verifies

\[
H^{DG}(0)=T_{\rm broken}+V_{\rm NL}+V_{\rm local}[\rho^*]+H_{\rm SIPG}
\]

from serialized component values before publication and after reading.

## Verification strategy

Each correction starts with a focused failing test.  Required fixtures cover
nonunit `hvol`, rank-local projector failure without deadlock, distinct volume
and SIPG face-action residual failures, zero-gap and fractional-occupation
continuation, final refresh at the former iteration boundary, batched
point-to-point materialization, distributed Hartree use without replicated
full-grid arrays, candidate-only expensive checks, and independent checkpoint
graphs with no duplicate metric payload.

All legacy routes and existing MPI decompositions remain protected.  No Si64
completion claim is allowed until the later checkpoint, GS-to-RT identity,
zero-field stationarity, and eight-rank acceptance tasks pass.
