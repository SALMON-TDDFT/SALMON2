# WF+PW DG Continuation Ground-State Design

## Purpose

The hybrid WF+PW real-time route must start from the self-consistent ground
state of the same complete discontinuous-Galerkin Hamiltonian that is used at
zero external field in real time.  A one-shot projection of the ordinary
real-space Kohn--Sham `hpsi` into the WF+PW basis does not establish that
condition.  The hybrid route therefore uses the converged divide-and-conquer
(DC) state as a seed and introduces the complete DG interface operator by an
adaptive continuation from zero to full strength.

This is a new, explicitly selected hybrid WF+PW route.  It does not change the
existing DC+LCFO, Wannier90, overlapping-Wannier, ordinary ground-state, or
ordinary real-time routes.

## Initial state and immutable continuation catalog

The initial density is exactly the converged total density produced by the DC
calculation:

\[
  \rho^{(0)}_{\lambda=0}=\rho_{\mathrm{DC}}^{\mathrm{conv}}.
\]

It is not replaced by a density reconstructed from a preliminary WF+PW
diagonalization.  The DC density payload and its fingerprint are copied into
the first continuation checkpoint.  The seed gate is

\[
 R_{\rho,\mathrm{seed}}=
 \frac{\lVert\rho^{(0)}_{\lambda=0}-
 \rho_{\mathrm{DC}}^{\mathrm{conv}}\rVert}
 {\max(1,\lVert\rho_{\mathrm{DC}}^{\mathrm{conv}}\rVert)}
 \leq \tau_{\mathrm{seed}}.
\]

The converged DC occupied space is projected into the symmetry-closed WF+PW
basis and orthonormalized in the DG metric.  Only its occupied projector is
tracked between iterations; raw eigenvector coefficients are not mixed.
Interface traces for the initial state are reconstructed from this occupied
projector.

Two different gauge-invariant matrices are used.  For an integer-occupied
subspace,

\[
 Q_{\mathrm{occ}}=C_{\mathrm{occ}}C_{\mathrm{occ}}^\dagger S^{DG}
\]

tracks the occupied subspace.  Physical density, electron number, energy, and
occupied interface density matrices use the occupation-weighted density
matrix

\[
 \Gamma_{\mathrm{occ}}=C f C^\dagger.
\]

States within a degenerate cluster may rotate among themselves.  Symmetry is
therefore never required of an individual eigenvector.

For one complete continuation attempt, the following catalog is immutable:

- the symmetry-closed WF+PW basis and its distributed row ownership;
- fragment geometry, physical face topology, normals, and quadrature;
- cutoff, selection, window, packet, and complement definitions;
- the DG metric and its numerical-rank decision;
- the symmetry representation and interface orbits;
- the SIPG penalty convention.

The retained WF+PW space is selected in complete symmetry blocks.  If a
cutoff intersects a degenerate multiplet, reciprocal star, or another
symmetry orbit, the whole block is retained or the whole block is omitted.
The retained-space projector must have negligible symmetry leakage,

\[
 R_{\mathrm{leak}}(g)=
 \lVert(1-Q_{\mathrm{ret}})D(g)Q_{\mathrm{ret}}\rVert.
\]

This closure requirement is distinct from convergence with respect to the
size of the omitted excited-state space.  A finite excitation cutoff may bias
a response even when the retained space is exactly symmetry closed; that
physical, observable-dependent cutoff convergence is outside the present
ground-state acceptance and must not be mislabeled as a symmetry violation.

If the final real-space DG residual proves that this finite basis is
insufficient, the route expands or revises the basis outside the continuation
loop and restarts from the DC seed.  It never changes basis silently inside a
lambda stage.

## Complete DG operator

For the fixed catalog, the zero-field Hamiltonian is

\[
 H_\lambda[\rho]=H_{\mathrm{volume}}[\rho]
                 +\lambda H_{\mathrm{interface}},
 \qquad 0\leq\lambda\leq1.
\]

The interface bilinear form contains every symmetric interior-penalty term:

\[
 a_\Gamma(u,v)=
 -\int_\Gamma \{\partial_nu\}[v]\,dS
 -\int_\Gamma [u]\{\partial_nv\}\,dS
 +\int_\Gamma \frac{\eta}{h}[u][v]\,dS.
\]

Thus the projected interface operator contains the numerical/consistency
flux, the adjoint-consistency flux, the penalty term, and both directions of
every complete neighboring-fragment coupling block.  Every physical face has
one canonical owner for assembly, but both matrix blocks are published.  The
assembled operator must pass Hermiticity, reciprocal-face, topology, and
internal-cancellation gates.

With the catalog fixed, the basis traces and the projected
`H_interface` are independent of the current eigenvector coefficients and
density.  Lambda scales this one operator uniformly.  The state-dependent
quantities are density, potential, occupied projector, occupied interface
observables, occupations, eigenvalues, and residuals.  Stale neighboring
eigenvectors are never inserted as an external boundary condition.

If the selected discretization requires a nontrivial DG overlap, the same
catalog constructs `S_DG`.  It remains fixed during the continuation.  A
lambda-dependent metric is permitted only if derived explicitly from the same
bilinear discretization; in that case it is scaled uniformly over every
symmetry-related interface and included in every residual and fingerprint.

## Coupled fixed point at one continuation stage

Let

\[
 0=\lambda_0<\lambda_1<\cdots<\lambda_N=1.
\]

At a trial stage, each inner iteration performs the following operations in
this order:

1. Build the full volume Hamiltonian from the current mixed density.
2. Add the uniformly scaled, preassembled complete interface operator.
3. Solve the distributed generalized eigenproblem.
4. Form the occupied `S`-projector and align the occupied subspace to the
   accepted or preceding projector in the `S` metric.
5. Reconstruct the output density and gauge-invariant occupied interface
   observables from that projector.
6. Evaluate all unmixed output-minus-input residuals.
7. If not converged, update density and interface observables with independent
   damping and repeat without advancing lambda.

The projector is

\[
 P=C_{\mathrm{occ}}C_{\mathrm{occ}}^\dagger S^{DG}.
\]

The principal subspace diagnostic is an `S`-metric projector difference or
equivalent principal-angle measure.  Eigenvectors may be Procrustes-aligned
for deterministic traces and output, and degenerate occupied clusters may be
diagonalized with symmetry operators, but aligned vectors are not themselves
mixed.

The independent updates are

\[
 \rho^{(m+1)}=\rho^{(m)}+
 \alpha_\rho\bigl(\rho[C^{(m)}_{\rm occ}]-\rho^{(m)}\bigr),
\]

\[
 T^{(m+1)}=T^{(m)}+
 \alpha_T\bigl(T[C^{(m)}_{\rm occ}]-T^{(m)}\bigr).
\]

`T` contains gauge-invariant occupied face density matrices and the derived
wavefunction-value and normal-derivative observables needed to check the DG
interface fixed point.  Density and trace residual histories remain separate.
The default trace damping is no larger than density damping.  Anderson, Pulay,
or Broyden acceleration may be added independently to either channel, but the
first implementation uses bounded damped updates and does not combine their
residual vectors.

## Residuals and stage acceptance

Every inner iteration measures

\[
 R_H=\frac{\lVert HC-SC\varepsilon\rVert}
 {\max(1,\lVert HC\rVert,\lVert SC\varepsilon\rVert)},
\]

\[
 R_\rho=\frac{\lVert\rho[C_{\rm occ}]-\rho\rVert}
 {\max(1,\lVert\rho\rVert)},\qquad
 R_T=\frac{\lVert T[C_{\rm occ}]-T\rVert}
 {\max(1,\lVert T\rVert)},
\]

\[
 R_S=\lVert C^\dagger S^{DG}C-I\rVert.
\]

A lambda stage is accepted only when all of the following pass the tolerance
assigned to that stage:

- generalized eigensystem residual;
- density residual;
- interface residual;
- occupied-projector change;
- `S`-orthogonality residual;
- electron-number error;
- Hamiltonian and metric Hermiticity;
- full-system symmetry covariance;
- finite values and valid occupied--unoccupied gap handling.

Intermediate stages may use inexact tolerances that tighten monotonically as
lambda approaches one.  Lambda one always uses the normal final acceptance
tolerances.

In addition to the coefficient-space residual, the code reconstructs selected
occupied and near-gap Ritz orbitals in real space, applies the complete DG
volume-plus-interface action there, and measures the normalized residual.  A
small coefficient residual cannot compensate for a failed real-space residual.
This is the finite-basis adequacy gate.

## Adaptive continuation and rollback

The controller starts from the accepted lambda-zero DC-seeded checkpoint.  A
trial step is selected from the previous accepted step and is constrained by
minimum and maximum bounds.  It considers:

- growth of `R_H`, `R_rho`, and `R_T`;
- occupied-projector displacement;
- the occupied--unoccupied gap;
- eigenvalue crossings and occupied-subspace discontinuity;
- iteration count and convergence rate;
- symmetry and real-space residuals.

Fast, well-separated convergence permits a bounded increase in the next
lambda step.  Residual growth, gap collapse, subspace discontinuity, or failure
within the iteration limit rejects the complete trial stage.  Rejection
atomically restores density, potential, projector, trace state, occupations,
eigenvalues, mixing histories, operator epochs, and all derived caches from the
last accepted checkpoint, then reduces the step.  No object with the rejected
epoch may survive rollback.

Lambda is a single scalar for the whole system.  It is applied simultaneously
to every interface in the same and different symmetry orbits.  Fragment-local
or face-local continuation parameters are forbidden.

## Ground-state symmetry contract

Before continuation, the basis must be closed under every required full-system
symmetry operation.  At every accepted stage the code evaluates

\[
 R_{\mathrm{sym},H}(g)=
 \frac{\lVert D(g)^\dagger H_\lambda D(g)-H_\lambda\rVert}
      {\max(1,\lVert H_\lambda\rVert)},
\]

and the analogous metric, retained-space, and occupied-subspace residuals.
For equal integer occupations this means

\[
 D(g)^\dagger Q_{\mathrm{occ}}D(g)=Q_{\mathrm{occ}}.
\]

For partial occupations, symmetry is evaluated using
`Gamma_occ`; symmetry-related degenerate states must receive compatible
occupations.  The face topology and interface blocks must map covariantly
under the same operation.  Diagonalization is not used as a symmetry repair.
Failure of the basis, retained-space closure, occupied-subspace covariance,
topology, or zero-field operator covariance closes the hybrid route with an
error.  No individual occupied or empty eigenvector is required to transform
as a one-dimensional invariant state.

## Fully refreshed lambda-one state

After the inner loop first satisfies all lambda-one conditions, the route
performs a final unmixed refresh:

1. Reconstruct density and all interface observables from the accepted
   occupied projector without mixing.
2. Rebuild the complete density-dependent volume operator.
3. Combine it with the full unscaled interface operator.
4. Solve or re-evaluate the generalized eigenproblem against this refreshed
   operator.
5. Re-evaluate coefficient-space and real-space `R_H`, `R_rho`, `R_T`, `R_S`,
   electron count, projector, Hermiticity, gap, and symmetry gates.

Only this refreshed fixed point may be published as the hybrid ground state.
No delayed density, old neighboring coefficient, old trace, or rejected-stage
cache may enter the final payload.

## Atomic GS checkpoint and RT handoff

The ground-state writer publishes one atomic checkpoint containing:

- the complete zero-field `H_DG(0)` payload;
- the complete `S_DG` payload;
- the distributed WF+PW basis catalog and ownership;
- cutoff, selection, window, packet, and complement metadata;
- occupied coefficients, occupations, and eigenvalues;
- final density and interface observables;
- symmetry and face topology metadata;
- DC seed-density fingerprint;
- basis, metric, operator, state, and complete-payload fingerprints;
- lambda history, rollback history, tolerances, and final residual receipt.

The matrix payload, state, and provenance are hashed together.  It is invalid
to reconstruct another operator independently and assign it the stored
fingerprint.

`main_tddft.f90` gains a separate hybrid WF+PW RT branch.  It reads the actual
metric, zero-field Hamiltonian, state, and catalog payloads from this file and
does not rebuild the zero-field operator independently.  Before the first time
step it rechecks generalized residual, `S` orthogonality, electron number,
Hermiticity, zero-field basis/operator/occupied-subspace symmetry, and payload
identity.  Any mismatch fails closed.

## Zero-field real-time acceptance

An end-to-end test starts RT from the accepted lambda-one checkpoint with no
external field.  At every sampled time it measures drift in:

- density;
- total energy;
- occupied `S`-projector;
- electron number;
- the complete DG Hamiltonian residual.

The stationary-state check compares occupied projectors or phase-aligned
subspaces, never raw coefficient differences.  Because the initial
occupied projector is a required symmetric ground-state object, its
stationarity already preserves that symmetry indirectly; a separate symmetry
drift gate is not required for the zero-field RT acceptance.  Symmetry of an
externally driven RT state is observable-, polarization-, gauge-, and
retained-excitation-space dependent and is outside this design.  The Si64
production acceptance uses eight MPI ranks, `OMP_NUM_THREADS=1`, and no
time-based termination.

## Failure handling and diagnostics

All stage decisions are collective.  Invalid topology, nonfinite data,
rank-disagreeing provenance, non-Hermiticity, symmetry loss, electron-number
failure, insufficient basis residual, or exhausted rollback limits aborts the
hybrid route without publishing a checkpoint.  Existing accepted checkpoints
remain intact.  Diagnostic output records every lambda proposal, acceptance or
rejection reason, residual channel, damping rate, gap, projector change,
symmetry maximum, real-space residual, and restored checkpoint epoch.

## Test strategy

Implementation follows strict TDD.  Each behavior is introduced by a focused
failing Python-runner or MPI fixture, the failure is observed, the minimal
implementation is added, and the same runner is observed passing before the
task commit.

The test layers are:

1. dense algebra tests for SIPG blocks, projector invariance, residuals, and
   symmetry covariance;
2. MPI tests for canonical face ownership, complete cross-fragment blocks,
   uniform lambda, trace exchange, and atomic rollback;
3. continuation tests for adaptive steps, gap/crossing rejection, independent
   damping, inexact tolerances, and lambda-one refresh;
4. checkpoint corruption, provenance, and exact-payload GS-to-RT tests;
5. zero-field RT stationarity tests using projector comparison, without an
   independent driven-state symmetry gate;
6. unchanged-route regression tests for DC+LCFO/Wannier90 and existing RT;
7. the final eight-rank Si64 run with no timeout.

No prior one-shot-LCFO result is accepted as evidence for this design.  The
interrupted calculation and its logs remain preserved as historical diagnostic
evidence only.
