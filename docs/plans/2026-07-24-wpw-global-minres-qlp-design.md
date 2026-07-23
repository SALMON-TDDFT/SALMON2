# WPW Global MINRES-QLP Correction Design

**Date:** 2026-07-24

**Status:** Approved

## Goal

Test whether a global matrix-free MINRES-QLP solution of

\[
P_L(H-\varepsilon_i S)P_R z_i=-P_Lr_i,\qquad Q^\dagger S z_i=0
\]

can produce an accepted correction for every retained state and improve both
occupied and extra-state fixed-H WPW residuals.

This is a separate default-off diagnostic comparison. It does not change the
existing GMRES route and does not promote either correction into normal outer
SCF, checkpoint schema, publication, or RT handoff.

## Evidence and decision

The global restarted-GMRES B=6 gate failed closed at the first fixed-H
correction. Seven of 160 state solves did not reach relative tolerance
`1E-2` within 32 iterations. The maximum explicit projected-equation residual
and equation defect were both `4.6571063211358618E-04`; no failed Arnoldi
breakdown or stale operator was reported.

The failure is therefore a convergence limitation rather than a callback,
provenance, or arithmetic-breakdown failure. The projected operator is
Hermitian but may be indefinite or singular, so MINRES-QLP is the next
comparison. A state-partitioned policy is deferred because excluding or
damping only the difficult seven states could hide the physical obstruction
before it is understood.

## Alternatives

### Independent MINRES-QLP route

Add a mutually exclusive default-off route that uses the same global
projected operator and explicit acceptance test as GMRES. This isolates the
solver comparison and is the selected approach.

### Automatic GMRES-to-MINRES-QLP fallback

Retry failed GMRES states with MINRES-QLP. This could reduce average work, but
it would confound the B=6 comparison and introduce solver-history semantics
before either route qualifies. It is rejected for the first gate.

### State-partitioned correction

Apply different policies to occupied, extra, or empirically difficult states.
This is cheaper but introduces state labels and acceptance exceptions. It is
reserved for a later design only if MINRES-QLP identifies a stable converged
subset with useful outer-residual improvement.

## Operator and Hermitian contract

For an S-orthonormal retained Ritz block `Q`, define

\[
P_R=I-Q Q^\dagger S,\qquad P_L=I-SQ Q^\dagger=P_R^\dagger .
\]

With Hermitian `H` and `S`,

\[
\mathcal A_i=P_L(H-\varepsilon_iS)P_R
\]

is Hermitian in the coefficient-space Euclidean inner product:

\[
\mathcal A_i^\dagger
=P_R^\dagger(H-\varepsilon_iS)^\dagger P_L^\dagger
=P_L(H-\varepsilon_iS)P_R.
\]

The implementation reuses the completed bounded H/S callbacks, global Gram
callback, current post-Ritz `Q`, epoch, and layout fingerprint from the GMRES
route. It may not assemble a global dense H or S matrix.

Before iteration, evaluate a deterministic Hermitian probe for each active
state batch. For finite probe blocks `X` and `Y`, require the collectively
reduced defect in

\[
X^\dagger\mathcal A_iY-(\mathcal A_iX)^\dagger Y
\]

to be below a scale-aware tolerance. A failed probe rejects the route before
Lanczos iteration. The probe is diagnostic validation, not a replacement for
the algebraic source-contract tests and dense-oracle fixture.

## MINRES-QLP algorithm

Implement complex Hermitian MINRES-QLP using a Lanczos three-term recurrence,
stable symmetric orthogonal transformations, and the QLP transition described
by the standard MINRES-QLP recurrence. The solver starts every state from zero
and computes the minimum-length solution when the projected system is
singular or numerically rank deficient.

Use:

- one maximum-iteration bound;
- one relative explicit-residual tolerance with the existing absolute floor;
- one bounded state-batch size;
- independent per-state active masks;
- identical collective callback ordering on all ranks;
- no inner preconditioner, damping, GMRES fallback, or state exclusion.

The first version uses scalar Lanczos recurrences per state while batching all
H/S applications for the active columns. Inactive columns are explicitly
zeroed but remain in each collective operator call. Store only the bounded
Lanczos and solution recurrence vectors required by MINRES-QLP; memory must
not grow with the iteration limit.

Transition from the MINRES phase to QLP when the estimated condition number
exceeds a configured transition threshold or a diagonal factor becomes small
relative to the running operator norm. Once a state enters QLP it may not
return to the MINRES phase during that solve.

## Configuration and routing

Add a default-`n` control:

```text
yn_dg_wpw_global_minres_qlp_correction
```

It is mutually exclusive with:

- `yn_dg_wpw_preconditioner`;
- `yn_dg_wpw_metric_preconditioner`;
- `yn_dg_wpw_h_epsilon_s_correction`;
- `yn_dg_wpw_global_projected_correction`;
- `yn_dg_wpw_s_orthogonal_pw`.

Expose bounded diagnostic controls:

```text
dg_wpw_global_minres_qlp_max_iterations
dg_wpw_global_minres_qlp_tolerance
dg_wpw_global_minres_qlp_state_batch
dg_wpw_global_minres_qlp_transition_condition
```

Use first-gate values maximum iterations `64`, tolerance `1E-2`, state batch
`8`, and transition condition `1E7`. Validate all controls even when the route
is disabled. The maximum iteration range is `[1,128]`, batch range `[1,16]`,
tolerance must be finite in `(0,1)`, and transition condition must be finite
and greater than one.

Only fixed-H and every fixed-H continuation trial/rollback may select the
callback. Normal production algebra remains callback-free.

## Acceptance and failure handling

Use the same projected right-hand side, zero-RHS classification, absolute
floor, final right projection, S-orthogonality test, snapshot validation, and
explicit projected-equation acceptance inequality as the GMRES route:

\[
\|b-\mathcal A_i z\|_2
\le \max(\tau_{\rm abs}, {\tt tolerance}\|b\|_2).
\]

The recursive MINRES residual and QLP condition/rank estimates are diagnostic
only. Every provisionally converged state receives a fresh explicit operator
application before acceptance.

Fail closed for:

- non-Hermitian probe failure;
- rank-disagreeing shapes, masks, or phase transitions;
- stale epoch or layout fingerprint;
- nonfinite Lanczos or QLP scalars/vectors;
- inconsistent Lanczos norm or loss of recurrence beyond a scale-aware bound;
- incompatible singular system detected by a nonzero residual floor;
- iteration exhaustion;
- final residual, S orthogonality, or projection defect above its bound.

A single failed state rejects the complete comparison callback. Do not fall
back to GMRES, a local inverse, diagonal preconditioning, or an uncorrected
direction after route selection.

## Diagnostics

Preserve the existing common correction diagnostics and add MINRES-QLP
specific per-state fields:

- phase at termination (`MINRES` or `QLP`);
- transition iteration;
- Lanczos iterations;
- compatible singular and incompatible classifications;
- estimated operator norm and condition number;
- estimated numerical rank;
- recursive and explicit final residual;
- solution norm and correction/raw amplification;
- maximum Hermitian-probe and Lanczos-recurrence defects.

At accepted correction counts 32/96/160, log occupied and extra maxima/counts
for these fields together with the existing outer residual, search-metric
rank/discarded weights, Ritz post/direct defects, and metric orthogonality.
On failure, log the same available diagnostics before returning collectively.

## Testing

Start with RED tests.

Extend the two-rank dense-oracle fixture with nonidentity SPD `S`,
noncommuting Hermitian `H`, and cases covering:

- nonsingular indefinite projected systems;
- compatible singular systems with a known minimum-length solution;
- incompatible singular systems that must fail closed;
- MINRES-to-QLP transition;
- zero projected RHS;
- independent state convergence and phase masks;
- iteration exhaustion;
- non-Hermitian callback corruption;
- stale operator and nonfinite callback results;
- scale invariance of acceptance and transition decisions.

Compare corrections with a dense constrained S-complement pseudoinverse
oracle. Verify explicit residual, minimum-length property, S orthogonality,
bounded iteration-storage source contracts, collective ordering, and
transactional diagnostics.

Run the existing GMRES, matrix-free SCF, generalized algebra, production
operator, occupied-W, and W-row-layout regressions. Perform the same temporary
parent-prerequisite MPI/EigenExa overlay build and resolve all
Critical/Important review findings before each task commit.

## Physical gate

Clone the Task 16 B=6 restart into a fresh `/tmp` directory. Preserve its
basis, normalized seed, cutoff, outer tolerance, continuation, search history,
and publication settings. Disable every existing correction route and enable
only global MINRES-QLP with the first-gate controls above.

The first acceptance boundary is stricter than the GMRES run: all 160 states
must pass the explicit correction-equation criterion and fixed-H must reach
inner 32. Only then compare inner 32/96/160 occupied and extra residuals,
search-metric effective rank/discarded weights, Ritz defects, solver phase,
condition/rank estimates, final `info`, and publication state against Task 16,
Task 19, the S-orthogonal comparison, fragment-local H-epsilon-S, and global
GMRES.

Improvement requires both occupied and extra residuals to improve without
solver nonconvergence, incompatible singular states, rank collapse, Ritz
inconsistency, or publication regression. No result automatically authorizes
promotion.
