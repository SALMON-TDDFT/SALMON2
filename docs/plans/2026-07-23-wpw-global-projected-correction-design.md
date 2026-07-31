# WPW Global Projected Correction Design

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

**Date:** 2026-07-23

**Status:** Approved

## Goal

Test whether a global matrix-free solution of the projected correction equation

\[
P_L(H-\varepsilon_i S)P_R z_i=-r_i,\qquad Q^\dagger S z_i=0
\]

can reduce both occupied and extra-state fixed-H WPW residuals without the
over-amplification and search-space rank collapse produced by fragment-local
inverse approximations.

This remains a default-off diagnostic comparison. It does not promote the
correction into normal outer SCF, checkpoint schema, publication, or RT
handoff.

## Evidence and decision

The fragment-local metric-block correction improved occupied residuals but
worsened extra states. The fragment-local generalized `H-epsilon S` spectral
correction then failed more strongly: at inner 160 it worsened the Task 16
occupied residual by 8.77 times and the extra residual by 85.4 times.
Correction/residual amplification was approximately 176 for occupied and 370
for extra states, even though the denominator floor never fired. The
search-metric effective rank fell to 206 and discarded residual weights rose
to 0.643 for occupied and 0.297 for extra states.

The next comparison therefore removes the fragment-local inverse
approximation. It solves the distributed correction equation through the
completed bounded H/S actions and global Gram callback.

## Alternatives

### Global projected restarted GMRES

Use a fixed-memory restarted GMRES solve for every retained state, batching H/S
applications across active states. GMRES does not require the projected
operator to be positive definite or exactly self-adjoint in the Euclidean
metric. This is the selected approach.

### MINRES-QLP

MINRES-QLP is attractive for Hermitian indefinite and nearly singular
problems. However, the S-metric left/right projectors require a careful
self-adjointness proof and a new robust QLP implementation. This is deferred
until the global GMRES gate establishes whether solving the global equation is
physically useful.

### State-partitioned local regularization

Use separate occupied and extra-state damping or inverse limits. This is less
expensive but adds an empirical state-label policy while retaining the failed
fragment-local approximation. It remains a fallback only if the global solve
shows useful occupied directions but unacceptable extra-state cost.

## Projected equation

At each fixed-H algebra step, the retained Ritz vectors satisfy

\[
Q^\dagger S Q\approx I.
\]

Define

\[
P_R x=x-Q(Q^\dagger Sx)
\]

and the compatible left projection

\[
P_L y=y-SQ(Q^\dagger y).
\]

The correction operator for state `i` is

\[
\mathcal A_i x=P_L(H-\varepsilon_iS)P_Rx.
\]

The right-hand side is `-P_L r_i`. The returned correction is projected once
more with `P_R` before entering the existing three-block Rayleigh--Ritz search.
All projection coefficients are computed through the existing global Gram
callback. No fragment-local dense inverse or global dense H/S assembly is
allowed.

The solver must measure the initial projected residual explicitly. A
zero-residual state returns a zero correction without entering Arnoldi.

## Restarted batched GMRES

Use classical restarted GMRES with modified Gram--Schmidt and one reorthogonal
pass. The first implementation has:

- one shared restart length and maximum iteration count;
- one relative residual tolerance with an absolute finite floor;
- one bounded state-batch size so Krylov storage does not scale as
  `restart*nretain`;
- independent state convergence masks and residual histories;
- batched H/S application for all currently active states;
- deterministic collective ordering even after individual states converge;
- no inner preconditioner and no state-dependent damping.

Each state owns its small Hessenberg matrix, Givens rotations, and least-square
right-hand side. Process retained states in deterministic contiguous batches
of at most eight states. Distributed Krylov vectors remain in the existing
local W/P row layout, so storage is bounded by
`(restart+1)*state_batch_size` vectors rather than
`(restart+1)*nretain`. Converged states contribute zero columns to later
batched operator calls so every rank executes the same callback sequence.

The correction callback remains compatible with
`dg_wpw_preconditioner(context,eigenvalues,rw,rp,zw,zp,info)`. Solver controls
and diagnostics are stored in the production context rather than changing the
matrix-free outer algebra interface.

## Configuration and routing

Add a default-`n` diagnostic control named
`yn_dg_wpw_global_projected_correction`. It is mutually exclusive with:

- `yn_dg_wpw_preconditioner`;
- `yn_dg_wpw_metric_preconditioner`;
- `yn_dg_wpw_h_epsilon_s_correction`.

The S-orthogonal coordinate comparison remains an independent fixed-H control
and stays `n` for the first physical gate. The new correction control requires
`yn_dg_wpw_fixed_h_relaxation='y'`.

Expose bounded diagnostic solver parameters only:

- `dg_wpw_global_correction_restart`;
- `dg_wpw_global_correction_max_iter`;
- `dg_wpw_global_correction_tolerance`;
- `dg_wpw_global_correction_state_batch`.

Use comparison defaults restart `8`, maximum iterations `32`, relative
tolerance `1.0E-2`, and state batch `8`. Collectively reject restart outside
`[1,16]`, maximum iterations outside `[1,64]`, nonfinite tolerance outside
`(0,1)`, or state batch outside `[1,16]`. A maximum iteration count smaller
than the restart length, or a nonmultiple of it, defines a shortened final
cycle and is valid. Before allocation, validate all dimension products against
integer overflow and the explicit restart/state-batch caps.

Repository defaults keep the route disabled. Normal production algebra stays
callback-free. Only fixed-H and its interface continuation may select the new
callback.

## Lifecycle

The global solver does not cache an H/S factor. Each application binds to the
current completed bounded operator through its epoch and layout fingerprint.
The callback snapshots those identifiers before GMRES and validates them
collectively after every batched H/S application and before returning.

Continuation may reuse solver work arrays, but not Krylov vectors, across
lambda changes. Every correction call starts from zero. A lambda trial or
rollback must use the newly installed operator epoch and may not continue an
old Arnoldi cycle.

## Collective safety and failure handling

Before the first operator call, collectively validate:

- communicator agreement;
- W/P row shapes and retained-state count;
- finite Ritz values, residuals, Q vectors, and solver controls;
- current bounded-operator epoch and layout fingerprint;
- finite and sufficiently orthonormal `Q^dagger S Q`.

After every projection, H/S action, Gram operation, Arnoldi update, and
least-square solve, reduce a failure flag before any rank can return.

Fail closed for:

- rank-disagreeing shapes or active-state masks;
- nonfinite vectors, Hessenberg entries, rotations, or residual estimates;
- stale epoch or fingerprint;
- Arnoldi breakdown that does not satisfy the residual tolerance;
- failure to reach tolerance within the maximum iteration count;
- final `Q^dagger S z` or explicit correction-equation residual above its
  acceptance bound.

Do not fall back to fragment-local, diagonal, or uncorrected search directions
after the route has been selected.

For each state let `b=-P_L r` and

\[
\tau_{\rm abs}=100\,\epsilon_{\rm mach}\max(1,\lVert b\rVert_2).
\]

If `||b||_2 <= tau_abs`, classify the state as a zero projected residual and
return a zero correction without Arnoldi. Otherwise accept only when the
explicitly recomputed final residual satisfies

\[
\lVert b-\mathcal A_i z\rVert_2
\le \max(\tau_{\rm abs},
         {\tt dg\_wpw\_global\_correction\_tolerance}\lVert b\rVert_2).
\]

A happy breakdown is success only when this same explicit criterion passes.
The recursive Hessenberg residual estimate is diagnostic and may not replace
the final explicit check.

## Diagnostics

Log the route and solver controls once at the fixed-H boundary. At inner
32/96/160 record occupied and extra-state maxima for:

- raw outer residual;
- initial and final projected correction-equation residual;
- relative residual reduction;
- GMRES iteration and restart counts;
- happy and failed breakdown counts;
- correction norm and correction/raw amplification;
- final `Q^dagger S z`;
- final explicit projected-equation defect;
- correction fraction removed by the final S projection;
- existing search-metric effective rank and discarded weights;
- existing Ritz post/direct defects and metric orthogonality.

Diagnostics must distinguish converged occupied and extra state counts. A
single failed state fails the comparison route, but its evidence is logged
before the collective failure is returned.

## Testing

Add a two-rank dense-oracle fixture with nonidentity SPD `S`, noncommuting
Hermitian `H`, multiple Ritz values, and a retained Q subspace. Verify:

- left and right projector identities;
- equality of matrix-free and dense projected operator applications;
- GMRES correction agreement with the dense constrained solve;
- zero-residual and happy-breakdown behavior;
- independent convergence masks with identical collective ordering;
- final S orthogonality and explicit projected-equation residual;
- collective rejection of nonfinite, malformed, stale, and nonconverged
  inputs.

Source contracts verify default-off input plumbing, mutual exclusion,
fixed-H/continuation-only routing, unchanged search-history propagation, and
the absence of a callback in normal production algebra.

Run focused two-rank tests, generalized-algebra and row-layout regressions,
the temporary parent-prerequisite overlay MPI/EigenExa build, and code review.
Do not import the parent's uncommitted prerequisite implementation into this
branch.

## Physical gate

Clone the Task 16 B=6 restart into a fresh directory. Preserve basis, seed,
cutoff, tolerance, continuation, search history, and publication settings.
Set diagonal, metric-block, H-epsilon-S, and S-orthogonal comparison controls
to `n`; enable only global projected correction.

Run eight MPI ranks through the existing fixed-H boundary. Compare inner
32/96/160 against Task 16, Task 19, the S-orthogonal gate, and the rejected
H-epsilon-S gate.

Improvement requires both occupied and extra residuals to decrease without
GMRES nonconvergence, severe search-rank loss, Ritz inconsistency, or
publication-safety regression. Do not promote the route automatically after
the gate.
