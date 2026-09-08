# WPW H-epsilon S Correction Design

**Date:** 2026-07-23

**Status:** Approved

## Goal

Test whether a state-dependent approximation to the correction equation

\[
(H-\varepsilon_i S)z_i=-r_i
\]

can remove the fixed-H WPW residual plateau without the extra-state amplification
observed for the fragment-local `S_block^{-1}r` correction.

This is a default-off diagnostic comparison. It does not promote S-orthogonal
coordinates into normal outer SCF, checkpoint, publication, or RT handoff.

## Evidence and decision

Task 19 showed that `S_block^{-1}r` improves occupied residuals by roughly
26--150 times but amplifies extra-state corrections and worsens the final
maximum residual. The S-orthogonal PW complement gate then removed W--P metric
overlap to `6.47E-13` but improved residuals by only 1.05--1.27 times.

The next comparison must therefore include the local Hamiltonian energy scale,
not merely rescale the metric inverse. The first implementation uses a
fragment-local generalized spectral approximation. A global iterative
correction solve is deferred because it would introduce a second distributed
iteration and stopping policy before the local hypothesis is tested.

## Architecture

For each fragment, assemble the existing bounded local Hermitian blocks
`H_f` and `S_f` in the owned W+P ordering. At the fixed-H boundary, solve

\[
H_f U_f=S_f U_f\Theta_f,\qquad U_f^\dagger S_f U_f=I.
\]

Store `Theta_f`, `U_f`, the operator epoch, layout fingerprint, dimensions,
conditioning diagnostics, and regularization policy in a new correction
factor object.

For Ritz state `i`, form the local spectral residual

\[
\widehat r_i=U_f^\dagger r_i
\]

and apply

\[
\widehat z_{ji}=-g(\theta_j-\varepsilon_i)\widehat r_{ji},\qquad
z_i=U_f\widehat z_i.
\]

The sign convention is tested against the dense correction-equation oracle.
The correction callback receives the current Ritz eigenvalues, so no solver
interface change is required.

## Regularization

Use one state-independent rule for occupied and extra states in the first
gate. Let

\[
s_f=\max(1,\max_j|\theta_j|,|\varepsilon_i|),\qquad
\delta_i=\delta_{\rm rel}s_f.
\]

For denominator `d=theta_j-epsilon_i`, use the signed floor

\[
\widetilde d=\operatorname{sign}(d)\max(|d|,\delta_i).
\]

At exact zero, choose the positive floor deterministically. Also apply an
explicit maximum inverse magnitude `1/delta_i`; this is equivalent to the
floor but is retained as a logged invariant. Do not add separate
occupied/extra mixing constants in the first implementation.

The default comparison value for `delta_rel` is selected by a dense fixture
and logged in the B=6 run. It must be exposed as a bounded finite input only if
the existing input layer has a suitable real-valued diagnostic-control
pattern; otherwise the first implementation uses one named module constant
and records it explicitly.

## Data flow and lifecycle

1. Build the factor only after the bounded real-space H/S operator is complete
   and the fixed-H interface lambda has reached its comparison value.
2. Bind validity to both operator epoch and layout fingerprint because this
   correction depends on H as well as S.
3. Before every application, collectively validate factor validity, Ritz
   eigenvalue count and finiteness, RHS shapes, and factor metadata on the
   bounded-operator communicator.
4. Apply the fragment-local spectral correction independently on every owner.
5. Project the assembled correction against the current Ritz vectors in the
   global S metric before it enters the three-block Rayleigh--Ritz search.
   This removes components already represented by the current Ritz subspace.
6. Preserve the existing retained-search update, convergence test, Ritz
   consistency diagnostics, and publication boundary.

The new route is mutually exclusive with the diagonal and metric-block
preconditioners. The normal production algebra remains callback-free. The
S-orthogonal complement control remains independent and off during the first
physical correction gate.

## Failure handling

Initialization and application fail closed and collectively for:

- nonfinite or non-Hermitian local H/S blocks;
- non-positive-definite or cutoff-failing `S_f`;
- generalized eigensolver failure;
- invalid regularization;
- rank-disagreeing dimensions or state counts;
- nonfinite Ritz eigenvalues, residuals, or corrections;
- stale operator epoch or layout fingerprint.

No partial-rank fallback, silent mode deletion, or unpreconditioned fallback is
allowed after the new mode has been selected.

## Diagnostics

Log once per fragment root:

- correction mode and regularization;
- factor epoch/fingerprint and dimension;
- S condition estimate;
- minimum and maximum generalized local eigenvalues.

At selected inner iterations 32/96/160, log occupied and extra:

- raw and corrected residual norms and their ratios;
- minimum denominator magnitude before and after flooring;
- number and residual weight of floored spectral components;
- correction component removed by the global S-orthogonal projection;
- search-metric effective rank and existing Ritz boundary defects.

## Testing

Add a two-rank dense oracle with nonidentity SPD `S_f`, noncommuting Hermitian
`H_f`, Ritz values on both sides of the local spectrum, and one near-pole
denominator. Verify:

- generalized spectral factors and S-orthonormal eigenvectors;
- equality to the dense regularized correction;
- signed-floor behavior and finite norm bounds;
- identical regularization for occupied and extra labels;
- global S-projection orthogonality;
- continued validity for no-op reads and collective stale rejection after
  H/S epoch or layout changes;
- collective rejection of malformed/nonfinite inputs.

Source contracts verify default-off plumbing, mutual exclusion, fixed-H and
continuation-only routing, no normal-production callback, and retained-search
propagation. Run focused MPI tests, the full MPI/EigenExa build, code review,
then a fresh eight-rank B=6 comparison against Tasks 16, 19, and the
S-orthogonal gate.

## Physical gate and interpretation

Keep the Task 16 basis, seed, cutoff, tolerance, history, continuation, and
publication settings. Disable the diagonal, metric-block, and S-orthogonal
comparison modes; enable only the new correction.

Promotion requires improvement of both occupied and extra residuals without
loss of Ritz consistency, search-metric rank collapse, or publication-safety
regression. Occupied-only improvement with extra degradation rejects the local
approximation, not the general `H-epsilon S` correction equation. If rejected,
the next design may use a global iterative correction solve or an explicitly
state-partitioned policy, based on the spectral-weight diagnostics.
