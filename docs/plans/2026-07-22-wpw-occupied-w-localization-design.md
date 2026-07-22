# Occupied-W Localization Diagnostic Design

## Purpose and scope

Determine whether the 128 runtime projected occupied Wannier functions are
broader than the earlier Wannier90 result (`sqrt(Omega) ~= 1.2 Angstrom`).
The first implementation is diagnostic only.  It does not rotate W, alter the
buffer gate, relax tolerances, or publish a production checkpoint.

The physical preflight is the current orthorhombic Si64 Gamma-point supercell.
Unsupported cell geometry is rejected explicitly rather than evaluated with
an orthorhombic approximation.

## Wannier90 spread contract

Use the discrete Marzari-Vanderbilt/Wannier90 Gamma-point per-W spread
functional, not a minimum-image RMS surrogate.  For each fragment-owned W,
assemble the diagonal projected neighbor links

`M_nn(b) = <W_n | exp(-i b.r) | W_n>`.

For the orthorhombic Gamma-only Si64 supercell use the six nearest reciprocal
neighbors `b=+/- (2 pi/L_a)e_a`.  Determine weights from
`sum_b w_b b_alpha b_beta = delta_alpha_beta`, giving
`w_(+/-a)=1/(2 |b_a|^2)`, and validate this identity numerically before use.
This is the shell written by the existing local Wannier90 `nnkp` path; check
its `G=+/-e_a` ordering when available.  Compute total per-W
`Omega [Angstrom^2]`, center, and `sqrt(Omega) [Angstrom]` from the diagonal
links.  Do not report `Omega_I`/`Omega_OD`: their full 128-W decomposition
would require cross-fragment link blocks.

Let `q_n(b)=Im log(M_nn(b)/N_n)` on the principal branch and
`N_n=<W_n|W_n>`.  With the +/- shell above, use the Wannier90 sign convention

`rbar_n = -sum_b w_b b q_n(b)`

and

`Omega_n = sum_b w_b {1-|M_nn(b)/N_n|^2 + [q_n(b)+b.rbar_n]^2}`.

Report `N_n`, require `abs(N_n-1)<=1e-8`, and use the normalized link in the
formula so harmless common scaling cannot change the measured width.  A norm
outside the tolerance is a diagnostic failure, not silently repaired input.
The +/-b links must satisfy conjugate-pair consistency; a single fixed-b link
is not Hermitian in general.

A minimum-image real-space RMS and radial cumulative norm may be printed as
separately named tail diagnostics, but neither may be labelled Wannier90
spread or used for the 1.2 Angstrom comparison.

## Data availability and periodic ownership

Each fragment owns 16 W rows and P contains their values.  For the Si64 cases,
the core has 12 points per direction and B>=6, so P spans at least one complete
24-point physical cell.  For each fragment-owned W and each canonical physical
grid point, choose exactly one P preimage in a deterministic half-open
24-point window.  If P contains multiple images, compare their values within a
named periodic-consistency tolerance and never sum them.  If P does not cover
every canonical point exactly once, fail collectively.

Exactly the fragment representative (`dc%id_frag==0`) evaluates the 16
diagonal link sets from its replicated P; other ranks contribute zero and no
fragment-communicator sum of the replicated P is performed.  Only the 16
per-W results per fragment and scalar diagnostics are collected globally, so
no global `128 W x grid` gather is introduced.  D is not the integration
source because spread needs the whole canonical cell and P already contains
the required values.

## Output and stopping behavior

Report per W: stable global ID, fragment owner, center, total `Omega_A2`, and
`width_A`.  Report global minimum,
mean, maximum, median, 90th percentile, count above 1.2 Angstrom, and the
widest W.  A vanishing first harmonic is reported as an undefined center and
a localization failure diagnostic; it does not erase the observable data
already accumulated.

The diagnostic runs before `buffer_sufficiency`.  Existing outer-shell failure
still stops fixed-H/publication/RT.  The measured width distribution is a
review checkpoint, not an automatic trigger based on an invented threshold.

## Conditional localization design boundary

If the exact diagnostic confirms excessive spread, create a second reviewed
plan for a deterministic 16-by-16 fragment-local unitary localization.  The
second design must define the exact MV objective, reciprocal weights, complex
Jacobi/gradient update, stationarity residual, sweep bound, tie breaking, and
rollback.  It must restrict mixing to equal-occupation blocks and prove
orthonormality, density/projector, owner identity, P values, and analytic
gradient invariance.

A global 128-by-128 rotation is out of scope: it conflicts with fragment row
ownership and would require all W values on every fragment P.

## Verification

Tests first establish the reciprocal-shell spread formula on analytic link
matrices with known centers/spreads.  A periodic-grid fixture crossing a seam
then verifies P-preimage uniqueness, invariance when B grows from 6 to 10, and
rejection of inconsistent duplicate images.  An MPI fixture distributes
fragment owners and must reproduce the serial 16-W block result.  Only after
these tests pass is a fresh Si64 B=6 diagnostic run allowed.
