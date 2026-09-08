# Buffer-Periodic Complete-s+p Wannier Design

## Scope

Replace the fixed retained-Wannier count used by the overlapping-Wannier
ground-state route with a complete atomic `s+p` projection manifest.  The
first production gate remains Si64, Gamma-only, non-SOI PZ, with the existing
DC candidate solve and with all direct-DG, LCFO, EigenExa, WPW, conventional
checkpoint, and real-time fallbacks disabled.

This change does not alter the normal DC LCFO plus EigenExa route.

## Basis dimension

Each core-owned atom contributes one real `s` channel and three real `p`
channels.  A shell may not be truncated:

\[
N_{\mathrm{shell}} = 1 + 3 = 4.
\]

The target dimension is derived from core ownership,

\[
N_{\mathrm{Wann},F}=4N_{\mathrm{atom},F}^{\mathrm{core}},
\qquad
N_{\mathrm{Wann}}=\sum_F N_{\mathrm{Wann},F}.
\]

For Si64 this gives 256 global functions.  For the 2x2x2 decomposition it
gives 32 functions per fragment.  The DC candidate window remains a larger,
independently converged outer space; candidate-window size is not the
published Wannier count.

## Projection manifest

Reuse the established `pseudo_channels` implementation:

- atom-major `(atom,l,m)` ordering;
- complete-shell validation;
- real spherical-harmonic projectors;
- pseudopotential radial/projector information where available;
- projectability and singular-value diagnostics.

The first tier is explicitly limited to `lmax=1`.  The manifest and
fingerprint record `s+p`; a future `s+p+d` tier must be a separate accepted
manifest and may not be reached by silently increasing a numeric target
count.

## Subspace construction

The DC core electron count defines the frozen occupied rank.  For Si64 this
is 16 occupied functions per fragment and 128 globally.  The target space
must contain that occupied projector exactly.

Project the complete `s+p` manifest into the fragment's buffer-supported DC
candidate space.  Remove the occupied component, orthonormalize the remaining
projected directions in the candidate metric, and retain the full-rank
complement required to complete the manifest.  Reject:

- a missing `s` or `p_m` member;
- manifest rank loss;
- a minimum singular value below the configured scale-aware cutoff;
- a target rank smaller than the frozen occupied rank;
- an occupied-inclusion failure.

Do not fill the target by taking an arbitrary prefix of periodic-localizer
eigenvectors.

## Localization, centers, and symmetry

Atomic `s+p` functions are seeds and a coverage contract, not fixed final
centers.  Localize the accepted target subspace with cell-periodic phase
operators on the core-plus-buffer box.  The final functions may form
symmetry-related `sp3` combinations centered on Si--Si bonds.

Construct one representative fragment basis and materialize equivalent
fragments by the authoritative periodic translation map.  Apply the same
mapping to values and gradients.  Verify:

- bond-center orbit closure;
- retained value and gradient streaming closure;
- complete `s+p` representation closure;
- positive global metric with no rejected rank.

The buffer exterior is not folded onto the same fragment core.  Retained
tails are assessed through buffer-width convergence and translated-fragment
closure.

## Ground-state acceptance

The initial comparison uses the same 256-function complete-`s+p` global
basis for buffer widths 6 and 8 and for both required decompositions.
Candidate-window convergence remains a separate outer-space test.

The gate requires converged `S`, `H`, `T`, `V_NL`, density, charge, energy,
generalized eigen residual, `S`-orthonormality, unmixed fixed point, and
route-checkpoint round trip.  Failure blocks coefficient real time; it does
not increase the shell tier automatically and does not invoke another route.

