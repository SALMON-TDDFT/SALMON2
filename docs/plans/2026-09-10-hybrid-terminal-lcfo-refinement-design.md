# Hybrid Terminal-LCFO Refinement Design

## Goal

Reduce formal divided Hybrid GS global generalized eigensolves from one or
more solves at every continuation stage to one terminal LCFO solve in the
normal case, with at most three additional global density refinements.

## Constraints

- MPI rank count remains equal to fragment count, with one fragment per rank.
- The converged conventional DC density is the frozen density throughout the
  fragment-local continuation.
- Kinetic, nonlocal pseudopotential, metric and full SIPG/DG interface terms
  are present from the start and remain fixed.
- Fragment solvers use only a few CG steps per local update and exchange only
  the existing neighboring-fragment interface data.
- The localized WF+PW construction basis remains the RT basis.  Terminal LCFO
  eigenvectors provide the initial coefficient matrix; they do not rotate the
  stored basis into delocalized eigenfunctions.
- Occupations are evaluated at 300 K by the existing occupation policy.
- No spatial-symmetry adaptation is imposed on individual localized basis
  functions.  No unmeasured symmetry receipt is introduced.
- Existing distributed-v5 checkpoint identity, full SHA-256 authentication,
  same-rank/same-fragment reuse policy and Exp RT contracts remain unchanged.

## Selected approach

Use adaptive terminal refinement.  The local continuation performs no global
ScaLAPACK eigensolve and no density update.  After the local basis and fixed
DG operator are ready, solve the complete generalized LCFO problem once,
construct occupations and density, and evaluate convergence.  If needed,
perform at most three additional Hartree+XC/local-potential refinements.  Thus
the total global eigensolve count is one to four.

Alternatives rejected were always performing all three refinements, which
wastes cubic work for a good DC seed, and never refining, which provides no
recovery for systems whose DG density differs materially from the DC density.

## Local continuation

The local phase starts from the converged DC density and immediately uses the
complete DG/SIPG coupling (`lambda=1`).  Each fragment performs the configured
small fixed number of CG steps, nominally three, then exchanges its interface
state with neighboring fragments.  Acceptance is based on local orbital,
metric, interface and electron-count diagnostics.  It does not require a
whole-system eigensystem, reconstructed global density or global spectral
symmetry diagnostic.

The local phase must preserve the DC density bit-for-bit.  Potential updates
needed by the local fragment solver are therefore made from that immutable
density.  The existing continuation rollback and finite-value guards remain,
but a rollback never triggers a global diagonalization.

## Terminal LCFO and adaptive refinement

After local continuation, assemble the distributed complete-system metric and
Hamiltonian and invoke the ScaLAPACK generalized eigensolver once.  Apply the
existing 300 K occupation policy and reconstruct the distributed density.
Record orbital residual, metric defect, projector defect, electron defect,
density change and total-energy change.

If the density and energy criteria are not satisfied, update only the
multiplicative local Hartree+XC potential using the existing history-aware
density mixer.  The metric, kinetic, nonlocal pseudopotential and SIPG
components remain immutable.  Recompose the Hamiltonian, solve LCFO again and
repeat for at most three additional solves.  A converged result exits early.

The density criterion uses `dg_dc_gs_final_density_tolerance`; orbital and
electron gates use `dg_dc_gs_final_orbital_tolerance` and
`dg_dc_gs_electron_count_tolerance`.  The energy change is recorded and uses
the existing final-energy tolerance where available; the implementation must
not invent a second independent user tolerance if an established one exists.

## Exhausted refinement

Failure to satisfy the density or energy criterion after three additional
solves is nonfatal by explicit policy.  The last finite, contract-valid LCFO
state is published.  The terminal acceptance receipt records:

- total global eigensolve count;
- additional refinement count;
- final density and energy changes;
- orbital, metric, projector and electron defects; and
- `refinement_converged=false`.

The distributed-v5 shard format itself remains byte-compatible: its fixed
eight-element `acceptance_receipts` array is already fully assigned to the
orbital/metric/projector/electron and spectral-rank contracts.  Refinement
metadata is therefore written to a separately authenticated versioned receipt
whose digest is bound to the v5 publication fingerprint.  Changing the v5
array extent or overloading one of its existing fields is forbidden.

A named warning is printed once by rank zero.  RT does not reject a checkpoint
solely because this flag is false, but all existing v5 validation, metric,
finite-value, energy-identity and zero-field stationarity checks remain active.

## Complexity target

Let `R` be the total construction rank and `P` the fragment/rank count.  The
current continuation performs `D` complete solves, each with leading cost
`O(R^3/P)`.  The new route uses one to four solves, normally one.  Local work
remains fragment-sized.  This change does not yet remove the `O(R^2/P)` dense
ScaLAPACK storage, the GS publication `O(R*g)` basis/gradient storage, or the
RT `O(r*Nocc)` occupied-column storage; those are separate optimizations.

## Validation

TDD must demonstrate:

1. a static and executable gate fails while a global solve remains in the
   local continuation loop;
2. a good DC seed produces exactly one global solve;
3. synthetic density histories produce two, three and four total solves;
4. no path exceeds three additional solves;
5. exhausted refinement publishes a warning receipt rather than aborting;
6. the local phase preserves the DC density and all fixed operator components;
7. only the local-potential component changes during global refinement;
8. 300 K occupations and electron count remain correct;
9. distributed-v5 publication and Exp RT startup remain valid on 1/2/4/8-rank
   focused tests; and
10. the current H4 production handoff, fresh Si8 reuse matrix, preserved Si64
    analyzer and conventional bulk-Si GS/RT remain valid.

Performance receipts must report local iterations, global solve count,
additional refinements, early-exit reason and the final convergence metrics so
that the expected cubic-work reduction is directly auditable.
