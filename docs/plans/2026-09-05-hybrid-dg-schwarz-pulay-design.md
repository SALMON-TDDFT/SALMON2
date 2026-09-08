# Hybrid DG Schwarz + Pulay Design

## Goal

Run the divided Hybrid ground-state iteration in the complete DG operator
from its first density epoch without repeatedly diagonalizing the complete
LCFO Hamiltonian.  Keep one MPI rank per fragment, exchange only neighboring
DG interface data during Hamiltonian and metric applications, use a bounded
number of CG updates, and perform exactly one complete LCFO generalized
eigensolve after the divided density iteration.

The divided Hybrid density mixer uses SALMON's existing history-capable
mixing implementation without changing the compatibility fingerprint of the
reused ordinary-DC seed.

## Current gap

The current production payload contains the fragment interior, nonlocal, and
SIPG interface terms, but the divided CG callback extracts only each
fragment's self block.  SIPG self contributions are present; off-diagonal
neighbor-fragment contributions are deferred until the terminal LCFO solve.
Consequently the initial divided iteration is block local rather than a full
DG iteration.

The current Si64 fixture also selects `method_mixing='simple'`.  Although the
adapter can call Pulay or Broyden, changing the ordinary method changes the DC
seed convergence fingerprint and prevents reuse of the already-published DC
seed.

## Distributed state representation

After the existing fragment-local WF+PW admission, form a common global trial
column inventory.  Its size is derived from the electron count, 300 K
occupations, degeneracy completion, and the configured guard policy.  It is
never a material-specific fixed constant.  A fragment that has fewer admitted
directions extends from its ordered DC/PW candidate catalog before the common
inventory is published.

Rank `f` owns only the coefficient rows belonging to fragment `f`:

```text
C = [C_1; C_2; ...; C_N],   C_f(nb_f, ntrial)
```

Every rank has the same `ntrial` and column IDs.  The common IDs bind a state
column across fragments; local fragment eigenstate ordinal alone is never
used as that binding.  Publication is collective and transactional.

## Schwarz DG operator application

The immutable production directory already identifies the owner and fragment
of every basis row and every SIPG face coupling.  Build a neighbor schedule
from that directory once per basis generation.

For each H or S application:

1. apply the local interior, local/nonlocal, and SIPG self blocks to `C_f`;
2. exchange only coefficient blocks needed by faces adjacent to fragment `f`;
3. apply the off-diagonal face rows to the received neighbor blocks;
4. return the result for the locally owned rows.

No complete coefficient array is all-gathered.  Non-neighbor communication in
H/S application is an error.  The neighbor manifest, basis generation,
row-owner fingerprint, face fingerprint, state-column fingerprint, and MPI
rank-to-fragment mapping are checked before every published epoch.

The update is block-Jacobi/Schwarz: all ranks read neighbor coefficients from
the accepted start of the CG step and publish their updated local row blocks
only after collective validation.  A failed rank rolls back the whole step.

## Local iterative solve

Use the existing preconditioned fragment CG machinery on the distributed
row blocks, capped by `dg_hybrid_fragment_cg_steps` per density epoch.  The
default remains three and convergence may stop after one or two steps.  The
algorithm uses global reductions only for small state-space quantities such
as `C^dagger S C`, Rayleigh energies, residual norms, and occupation totals.
It does not diagonalize the complete DG basis Hamiltonian during the divided
SCF loop.

Columns are globally S-orthonormalized after each accepted step.  At 300 K,
global Rayleigh energies determine one common chemical potential and the
number of occupied/guard columns.  If the upper occupied tail is unresolved,
all ranks extend the common inventory transactionally from their local
DC/PW candidates.  Capacity exhaustion is a hard diagnostic, not a hidden
fallback.

## Density and mixing

Add a divided-Hybrid-specific input control:

```text
dg_hybrid_divided_mixing = 'pulay'
```

Allowed values are `inherit`, `simple`, `pulay`, and `broyden`; the divided
Hybrid default is `pulay`.  `inherit` reproduces `method_mixing`.  This control
belongs to the post-seed Hybrid algorithm and therefore does not alter the
ordinary DC seed compatibility contract.  The existing exact MPI count and
rank-fragment mapping requirements remain unchanged.

Reuse SALMON's existing `copy_density`, `pulay`, `wrapper_broyden`, and simple
mixing implementations and their existing parameters.  History persists
between density epochs.  It is reset only when the basis generation or common
state inventory changes, or after an explicitly diagnosed rollback.  The log
records the selected method, history length/reset reason, mix rate, density
convergence value, electron defect, CG steps, and basis extensions.

## Terminal LCFO

After divided-density convergence, refresh the total potential once, assemble
the final complete H and S rows, and perform exactly one generalized LCFO
eigensolve.  The LCFO occupied subspace remains authoritative for symmetry,
final density/state acceptance, and the real-time initial state.  No density
update follows LCFO unless the user explicitly selects the already-defined
optional short global refinement mode.

## Error handling

Reject collectively and preserve the last accepted state for:

- unequal MPI rank and fragment counts or changed rank-fragment mapping;
- inconsistent common state-column inventory;
- missing, duplicate, or non-neighbor face communication;
- stale basis, metric, face, operator, or mixing-history fingerprints;
- non-finite H/S action, loss of S rank, or failed global orthogonality;
- unresolved 300 K occupation tail with no admissible extension;
- rank-local update or history-mixing failure.

There is no automatic switch back to self-block-only CG, simple mixing, or a
repeated complete diagonalization.

## Verification

Tests must cover:

1. two unequal fragments with different local basis sizes but one common
   state-column inventory;
2. H/S action equality between neighbor-only Schwarz exchange and an explicit
   assembled DG matrix on 2, 4, and 8 ranks;
3. rejection of missing, duplicate, non-neighbor, and stale-generation data;
4. one-to-three CG steps per density epoch and transactional rollback;
5. Pulay history persistence and reset on basis extension, with the ordinary
   DC seed fingerprint unchanged;
6. 300 K global electron count and dynamic tail extension;
7. exactly one terminal LCFO solve and no post-LCFO density update by default;
8. Si64 comparison of simple and Pulay density-convergence histories while
   reusing the same compatible eight-rank DC seed.

