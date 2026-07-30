# Direct pseudo-atomic s+p complement design

## Problem

Projecting the complete pseudo-atomic s+p shell into a finite fragment
eigenstate window leaves the retained Wannier space strongly dependent on
the empty-state cutoff.  This remains true at the computed
`nstate_frag=400` limit and is insensitive to whether the radial seed is
a pseudo-atomic orbital or a nonlocal projector.

## Design

For each buffer-periodic fragment, append the complete pseudo-atomic s+p
values and their periodic finite-difference gradients to the fragment
candidate block.  Extend the core-owned occupied coefficient matrix with
zero rows for these appended functions.  The existing construction then
forms the direct sum of:

- the occupied fragment subspace represented by fragment eigenstates;
- the complete raw pseudo-atomic s+p complement represented by the
  appended buffer-periodic functions.

Passing the same raw s+p values as `projection_seed_values` remains useful:
the existing rank, direct-sum, localization, symmetry, center, and
inclusion gates verify the augmented representation.  Because the raw
shell is contained exactly in the augmented candidate block, adding more
empty fragment eigenstates cannot change the retained occupied-plus-s+p
space.

The configured candidate-window count continues to mean the number of
fragment eigenstates.  The internal construction rank is that count plus
the complete-shell channel count.  Evidence records both quantities.

No change is made to normal DC, LCFO, EigenExa, WPW, normal checkpoint
publication, conventional RT, or Hamiltonian pseudopotential operators.

## Validation

The RED fixture constructs the same occupied-plus-s+p problem with two
different numbers of irrelevant empty candidates and requires identical
retained projectors.  Production source contracts require raw seed values
and periodic gradients to be appended and occupied coefficients to be
zero-extended.

After focused 1/2/4/8 MPI verification and clean overlay build, fresh
Si64 c192 and c256 rows must pass their individual gates and checkpoint
reuse.  The unmodified public `1e-4` window gate decides Task 9.
