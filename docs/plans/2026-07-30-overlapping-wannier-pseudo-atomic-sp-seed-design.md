# Pseudo-atomic s+p seed design

## Context

The buffer-periodic overlapping-Wannier route retains 16 occupied
functions and the 32-dimensional complete s+p complement in every Si64
candidate-window diagnostic.  Candidate windows from 192 through the
fragment-state limit of 400 nevertheless change the metric and total
energy substantially.  The current seed adapter uses the pseudo-atomic
orbital table only for a local-reference channel and substitutes the
short-ranged nonlocal projector table for other channels.

## Decision

Every complete s+p seed channel uses the pseudo-atomic orbital
`pp%upptbl_ao(:,l,species)` on `pp%rad`.  The nonlocal projector
`pp%udvtbl` remains part of the Hamiltonian and operator fingerprint, but
is not used as a Wannier projection seed.

The manifest remains one s and three p functions per core-owned atom.
Periodic-image evaluation, symmetry closure, occupied inclusion, buffer
support, and the fixed 48-function fragment target are unchanged.  The
normal DC, LCFO, EigenExa, WPW, checkpoint-publication, and conventional
RT routes are unchanged.

## Validation

A focused numerical fixture distinguishes pseudo-atomic orbitals from
nonlocal projectors and requires the seed adapter to select the former
for both local and nonlocal s+p channels.  Existing projection,
construction, solver, SCF, row-storage, and forbidden-route contracts
must pass.

The clean parent-prerequisite overlay then reruns fresh Si64 c192 and
c256, buffer-5, 2x2x2 rows.  Task 9 passes only if the public raw-evidence
checker accepts the required window gate.  Otherwise the measured
failure is recorded and Task 10 remains blocked.
