# Si64 Hybrid Self-Consistent Ground-State Design

## Objective

Construct the physical ground state after the Wannier plus windowed-plane-wave
hybrid basis has been fixed.  The result must be self-consistent with the density-
dependent SALMON Hamiltonian and must not be inferred from a single generalized
diagonalization.  Si64, not Si8, is the first material gate because the smaller
cell can make otherwise valid windows, Wannier tails, and fragment boundaries
interact artificially.

## Decision

Use two solver stages with one common outer density loop.

1. Establish a reference Si64 calculation with generalized EigenExa inside the
   existing Pulay density-mixing loop.
2. Replace only the inner eigensolver with an adaptive block-CG implementation
   and compare it to the converged reference.

Starting directly with block-CG is rejected because a failure would not identify
whether the inner subspace iteration or the outer density mixing caused the
oscillation.  A one-shot generalized eigensolve is also rejected because the
hybrid density changes the Hartree and exchange-correlation potentials.

## Fixed-Basis SCF Data Flow

The hybrid basis catalog, window functions, selected Wannier blocks, projected
PW packets, sparse graph, and overlap matrix are fixed before this SCF starts.
They must not be reselected in the density loop.

Each outer iteration performs:

1. start from the mixed full-cell density;
2. update the SALMON Hartree, XC, local, and nonlocal Hamiltonian state;
3. rebuild the sparse hybrid Hamiltonian in the unchanged basis;
4. solve `H C = S C epsilon` for the occupied subspace;
5. assign the authoritative LCFO-derived occupations and check total electrons;
6. reconstruct the full-cell output density from all occupied hybrid states;
7. measure density, energy, eigensystem, metric-orthogonality, and symmetry
   residuals;
8. update the input density through the existing SALMON Pulay machinery;
9. publish the RT checkpoint only after all convergence gates pass.

The initial density is the converged DC+LCFO density.  The propagated state is a
distributed coefficient matrix `C_owned(nowned,noccupied)` plus occupations,
not the single arbitrary coefficient vector used by the current primitive RT
fixture.

## Reference EigenExa Stage

The first Si64 gate uses generalized EigenExa to remove inner-iteration ambiguity.
It records the generalized spectrum, occupied projector, density, total energy,
electron count, raw and mixed density residuals, Hamiltonian symmetry residual,
and memory receipts on every SCF iteration.  Degenerate eigenvector columns are
never compared directly; comparisons use the occupied `S`-metric projector.

This stage is diagnostic and establishes the target SCF fixed point.  It does not
authorize a dense replicated production Hamiltonian or a new basis selection.

## Adaptive Block-CG Stage

The scalable solver reuses the preceding SCF occupied subspace as its trial block.
It does not use a large fixed iteration count.

The inner tolerance is tied to the current outer density residual.  Early outer
iterations use a loose eigensystem target and a small iteration cap.  The target
tightens only as the density converges.  An inner cycle stops when any of these
conditions holds:

- the occupied generalized residual reaches its adaptive target;
- the residual reduction stalls for a bounded number of iterations;
- residual growth or alternating growth indicates oscillation;
- the small adaptive iteration cap is reached.

Inactive metric rows remain exactly zero.  Each accepted block is reorthogonalized
in the `S` metric, and the occupied projector rather than individual degenerate
vectors is tracked between iterations.

Over-solving is prohibited: the inner residual only needs to be safely below the
error scale currently allowed by the outer density residual.  The implementation
must emit the requested and achieved inner tolerances, iteration count, residual
history summary, and stop reason.

## Pulay Stability

Reuse the existing SALMON density mixing state and history format.  Do not create
a second independent mixing algorithm for the hybrid route.

If the density residual grows excessively or alternates, reject the newest Pulay
history contribution, reset the affected history, and reduce the mixing factor.
An inner solve that fails its acceptance gates must not enter Pulay history.  A
later stable decrease may restore the configured mixing factor gradually, but no
automatic increase occurs without bounded evidence.

## Acceptance Gates

Checkpoint publication requires all of the following:

- density residual and total-energy change below their configured tolerances;
- generalized residual `||H C - S C epsilon||` below its final target;
- `||C^H S C - I||` below tolerance;
- reconstructed electron count equal to the requested count;
- finite occupations, eigenvalues, density, and operator entries;
- Hamiltonian Hermiticity and required symmetry receipts below tolerance;
- unchanged basis, window, packet, complement, metric, and position provenance;
- agreement of the occupied projector, density, and energy with the Si64
  EigenExa reference within the declared diagnostic tolerances.

Failure is collective and leaves no publishable RT checkpoint.

## Testing and Rollout

Focused tests first cover adaptive stopping, residual stagnation, alternating
residuals, Pulay-history rejection, metric-rank deficiency, degenerate occupied
blocks, electron-count rejection, and MPI decomposition invariance.  A synthetic
nonlinear Hamiltonian fixture separates inner-solver and outer-mixing failures.

The material gate is then Si64 with the previously agreed MPI rank count and
`OMP_NUM_THREADS=1`.  Run the EigenExa reference to convergence before enabling
block-CG.  Compare full convergence histories and memory receipts, not only the
last eigenvalues.  Si8 may remain a fast algebraic test but is not the physical
acceptance system.

RT integration remains disabled until the multi-orbital checkpoint and the
self-consistent Si64 ground state pass this gate.
