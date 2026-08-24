# WF+PW LCFO Divided-SCF Design

## Objective

Reduce the repeated cost of the Hybrid ground-state calculation without
inventing a second SCF model.  Extend the existing divide-and-conquer (DC)
ground-state method from its fragment orbital space to a Wannier-function plus
windowed-plane-wave (WF+PW) local space, then perform the LCFO reconstruction
and one distributed full-system diagonalization after the DC density has
converged.

The production scope is restricted to materials with a finite band gap.  Metals
are out of scope.  The first physical acceptance system is Si64 with eight MPI
ranks and `OMP_NUM_THREADS=1`.

## Governing Decision

Preserve the established DC/LCFO separation.

1. The DC stage determines the self-consistent density from fragment solves on
   core-plus-buffer domains.
2. The LCFO stage uses the converged DC potential to construct one coherent
   full-system coefficient-space eigenproblem.
3. LCFO is diagonalized once.  Its density is not returned to the SCF loop and
   no new post-LCFO density-convergence gate is introduced.

This is the same accuracy model as conventional DC+LCFO.  The new work changes
the local basis from the existing fragment orbital basis to WF+PW; it does not
replace DC convergence by a new fragment-stitching SCF.

## Alternatives Considered

### Repeated full-system Hybrid diagonalization

This gives a direct reference fixed point, but repeats an approximately cubic
global eigensolve on every SCF iteration and retains quadratic distributed
matrix storage.  Keep the converged Si64 result as a validation oracle, not as
the scalable production route.

### Sparse density-matrix purification or Fermi-operator expansion

For gapped systems this can eventually give linear scaling.  In the present
nonorthogonal WF+PW metric it adds chemical-potential, truncation, and metric
purification machinery before the simpler DC extension has been validated.
It is deferred.

### Selected design: DC-SCF followed by WF+PW LCFO

This reuses the mature DC density construction, occupation handling, mixing,
buffer semantics, and convergence controls.  Only the local retained space and
the subsequent LCFO basis are generalized to WF+PW.  One final distributed
eigensolve is accepted, as it is in the existing real-space LCFO method.

## SCF Architecture

Each fragment communicator owns its existing real-space core and buffer.  The
buffer must continue to cover the finite-difference stencil and pseudopotential
projector support required by the DC calculation.  It is an operator halo, not
an additional owner of physical density.

Within that domain, construct a bounded WF+PW local basis:

- retained fragment states supply the WF sector;
- windowed PW packets supply the complementary sector;
- projection against the WF sector and metric rank filtering remove linear
  dependence;
- basis selection is fixed before the density iteration and is not changed by
  SCF history.

Every SCF iteration then follows the existing DC route:

1. apply the current total-system potential restricted to each fragment;
2. solve the fragment problem in its fixed WF+PW local space using the existing
   DC eigensolver/CG policy;
3. determine occupations and the common chemical potential through the existing
   DC electron-count procedure;
4. contribute only the authoritative core density from each fragment;
5. assemble `dc%rho_tot` and update Hartree, XC, local, and nonlocal terms using
   the existing DC infrastructure;
6. mix and test convergence using the existing DC settings.

No full-system Hybrid coefficient matrix or full-system real-space orbital
array is materialized during this loop.

## Convergence and Accuracy

Use the existing SALMON `convergence` and `threshold` inputs without adding a
Hybrid-specific density tolerance.  In the standard DC configuration,
`convergence='rho_dne'` tests

\[
  N_e^{-1}\int |\rho_i(\mathbf r)-\rho_{i-1}(\mathbf r)|\,d\mathbf r.
\]

The existing `norm_rho`, `norm_rho_dng`, and other supported choices retain
their current meanings.  Existing DC mixing, chemical-potential handling, and
bounded local CG behavior remain authoritative.

Accuracy is controlled in the same way as DC+LCFO:

- converge fragment and buffer dimensions;
- converge the retained fragment-state count;
- converge the windowed-PW cutoff or packet count;
- compare energy, gap, density, and occupied projector against the existing
  full-system Si64 Hybrid reference during development.

The reference comparison is a validation study.  It does not become a new
runtime post-LCFO density gate.

## WF+PW LCFO Stage

After DC convergence, freeze the converged density and potential.  Construct
the LCFO overlap and Hamiltonian matrices in the combined WF+PW basis.  Basis
indices have deterministic center/fragment ownership.  Matrix elements are
evaluated from the owning fragment and the required buffer or neighboring
support, including kinetic, local, nonlocal pseudopotential, and DG boundary
terms exactly once.

Send owned matrix blocks directly into the distributed eigensolver layout.
Do not replicate complete `H` or `S` on every rank and do not generate a
full-system real-space orbital array.  Solve the generalized eigenproblem once:

\[
  H C = S C \varepsilon.
\]

Use EigenExa when the accepted Gamma route is represented as a real problem.
Use the existing complex distributed ScaLAPACK backend when the retained WF+PW
gauge is genuinely complex.  This backend choice does not change the one-shot
LCFO semantics.

Keep coefficients distributed.  Reconstruct density diagnostics and TDDFT
initial states in bounded orbital and spatial tiles.

## Validation and Failure Semantics

The WF+PW-specific runtime checks are limited to conditions needed for a valid
generalized eigenproblem and state publication:

- finite basis and operator values;
- acceptable local and global overlap rank;
- Hermitian `H` and `S` within the existing numerical tolerance policy;
- generalized residual `HC-SC epsilon` and metric orthogonality `C^HSC-I`;
- requested electron and state counts;
- unchanged basis, window, packet, metric, and ownership provenance;
- collective agreement before checkpoint publication.

A failed check produces no TDDFT checkpoint.  Do not hide a failure by adding
material-specific mixing, an empirical boundary correction, or a repeated
global SCF loop.

## Scaling

For fixed physical fragment and buffer sizes, the divided WF+PW SCF has linear
aggregate fragment work; the total-system Hartree operation may retain its
`O(N log N)` cost.  Persistent local-basis and density storage is linear.

The final distributed LCFO diagonalization remains approximately `O(M^3)` in
work and `O(M^2)` in distributed matrix storage for total WF+PW basis size
`M`.  This is explicitly accepted for the first production implementation.
The operator assembly and ownership interfaces must nevertheless remain block
sparse so that a later sparse iterative eigensolver can replace only the final
backend without changing DC-SCF physics.

## Test Strategy

Development follows TDD in this order:

1. synthetic fragment tests for WF/PW complement rank and deterministic basis
   ownership;
2. tests showing the new driver reproduces the existing DC density/convergence
   behavior when the PW sector is empty;
3. local WF+PW fragment Hamiltonian and density tests against direct bounded
   references;
4. MPI tests for exactly-once core density, shared electron count, and
   decomposition invariance;
5. LCFO `H/S` block assembly tests, including neighboring projector and buffer
   support;
6. one-shot distributed eigensolver tests for residual and metric
   orthogonality;
7. Si64 comparison with the converged full-system Hybrid-SCF oracle using eight
   MPI ranks and one OpenMP thread;
8. buffer and PW-space convergence studies;
9. zero-field TDDFT propagation from the accepted WF+PW LCFO state.

Si8 may be used for algebraic fixtures but is not a physical acceptance gate.
No time-based cutoff is used for the Si64 validation.
