# Hybrid Divided-SCF and One-Shot LCFO Production Design

## Objective

Make the scalable Hybrid ground-state route follow the established real-space
DC+LCFO accuracy model:

1. converge the density with fragment-local WF+PW solves;
2. include only bounded neighboring interface/projector communication during
   that local stage;
3. freeze the converged density and solve the complete distributed LCFO
   generalized eigenproblem once;
4. optionally perform a small, user-requested number of global LCFO density
   refinements when higher self-consistency is wanted.

The repeated complete-Hybrid continuation remains available as a validation
oracle.  It is not the default production algorithm.

This design refines the earlier
`2026-08-24-wpw-lcfo-divided-scf-design.md`.  It incorporates the later
localization-first, energy-window certification, complete-v3 checkpoint, and
strict DC-seed compatibility work.

## Governing Decisions

### Default route

`yn_dg_hybrid_divided_scf='y'` selects the production route.  Its default
post-fragment behavior is one complete LCFO eigensolve.  The LCFO density is
not fed back automatically.

### Optional refinement

Add an integer control `dg_hybrid_lcfo_refine_steps`, defaulting to zero.
After the first LCFO solve, each requested refinement step:

1. reconstructs the physical LCFO density on uniquely owned core points;
2. records its difference from the potential-generating density;
3. mixes it using the existing DC density-mixing policy;
4. updates the total Hartree/XC/local potential;
5. reassembles the density-dependent Hamiltonian rows; and
6. performs one new complete LCFO solve.

Thus zero requests exactly one complete eigensolve, while `N>0` requests
exactly `N+1`.  There is no implicit automatic fallback from the one-shot mode
to repeated global SCF.  A user who requests refinement accepts its extra
cost.  The final published eigensystem must correspond to the final updated
potential epoch.

### Reference route

`yn_dg_hybrid_continuation_scf='y'` keeps the current repeated complete
eigensolve algorithm as an oracle for accuracy and convergence studies.  It
must not be entered by the divided route unless explicitly selected by the
user.

## Complete Production Sequence

1. Load a conventional DC seed only when its MPI rank count and exact
   rank--fragment mapping match the current calculation.  Preserve this rule
   in the production format.  If no compatible seed exists, run conventional
   DC once and publish a new seed.
2. Run conventional LCFO and unconstrained localization to obtain the WF
   sector.  Individual WFs are optimized for locality and are not required to
   transform symmetrically.
3. Build the user-cutoff-controlled windowed-PW complement.  The retained
   basis count is derived from the material, fragment size, requested energy
   window, metric rank, and PW cutoff; no material-specific count such as 384
   is permitted.
4. Freeze the fragment WF+PW catalogs, ownership, overlap metric, broken-volume
   kinetic/nonlocal data, and DG interface payload.
5. Run the divided density SCF described below.
6. Freeze its converged density and potential, assemble the complete
   distributed LCFO `H/S`, and solve once.
7. If requested, perform exactly the specified number of LCFO refinement
   steps.
8. Certify the physical occupied/energy-window eigenspace, localize only the
   certified RT space by a unitary gauge change, and publish the complete-v3
   ground-state checkpoint.

## Divided Density SCF

### Fragment operator

Each fragment solves only in its fixed local WF+PW space.  The local operator
contains the fragment self block of

\[
H[\rho] = T_{\rm broken} + V_{\rm NL} + V_{\rm local}[\rho]
          + H_{\rm SIPG}.
\]

The finite-difference/buffer action, DG face data, and nonlocal-projector
support may require data owned by neighboring fragments.  Communication is
limited by the accepted face-neighbor and projector-support graphs.  A
nonlocal projector is not assumed to stop at the nearest Cartesian face.

Off-diagonal fragment hybridization is not diagonalized globally inside this
loop.  Its complete matrix elements are retained for the final LCFO solve.
This is the deliberate DC+LCFO approximation.  Turning the local loop into a
globally coupled occupied-subspace iteration would define a different
algorithm and is out of scope for this production mode.

### Global operations that remain global

The following operations are not made fragment-private:

- Hartree/Poisson update;
- common chemical potential and occupation determination;
- total electron-count reduction;
- unique-core density assembly and density mixing;
- density convergence reduction; and
- final LCFO solve and physical-space certification.

Each fragment solve publishes its eigenvalues and core weights to the common
occupation routine.  A separate chemical potential per fragment is forbidden.
At finite temperature the state inventory must include enough unoccupied tail
states to make the common occupation and electron count well defined.

### Density ownership and electron count

Buffers are operator halos and never own physical density.  Every physical
core point has exactly one owner.  The assembled density must be finite and
must satisfy the total electron count at every iteration within the existing
DC electron tolerance.  Failure is collective and fail-closed; the driver must
not ignore an electron count returned by a density callback.

### Convergence semantics

The divided route uses the exact same formulas, normalization, grid-volume
factors, and thresholds as the authoritative conventional DC implementation.
The formulas must be extracted into a shared helper rather than copied with
similar names.  In particular, the present divided-driver interpretations of
`rho_dne`, `norm_rho`, and `norm_rho_dng` are not authoritative and must be
replaced by the shared DC calculation.

Failure to converge within `nscf` is fatal for this route.  No final LCFO state
or RT checkpoint is published from an unconverged divided density.

## Final LCFO and Optional Refinement

The converged divided density defines a frozen potential.  Assemble all
row-distributed WF+PW matrix elements, including broken-volume kinetic, local,
nonlocal-projector, and SIPG self/cross-fragment contributions, exactly once
per potential epoch.  Do not replicate full `H` or `S` on each MPI rank and do
not materialize full-cell orbitals for all states.

Solve

\[
H C = S C \varepsilon
\]

with the existing distributed generalized eigensolver.  Derive occupations
and the HOMO from this complete spectrum.  Starting from the user energy
cutoff, extend the retained physical window through the first higher
degenerate/symmetry-connected cluster required to close the LCFO eigenspace.
This retained rank is dynamic.

For each optional refinement epoch, the fixed kinetic, nonlocal, metric, face,
and ownership payloads remain immutable.  Only density-dependent potential
rows and derived receipts are refreshed.  The implementation records

\[
\Delta\rho_i = \|\rho_{\rm LCFO}^{(i)}-\rho_{\rm potential}^{(i)}\|
\]

using the shared DC density norm.  `Delta rho` is a diagnostic in the default
DC+LCFO accuracy model, not an automatic gate or hidden iteration trigger.

## Symmetry and RT Semantics

Neither individual localized WFs nor the complete construction basis need be
symmetry closed.  Those measurements remain diagnostics only.

Acceptance applies to physical quantities reconstructed from the final LCFO
eigensystem:

- occupied projector;
- the requested energy-window subspace after cluster extension;
- density;
- projected Hamiltonian covariance; and
- final localized RT basis as a unitary gauge of the certified eigenspace.

The final checkpoint contains coefficients, eigenvalues, occupations,
certified rank, symmetry receipts, density/potential epoch, operator and basis
fingerprints, and exact MPI ownership provenance.  RT evolves only inside the
certified space.  A one-shot DC+LCFO state is an approximate self-consistent
initial state in the same sense as conventional DC+LCFO; zero-field drift and
`Delta rho` must be reported.  Users needing lower drift select explicit
refinement steps.

The checkpoint distinguishes the potential-generating density from the
LCFO-reconstructed orbital density.  The final coefficients and eigenpair
residual belong to the former; the latter supplies the reported physical
density and the next optional refinement input.  One-shot mode must not label
the two densities as equal or silently recompute a different Hamiltonian while
retaining the old eigenpair receipt.

## Failure and Publication Rules

No checkpoint is published unless all of the following hold collectively:

- fragment catalogs, ownership, and fixed payload are finite and immutable;
- divided SCF converged under the shared DC criterion;
- every divided density has the requested electron count;
- final `H/S` are finite, Hermitian, and have an acceptable metric rank;
- `HC-SC epsilon` and `C^HSC-I` pass their tolerances;
- occupations are physical and conserve electron count;
- the certified physical occupied/window space passes its symmetry receipts;
- potential-generating density, LCFO density, coefficients, and operator are
  separately epoch-tagged with no cross-epoch receipt reuse; and
- every complete-v3 fingerprint and ownership field is present.

Construction-basis nonclosure alone is not a failure.  Changing rank count or
rank--fragment mapping makes a conventional DC seed ineligible for reuse; it
must never be silently migrated.

## Existing Code to Reuse

- `dg_hybrid_divided_scf.f90`: callback driver, after convergence and electron
  semantics are corrected.
- `dcdft.f90` / `scf_iteration_dft.f90`: authoritative occupation, chemical
  potential, density assembly, mixing, and convergence behavior; extract
  narrow shared helpers rather than duplicate their algebra.
- `dg_hybrid_fragment_solver.f90`: fragment generalized solve, after it accepts
  the production DG self block and the occupations obtained for the current
  iteration.
- `s_dg_hybrid_fixed_payload` and production face/projector graphs: immutable
  metric, volume, nonlocal, and SIPG data.
- `dg_hybrid_generalized_eigensystem.f90`: complete distributed LCFO solve.
- the localization-first LCFO certification and complete-v3 checkpoint path:
  extract a terminal finalization routine shared by one-shot, refined, and
  reference routes.

The current divided implementation cannot be promoted unchanged: it currently
uses ordinary local `hpsi` with an identity metric, copies stale
`system%rocc`, ignores the returned divided electron count, and implements
three convergence labels with formulas different from conventional DC.

## Test Strategy

Implementation follows TDD in this order:

1. shared pure tests proving bitwise-identical conventional/divided convergence
   metrics for all supported convergence modes;
2. MPI tests for a common chemical potential, current-iteration occupations,
   finite-temperature tails, and exact electron count;
3. fragment-operator fixtures comparing volume, nonlocal-projector, SIPG self,
   and neighbor-support contributions with a direct small reference;
4. unique-core density and rank-decomposition invariance tests;
5. route tests requiring exactly one complete eigensolve for refinement count
   zero and exactly `N+1` for count `N`;
6. final LCFO residual, metric orthogonality, dynamic window/cluster extension,
   physical symmetry, and complete-v3 checkpoint tests;
7. DC-seed tests retaining the exact rank-count and rank--fragment reuse rule;
8. Si64 eight-rank validation against the saved repeated-continuation oracle,
   recording energy, gap, density difference, occupied projector, symmetry,
   time, and memory; and
9. separate-directory zero-field RT runs for one-shot and a small refinement
   count, reporting stationarity rather than silently adding iterations.

All existing dirty changes and verification directories remain untouched.
