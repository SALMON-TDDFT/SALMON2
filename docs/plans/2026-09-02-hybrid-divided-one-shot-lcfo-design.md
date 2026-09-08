# Hybrid Divided-SCF and One-Shot LCFO Production Design

## Objective

Make the scalable Hybrid ground-state route follow the established real-space
DC+LCFO accuracy model:

1. construct unconstrained WFs independently from the reusable DC orbitals on
   each fragment, without a preliminary complete LCFO solve or a complete-cell
   construction-WF Wannier90 solve;
2. converge the density with deliberately short, warm-started fragment-local
   WF+PW subspace updates;
3. include only bounded neighboring interface/projector communication during
   that local stage;
4. freeze the converged density and solve the complete distributed LCFO
   generalized eigenproblem once;
5. optionally perform a small, user-requested number of global LCFO density
   refinements when higher self-consistency is wanted.

The repeated complete-Hybrid continuation remains available as a validation
oracle.  It is not the default production algorithm.

This design refines the earlier
`2026-08-24-wpw-lcfo-divided-scf-design.md`.  It incorporates the later
localization-first, energy-window certification, complete-v3 checkpoint, and
strict DC-seed compatibility work.

## Governing Decisions

### Core-centered construction WF selection (user decision, 2026-09-04)

The authoritative amendment is
[`2026-09-04-hybrid-core-center-selection-design.md`](2026-09-04-hybrid-core-center-selection-design.md).
Preserve every generated WF in the raw cache, but admit to a fragment's
production catalog only WFs whose actual centers belong to its half-open
periodic core. Keep the selected columns' buffer data. Project PW against the
selected WF union, certify the actual core metric and seed/support
reconstruction, and use a selected-space fixed-reference preconditioner.
Center selection alone is neither a rank nor an accuracy certificate.

This supersedes the former unconditional all-WF production retention and
zero-padded raw seed-map requirements. Gauge covariance applies inside the
frozen selected space, not to rotations mixing selected and excluded columns.
The raw-cache all-column/reconstruction guarantees remain unchanged. Task 8
main integration waits for the amendment's tests; no production switch has
been made merely by accepting this design.

### Production MPI scope (user decision, 2026-09-04)

The conventional-DC to divided-Hybrid production route uses exactly one MPI
rank per fragment: `MPI size == fragment count`, with a bijective rank--fragment
mapping. Each rank retains all WF+PW columns, local H/S blocks, and CG state
for its fragment. Fragment construction and local updates are parallel across
fragments; there is no intra-fragment MPI orbital or coefficient-column split.
Keep the required inter-fragment DG interface/projector communication and the
established global DC reductions/potential and final LCFO operations.

This scope is chosen to keep the ordinary DC handoff simple. Multi-rank
fragment distribution is not a prerequisite or a planned production feature
in this implementation. Existing generic distributed kernels/tests may remain,
but must not add redistribution requirements to the one-rank handoff. DC seed
reuse still requires the exact MPI rank count and exact rank--fragment mapping;
never silently remap or repartition an incompatible seed.

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

The complete-cell construction-WF Wannier90 path and the dense fragment
generalized eigensolver remain explicit small-system/reference backends.  They
are not automatic recovery paths for a failed production fragment update.

### Construction localization and fragment eigensolver budgets

The production construction WFs are generated once per basis epoch on
`dc%icomm_frag`, not on `dc%icomm_tot`.  Each fragment uses its restored or
freshly converged DC orbitals, buffer-support candidates, and only those local
atomic projector directions that add metric rank.  The retained rank is
dynamic.  Wannier90 performs an unconstrained square unitary rotation of the
whole accepted raw fragment subspace. All columns remain in the raw cache;
the separate production catalog is selected by core-centered ownership under
the September 4 amendment. No gauge matching is imposed between fragments.
Every fragment uses a fragment-specific seed and artifact directory.

The density loop uses a separate control
`dg_hybrid_fragment_cg_steps`, defaulting to three.  This is an update budget,
not an eigenpair-convergence requirement.  A fragment update may finish before
the budget when its intermediate tolerance is reached; reaching the budget is
also a successful outcome when the Ritz values and coefficients are finite,
the metric rank is retained, `C^H S C-I` passes, and the residual has not grown
beyond `dg_dc_gs_allowed_residual_growth` times the residual of the safe
entry subspace under the current operator.  Here the reported residual is the
largest distributed two-norm of `H C_j-epsilon_j S C_j`, divided by
`max(1,abs(epsilon_j),norm2(H C_j))`.  Early success uses
`dg_dc_gs_intermediate_orbital_tolerance`, and metric orthogonality uses
`dg_dc_gs_orthogonality_tolerance`.  A rejected trial restores the best safe
subspace.  Returning that entry subspace without advancement is successful if
it still passes the finite/rank/orthogonality gates; failure to improve is then
handled by the outer density-convergence limit.  The strict block-CG API keeps
its existing converged-or-fail semantics for reference tests.

## Complete Production Sequence

1. Load a conventional DC seed only when its MPI rank count and exact
   rank--fragment mapping match the current calculation.  Preserve this rule
   in the production format.  If no compatible seed exists, run conventional
   DC once and publish a new seed.
2. On every fragment communicator, form a fragment-local candidate
   space directly from the DC `rwf` payload plus accepted buffer/projector
   directions, remove only metric-null directions, and run unconstrained
   fragment-local Wannier90 once.  There is no preliminary complete LCFO
   diagonalization.  Individual WFs are optimized for locality and are not
   required to transform symmetrically. Preserve the raw result and select
   core-centered WFs into a separately fingerprinted production catalog.
3. Build the user-cutoff-controlled windowed-PW complement against the selected
   WF union, verify the actual core metric and project/certify DC seeds in the
   selected WF+PW space. The retained
   basis count is derived from the material, fragment size, requested energy
   window, metric rank, and PW cutoff; no material-specific count such as 384
   is permitted.
4. Freeze the fragment WF+PW catalogs, ownership, overlap metric, broken-volume
   kinetic/nonlocal data, and DG interface payload.
5. Run the divided density SCF described below, retaining each fragment's
   coefficient cache and applying at most the requested number of warm-started
   LOBPCG updates per density iteration.
6. Freeze its converged density and potential, assemble the complete
   distributed LCFO `H/S`, and solve once.
7. If requested, perform exactly the specified number of LCFO refinement
   steps.
8. Certify the physical occupied/energy-window eigenspace.  When an RT
   checkpoint is requested, localize the complete certified RT space by a
   unitary gauge change and publish complete-v3.  This required RT-publication
   step is a distinct whole-system operation; it is not the construction-WF
   path.

## Fragment-Local Construction WFs

The conventional DC seed already stores each rank's real fragment orbitals,
fragment eigenvalues, occupations, density, and potential.  The production
route consumes those local orbitals directly.  It must not call `dc_lcfo` to
manufacture complete-system orbitals merely to localize them again.

For fragment `f`, construct a full-rank local candidate matrix on its core plus
accepted operator buffer.  Candidate directions may include DC occupied and
guard orbitals and fragment-local atomic projectors, but a fixed material count
is forbidden.  A distributed metric decomposition removes null directions and
records the retained rank.  The resulting square subspace is passed to
Wannier90 on `dc%icomm_frag` with a seed namespace containing the fragment ID.
The namespace also contains the basis generation so artifacts from different
epochs cannot collide.  Setup and run occur exactly once per fragment and
basis epoch, outside the SCF loop. All accepted columns survive in the raw
cache; final centers, values and transform ordering must remain aligned.
Production center selection is a separate, loss-checked operation.

Wannier90 receives a square `num_bands=num_wann` space with no disentanglement.
Its required eigenvalue array is therefore an auxiliary finite zero label for
every retained candidate, including appended buffer/projector directions; it
is not a fragment spectrum.  The saved physical DC eigenvalues remain in a
separate payload and only those, followed by later `H/S` Ritz values, may drive
occupation or state-extension decisions.

Buffers from different fragments may overlap and independent fragment gauges
may differ by phases, permutations, or general unitary rotations.  The union
is therefore certified through its complete overlap matrix.  Near-null global
directions are handled by a gauge-invariant metric eigenspace compression, not
by deleting named WFs.  Interface tails, periodic-wrap tails, and nonlocal
projector support must remain present.  Missing support, duplicate ownership,
or insufficient metric rank is a collective failure.

Keep the full raw construction cache in addition to two explicitly related
production catalogs. The divided-density catalog keeps every admitted
core-centered fragment WF and its fragment-assigned PW
complement in the original block layout; local `H_ff/S_ff`, warm-start
coefficients, and rank--fragment ownership always use this uncompressed
catalog.  Separately, assemble the complete union metric.  If cross-fragment
overlap creates numerical null directions, build a rectangular,
gauge-covariant union transform that removes only that null eigenspace and use
it for the terminal complete-LCFO catalog.  When the union has full rank this
transform is the identity and `S` remains a valid nonidentity generalized
metric.  The complete transform may mix fragment columns and therefore is
never fed back into fragment ownership or divided SCF.

Store the original raw seed map, the recomputed DC seed projection into the
selected WF+PW catalog, and the
uncompressed union to the terminal complete catalog.  Composing them must
reconstruct every retained seed orbital through the terminal catalog whenever
null compression is applied.  Projectors, occupied density, and eventual LCFO
observables must be invariant under independent unitary gauges in every
selected fragment block. No invariance is asserted for rotations that change
which raw columns meet the center-selection policy.

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

The fragment eigensolver advances only the occupied-plus-guard subspace, not
all local basis vectors.  Its initial state count comes from the DC spectrum,
occupation tail, requested energy window, and degenerate-shell closure.  If
the common-chemical-potential calculation reports an insufficient high-energy
tail, that fragment state inventory grows monotonically; it never shrinks
during the SCF.  The occupation adapter returns a per-fragment extension mask.
For fragment `f`, its terminal numerically degenerate shell contributes
`q_tail(f)=sum(occupation*core_norm)`; it is marked when
`q_tail(f)>electron_tolerance/n_frag`.  At zero temperature, a terminal shell
that is occupied or intersects the common chemical potential is also marked.
Each marked fragment first adds the complete next unused DC-seed eigenspace,
ordered by the saved seed eigenvalue and closed over numerical degeneracy.  The
seed-to-WF coefficient map makes this selection invariant under the later
fragment Wannier90 gauge.  After the seed eigenspaces are exhausted, form the
next candidate pool from a complete kinetic-energy shell of the windowed PW
catalog plus any still-unused projector-support directions.  S-project the
whole pool against the accepted state space and diagonalize `H/S` only inside
that new pool; never rank individual localized basis columns by their diagonal
`H_ii/S_ii`.  Stable global candidate IDs break ties between otherwise equal
shells, but do not select a direction inside a degenerate pool.  Occupations
are then solved again.  Expansion repeats until the aggregate tail gate
passes; exhausting a fragment's metric rank before that point is a collective
failure.

The previous safe coefficients are reused whenever fragment ID, basis
generation, and metric fingerprint agree.  Equal state counts reuse all
columns.  On expansion, the old columns are embedded unchanged and only new
directions are deterministically S-orthogonalized and appended; the new search
directions start with zero history.  A changed potential epoch does not
invalidate the warm start.

The production updater uses an `[X,R,P]` LOBPCG trial space and performs at
most `dg_hybrid_fragment_cg_steps` updates.  It does not diagonalize the full
fragment basis as a hidden cold start.  Its residual is an intermediate SCF
diagnostic, not a published physical eigenpair receipt.  Exact eigenpair and
symmetry acceptance is deferred to the one complete LCFO solve.

### Fixed-frame fragment preconditioner (2026-09-03 amendment)

The square-map construction below describes the all-retained case. When center
selection removes columns, the September 4 amendment's explicit rectangular
projected-frame API replaces this construction. Do not pass a sliced
`U^dagger` to the existing square-unitary API. The denominator safety, epoch
validation and selected-space physical covariance requirements still apply.

The approved short-update gauge-invariance requirement rules out rebuilding a
coordinate-diagonal preconditioner from the currently localized WF columns.
Freeze the pre-Wannier90 metric-compressed candidate frame for a basis
generation, together with the invariant projected-PW catalog. Do not rerun
metric compression or choose new reference directions after changing the WF
gauge. This reference choice does not constrain Wannier90 localization.

Write the current uncompressed fragment basis as `B` and the fixed physical
reference as `F=B Q`. For the stored Wannier90 rotation `U`, use
`Q=block_diag(U^dagger,I_PW)`, including its stored phase/order corrections.
The buffer/projector directions already retained before Wannier90 belong to
that same reference block. Keep Q separate from the union-to-complete map.
Q is unitary in coefficient coordinates; F need not be orthonormal in the
current broken-volume fragment metric.

At each local operator epoch compute only
`h_a=diag(Q^dagger H_ff Q)` and `s_a=diag(Q^dagger S_ff Q)`. On a raw residual
column r_j apply `z_j=Q D_j^-1 Q^dagger r_j`, where
`D_j(a)=h_a-epsilon_j*s_a` is regularized at a collective scale using the
numerical tolerance. Preserve the sign of resolvable nonzero denominators;
values indistinguishable from zero at roundoff use a deterministic positive
floor. Reject nonfinite arithmetic and nonpositive reference metric diagonals.
This is an approximate inverse, not a local eigenproblem or a dense-solve
fallback. Identity preconditioning is not a production substitute.

For a WF gauge V, `B'=B V`, `Q'=V^dagger Q`, `H'=V^dagger H V`,
`S'=V^dagger S V`, and `r'=V^dagger r` give `z'=V^dagger z`. A diagonal
formed directly from H' lacks this property. The production fixed frame must
be built from the stored pre-localization transform, not replaced by identity
after a gauge change. A physical reference fingerprint stays fixed while the
coordinate-map fingerprint can change. No claim is made about rebuilding the
reference from a different DC seed or a different metric-compression gauge.

The bounded updater passes its current Rayleigh values to an explicit shifted
preconditioner callback. The original residual-only callback remains available
for existing fixtures, with exactly one callback selected and no silent
fallback. Row layout, fragment, generation, basis/metric/reference identity,
operator fingerprint and potential epoch bind each prepared preconditioner.
Stale, rank-disagreeing or invalid requests fail collectively and publish no
partial output. All communication stays inside the fragment communicator.

The input-controlled update budget, dynamically extended state inventory,
one-time construction-WF localization, final one-shot LCFO, and exact MPI
rank-count/rank--fragment restart rule remain unchanged.

### Global operations that remain global

The following operations are not made fragment-private:

- Hartree/Poisson update;
- common chemical potential and occupation determination;
- total electron-count reduction;
- unique-core density assembly and density mixing;
- density convergence reduction;
- final LCFO solve and physical-space certification; and
- unitary localization of the final certified RT eigenspace when an RT
  checkpoint is requested.

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
per potential epoch in the uncompressed union, then congruence-transform them
with the immutable terminal complete-catalog map before solving.  Do not
replicate full `H` or `S` on each MPI rank and do not materialize full-cell
orbitals for all states.

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
Intermediate LCFO epochs do not invoke the certified-RT localizer.  When RT
checkpoint output is requested, that whole-system unitary localization occurs
exactly once, after the final requested refinement solve.

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
- each construction-WF Wannier90 call occurs exactly once on its fragment
  communicator with a collision-free seed, and no fragment publishes a
  partial basis when another fragment fails;
- fragment buffers cover every accepted face and nonlocal-projector support,
  and the union metric retains the required rank under independent fragment
  gauge rotations;
- every short LOBPCG update is finite, metric-orthonormal, rank preserving, and
  within the residual-growth safety bound;
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
- `dg_dc_seed_checkpoint.f90`: exact rank-count/rank--fragment provenance and
  the saved fragment `rwf`, eigenvalue, occupation, density, and potential
  payload used to avoid a preliminary complete LCFO solve.
- `dg_overlapping_wannier_w90.f90`: communicator-parametric Gamma setup,
  matrix assembly, run, and transform.  Add a fragment orchestration layer and
  unique artifact namespaces; do not fork the numerical wrapper.
- `dg_hybrid_block_cg.f90`: preserve the strict solver and add a separate
  bounded fragment-subspace advancement API with successful cap semantics.
- `dg_hybrid_fragment_solver.f90`: fragment generalized solve, after it accepts
  the production DG self block and the occupations obtained for the current
  iteration.
- `s_dg_hybrid_fixed_payload` and production face/projector graphs: immutable
  metric, volume, nonlocal, and SIPG data.
- `dg_hybrid_generalized_eigensystem.f90`: complete distributed LCFO solve.
- the localization-first LCFO certification and complete-v3 checkpoint path:
  extract a terminal finalization routine shared by one-shot, refined, and
  reference routes.

The remaining divided implementation cannot be promoted unchanged: its
production callback still uses ordinary local `hpsi` with an identity metric,
the construction path still performs a preliminary complete LCFO/global
Wannier90 operation, and the fragment spectrum path still diagonalizes the
whole local basis.  The shared convergence, electron-count, and common-μ
semantics established in Tasks 1--4 remain authoritative.

## Test Strategy

Implementation follows TDD in this order:

1. shared pure tests proving bitwise-identical conventional/divided convergence
   metrics for all supported convergence modes;
2. MPI tests for a common chemical potential, current-iteration occupations,
   finite-temperature tails, and exact electron count;
3. fragment-Wannier MPI fixtures proving exactly one setup/run per fragment
   communicator and basis epoch, variable ranks, unique seed namespaces,
   collective failure, buffer-tail coverage, and independence of the final
   LCFO observables under separate fragment unitary gauges;
4. fragment-operator fixtures comparing volume, nonlocal-projector, SIPG self,
   and neighbor-support contributions with a direct small reference;
5. bounded LOBPCG tests against those fixed fragment self blocks, proving a
   three-step safe nonconverged update, warm starts, monotone state extension,
   retained metric rank, and unchanged strict solver semantics;
6. unique-core density and rank-decomposition invariance tests;
7. route tests requiring no preliminary complete LCFO solve, no construction
   Wannier90 call inside SCF, exactly one final complete eigensolve for zero
   refinement steps, and exactly `N+1` final/refinement solves for count `N`;
8. final LCFO residual, metric orthogonality, dynamic window/cluster extension,
   physical symmetry, and complete-v3 checkpoint tests;
9. DC-seed tests retaining the exact rank-count and rank--fragment reuse rule;
10. Si64 eight-rank validation against the saved repeated-continuation oracle,
   recording construction-WF time, energy, gap, density difference, occupied
   projector, symmetry, total time, and memory; and
11. separate-directory zero-field RT runs for one-shot and a small refinement
   count, reporting stationarity rather than silently adding iterations.

All existing dirty changes and verification directories remain untouched.

## Superseding route-removal decision, 2026-09-05

After the separated fragment-local production route is implemented and its
focused numerical and source-contract tests are GREEN, remove the obsolete
divided-only implementation rather than retaining a hidden fallback.  The
boundary is behavioral, not merely a list of routine names:

- remove the divided callback that applies ordinary grid `hpsi` to every
  basis tile, its identity-metric companion, and the whole-local-basis dense
  fragment eigensolver/occupation adapter;
- remove the preliminary-complete-LCFO-dependent divided branch and any state
  padding or raw-retained-rank initialization used only by that branch;
- retain conventional overlapping-Wannier ground state, continuation, shared
  DC potential/mixing helpers, production operator assembly, and the single
  terminal complete LCFO path;
- do not delete shared routines solely because the new divided route no longer
  calls them; another production or reference route must be checked first; and
- perform deletion only after the new path passes direct-construction,
  admission, production-H/S, thermal occupation, bounded-update and divided-SCF
  tests.  Rerun those tests after deletion so no result depends on the legacy
  implementation remaining in the source.

This staged replacement is preferred to immediate deletion because it gives a
testable behavioral handoff while still leaving one supported divided route at
the end.  There is no runtime option or automatic fallback to the removed
algorithm.
