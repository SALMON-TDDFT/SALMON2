# Hybrid Localization-First LCFO Symmetry Design

**Status:** Approved on 2026-09-01

## Goal

Construct the Hybrid DG basis for locality first, and require crystallographic
symmetry only of the physically used full-system LCFO eigenspaces.  The
localized WF+PW basis may be nonclosed as a whole provided that it reproduces
a symmetric occupied ground state and a user-selected low-energy excitation
window.  Reuse a converged conventional DC state while iterating on
localization, PW, LCFO, and symmetry-window controls.

## Motivation

The current Hybrid continuation route symmetry-adapts the occupied seeds,
constructs fixed-center and translation-character sectors, supplies a `.dmn`
file to Wannier90, and validates several pre-Wannier covariance receipts.  This
protects symmetry of the individual WFs, but it restricts the gauge before
localization and can make the WFs less local than necessary.

Individual WFs do not need to transform as a closed representation.  The
physical requirement is that the full-system Hamiltonian solved in the
retained WF+PW variational space produce a symmetric occupied projector,
density, and required low-energy eigenspace.  Plane waves can supply symmetry
content absent from the localized WFs.  A fixed target such as 384 states is
also not a physical input: the required rank varies with material spectrum and
must be derived from an energy window.

Repeated conventional DC solves are unnecessary when only downstream
localization or LCFO controls change.  A dedicated post-DC seed checkpoint is
therefore part of this design.

## Chosen Architecture

The `yn_dg_hybrid_continuation_scf='y'` route uses a localization-first arm:

1. Obtain a converged conventional DC state, either by solving it or by loading
   an exact compatible DC seed.
2. Build raw seeds from the DC occupied states and the complete all-atom s+p
   projection manifest.
3. Apply only metric rank detection and metric orthonormalization to the raw
   seeds.
4. Run ordinary Wannier90 localization with `site_symmetry=.false.` and without
   a `.dmn` file.
5. Add the PW orthogonal complement selected by the user PW cutoff, preserving
   complete boundary shells and authoritative nonidentity reciprocal-symmetry
   orbits.  The current identity-only Hybrid LCFO selector is replaced; WF
   symmetry deferral must not discard the physical reciprocal action.
6. Assemble the spatially divided DG metric and Hamiltonian in the localized
   WF+PW basis.
7. Solve the full-system LCFO generalized eigensystem and certify the occupied
   space plus a user-selected excitation-energy window.
8. Publish the localized basis, operators, occupied coefficients, and physical
   symmetry receipts for real-time initialization.

Other overlapping-Wannier routes retain their existing behavior.  The shared
Wannier90 helper receives an explicit constrained/unconstrained mode so that an
unconstrained Hybrid run cannot accidentally reuse a stale `.dmn` file or
publish constrained-WF provenance.

## Localization Boundary

The localization-first arm retains:

- raw occupied LCFO rows;
- the complete s+p projection catalog;
- ordinary metric rank and orthonormality checks;
- crystallographic point-operation maps, products, translations, rotations,
  generators, and identity/integrity checks;
- Wannier90 overlap/projection construction;
- finite, dimensionally valid, unitary Wannier90 transforms;
- Wannier centers, spreads, localization convergence, and projected spatial
  diagnostics.

It bypasses:

- fixed-center group and character construction;
- occupied translation/point averaging;
- spectral basin/orbit symmetry adaptation;
- translation-character sector reconstruction;
- fixed-center `.dmn` generation;
- `site_symmetry=.true.`;
- pre- and post-Wannier affine-closure acceptance gates.

The crystallographic maps are still authoritative physical data.  They are
kept for projected-action diagnostics and the post-LCFO physical tests, not to
constrain the WF gauge.

Projected operations on the complete localized basis need not be unitary or
closed.  Routines on this arm must not pass them into helpers whose contracts
assume an exact representation.

## Localization Quality and Rank Policy

Localization-first means that the fixed-rank raw seed span is optimized by the
ordinary unconstrained Wannier90 spread functional.  It does not introduce a
post-Wannier per-WF pruning rule.  All `ntarget` WFs are retained, and their
unitary gauge therefore preserves the input span and occupied representability.
Discarding a WF by a material-independent spread cutoff could remove a required
symmetry partner or occupied direction and is outside this implementation.

Acceptance requires a converged Wannier90 run, finite centers and per-WF
spreads in `bohr^2`, a finite total spread, and a unitary transform.  The run
publishes the minimum, maximum, mean, and total spreads and the iteration count.
There is no universal maximum-spread pass/fail threshold; locality is compared
through these receipts when the user changes localization inputs.  Existing
localization tolerance keywords must either be wired to a documented
Wannier90 convergence quantity or explicitly retained as diagnostic-only
compatibility inputs; they are never silently treated as WF-retention cutoffs.

## User-Controlled Cutoffs

Four existing/new controls affect different stages and must not be aliased:

- `energy_cut` is the existing fragment-state energy cutoff used while forming
  the upstream DC-LCFO fragment basis.
- `lambda_cut` is the existing overlap-eigenvalue cutoff used for fragment
  metric cleanup.  Removing a mode here cannot be repaired by a later PW or
  symmetry-window setting.

- `wannier_pw_cutoff` is the existing kinetic-energy cutoff for PW augmentation.
  It remains user supplied in the selected SALMON input energy unit and is
  converted with `uenergy_to_au`.  A cutoff boundary never splits a degenerate
  PW shell or a required reciprocal-symmetry orbit.
- `dg_hybrid_symmetry_energy_window` is a new excitation window measured from
  the HOMO.  It is supplied in the same input energy unit and converted
  independently.  It does not choose or truncate the RT propagation basis.

The Hybrid continuation route requires finite `energy_cut`, positive finite
`lambda_cut`, and positive finite `wannier_pw_cutoff`.  In this route a zero PW
cutoff does not mean “G=0 only” or “disabled”; it is rejected.  The new symmetry
window defaults to exactly `-1d0`, which is the only legacy-rank sentinel.
Other negative or nonfinite values are invalid.  Nonnegative values select the
energy-window mode.

`wannier_pw_max` is not a physical rank selector on this route.  If retained as
an allocation safety limit, it may never truncate a completed energy shell or
reciprocal orbit; such a conflict is reported as a capacity failure rather than
silently clipping the basis.

PW shell membership uses an energy-dimensional machine-precision tolerance,
not `dg_ow_symmetry_tolerance`.  For a PW energy `E_g`, use

```text
tau_pw = 64*epsilon(1d0)*max(1 hartree, abs(wannier_pw_cutoff), abs(E_g))
```

Include every vector below `wannier_pw_cutoff + tau_pw`, complete the numerical
boundary shell, and then close its authoritative reciprocal orbit.  The
requested cutoff, maximum retained PW energy, shell-extension size, and
reciprocal-orbit extension are published.  The physical reciprocal operation
catalog is required to be nonidentity when the crystal catalog is nonidentity;
the Hybrid selector may not replace it by an identity-only catalog.

For `dg_hybrid_symmetry_energy_window >= 0`, let

```text
E_cut = E_HOMO + dg_hybrid_symmetry_energy_window
```

The target contains every LCFO state with energy at or below `E_cut`, never
fewer than the occupied rank, and is extended through the complete numerical
degeneracy cluster at its upper boundary.  The HOMO is the last state whose
authoritative occupation exceeds `64*epsilon(1d0)`.

Energy-degeneracy comparisons use an energy-specific numerical tolerance,
separate from `dg_ow_symmetry_tolerance`.  For adjacent eigenvalues `E_i` and
`E_j`, use

```text
tau_deg = max(dg_dc_gs_final_orbital_tolerance, 64*epsilon(1d0))
          * max(1 hartree, abs(E_i), abs(E_j), abs(E_cut))
```

Window membership is `E_i <= E_cut + tau_deg`; a boundary cluster continues
while adjacent gaps are at most `tau_deg`.  All energies and tolerances in this
comparison are already in atomic units.

The current distributed backend uses full `PZHEEVD` and computes the complete
retained spectrum even when its interface returns only a prefix.  The initial
implementation therefore performs one full generalized eigensolve per
continuation candidate and selects the energy window from that result.  It
must not repeat the same full solve while geometrically increasing a requested
prefix.

If a future backend genuinely supports a partial spectrum, it may start with
the occupied requirement plus a proof state and grow monotonically:

```text
k_new = min(N_basis, max(k + 16, ceil(1.5*k)))
```

The selected window is certified only after one state above the cutoff and
boundary cluster has been obtained.  Reaching the full retained basis without
that proof state is a basis-ceiling failure, not a successful full-basis
special case.  The message identifies the PW/basis capacity as insufficient;
the initial implementation never changes the user cutoff automatically.  This
coverage rule is intentionally stricter than merely enumerating every Ritz
state available in a finite retained space: a basis whose highest Ritz value
does not bracket the user-requested energy cannot certify physical coverage of
that energy window.

A negative energy window retains the previous dynamic-rank behavior with an
explicit compatibility warning.  Production acceptance inputs, including
Si64, set the energy window explicitly.  No production code contains a
material-specific 384-state target.

## Occupation and HOMO Policy

The authoritative SALMON occupation policy, electronic temperature, electron
count, and spin convention are reapplied to the ascending final LCFO spectrum;
occupations are not copied by state identity from pre-localization orbitals.
The chemical potential is solved consistently from the complete retained
spectrum.  Insufficient state capacity for the requested electron count is
fatal.

`noccupied` means the number of leading final LCFO columns whose occupation is
greater than `64*epsilon(1d0)`.  The checkpoint stores and propagates exactly
those columns, while the omitted occupation tail must be below the configured
electron-count tolerance.  `E_HOMO` is the eigenvalue of the last such column.
The density, projector, energy-window anchor, and RT occupation kernel all use
this same definition.

## Physical Symmetry Acceptance

For the metric-orthonormal LCFO coefficients `C`, solve

```text
H C = S C epsilon
C^H S C = I
```

For every retained physical symmetry operation, acceptance requires all four
finite defects to remain within `dg_ow_symmetry_tolerance`:

1. occupied eigenspace and occupied-projector covariance;
2. degeneracy-complete energy-window subspace closure;
3. Hamiltonian energy covariance within that target subspace;
4. input/output real-space density covariance.

Individual eigenvectors inside a degenerate space are gauge dependent, so the
test acts on subspaces and projectors, not state labels.

The following remain finite, nonnegative diagnostics and are not physical
acceptance gates:

- individual-WF covariance;
- Wannier-center orbit closure;
- full WF+PW retained-basis closure;
- full retained-space metric and Hamiltonian covariance.

The code does not explicitly symmetrize the Hamiltonian, coefficients,
occupations, or density.  A failure therefore demonstrates insufficient
variational content or another physical inconsistency rather than being hidden
by projection.

## Conventional DC Seed Reuse

Add dedicated controls in `&dc`:

```text
dg_dc_seed_mode = 'off' | 'write' | 'read' | 'auto'
dg_dc_seed_directory = '<path>'
```

The default is `off` for compatibility.

- `write`: run conventional DC, then publish a new seed before WF generation.
- `read`: require and load an exact compatible seed; never fall back to a DC
  solve.
- `auto`: load a valid seed if present; if no committed seed exists, run DC and
  publish it.  A present but incompatible or corrupt seed is fatal rather than
  silently recomputed or overwritten.
- `off`: retain the existing conventional DC behavior.

The minimal rank-local payload contains:

- `spsi%rwf` with exact allocated bounds;
- the owned `dc%rho_tot_s(1)%f` slab;
- the owned `dc%vloc_tot(1)%f` slab;
- `energy%esp`, `system%rocc`, and `system%mu`;
- the accepted convergence residual and iteration provenance.

Derived fragment potentials and work arrays are rebuilt after load.  Mixing
history, communicators, grids, stencil workspaces, WFs, PWs, and LCFO data are
not serialized.

The immutable manifest covers cell and grid, atom/species order, fragment and
buffer topology, electron/state/spin/Gamma-real configuration, occupation
policy, XC, stencil, pseudopotential content, rank-to-fragment mapping, local
array bounds, and the convergence contract.  Downstream localization,
Wannier90, PW cutoff, LCFO, and symmetry-window controls are deliberately
excluded.

The formal format permanently requires the same MPI rank count, the same
rank-to-fragment mapping, and the same local array bounds.  DC seed
redistribution is not a future compatibility goal.  This restriction applies
only to the DC seed; the separate GS-to-RT checkpoint retains its existing rank
redistribution capability.  Compatibility fingerprints include the exact
global-to-local grid, orbital, and fragment ownership maps; equal counts and
bounds alone are not considered proof of identical ownership.

Each rank writes a versioned shard with a digest.  Publication uses temporary
shards followed by an atomic manifest-last commit.  The reader validates the
complete shard set, publication identifier, ordered global digest, shapes and
bounds, finite data, electron count, convergence residual under the current
threshold, and every immutable fingerprint.

## GS-to-RT Data Flow

The localized DG basis functions remain the propagation basis.  They are not
rotated into delocalized full-system eigenstates.  Full-system occupied states
are represented by row-owned coefficient columns:

```text
psi_n(r) = sum_mu phi_mu(r) C(mu,n)
```

The GS checkpoint stores:

- the complete retained localized basis description;
- sparse metric and Hamiltonian components;
- row ownership and basis ordering;
- occupied coefficient columns, eigenvalues, occupations, and density;
- position, nonlocal, face, pseudopotential, and provenance receipts;
- the energy-window mode, window size, HOMO and cutoff energies, solved rank,
  selected rank, boundary-cluster rank, proof-state status and energy, and four
  physical symmetry defects;
- canonical basis catalog IDs, generations, ordering, and their fingerprint.

These fields form a named, versioned checkpoint receipt rather than extending
the current positional real array.  They participate in the payload digest and
the basis/operator provenance chain.  This changes the Hybrid ground-state
checkpoint to version 3.  The new Hybrid continuation RT route requires a
version-3 physical receipt; it never infers one from a version-2 positional
array.  Any retained legacy reader remains confined to its legacy route.

Nonoccupied eigenvectors used only to certify the GS window are not required
in the RT checkpoint.  The complete WF+PW basis remains available to each
occupied state during propagation:

```text
i S dC_n(t)/dt = H(t) C_n(t)
```

Consequently, the external field can move occupied amplitudes into every
retained basis direction even though only occupied initial coefficient columns
are stored.

RT startup revalidates `H*C-S*C*epsilon`, `C^H*S*C-I`, electron count,
reconstructed density, occupied projector and density symmetry, and operator
and basis fingerprints.  It authenticates the versioned GS energy-window
receipt through the payload digest and matching provenance, but does not
recompute target-subspace closure because nonoccupied target vectors are not
stored.  Full-operator covariance is finite/nonnegative diagnostic data only.
Zero-field RT retains the stationarity gates for density, energy, projector,
charge, and Hamiltonian residual.

## Provenance and Diagnostics

Map provenance and localization provenance are distinct.  An unconstrained WF
set is never labeled symmetry certified.  Replay export is mode aware and does
not copy or consume `.dmn` in unconstrained mode.

Rank-zero receipts include at least:

```text
[DG-DC-SEED] mode=... publication_id=... scf_skipped=...
[HYBRID-WF-LOCALIZATION] symmetry_constraint=off ...
[HYBRID-RETAINED-BASIS-SYMMETRY] closed=... defect=...
[HYBRID-PW-CUTOFF] requested=... effective=... shell_added=... orbit_added=...
[HYBRID-LCFO-WINDOW] delta_e=... homo=... cutoff=... solved_rank=... window_rank=... proof=... proof_energy=...
[HYBRID-LCFO-SYMMETRY] occupied=... target=... energy=... density=... worst_operation=...
[HYBRID-RT-HANDOFF] payload_fingerprint=... projector_symmetry=...
```

The receipts must distinguish basis construction rank, full retained rank,
occupied rank, and energy-selected physical target rank.

## Error Classification

Structural corruption remains immediately fatal:

- invalid dimensions, group maps, products, identities, or ownership;
- metric rank loss, singular/indefinite metric, failed collective, or failed
  distributed solve;
- nonfinite data, transform nonunitarity, Hamiltonian non-Hermiticity, or
  fingerprint mismatch;
- incomplete, corrupt, stale, rank-incompatible, or mapping-incompatible DC
  seed data.

Physical insufficiency is also an unsuccessful calculation, but it is reported
as a capacity or symmetry failure rather than an internal structural error:

- no proof state above the requested energy window before the basis ceiling;
- occupied, target, target-energy, or density symmetry defect above tolerance.

These messages report the WF rank, PW rank, full rank, PW cutoff, HOMO, window
cutoff, solved and selected ranks, proof-state status, worst operation, and all
four physical defects.  The user can then increase the PW cutoff, enlarge the
retained basis, or reduce the requested excitation window.

## Verification Strategy

Focused TDD coverage includes:

- unconstrained Wannier90 input has no `.dmn` and sets
  `site_symmetry=.false.`;
- raw seeds undergo metric orthonormalization without symmetry adaptation;
- constrained legacy routes retain their existing behavior;
- energy-window selection at zero window, between levels, exactly on a level,
  and through split or degenerate multiplets;
- the same energy window selects different ranks for different spectra;
- one full solve with prefix-free energy-window analysis for the current
  backend, plus monotone growth tests only for a genuine partial backend;
- distinct `energy_cut`, `lambda_cut`, PW-cutoff, and symmetry-window semantics,
  including invalid, zero, nonfinite, requested, and effective boundaries;
- authoritative nonidentity reciprocal PW orbit closure;
- final-spectrum occupation/chemical-potential reconstruction and a single
  HOMO definition shared by GS, checkpoint, and RT;
- checkpoint-v3 receipt authentication, version-2 rejection on the new route,
  canonical basis ordering, and proof-state energy provenance;
- WF nonclosure with successful low-energy WF+PW symmetry recovery;
- removal of a required PW partner causes a clear physical insufficiency;
- no production dependence on literal Si64 ranks such as 128 or 384;
- DC seed write/read round trip, strict and auto behavior, corruption,
  truncation, missing shards, interrupted publication, immutable mismatch, and
  downstream-control compatibility;
- permanent rejection of DC seed MPI-rank, rank-to-fragment, or local-bound
  mismatch;
- a cold Si64 run writes a seed and a subsequent run contains no conventional
  `DC #SCF =` iterations while reproducing downstream results;
- localization, PW, and symmetry-window changes reuse the same valid DC seed;
- GS checkpoint to zero-field RT stationarity and finite-field propagation in
  the complete retained basis;
- rank-distribution-invariant occupied-projector and reconstructed-density
  checks across supported GS and RT rank layouts;
- supported focused MPI tests on 1, 2, 4, and 8 ranks, followed by the full
  build and protected regression suite.

Large Si64 outputs and existing diagnostic logs remain untracked evidence and
are never included in implementation commits.
