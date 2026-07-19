# WPW GS Bounded Adapter Design

## Scope

This supplement completes Task 5C without changing the approved `windowed_kg`
basis, fragment-root ownership, canonical-face convention, or DG-DC convergence
gates.  It replaces the provisional full-vector context apply with a bounded
GS-neutral operator adapter and places self-consistent operator publication
after the existing LCFO WW blocks are available.

The dense fixed-basis solver remains a mathematical oracle only.  Production
stores owned W/P coefficient slices and bounded support metadata; it never
allocates a global H, S, density matrix, coefficient vector, WP block, or PP
block.

## Lifecycle and ordering

The current geometry and trace stages remain before `hpsi_basis`, because they
consume the explicit basis and its gradients.  They publish only unpublished
quadrature candidates and prepared halo state.  They do not bind callbacks.

After `hpsi_basis` and `calc_hamiltonian_matrix` have produced the LCFO WW
local and canonical cross-face blocks, the fragment roots perform one
transactional publication:

```text
prepared geometry/trace epoch
  -> reduced WP/PP candidates
  -> current WW local/halo blocks
  -> sparse H/S blocks
  -> coefficient schedule
  -> bounded callbacks
  -> SCF iteration state
```

Any failure leaves the previous valid operator epoch untouched.  No callback
is visible while only some WW/WP/PP components have been replaced.

## GS-neutral bounded operator adapter

A common module owns only plain production data:

- production communicator and epoch/fingerprint;
- owned and required-input W/P stable IDs;
- compact WW local/cross blocks;
- sparse WP and PP blocks;
- a fixed coefficient fetch/reverse-reduction schedule.

The module exposes batched `apply_h`, `apply_s`, and `global_gram` callbacks
using the existing `dg_wpw_matrix_free_operator` interfaces.  It has no
dependency on `s_dg_fragment_rt`, LCFO host variables, or RT integrator state.
GS and RT construct this same common context from their own storage adapters.
The owner-targeted coefficient exchange is also moved to the common layer.
No module under `src/common` imports `src/rt`; the RT adapter consumes the
common exchange implementation.

Each apply validates all local extents and finiteness collectively before
communication.  It fetches only stable IDs listed by the compact blocks, applies
WW/WP/PP locally, reverse-reduces W contributions to their owners, and publishes
owned output only after the complete exchange succeeds.  Ownership is computed
from bounded fragment-root layout data; no global owner array is stored.

## WW import

The LCFO adapter does not import monolithic `mat_H_local` as one mutable WW
block.  After `calc_hamiltonian_matrix` it publishes separately fingerprinted
fixed kinetic, fixed nonlocal, fixed self/cross-face, and local-potential WW
components.  Existing `mat_H_weak_kinetic`, `mat_H_weak_nonlocal`,
`mat_H_weak_potential`, `mat_H_surface_self`, and canonical entries of
`halo(:)%mat_H_surface_cross` provide their provenance.  The iterative map
replaces only the local-potential component and recomposes total WW once,
preventing stale potential and double counting.

The import uses `n_basis`, `index_basis`, and existing halo neighbor identity.
Each physical cross face is imported once; the two periodic faces of a
two-fragment axis remain distinct records.  It does not traverse `1:n_frag`
and does not assemble `h_full`.

The overlap WW action follows the accepted LCFO basis convention.  If the
current basis is explicitly S-orthonormal, the adapter records WW identity;
otherwise it imports the existing compact overlap blocks.  It must never infer
identity from a tolerance without recording and fingerprinting that policy.

## Self-consistent potential map

The SCF state stores only retained owned W/P coefficients, occupations, and a
rank-local real-space density over uniquely owned fragment-core points.  The
potential map:

1. fetches only coefficient support required on the local core;
2. reconstructs occupied wavefunctions and density on uniquely owned points;
3. verifies charge and finiteness collectively;
4. forms `rho_candidate=(1-mix)*rho_in+mix*rho_raw`;
5. calls the existing distributed Hartree and LDA paths exactly once for
   `rho_candidate` and does not separately mix the resulting potential;
6. rebuilds potential-dependent WW/WP/PP volume contributions from that
   candidate potential;
7. combines them with fixed kinetic, nonlocal, and canonical-face components;
8. publishes a new operator epoch only after all ranks validate it.

Fixed face terms are not reintegrated during an SCF iteration.  PP remains
periodic H1 and receives no face term.  A failed map returns no density or
operator update.

### Communicator choreography

The eigensolver algebra runs on fragment roots in the production communicator,
but density reconstruction and existing Hartree/LDA collectives require every
rank in `dc%icomm_frag` and `dc%icomm_tot`.  Root-only callbacks therefore
never enter total-rank collectives.

An outer all-rank driver advances a fixed command/epoch state machine.  During
an algebra command, production roots apply H/S while fragment nonroots wait at
the next fragment command broadcast.  For a potential-map command, each root
broadcasts only its fragment's required occupied coefficient slices; all
fragment ranks reconstruct their assigned core points and enter Hartree/LDA in
the same order.  The next command is not issued until fragment and production
status reductions complete.  Failure broadcasts a terminal command so no
nonroot remains in a service loop.

`run_dg_wpw_matrix_free_scf` is split into stepwise algebra and convergence
operations.  It does not hide a variable number of all-rank potential-map
calls inside a root-only solver invocation.

### Fragment-core and Hartree/LDA layout bridge

The density bridge uses the existing global real-grid ownership map to route
each uniquely owned fragment-core point to exactly one global-grid owner.
Messages contain the unwrapped global grid ID, density value, integration
weight, epoch, and source fragment.  The bridge rejects duplicate or missing
grid IDs and verifies that the redistributed charge equals the fragment-core
charge before invoking Hartree/LDA.  The returned potential is routed back only
to fragment-core and bounded support points required by WW/WP/PP quadrature;
no rank stores or scans the complete fragment set.

Boundary points retain the existing SALMON grid integration weights.  Periodic
images map to one global grid ID for density/Hartree ownership while their
distinct image identity remains available to canonical-face geometry.

### Energy functional

The map has one candidate state tuple
`(rho_in,rho_raw,rho_candidate,V_candidate,E_candidate,H_candidate)`.  Only
density is mixed.  `V_candidate`, `E_candidate`, and every potential-dependent
operator block are evaluated from `rho_candidate` in the same epoch.  A
cross-epoch combination or separately mixed potential is invalid.

The production energy gate uses the existing SALMON LDA/Hartree energy
components evaluated at `rho_candidate`: band/kinetic, local and nonlocal ionic, Hartree with its
double-counting correction, XC energy and potential correction, and ion-ion
energy.  Occupations are those passed to the retained subspace.  Energy and
operator publication share one density/potential epoch; a scalar computed from
`rho_in`, `rho_raw`, or another epoch is rejected.  The fixed-point map repeats
this evaluation with `mix=1`, requires `rho_raw=rho_in` within tolerance, and
requires both density and total-energy reproduction.

Mixer history, previous residuals, iteration counter, candidate density,
candidate potential, candidate energy, and candidate operator form one
transaction object.  Mixing never mutates live history in place.  A failure
after mixing, redistribution, Hartree/LDA, energy evaluation, or operator
rebuild discards the entire candidate and leaves the prior transaction epoch
and values unchanged.

### Shared GS/RT provenance

The common context records ownership kind/version, communicator layout hash,
owned/support ID hashes, deterministic block ordering, metric policy/cutoff,
WW component provenance, WP/PP operator convention, canonical-face convention,
coefficient-schedule hash, and operator epoch.  Block entries are sorted by
stable IDs and component kind.  Task 6 serializes these fields without
reinterpretation, and Task 7 reconstructs the same context and compares
rank-local H/S action and fingerprints before RT starts.

## Initialization and convergence

The retained dimension is exactly `n_occ + dg_wpw_extra_states`, bounded by
the distributed basis dimension.  Initial vectors are generated
deterministically from stable IDs and then S-orthonormalized through the common
callbacks; conventional orbital projection is not used.

`run_dg_wpw_matrix_free_scf` retains the approved gates: density, potential,
energy, occupied projector, generalized residual, metric orthonormality,
charge, gap, and one-map fixed-point reproduction.  A nonpositive retained
dimension, insufficient metric rank, missing gap state, or any failed gate is
fatal for the explicit production route.

## Reachability

`main_dft` keeps the existing explicit `yn_dc_lcfo_flux` entry and additionally
requires `yn_dg_wpw_production='y'` for this consumer.  `lcfo_flux` owns the
ordered construction because it owns the explicit basis, fragment WW blocks,
and real-space potential.  The source contract verifies an actual call chain,
not the presence of an unused import or comment.

## Tests and failure cases

TDD fixtures cover:

1. bounded two-rank WW/WP/PP H/S action against a dense oracle;
2. differing W ownership and neighbor support cardinalities;
3. stale epoch, missing coefficient, nonfinite input, invalid fingerprint,
   truncated header, and reverse-reduction failure;
4. LCFO WW import with separate two-fragment periodic faces;
5. one-map density/potential/operator rebuild with transactional rollback;
6. retained dimension, S-orthonormality, charge, gap, all convergence gates,
   and fixed-point reproduction;
7. fragment-count scaling with constant bounded support degree;
8. source rejection of full vectors, global matrices, global fragment scans,
   dense production SCF, and callback binding before WW publication.

Task 5C is complete only when focused fixtures, all WPW contracts, the full
MPI/EigenExa build, and `git diff --check` pass, followed by a findings-first
review with no Critical or Important findings.
