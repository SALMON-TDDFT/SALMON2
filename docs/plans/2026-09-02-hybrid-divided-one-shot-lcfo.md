# Hybrid Divided-SCF and One-Shot LCFO Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build construction WFs directly and independently from reusable DC fragment orbitals, converge the WF+PW density with bounded warm-started fragment updates, and perform exactly one complete LCFO solve by default before optional refinement and complete-v3 RT publication.

**Architecture:** Tasks 1--4 establish shared conventional-DC convergence, electron-count, fragment spectrum/density, and common-chemical-potential semantics.  The remaining production path consumes the exact-rank-compatible DC `rwf` seed directly, performs one unconstrained Wannier90 rotation per fragment communicator and basis generation, and augments it with the dynamic windowed-PW complement.  Divided SCF retains an uncompressed fragment-owned catalog and advances only occupied-plus-guard states with a default three-step warm-started `[X,R,P]` LOBPCG update in the correct broken-volume/SIPG self block; a separate immutable union map removes only complete-system metric-null directions for terminal LCFO.  No preliminary complete LCFO or complete-cell construction-WF Wannier90 operation is allowed; the complete distributed LCFO eigensolve occurs once after divided density convergence, while dense fragment solves and repeated complete-Hybrid continuation remain explicit reference backends.

**Tech Stack:** Fortran 2008, SALMON DC/LCFO, MPI, ScaLAPACK, Wannier90, CMake, standalone Python source-contract tests, linked Fortran MPI fixtures.

---

The worktree is intentionally dirty.  Preserve every existing modification and
all verification directories.  Never use `git add -A`, `git commit -a`,
`git reset`, `git checkout --`, or `git clean`.  Use `git add -p` for every
already modified file, inspect `git diff --cached`, and commit only the current
task's hunks.  Do not create another worktree.

The earlier eight-rank continuation attempt has ended and its output remains
reference evidence.  Do not restart a heavy Si64 calculation before Task 11
and do not remove or overwrite any existing verification output.  Source-level
and small MPI tests may proceed.

The production scope remains the accepted Gamma, non-SOI, PZ-LDA, gapped
Hybrid route.  Do not broaden theory scope in this plan.

### Task 1: Share the authoritative conventional-DC density convergence metric

**Files:**

- Create: `src/gs/dc/dc_scf_convergence.f90`
- Create: `tests/dg/test_dc_scf_convergence_mpi.f90`
- Create: `tests/dg/run_dc_scf_convergence_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/scf_iteration_dft.f90:497-537`
- Modify: `src/gs/dc/dg_hybrid_divided_scf.f90:82-101`

**Step 1: Write the failing MPI fixture**

For a distributed density difference with known absolute sum and square sum,
require the conventional definitions:

```fortran
call reduce_dc_density_convergence(comm,'rho_dne',local_abs,local_square,&
  hvol,nelec,ngrid,value,ok,message)
call require(ok .and. abs(value-hvol*global_abs/nelec)<1d-14,&
  'rho_dne normalization changed')

call reduce_dc_density_convergence(comm,'norm_rho',local_abs,local_square,&
  hvol,nelec,ngrid,value,ok,message)
call require(ok .and. abs(value-global_square)<1d-14,&
  'norm_rho normalization changed')

call reduce_dc_density_convergence(comm,'norm_rho_dng',local_abs,local_square,&
  hvol,nelec,ngrid,value,ok,message)
call require(ok .and. abs(value-global_square/real(ngrid,real64))<1d-14,&
  'norm_rho_dng normalization changed')
```

Also require collective rejection of unsupported modes, non-finite
accumulators, non-positive volume/electron/grid counts, and rank-disagreeing
scalar controls.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dc_scf_convergence_mpi.py`

Expected: compile failure because `dc_scf_convergence` does not exist.

**Step 3: Implement the shared reducer**

Expose only this narrow API:

```fortran
module dc_scf_convergence
  use,intrinsic::iso_fortran_env,only:real64
  implicit none
  private
  public::reduce_dc_density_convergence
contains
  subroutine reduce_dc_density_convergence(comm,mode,local_abs,local_square,&
      hvol,electron_count,global_point_count,value,ok,message)
    integer,intent(in)::comm,global_point_count
    character(*),intent(in)::mode
    real(real64),intent(in)::local_abs,local_square,hvol,electron_count
    real(real64),intent(out)::value
    logical,intent(out)::ok
    character(*),intent(out)::message
```

Use one collective reduction of `[local_abs,local_square]`, then apply exactly
the formulas in `scf_iteration_dft.f90`.  Do not take square roots and do not
replace `rho_dne` by a maximum norm.

**Step 4: Route both implementations through the helper**

Keep each caller's existing local grid loop, but send its local absolute and
square sums to the shared reducer.  `norm_pot` and `norm_pot_dng` remain in the
ordinary SCF implementation and are outside the divided-route scope.

**Step 5: Run GREEN and protected checks**

Run:

```text
python3 tests/dg/run_dc_scf_convergence_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
cmake --build build-hybrid-commit -j2
```

Expected: MPI 1/2/4 ranks PASS, divided-driver tests PASS, and the build
completes.

**Step 6: Commit only Task 1**

```text
git add src/gs/dc/dc_scf_convergence.f90 tests/dg/test_dc_scf_convergence_mpi.f90 tests/dg/run_dc_scf_convergence_mpi.py
git add -p src/gs/dc/CMakeLists.txt src/gs/scf_iteration_dft.f90 src/gs/dc/dg_hybrid_divided_scf.f90
git diff --cached --check
git diff --cached
git commit -m "fix(dc): share divided density convergence semantics"
```

Checkpoint after this task: report the exact RED and GREEN evidence before
continuing.

### Task 2: Enforce divided-SCF electron count and terminal potential refresh

**Files:**

- Modify: `src/gs/dc/dg_hybrid_divided_scf.f90`
- Modify: `tests/dg/test_dg_hybrid_divided_scf_mpi.f90`
- Modify: `src/gs/main_dft.f90:3613-3622`

**Step 1: Extend the fixture RED**

Add `core_weights`, `expected_electron_count`, and `electron_tolerance` to the driver fixture.
Require collective rejection when any density callback returns NaN or when

```fortran
abs(electron_count-expected_electron_count) > electron_tolerance
```

Require the rejected path not to call the mixer or publish a converged density.
Independently integrate `sum(core_weights*new_density)` across the total
communicator and reject disagreement with either the callback count or the
target count.  Apply the same gate to the mixed density.
On convergence, require one terminal `update_total_potential(new_density)`
after the accepted density callback, so the returned density and the potential
used by final LCFO have the same potential epoch.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py`

Expected: fixture compile failure on the new driver arguments or an assertion
showing that the bad electron count was accepted.

**Step 3: Extend the driver contract**

Use this signature tail:

```fortran
...,mix_dc_density,maximum_iterations,core_weights,expected_electron_count,electron_tolerance,&
converged_density,iterations,convergence_value,electron_defect,ok,message)
```

Validate all scalars collectively before iteration.  After every core-density
callback, check finiteness and electron count before calculating convergence.
When convergence passes, refresh the potential from `new_density`; publish
nothing if that refresh fails.

**Step 4: Update the production call and run GREEN**

Pass `dc%elec_num_tot` and `dg_dc_gs_electron_count_tolerance` from
`main_dft.f90`.  Emit the final electron defect in the divided-SCF receipt.

Run:

```text
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/run_dg_dc_seed_state_mpi.py
```

Expected: all PASS for their configured rank sets.

**Step 5: Commit**

```text
git add -p src/gs/dc/dg_hybrid_divided_scf.f90 tests/dg/test_dg_hybrid_divided_scf_mpi.f90 src/gs/main_dft.f90
git diff --cached --check
git diff --cached
git commit -m "fix(dg): gate divided density electron count"
```

### Task 3: Separate fragment eigensolving from occupation-owned density

**Files:**

- Modify: `src/gs/dc/dg_hybrid_fragment_solver.f90`
- Modify: `tests/dg/test_dg_hybrid_fragment_solver_mpi.f90`
- Modify: `tests/dg/run_dg_hybrid_fragment_solver_mpi.py`

**Step 1: Write RED split-phase tests**

Require a spectrum phase that returns coefficients/eigenvalues/core norms
without accepting occupations, followed by a density phase that accepts the
current iteration's occupations.  Verify that the compatibility wrapper gives
the same density as the two explicit phases.

```fortran
call solve_dg_hybrid_fragment_spectrum(...,coefficients,eigenvalues,&
  core_norms,residual,orthogonality,ok,message)
call reconstruct_dg_hybrid_fragment_density(basis,coefficients,occupations,&
  core_mask,point_weights,density,electron_count,ok,message)
```

Reject negative/greater-than-spin-degeneracy occupations, extent mismatches,
non-finite coefficients, and duplicated/non-unique core ownership.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dg_hybrid_fragment_solver_mpi.py`

Expected: compile failure because the two new public routines are absent.

**Step 3: Extract without changing algebra**

Move the existing generalized solve into
`solve_dg_hybrid_fragment_spectrum`; move density reconstruction and electron
integration into `reconstruct_dg_hybrid_fragment_density`.  Keep
`solve_dg_hybrid_fragment_basis` as a thin compatibility wrapper until all
callers are migrated.

**Step 4: Run GREEN**

Run:

```text
python3 tests/dg/run_dg_hybrid_fragment_solver_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
```

Expected: PASS for 1/2/4 ranks.

**Step 5: Commit**

```text
git add -p src/gs/dc/dg_hybrid_fragment_solver.f90 tests/dg/test_dg_hybrid_fragment_solver_mpi.f90 tests/dg/run_dg_hybrid_fragment_solver_mpi.py
git diff --cached --check
git diff --cached
git commit -m "refactor(dg): separate fragment spectrum and density"
```

### Task 4: Reuse one common DC chemical potential and current occupations

**Files:**

- Create: `src/gs/dc/dc_fragment_occupation.f90`
- Create: `tests/dg/test_dc_fragment_occupation_mpi.f90`
- Create: `tests/dg/run_dc_fragment_occupation_mpi.py`
- Modify: `src/gs/occupation_kernel.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/dc/dcdft.f90:743-933`
- Modify: `src/gs/main_dft.f90:4911-4939`

**Step 1: Write the common-kernel RED test**

Build two gapped fragment spectra with unequal core norms.  Require one common
chemical potential, occupations in `[0,wspin]`, and

```fortran
sum(occupation(fragment,state)*core_norm(fragment,state)) == Ne
```

within the existing electron tolerance.  Test zero temperature, finite
temperature with an explicitly supplied guard-state tail, degeneracy at the
Fermi edge, decomposition invariance, and collective rejection when the
available spectrum cannot carry `Ne`.

**Step 2: Run RED**

Run: `python3 tests/dg/run_dc_fragment_occupation_mpi.py`

Expected: compile failure because `dc_fragment_occupation` is absent.

**Step 3: Generalize the existing authoritative numerical kernel**

Add a weighted-state entry point to `occupation_kernel.f90` and make the
existing `solve_spectrum_occupations` a shape adapter around it:

```fortran
subroutine solve_weighted_state_occupations(eigenvalues,state_weights,&
    electron_target,electronic_temperature,maximum_occupation,&
    occupations,chemical_potential,electron_count,ok,message)
```

Preserve the existing Fermi function, bracketing, exact zero-temperature
fallback, and stopping semantics.  A state weight is the unique-core norm of
that fragment state; it is not a fragment-wide k-point weight.

**Step 4: Add the distributed fragment adapter**

Expose a routine independent of SALMON structure types:

```fortran
subroutine determine_dc_fragment_occupations(comm,energies,core_norms,&
    representative_mask,temperature,wspin,expected_electrons,tolerance,&
    chemical_potential,occupations,electron_count,ok,message)
```

Only one representative per fragment contributes to the total communicator;
broadcast the result within each fragment communicator.  Do not introduce a
per-fragment chemical potential.

**Step 5: Route conventional and Hybrid callers through it**

Make `ne2mu_dcdft` use the extracted kernel after its existing spectrum/core
norm gathering.  Change `solve_dg_hybrid_divided_fragments` to:

1. solve all fragment spectra;
2. gather representative eigenvalues and core norms;
3. determine the common chemical potential and occupations;
4. reconstruct the local core density with those occupations.

Delete the stale `fragment_occupations=system%rocc(...)` assignment from the
production path.

**Step 6: Run GREEN and conventional protection**

Run:

```text
python3 tests/dg/run_dc_fragment_occupation_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_solver_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
cmake --build build-hybrid-commit -j2
```

Expected: all PASS and the build completes.

**Step 7: Commit**

```text
git add src/gs/dc/dc_fragment_occupation.f90 tests/dg/test_dc_fragment_occupation_mpi.f90 tests/dg/run_dc_fragment_occupation_mpi.py
git add -p src/gs/occupation_kernel.f90 src/gs/dc/CMakeLists.txt src/gs/dc/dcdft.f90 src/gs/main_dft.f90
git diff --cached --check
git diff --cached
git commit -m "fix(dc): share fragment occupation policy"
```

### Task 5: Build construction WFs once per fragment and a generalized PW complement

**Files:**

- Create: `src/gs/dc/dg_hybrid_fragment_wannier.f90`
- Create: `tests/dg/test_dg_hybrid_fragment_wannier_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_fragment_wannier_mpi.py`
- Create: `tests/dg/test_dg_hybrid_fragment_wannier_lcfo_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_fragment_wannier_lcfo_mpi.py`
- Modify: `src/common/dg_hybrid_wannier_complement.f90`
- Modify: `tests/dg/test_dg_hybrid_wannier_complement_mpi.f90`
- Modify: `src/gs/dc/dg_hybrid_projected_fragment_pipeline.f90`
- Modify: `tests/dg/test_dg_hybrid_projected_fragment_pipeline_mpi.f90`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write a RED communicator/epoch Wannier fixture**

Split 2, 4, and 8 total ranks into deterministic fragment communicators.  Feed
each fragment a different number of reusable DC `rwf` columns (including 2
and 3, so no test encodes the Si64 value 384).  Stub the Wannier90 setup and
run entry points and require:

- exactly one setup and one run on the root of each fragment communicator for
  one `(fragment_id,basis_generation)` pair;
- no setup or run on the non-root ranks of that fragment;
- one bitwise-identical transform broadcast to every rank in the fragment;
- distinct seed directories of the form
  `fragment-%06d/generation-%08d`, with the fragment ID and basis generation
  included in the receipt and fingerprint;
- every locally generated WF column retained, even when its center lies in an
  overlapping neighboring buffer;
- reconstruction of every input DC orbital from the returned seed-to-WF
  coefficients, with invariant occupied projector and occupation-weighted
  density before and after the unconstrained rotation;
- no `.dmn` symmetry constraint or cross-fragment gauge-matching input;
- a collective total-communicator failure if any fragment setup, run, output
  parse, or broadcast fails; and
- reusing a successful receipt for the same immutable seed/basis fingerprint,
  but a new single invocation after a genuine basis-generation change.

The fixture must call construction localization outside any density-iteration
callback and assert that callback invocation counts remain zero during a mock
SCF loop.

**Step 2: Run the fragment-Wannier RED test**

Run: `python3 tests/dg/run_dg_hybrid_fragment_wannier_mpi.py`

Expected: compile failure because `dg_hybrid_fragment_wannier` does not exist.

**Step 3: Implement the fragment construction-WF contract**

Expose SALMON-structure-independent metadata and numerical APIs.  The
production wrapper may call the existing Wannier90 adapter internally, but
the fragment module owns validation, namespacing, broadcasts, and receipts:

```fortran
type s_dg_hybrid_fragment_wannier_receipt
  integer :: fragment_id
  integer :: basis_generation
  integer :: candidate_rank
  integer :: retained_rank
  integer :: setup_count
  integer :: run_count
  integer(int64) :: seed_fingerprint
  integer(int64) :: basis_fingerprint
  integer(int64) :: transform_fingerprint
  real(real64) :: seed_reconstruction_defect
end type

type s_dg_hybrid_fragment_wannier_cache
  logical :: valid
  type(s_dg_hybrid_fragment_wannier_receipt) :: receipt
  complex(real64),allocatable :: wannier_values(:,:)
  complex(real64),allocatable :: candidate_compression(:,:)
  complex(real64),allocatable :: wannier_transform(:,:)
  complex(real64),allocatable :: dc_seed_coefficients_in_wannier(:,:)
  real(real64),allocatable :: physical_dc_seed_energies(:)
  real(real64),allocatable :: physical_dc_seed_occupations(:)
end type

subroutine build_dg_hybrid_fragment_wannier(comm_total,comm_fragment,&
    fragment_id,basis_generation,seed_directory,grid_ids,grid_weights,&
    dc_seed_values,dc_seed_energies,dc_seed_occupations,&
    buffer_candidate_values,projector_candidate_values,metric_tolerance,&
    fragment_real_lattice,fragment_reciprocal_lattice,atom_symbols,&
    atoms_cart,fractional_coordinates,num_iter,localization_tolerance,&
    coordinator_byte_limit,&
    cache,ok,message)
```

The input columns come directly from the exact-rank-compatible DC seed
checkpoint; do not call preliminary `dc_lcfo` and do not form complete-cell
LCFO eigenvectors.  Run unconstrained Wannier90 only on the fragment
communicator root, then broadcast the accepted transform.  Preserve all input
DC orbital columns and every independent accepted buffer/projector direction;
remove only candidate-space metric null modes before localization, then
preserve every retained column and its span.  Reject unexpected rank loss,
non-finite output, a transform that is not unitary within tolerance, a changed
fragment/basis fingerprint, or a seed directory collision.  Do not select WFs
by post-localization centers and do not compare or align gauges between
fragments.

Use the accepted fragment lattice, fragment atom list, fractional grid
coordinates, and retained candidate columns to call the existing
`setup_dg_w90_gamma_library`, `assemble_dg_w90_gamma_matrices`,
`run_dg_w90_gamma_library`, and `apply_dg_w90_gamma_transform` in that order on
`comm_fragment`.  Pass `DG_W90_UNCONSTRAINED`; do not fork Wannier90 numerical
code into the new orchestrator.

Use square `num_bands=num_wann` localization with no disentanglement.  Supply
Wannier90's required eigenvalue array as finite zero-valued auxiliary labels
for every retained candidate, including added buffer/projector directions.
Keep the physical `dc_seed_energies` separately in the cache and never return
or consume those auxiliary labels as fragment eigenvalues, occupation input,
or extension energies.  The MPI fixture must include more retained candidates
than DC eigenvalues and prove that physical seed energies remain unchanged.

`cache` is `intent(inout)`.  An invalid cache performs setup/run once and then
stores every returned array and receipt in memory.  A valid cache with the
same fragment ID, basis generation, seed fingerprint, and basis fingerprint
returns those stored arrays with zero Wannier90 calls.  A key mismatch within
the same generation is a hard stale-cache error; a caller starting a new basis
generation must first create a fresh cache, which writes to the new generation
namespace rather than overwriting prior artifacts.

Return both the candidate-space metric compression and the square Wannier90
unitary, plus the coefficients that reconstruct every original DC seed
orbital in the final fragment-WF basis.  Certify the reconstruction in the
fragment metric and certify invariance of the DC occupied projector and
occupation-weighted density.  Task 8 uses this map, padded by zero PW
coefficients, to initialize the occupied-plus-guard `X`; it must not infer the
inverse transformation from WF ordering or centers.

Keep `candidate_rank` and `retained_rank` runtime-sized.  Occupied-plus-guard
selection belongs to the fragment eigensolver state policy in Task 7; this
routine must not contain a material-specific state count.

**Step 4: Write the RED generalized-complement and gauge-invariance fixture**

Construct two overlapping fragment WF blocks with a nonidentity union Gram
matrix, raw plane waves, and nonzero cross-fragment `H` and `S` rows.  Require
the PW complement to satisfy the retained-union condition

```text
W^dagger S (P - W G^+ W^dagger S P) = 0,
G = W^dagger S W,
```

within tolerance.  Apply independent phase, permutation, and full unitary
rotations inside every complete fragment WF block and require invariant:

- union metric rank and accepted basis span;
- unchanged fragment ownership and seed coefficients in the uncompressed
  divided-SCF catalog, plus a separately fingerprinted union-to-complete map;
- complete generalized eigenvalues;
- occupied projector and reconstructed density;
- generalized residual and orthogonality receipts.

Add negative cases for a missing local WF column, duplicate global basis ID,
an indefinite Gram matrix, metric rank loss beyond the declared tolerance,
and insufficient PW/projector/tail coverage.  The negative missing-column case
must fail rather than silently treating individual WF symmetry or centers as
an acceptance gate.

Include a cross-fragment near-null direction.  Require its removal only from
the terminal complete catalog, require all local fragment columns to remain,
and verify seed reconstruction after composing the fragment embedding with the
union-to-complete map.  Reject compression if a seed projector loses rank.

**Step 5: Run the generalized-complement RED test**

Run:

```text
python3 tests/dg/run_dg_hybrid_fragment_wannier_lcfo_mpi.py
python3 tests/dg/run_dg_hybrid_wannier_complement_mpi.py
```

Expected: compile failure on the new generalized projection entry point or an
assertion showing that the old orthonormal-WF formula is invalid for the
fragment-union metric.

**Step 6: Add the generalized fragment-union projection**

Keep the existing exact orthonormal-WF routine for its current callers and add
a separate API:

```fortran
subroutine compute_dg_hybrid_generalized_wannier_projection_tile(comm,&
    global_row_count,row_ids,weights,wannier_values,pw_tile,&
    wannier_fingerprint,packet_fingerprint,first_column,metric_tolerance,&
    coefficients,projected_pw,metric_rank,metric_condition,&
    orthogonality_defect,workspace_peak_bytes,fingerprint,ok,message)

subroutine build_dg_hybrid_complete_union_map(comm,global_row_count,row_ids,&
    weights,uncompressed_basis_values,metric_tolerance,&
    complete_basis_transform,complete_basis_values,metric_rank,&
    metric_condition,projector_fingerprint,ok,message)
```

Assemble the distributed Hermitian Gram matrix, diagonalize it collectively,
reject negative modes, and form the Moore--Penrose inverse only over the
declared retained metric range.  Project with `G^+`; never assume `G=I` for a
union of independently localized fragment blocks.  This first routine leaves
the original WFs in their uncompressed fragment catalogs and returns only the
PW complement and metric receipts.

After every fragment WF+PW catalog is assembled, the second routine checks the
complete union metric.  If it is full rank, return the identity transform so
fragment columns remain unchanged.  If it has cross-fragment numerical null
directions, return a rectangular transform spanning the retained metric
eigenspace and the corresponding terminal complete-basis values.  This map is
used only to congruence-transform the final LCFO `H/S`; it is never used as a
fragment basis.  Near-null removal must be defined by the union metric
eigenspace, not by fragment order or individual WF labels, so the retained
projector is invariant under complete within-fragment unitary rotations.
Certify the retained rank, condition estimate, weighted orthogonality defect,
and unique global row IDs in these numerical routines.

Compare the retained metric projector, not individual complete-transform
columns or their raw fingerprints, across exactly degenerate metric clusters.
Such columns may differ by a harmless unitary gauge while the Hybrid span and
LCFO observables remain invariant.

Route the projected-fragment pipeline through the generalized entry point when
its WFs come from independently localized fragment blocks.  Keep the old
orthonormal fast path only when its caller supplies and fingerprints an
explicit globally orthonormal frame; never choose it merely because every
fragment block is internally orthonormal.  The pipeline, which owns fragment
catalogs, packet-neighbor metadata, buffer rows, and projector graphs, performs
the separate interface/periodic-wrap/nonlocal-projector/tail coverage gates
before it publishes two related products: uncompressed, fragment-owned bases
for divided SCF and the immutable union-to-complete map for terminal LCFO.

```fortran
type s_dg_hybrid_dual_basis_catalog
  logical :: valid
  type(s_dg_hybrid_fragment_basis),allocatable :: fragment_bases(:)
  complex(real64),allocatable :: union_to_complete(:,:)
  integer :: uncompressed_rank
  integer :: complete_rank
  integer(int64) :: fragment_catalog_fingerprint
  integer(int64) :: complete_map_fingerprint
end type
```

The projected-fragment pipeline finalization returns this catalog
collectively.  Fragment entries and their coefficient coordinates are never
rewritten by `union_to_complete`.

Compose each fragment's seed-to-WF map with its embedding into the
uncompressed union and with the union-to-complete map.  The Task 5 integration
fixture must reconstruct the same seed orbital, occupied projector, and
density through both the fragment-local coordinates and the terminal complete
coordinates.  A rank compression that removes any physical seed direction is
a collective failure.

**Step 7: Run GREEN and protected basis tests**

Run:

```text
python3 tests/dg/run_dg_hybrid_fragment_wannier_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_wannier_lcfo_mpi.py
python3 tests/dg/run_dg_hybrid_wannier_complement_mpi.py
python3 tests/dg/run_dg_hybrid_projected_fragment_pipeline_mpi.py
python3 tests/dg/run_dg_hybrid_production_fragment_basis_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
cmake --build build-hybrid-commit -j2
```

Expected: PASS for every configured decomposition; receipts show one
construction Wannier90 call per fragment and generation and no fixed global
state count, and the build completes.

**Step 8: Commit only Task 5**

```text
git add src/gs/dc/dg_hybrid_fragment_wannier.f90 tests/dg/test_dg_hybrid_fragment_wannier_mpi.f90 tests/dg/run_dg_hybrid_fragment_wannier_mpi.py tests/dg/test_dg_hybrid_fragment_wannier_lcfo_mpi.f90 tests/dg/run_dg_hybrid_fragment_wannier_lcfo_mpi.py
git add -p src/common/dg_hybrid_wannier_complement.f90 tests/dg/test_dg_hybrid_wannier_complement_mpi.f90 src/gs/dc/dg_hybrid_projected_fragment_pipeline.f90 tests/dg/test_dg_hybrid_projected_fragment_pipeline_mpi.f90 src/gs/dc/CMakeLists.txt
git diff --cached --check
git diff --cached
git commit -m "feat(dg): localize construction WFs per fragment"
```

### Task 6: Build fragment self blocks from the fixed broken-volume/SIPG payload

**Files:**

- Create: `src/gs/dc/dg_hybrid_divided_operator.f90`
- Create: `tests/dg/test_dg_hybrid_divided_operator_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_divided_operator_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/main_dft.f90:3440-3612`
- Modify: `tests/dg/check_dg_hybrid_continuation_route.py`

**Step 1: Write a RED two-fragment operator fixture**

Create a row-distributed synthetic payload with:

- broken-volume kinetic rows;
- a nonlocal projector whose support reaches a fragment not identified only by
  a Cartesian face;
- local-potential rows;
- SIPG self and cross-fragment rows; and
- a nonidentity metric.

Require `extract_dg_hybrid_fragment_self_block` to return the exact `H_ff` and
`S_ff`, including all diagonal/self pieces of the SIPG and nonlocal operators.
Require `compose_dg_hybrid_complete_rows` to reconstruct the direct full
reference in the uncompressed union with every cross-fragment contribution
exactly once, then apply the Task 5 union-to-complete congruence transform to
both `H` and `S`.  Use a rectangular transform in one case and require that it
does not alter any extracted local self block.  Check Hermiticity, positive
retained metric rank, and rank-decomposition invariance.

Extend the continuation route checker to require construction of the basis
directory, interior rows, projector support, faces, and frozen variational
payload before the divided/reference branch.  Require both branches to consume
the same immutable payload fingerprint and forbid a second payload freeze in
either branch.

**Step 2: Run RED**

Run:

```text
python3 tests/dg/run_dg_hybrid_divided_operator_mpi.py
python3 tests/dg/check_dg_hybrid_continuation_route.py --interface-only
```

Expected: compile failure because `dg_hybrid_divided_operator` is absent and
source-contract failure because payload construction is still branch-local.

**Step 3: Implement narrow matrix APIs**

```fortran
subroutine extract_dg_hybrid_fragment_self_block(comm_fragment,fragment_id,&
    row_ids,basis_fragment,fixed_payload,local_rows,hff,sff,ok,message)

subroutine compose_dg_hybrid_complete_rows(comm,row_ids,fixed_payload,&
    local_rows,union_to_complete,hamiltonian_rows,metric_rows,&
    operator_fingerprint,ok,message)
```

Use the frozen basis directory and row IDs; never infer a fragment from a
rank-local offset.  Fragment extraction always uses the uncompressed,
fragment-owned Task 5 catalog.  The fixed payload contributes kinetic,
nonlocal, metric, and SIPG rows in that uncompressed union.  Only `local_rows`
changes with density.  Complete composition applies
`T^H H_union T` and `T^H S_union T` with the immutable
`union_to_complete=T`; it never rewrites fragment basis IDs or caches.

**Step 4: Share payload construction between divided and reference routes**

Move the already accepted basis-directory, interior, nonlocal, face, and
`freeze_dg_hybrid_variational_payload` construction before the
`yn_dg_hybrid_continuation_scf`/`yn_dg_hybrid_divided_scf` branch.  Keep one
immutable payload and one fingerprint contract for both routes.

**Step 5: Run GREEN**

Run:

```text
python3 tests/dg/run_dg_hybrid_divided_operator_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_support_redistribution_mpi.py
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/check_dg_hybrid_continuation_route.py --interface-only
cmake --build build-hybrid-commit -j2
```

Expected: all PASS and the build completes.

**Step 6: Commit**

```text
git add src/gs/dc/dg_hybrid_divided_operator.f90 tests/dg/test_dg_hybrid_divided_operator_mpi.f90 tests/dg/run_dg_hybrid_divided_operator_mpi.py
git add -p src/gs/dc/CMakeLists.txt src/gs/main_dft.f90 tests/dg/check_dg_hybrid_continuation_route.py
git diff --cached --check
git diff --cached
git commit -m "feat(dg): project fixed DG fragment self blocks"
```

### Task 7: Add bounded warm-started fragment subspace updates and dynamic state extension

**Files:**

- Create: `src/gs/dc/dg_hybrid_fragment_subspace.f90`
- Create: `tests/dg/test_dg_hybrid_fragment_subspace_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_fragment_subspace_mpi.py`
- Modify: `src/gs/dc/dc_fragment_occupation.f90`
- Modify: `tests/dg/test_dc_fragment_occupation_mpi.f90`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write RED tests for a successful bounded update**

Keep the current strict `solve_dg_hybrid_block_cg` contract intact: reaching
its iteration cap without its eigensolver target must still return failure.
Add a separate fixture for the density-SCF update API using distributed
generalized Hermitian problems.  Require:

- `maximum_steps=3` performs at most three `[X,R,P]` Rayleigh--Ritz updates;
- a finite, metric-orthonormal, improved state is returned with
  `advanced=.true.` even when `eigensolver_converged=.false.`;
- the next density epoch reuses both accepted `X` and prior direction `P`;
- a supplied diagonal/preconditioner callback is applied to `R`; an identity
  callback is allowed only when explicitly supplied by a fixture;
- the trial dimension never exceeds `min(global_count,3*nstate)`, and the
  fixture starts from a genuine occupied-plus-guard subspace rather than
  setting `nstate=global_count` automatically;
- the reported residual is the largest statewise distributed two-norm
  `||H C_j-epsilon_j S C_j||_2 / max(1,|epsilon_j|,||H C_j||_2)`, not a local
  max norm;
- residual below `dg_dc_gs_intermediate_orbital_tolerance` reports early
  eigensolver convergence;
- a trial exceeding the entry residual by more than
  `dg_dc_gs_allowed_residual_growth` rolls back to the safe entry state and
  still reports a valid bounded update; and
- when every trial is rejected but the entry state remains finite,
  rank-preserving, metric-orthonormal, and inside the growth bound, return
  `ok=.true.`, `advanced=.false.`, and
  `stop_reason='safe_entry_retained'`; and
- non-finite algebra, an indefinite/rank-deficient metric, failed operator
  callbacks, or no safe entry state remain collective hard failures.

Compare cold and warm second epochs on the same perturbed operator and require
the warm path to use its saved direction and achieve no worse certified
residual.  Also assert that the density outer loop, not this bounded routine,
owns the final convergence decision.

**Step 2: Run the bounded-update RED test**

Run:

```text
python3 tests/dg/run_dg_hybrid_block_cg_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_subspace_mpi.py
```

Expected: the strict fixture remains PASS and the new fixture fails to compile
because `dg_hybrid_fragment_subspace` and the bounded API are absent.

**Step 3: Implement the separate `[X,R,P]` update API**

Expose persistent state keyed by fragment and basis generation:

```fortran
type s_dg_hybrid_fragment_subspace_state
  integer :: fragment_id
  integer :: basis_generation
  integer :: state_count
  integer(int64) :: basis_fingerprint
  integer(int64) :: metric_fingerprint
  complex(real64),allocatable :: vectors(:,:)
  complex(real64),allocatable :: directions(:,:)
end type

subroutine advance_dg_hybrid_fragment_subspace(comm,global_count,row_ids,&
    fragment_id,basis_generation,basis_fingerprint,metric_fingerprint,&
    apply_h,apply_s,apply_preconditioner,maximum_steps,intermediate_tolerance,&
    orthogonality_tolerance,allowed_residual_growth,state,eigenvalues,&
    iterations,relative_residual,eigensolver_converged,advanced,&
    stop_reason,workspace_peak_bytes,fingerprint,ok,message)
```

At each step call the explicit preconditioner on the raw residual, then form a
metric-orthogonalized trial block from accepted `X`, preconditioned `R`, and
saved direction `P`; solve only that small generalized problem.  Retain no
more than `nstate` vectors and the new conjugate directions.  Validate row
ownership, every callback result, all scalar controls, and the nonzero
basis/metric fingerprints collectively.  Compute the normalized distributed
two-norm receipt explicitly.
Return the best safe state at the step cap.  A cap is `stop_reason='step_cap'`,
not solver failure.  Returning the unchanged safe entry is also a successful
bounded density step with `advanced=.false.`; only `ok=.false.` prevents
density reconstruction.

Do not weaken or silently redirect the strict full-convergence entry point.
Dense fragment diagonalization remains available only when an explicit
reference fixture or reference backend calls it; it is not automatic recovery
for a production bounded-update failure.

**Step 4: Write RED tests for runtime-sized occupied-plus-guard extension**

Extend the occupation fixture and the new subspace fixture with unequal
fragment state counts and a terminal degenerate shell.  Require an extension
mask per fragment when:

- the total represented maximum occupation capacity is below the requested
  electron count, before invoking the occupation kernel; in that case mark
  every non-exhausted fragment to add its next complete shell; or
- at finite temperature, the charge in the final represented degenerate shell
  exceeds `electron_tolerance/nfragment`; or
- at zero temperature, an occupied or Fermi-intersecting terminal shell is not
  represented completely.

Extension first adds the complete next unused DC-seed eigenspace, ordered by
the saved seed eigenvalue and mapped through Task 5's seed-to-WF coefficients.
After those eigenspaces are exhausted, it adds a complete kinetic-energy shell
from the windowed-PW catalog together with still-unused projector-support
directions, S-projects the whole candidate pool against `X`, and solves the
small `H/S` problem only within that pool.  Stable global candidate IDs order
equal pools but never select a named localized basis column inside a degenerate
pool.  It must not use per-column `H_ii/S_ii`, add a fixed number such as 384,
stop at a material-specific energy, or use individual-WF symmetry as the
state-count gate.  If a requested fragment has no remaining candidate, all
ranks fail collectively with an insufficient-spectrum diagnostic.

Apply independent phase, permutation, and full unitary rotations to each
fragment construction-WF block and require the same extension shell, accepted
state projector, common-mu density, and residual receipts.  This test prevents
the short density SCF from depending on the arbitrary unconstrained Wannier90
gauge.

When the state count grows, require the old `X` columns to remain bitwise
unchanged before reorthogonalization, the new columns to be S-orthogonalized
against them, old valid `P` history to be retained for old columns, and new
direction columns to be initialized to zero.  A basis-generation or metric
fingerprint change must invalidate the whole cache instead of embedding stale
coefficients.

**Step 5: Run the extension RED tests**

Run:

```text
python3 tests/dg/run_dc_fragment_occupation_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_subspace_mpi.py
```

Expected: compile or assertion failure because capacity preflight, tail masks,
gauge-invariant whole-shell extension, and cache embedding are absent.

**Step 6: Implement the extension policy and cache transitions**

Add `assess_dc_fragment_occupation_capacity`, which collectively compares
`sum(wspin*core_norm)` with the target before calling the occupation kernel.
When capacity is insufficient and any fragment still has unused candidates,
return `needs_extension=.true.` for every such fragment and no occupations;
the caller adds one complete next shell per marked fragment and repeats.  If
all fragments are exhausted, fail collectively.  This prevents the existing
kernel's correct insufficient-capacity error from preempting a possible
dynamic extension.

```fortran
subroutine assess_dc_fragment_occupation_capacity(comm,core_norms,&
    representative_mask,can_extend,wspin,expected_electrons,tolerance,&
    capacity_sufficient,needs_extension,ok,message)
```

After capacity is sufficient, add an optional `needs_extension(:)` tail result
to `determine_dc_fragment_occupations` without changing existing caller
semantics when it is absent.  A successful common-mu solve may still request
the per-fragment finite-temperature or zero-temperature terminal-shell
extensions above; repeat spectrum/occupation determination until those masks
clear.

In `dg_hybrid_fragment_subspace`, provide a deterministic extension routine
that accepts the saved seed-eigenspace catalog and the globally identified
PW/projector candidate pools, chooses the whole next invariant shell within
the existing energy comparison tolerance, embeds the accepted cache, and
records source kind, old/new state counts, and shell bounds in the receipt.
The production caller will repeat extension until every fragment mask is false
or candidates are exhausted.

**Step 7: Run GREEN and protected occupation/solver tests**

Run:

```text
python3 tests/dg/run_dg_hybrid_block_cg_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_subspace_mpi.py
python3 tests/dg/run_dc_fragment_occupation_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_solver_mpi.py
python3 tests/dg/run_dg_hybrid_occupation_policy_mpi.py
cmake --build build-hybrid-commit -j2
```

Expected: PASS for 1/2/4 ranks where configured; the bounded path accepts a
safe three-step result, all state counts are derived at runtime, and the build
completes.  The existing strict block-CG source and fixture remain unchanged;
their fresh PASS protects the converged-or-fail reference contract.

**Step 8: Commit only Task 7**

```text
git add src/gs/dc/dg_hybrid_fragment_subspace.f90 tests/dg/test_dg_hybrid_fragment_subspace_mpi.f90 tests/dg/run_dg_hybrid_fragment_subspace_mpi.py
git add -p src/gs/dc/dc_fragment_occupation.f90 tests/dg/test_dc_fragment_occupation_mpi.f90 src/gs/dc/CMakeLists.txt
git diff --cached --check
git diff --cached
git commit -m "feat(dg): add bounded fragment subspace updates"
```

### Task 8a: Certify the fixed-frame preconditioner before production wiring

The user-approved 2026-09-03 amendment replaces the gauge-dependent current-WF
diagonal preconditioner. Execute
`docs/plans/2026-09-03-hybrid-fixed-frame-preconditioner.md` first, including its
checkpoint. Task 8 below remains unfinished until its production route is
actually connected and all its tests pass.

### Task 8: Connect the production divided loop to fragment-local construction and bounded updates

**Progress, 2026-09-03 (Task 8 remains incomplete):** Input-side
`dg_hybrid_fragment_cg_steps` (default 3, range 1--256, broadcast and log),
mutually exclusive Hybrid route flags, and divided-route positive finite PW
cutoff validation are implemented. A separate
`measure_dg_hybrid_fragment_subspace` API now supports the zero-remaining-budget
case after extension: it validates S orthonormality and computes current
Rayleigh values/residuals without CG, Ritz rotation or X/P changes. Stale
basis/metric provenance retains the existing cache-invalidation semantics.
Values remain in coefficient-column order, so production occupation packing
must sort energies and core norms together and undo that permutation for the
returned occupations. The old production main route is not yet replaced;
construction from DC seeds, remaining-budget accounting, dynamic extension,
fixed-frame callback wiring, dual-catalog integration and route tests below
remain required. No Si64 calculation has been started.

**Progress, 2026-09-04 (Task 8 remains incomplete):** The common occupation
adapter now accepts explicit `allow_unordered=.true.` for measured/bounded
states. It collectively agrees this option, stably sorts gathered energy/core
weight pairs for occupation and terminal-tail decisions, then returns
occupations in the original coefficient order. The divided production call
opts in; conventional callers retain the original sorted-spectrum requirement.
Tests cover crossing measured Rayleigh values with unchanged X/P, unequal core
weights, padding, degenerate terminal shells, finite temperature, invalid
representative data and rank-disagreeing policy. The pending construction and
bounded-update main-loop replacement described above is still required.

**Coordinate handoff checkpoint, 2026-09-04 (Task 8 remains incomplete):**
`export_dg_hybrid_fragment_coordinates` validates the saved fragment-Wannier
cache against the expected seed/basis/generation and rank-local grid layout,
then exports replicated WF-coordinate maps: the actual saved `U^dagger` for
the fixed reference frame and the saved physical DC seed coefficients.
It never reruns Wannier90 or changes the cache, and publishes neither output
on failure. The caller must still append PW identity/zero blocks, distribute
coefficient rows, and connect these maps to the production divided loop.
MPI fixtures verify complex frame recovery, DC orbital reconstruction,
corruption/generation/rank-layout rejection and unchanged W90 call counts.
This helper does not replace the old main-route construction or complete Task 8.
Verification: fragment-Wannier MPI fixture passed on 2/4/8 ranks; fixed-frame
preconditioner and bounded subspace fixtures passed on 1/2/4/8 ranks;
`cmake --build build-hybrid-release -j 4` completed successfully. Read-only
review found no Critical/Important issues; its suggested caller-side
fingerprint/grid-layout mismatch tests were added and passed. No Si64 run.

**Seed-to-CG checkpoint, 2026-09-04 (Task 8 remains incomplete):**
`initialize_dg_hybrid_fragment_subspace` now selects the occupied spectral
prefix plus a caller-specified guard count and optional physical DC energy
window, closing the boundary degenerate shell. It selects mapped physical DC
seed columns, not named WF columns, certifies/normalizes the current S metric,
and initializes zero search history. It returns original seed column indices
for the later extension catalog. The optional energy window is not the PW
kinetic cutoff, nor does this initial selection certify a tail or symmetry.
Insufficient saved guard inventory saturates at the saved seed count; later
capacity/tail-driven extension remains mandatory. Initial selection that
exhausts or exceeds the fragment basis is rejected without truncating a shell
or changing an accepted state. The global-count limit and workspace bounds
are checked collectively before publication.

The actual saved complex Wannier map is now exercised through the initializer
with distributed reversed coefficient rows, zero-row ranks, and zero padded
PW components; reconstructed initial orbitals match the original DC seeds.
Subspace MPI tests pass on 1/2/4/8 ranks, the Wannier handoff on 2/4/8,
preconditioner regression on 1/2/4/8, and the release build succeeds.
The full-basis rejection regression was observed failing before its fix.
Main-route construction, PW/catalog assembly, per-epoch budget accounting,
capacity/tail extension and the production callback replacement are still
pending. No production SCF or Si64 run was started at this checkpoint.

**Epoch-budget checkpoint, 2026-09-04 (Task 8 remains incomplete):**
`advance_dg_hybrid_fragment_epoch` owns a persistent private per-fragment
budget around the bounded updater and measurement-only path. It replenishes
the budget only on a newer density epoch; repeated calls, including after
state-count extension, share the remaining steps. Same-epoch limit changes,
older/rank-disagreeing epochs and fragment/basis/metric/communicator-layout
mismatches are rejected. Attempted updates are charged even when the safe
entry is retained. Once exhausted, calls only measure the current state and
leave X/P unchanged; a callback/metric failure during execution poisons the
budget instead of silently retrying. The shifted preconditioner callback is
mandatory. The caller must keep one budget object per fragment, use the true
outer-density epoch (not the extension pass index), and pass unchanged H/S
within that epoch. Operator-provenance checks remain the responsibility of
the existing fixed-frame preconditioner callback.

MPI tests exercise three consumed steps, repeated measurement, actual whole
shell extension from two to four states without replenishment, next-epoch
renewal, invalid epoch/limit agreement and failure without implicit retry.
This is still a reusable callback component, not replacement of the old
production main route. That route's construction/catalog/occupation/density
wiring and invocation of this budgeted updater remain pending.
Verification: subspace/epoch and preconditioner MPI fixtures passed on
1/2/4/8 ranks, saved-Wannier handoff passed on 2/4/8, and the release build
completed successfully. Review found no Critical/Important issues. Both minor
suggestions were addressed: measurement workspace is recorded, and an early
target leaving three unused steps followed by an extended-state update is
covered. The workspace regression failed before its fix. No Si64 run.

**Current-coefficient core-weight checkpoint, 2026-09-04 (Task 8 remains
incomplete):** `measure_dg_hybrid_fragment_core_norms` computes each current
state's weighted norm only on the core, in coefficient-column order, without
H/S application or diagonalization. It shares the private reconstruction
worker with the unchanged public density API, using zero occupations to
obtain the independent state norms. It retains the existing within-fragment
replicated point layout and distributed coefficient-row ownership contract;
the resulting norm vector is replicated, not summed again across those ranks.
Only the fragment representative contributes it to total-communicator capacity
and occupation decisions. Invalid/nonfinite or overflowing coefficients do not
publish a norm vector.

The MPI fixture now covers 1/2/4/8 ranks including zero coefficient-row ranks,
complex mixed non-eigenstates, buffer exclusion and unchanged coefficients.
It passes current norms through capacity preflight and unordered common
occupation determination into density reconstruction, recovering the same
electron count without another operator application. This closes the missing
coefficient-to-occupation interface, but the production main callback still
uses the old spectrum route. Direct fragment-Wannier construction, catalog
assembly and replacement of that callback remain pending; no Si64 run.
Verification: fragment solver/core-weight fixture and bounded subspace fixture
passed on 1/2/4/8 ranks; divided SCF passed on 1/2/4. Existing divided DC
controls and LCFO route source checks passed (these do not yet certify the
new Task 8 production route). The release build completed successfully.
Read-only review found no Critical/Important issues.

**Occupation-loop checkpoint, 2026-09-04 (Task 8 remains incomplete):**
`run_dc_fragment_occupation_epoch` connects refresh, represented-capacity
preflight, common occupations, terminal-tail masks and invariant-shell
extension. All passes receive the same outer epoch. Capacity shortage extends
every non-exhausted fragment; otherwise only tail-marked fragments extend.
Callbacks must grow the state count within the fixed basis rank, and the next
refresh must match the reported count. A persistent budget remains owned by
each fragment callback, so the controller cannot replenish it on an extension
pass. The final occupation vector is published in local coefficient-column
order only after all tails clear. Exhausted tails, non-progress, callback
failure and rank-disagreeing fragment inventories/data fail collectively.
This driver requires one fragment per rank, permits multiple ranks per
fragment, and requires exactly one representative per fragment.

Tests cover capacity then tail extension, selective fragment extension,
zero/finite temperature, failure/exhaustion, and representative-data mismatch.
The real bounded updater plus real two-to-four-state shell extension is
connected in the subspace fixture: two occupation passes consume exactly
three CG updates in total. Occupation tests pass on 1/2/4 ranks; bounded
subspace and core-weight fixtures pass on 1/2/4/8; divided SCF passes on
1/2/4; the release build passes. Review found no Critical/Important issues;
its metadata extent-bound suggestion was incorporated.
Production main still needs the direct seed/Wannier/catalog construction and
concrete H/S/preconditioner/density callback wiring. This checkpoint does not
claim that the old main route has been replaced, and does not start Si64.

**DC tensor packing checkpoint, 2026-09-04 (Task 8 remains incomplete):**
`pack_dg_hybrid_fragment_dc_seed` converts the saved real Gamma, single-spin
7D orbital tensor and physical spectra into state-major fragment WF inputs.
It retains the entire periodic core+buffer cell and excludes communication
halos only; it applies no taper, truncation, normalization or new boundary
condition. Cell-local x-fast IDs and fractional coordinates are explicit.
Unique owned points must cover the complete cell, with finite owned values
and rank-consistent spectra/volume. Outputs are published only on success.
The supplied spatial communicator must hold the complete replicated orbital
inventory: orbital-partitioned tensors are rejected pending explicit
redistribution, not silently interpreted as a complete fragment seed.
This helper has not yet been connected to main or the production catalog.

Tests cover 3D indexing, empty owners, NaN communication halos, missing or
duplicate owned points, malformed/unallocated tensors and rank-disagreeing
spectra. Review identified a NaN-volume ordered-comparison trap; a regression
reproduced it, and finite validation is now a separate collective before
comparisons. The complete fragment Wannier fixture passes on 2/4/8 ranks;
the release build and diff whitespace check pass after the fix. No material
calculation or new worktree was started.

**Orbital redistribution checkpoint, 2026-09-04 (Task 8 remains incomplete):**
The DC seed packer now accepts an explicit `orbital_comm` (the production
communicator to supply is `info%icomm_o`). Each orbital group must be a
subgroup of the fragment communicator, share the same owned spatial slab,
and own every seed orbital exactly once. Only owned spatial values are
assembled; NaN communication halos are never read. A complete slab tensor
is temporarily replicated within that orbital group, not across the full
system. Only the orbital-group root publishes its slab to the fragment WF
builder, so points are not duplicated. Empty orbital owners are supported;
without the optional communicator the previous full-orbital contract remains.
This does not relax checkpoint MPI-size/rank--fragment reuse restrictions.

The new API first failed its compile regression, then passed mixed orbital
and spatial decompositions on 2/4/8 total ranks. Missing/out-of-range orbitals,
inconsistent slabs and nonfinite owned values fail collectively. Review found
the foreign-member communicator hazard: the specific rejection regression
failed before local group-membership validation was added. All fragment WF
tests, the release build and diff check pass after that fix. Main-route
construction/callback wiring is still pending; no Si64 calculation was run.

**Direct construction entry checkpoint, 2026-09-04 (Task 8 remains incomplete):**
`build_dg_hybrid_fragment_wannier_from_dc_seed` now composes the orbital-aware
DC packer and the fragment Wannier builder. Candidate grid IDs must match the
unique packed layout exactly. Packing failures are synchronized across all
fragments before entering construction, so one bad fragment cannot leave its
peers waiting in the next total-level collective. Existing cache publication
remains transactional. The MPI fixture exercises the complete entry, physical
seed reconstruction from WF coefficients, one construction per generation,
cache reuse and one-rank input/layout failure without cache mutation or W90
re-entry. External W90 calls use the fixture's stubs; this is not a material
localization benchmark. The fragment WF fixture passes on 2/4/8 ranks and
the release build and whitespace check pass. Review found no Critical/Important issues.

`check_dg_hybrid_fragment_wannier_route.py` is now present and intentionally
RED: main still dispatches through the old complete-system construction.
It requires a separated divided entry using direct DC construction and the
bounded update/occupation kernels, and rejects legacy construction/solver
calls in that entry. The old shared main routine consumes preliminary LCFO
results throughout its construction block; inserting an isolated new call
would not safely replace that dependency. Production separation and the
uncompressed-catalog/payload/callback wiring remain the next work, including
mapping raw DC cell order with `dc%jxyz_tot` rather than assuming a centered
core in the legacy buffer layout. No main-route switch or Si64 run is claimed.

**Physical grid mapping checkpoint, 2026-09-04 (Task 8 remains incomplete):**
`map_dg_hybrid_fragment_dc_grid` maps packed cell IDs through the authoritative
`dc%jxyz_tot` and marks core ownership with raw indices `1:core_shape`, matching
the existing DC density aggregation. It preserves row order, periodic buffer
tails and WF values; it neither centers the core nor truncates/tapers a WF.
It validates rank-consistent geometry and maps, complete unique raw-cell
coverage, physical index ranges and nonduplicated core-axis coordinates.
Physical aliases in the buffer are permitted. This is a per-fragment map;
unique physical-core ownership across fragments still requires downstream
certification. Outputs remain unallocated on failure.

MPI tests include reversed/distributed rows, empty owners, periodic wrapping,
invalid mappings and IDs, and duplicate core coordinates. The direct DC-to-WF
fixture now maps reconstructed orbitals back to the physical core and compares
both density and electron count against the original DC seed. The 2/4/8-rank
fixture, release build and whitespace check pass; review found no
Critical/Important issues. Main and the old LCFO construction route have not
been modified; the production connection gate remains RED.

**Files:**

- Modify: `src/io/salmon_global.f90:480-545`
- Modify: `src/io/inputoutput.f90:630-700,1155-1220,1885-1960,2955-3000,3130-3170`
- Modify: `src/gs/main_dft.f90:3613-3668,4867-4995`
- Create: `tests/dg/check_dg_hybrid_fragment_wannier_route.py`
- Modify: `tests/dg/check_dg_hybrid_divided_dc_controls.py`
- Modify: `tests/dg/check_dg_hybrid_divided_lcfo_route.py`
- Modify: `tests/dg/check_dg_hybrid_localization_first_inputs.py`
- Modify: `tests/dg/test_dg_hybrid_divided_scf_mpi.f90`
- Modify existing untracked file: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in`

**Step 1: Make the route checks RED**

Require the divided branch to use
`build_dg_hybrid_fragment_wannier`, the generalized fragment-union PW
projection, `extract_dg_hybrid_fragment_self_block`, the bounded fragment
subspace update, and the split spectrum/occupation/density phases.  Forbid:

- preliminary `dc_lcfo` or another complete generalized solve before divided
  density convergence;
- construction-WF Wannier90 setup/run on `dc%icomm_tot`;
- any construction-WF Wannier90 setup/run from inside the density loop;
- production calls to `apply_dg_hybrid_divided_fragment_hpsi`, the identity
  metric callback, or dense fragment eigensolving; and
- a literal or derived material-specific retained-state count such as 384.

Require exactly one construction setup/run per fragment communicator and
basis generation outside the loop.  The later certified RT localizer remains
a separate complete-system operation and is required only when an RT
checkpoint is requested; the source check must not confuse it with the
construction-WF call.

Require one local-update budget per fragment communicator and outer density
iteration, shared across that fragment's state-extension passes, and require
the initial production state count to come from the occupied-plus-guard policy
rather than the full fragment basis.
Require local `H_ff/S_ff`, seed coefficients, and density reconstruction to use
only the uncompressed fragment catalog.  Require the union-to-complete map to
remain immutable and unused until complete row composition after divided
convergence.

Require a finite fixed-frame production preconditioner constructed from each
fragment self block and the saved pre-localization-to-WF map. An identity
callback remains fixture-only. Test the actual saved Wannier90 transform,
including complex phases/order, with `F=B Q`; passing identity as Q after
localization is forbidden. The PW complement uses the unchanged fixed PW
reference. Check physical preconditioned residuals under independent full WF
unitary rotations, not only phase/permutation changes.

Require `yn_dg_hybrid_divided_scf`, `yn_dg_hybrid_continuation_scf`, and
`yn_dg_hybrid_scf` to be mutually exclusive.

Extend the divided MPI fixture so 1/2/4-rank layouts reconstruct the same
unique-core density and electron count from the same fragment projectors.
Overlapping buffer points must contribute no physical density ownership, and
duplicated or missing core ownership must fail collectively.

Add `dg_hybrid_fragment_cg_steps` as an integer input with default `3`; require
collective broadcast, logged value, and rejection outside `[1,256]`.  Preserve
the existing user-provided positive finite `wannier_pw_cutoff`; do not replace
it by a hidden fixed cutoff or reinterpret it as a symmetry guarantee.

**Step 2: Run the source/input RED tests**

Run:

```text
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_fragment_wannier_route.py
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
```

Expected: FAIL because the old divided callback still calls ordinary `hpsi`
with an identity metric, performs complete-cell construction before the loop,
the bounded-step input is absent, and the current production density fixture
does not yet enforce decomposition-invariant unique-core ownership.

**Step 3: Wire construction directly from the reusable DC seed**

After the exact MPI-rank-count and rank--fragment mapping checks accept the DC
seed checkpoint:

1. take each fragment's saved `rwf`, spectrum, occupation, density, and
   potential directly from the seed payload;
2. run the Task 5 unconstrained construction localization once on that
   fragment communicator and retain every resulting column;
3. add plane waves selected dynamically by the user input
   `wannier_pw_cutoff`, project them with the generalized union metric, and
   certify metric rank plus buffer/projector/tail coverage;
4. freeze both the uncompressed fragment catalog and the separately
   fingerprinted union-to-complete transform, then build the common
   variational payload in uncompressed union coordinates; and
5. initialize each fragment's occupied-plus-guard `X` from
   `cache%dc_seed_coefficients_in_wannier`, append zero coefficients for the PW
   complement, and use the seed spectrum/occupations to choose its runtime
   columns without a complete-cell eigensolve.

There is no preliminary full LCFO and no complete-cell construction-WF
Wannier90 operation.  A fragment/basis-generation fingerprint mismatch is a
hard error and triggers neither silent gauge repair nor global fallback.

**Step 4: Wire the bounded divided-density callbacks**

At each divided iteration:

1. update the total potential from the current unique-core density;
2. assemble the density-dependent local rows;
3. extract each fragment `H_ff` and `S_ff` from the same fixed payload used by
   final LCFO, always in the uncompressed fragment coordinates;
4. advance the fragment cache by at most
   `dg_hybrid_fragment_cg_steps` using the warm `[X,R,P]` state and a
   Task 8a's shifted callback: `Q D_j^-1 Q^dagger r_j`, with diagonals
   `diag(Q^dagger H_ff Q)-epsilon_j*diag(Q^dagger S_ff Q)` in the fixed
   pre-localization reference frame, not the current WF coordinates;
5. preflight the total represented occupation capacity; when it is
   insufficient, extend every non-exhausted fragment by its next invariant
   shell and repeat steps 4--5 without resetting the remaining update budget;
6. only after capacity is sufficient, determine current occupations with one
   common DC chemical potential;
7. extend only fragments whose terminal occupation/tail mask requests the
   complete next shell, repeating steps 4--6 until every mask clears;
8. reconstruct the unique-core density from current uncompressed fragment
   coefficients and
   occupations; and
9. use the shared DC convergence, electron-count, potential-refresh, and
   mixing routines.

Initialize
`remaining_fragment_steps(fragment)=dg_hybrid_fragment_cg_steps` for every
fragment at the start of each outer density iteration.  Decrement a
fragment's own counter by every `[X,R,P]` update across that fragment's
state-extension passes; never reset it at either step 5 or step 7.  Budgets for
different fragment communicators are independent and may advance
concurrently.  Appending a
shell performs deterministic S-orthogonalization plus finite Rayleigh-quotient
evaluation for the new directions.  If the update budget is exhausted, those
safe estimates participate in the repeated occupation/tail decision and the
new directions receive their first LOBPCG update in the next density epoch.
Thus extension may repeat without permitting more than the user-requested
number of local updates in one outer iteration.

Use the Task 8a denominator floor with its roundoff-zero convention and reject
nonfinite inputs/outputs. Log the physical reference, coordinate map and
preconditioner fingerprints with the local operator epoch. Require the exact
epoch and H/S provenance on every application so a stale preconditioner cannot
be reused after a potential change. Pass the current Rayleigh values from the
bounded updater; do not freeze them at the beginning of the three-step call.

The local loop contains no complete generalized eigensolve and no Wannier90
call.  Reaching the local three-step cap is accepted when the updater returns a
finite safe state; the outer density criterion decides SCF convergence.  The
old ordinary-`hpsi`, strict-CG, and dense-solve paths may remain for isolated
fixtures or an explicitly selected small/reference backend, but are not
automatic production recovery.

The production initializer must choose `nstate` from the DC occupied count,
guard/tail policy, requested window, and degenerate-shell closure.  It must not
initialize `nstate=fragment_basis_count` merely because the basis is
available.  Equality is permitted only after recorded dynamic shell extension
has exhausted the fragment metric rank; at that point the tail must clear or
the calculation fails collectively.

Log per-iteration fragment state counts, extension events, common chemical
potential, electron defect, local steps, normalized residual, rollback/stop
reason, and basis/operator fingerprints.  These receipts must be finite and
decomposition-consistent.

**Step 5: Run GREEN and build**

Run:

```text
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_fragment_wannier_route.py
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
python3 tests/dg/run_dg_hybrid_fragment_wannier_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_wannier_lcfo_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_subspace_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/run_dg_hybrid_divided_operator_mpi.py
python3 tests/dg/run_dc_fragment_occupation_mpi.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
cmake --build build-hybrid-commit -j2
```

Expected: all PASS; source receipts show no preliminary complete solve, one
construction localization per fragment/generation, and no complete eigensolve
or localization call inside the divided SCF loop.

**Step 6: Commit**

```text
git add -N tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in
git add -p tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in
git add tests/dg/check_dg_hybrid_fragment_wannier_route.py
git add -p src/io/salmon_global.f90 src/io/inputoutput.f90 src/gs/main_dft.f90 tests/dg/check_dg_hybrid_divided_dc_controls.py tests/dg/check_dg_hybrid_divided_lcfo_route.py tests/dg/check_dg_hybrid_localization_first_inputs.py tests/dg/test_dg_hybrid_divided_scf_mpi.f90
git diff --cached --check
git diff --cached
git commit -m "feat(dg): run bounded fragment-local density SCF"
```

### Task 9: Extract one terminal LCFO certification and complete-v3 publication

**Files:**

- Create: `src/gs/dc/dg_hybrid_lcfo_finalization.f90`
- Create: `tests/dg/test_dg_hybrid_lcfo_finalization_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_lcfo_finalization_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/main_dft.f90:5060-5730`
- Modify: `tests/dg/check_dg_hybrid_divided_lcfo_route.py`
- Modify: `tests/dg/check_dg_hybrid_continuation_route.py`

**Step 1: Write RED finalization tests**

Use a small complete eigensystem with an occupied block and a degenerate
unoccupied boundary cluster.  Require one call to the complete generalized
solver, current-spectrum occupations, extension from the requested cutoff to
the first complete degenerate cluster that closes the requested physical
window, physical occupied/window symmetry checks against the symmetry group
of the actual complete system, the complete fragment-local construction span,
and a complete version-3 payload.  Individual construction WFs and individual
fragment blocks need not transform symmetrically; their nonclosure remains a
diagnostic and never rejects an otherwise complete LCFO subspace.  Reject
missing operator, density, selection, pseudopotential, ownership, fragment-WF,
or generalized-metric fingerprints.

Use the terminal complete catalog and the congruence-transformed `H/S` from
Task 6.  Require a mismatched union-to-complete fingerprint, an uncompressed
row extent, or accidental fragment-coordinate coefficients to fail before the
complete solve.  Run `prepare_rt_checkpoint=.false.` and `.true.` cases;
both solve once, but only the latter invokes the distinct complete-system RT
localizer and prepares a version-3 payload.

Extend the divided and continuation route checkers at the same time.  Require
both routes to call the extracted finalizer, require the zero-refinement
divided branch to reach it exactly once after divided convergence, forbid the
old version-2 publisher there, and keep the reference continuation scheduler
outside the extracted terminal routine.

**Step 2: Run RED**

Run:

```text
python3 tests/dg/run_dg_hybrid_lcfo_finalization_mpi.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
```

Expected: compile failure because the finalization module is absent and route
failures because terminal logic is still inline.

**Step 3: Extract the terminal logic**

Expose a routine that performs no density iteration:

```fortran
subroutine finalize_dg_hybrid_lcfo_once(comm,operator_epoch,row_ids,&
    hamiltonian_rows,metric_rows,fixed_payload,complete_basis_values,&
    union_to_complete_fingerprint,grid_ids,grid_weights,&
    physical_symmetry,selection,provenance,prepare_rt_checkpoint,solve_complete,&
    result,lcfo_density,checkpoint_payload,receipt,ok,message)
```

Move, without weakening, the current continuation logic for:

- `solve_dg_hybrid_generalized_complete_once`;
- `derive_dg_hybrid_occupation_policy`;
- occupied density/projector reconstruction;
- requested/certified/proof energy-window selection;
- physical occupied/window symmetry certification;
- certified RT basis localization, only when an RT checkpoint is requested;
- component and energy receipts; and
- complete-v3 payload authentication/publication preparation.

The requested energy is a user cutoff, not itself a promise of a closed
subspace.  Extend above it through the first complete boundary cluster needed
by the physical LCFO window policy.  This extension uses LCFO eigenvalues and
the complete-system symmetry action, not WF centers or a fixed orbital count.
If the input supplies no valid symmetry operations for the actual system,
retain spectral cluster closure and report symmetry as unavailable rather than
inventing a higher-symmetry fragment gate.

Do not move the continuation stage scheduler, lambda trials, density mixing, or
repeated candidate loop.  Publication remains a separate final call so an
optional refinement does not write intermediate checkpoints.

**Step 4: Route both algorithms through the terminal routine**

The divided route composes `H[rho_divided]` and calls the routine once.  The
reference continuation calls it only after its accepted final refresh.  The
default divided path must not call the old version-2 occupied-checkpoint
publisher.  It must not rebuild or relocalize the construction WFs.  Without
an RT checkpoint request, finalization returns the LCFO eigensystem and density
without invoking the complete-system certified-RT Wannier90 stage.

`prepare_rt_checkpoint` is true only for the final requested LCFO epoch when
the user requested checkpoint output.  It is false for diagnostic solves and
every nonterminal refinement epoch.

**Step 5: Run GREEN**

Run:

```text
python3 tests/dg/run_dg_hybrid_lcfo_finalization_mpi.py
python3 tests/dg/run_dg_hybrid_occupation_policy_mpi.py
python3 tests/dg/run_dg_hybrid_low_energy_symmetry_mpi.py
python3 tests/dg/run_dg_hybrid_certified_rt_basis_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
cmake --build build-hybrid-commit -j2
```

Expected: all PASS; zero-refinement divided fixtures report exactly one
complete eigensolve and version 3, and the build completes.

**Step 6: Commit**

```text
git add src/gs/dc/dg_hybrid_lcfo_finalization.f90 tests/dg/test_dg_hybrid_lcfo_finalization_mpi.f90 tests/dg/run_dg_hybrid_lcfo_finalization_mpi.py
git add -p src/gs/dc/CMakeLists.txt src/gs/main_dft.f90 tests/dg/check_dg_hybrid_divided_lcfo_route.py tests/dg/check_dg_hybrid_continuation_route.py
git diff --cached --check
git diff --cached
git commit -m "feat(dg): finalize one-shot LCFO into certified RT space"
```

### Task 10: Add explicit user-controlled LCFO refinement epochs

**Files:**

- Create: `src/gs/dc/dg_hybrid_lcfo_refinement.f90`
- Create: `tests/dg/test_dg_hybrid_lcfo_refinement_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_lcfo_refinement_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/io/salmon_global.f90:480-490`
- Modify: `src/io/inputoutput.f90:630-645,1155-1170,1885-1905,2955-2975,3130-3170`
- Modify: `src/gs/main_dft.f90:3613-3670`
- Modify: `src/rt/dg/rt_dg_hybrid_checkpoint.f90`
- Modify: `src/rt/dg/rt_dg_hybrid_initialization.f90:640-675`
- Modify: `src/rt/dg/rt_dg_hybrid_stationarity.f90`
- Modify: `src/rt/main_tddft.f90:330-405`
- Modify: `tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_initialization_mpi.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_stationarity_mpi.f90`
- Modify: `tests/dg/check_dg_hybrid_localization_first_inputs.py`
- Modify: `tests/dg/check_dg_hybrid_divided_lcfo_route.py`

**Step 1: Write RED controller and input tests**

Require `dg_hybrid_lcfo_refine_steps=0` by default and broadcast it on all
ranks.  Reject negative values, positive values unless the divided route is
selected, and simultaneous divided/continuation/full-Hybrid route flags.  A
callback-counting fixture must require:

- `0` refinements: one complete solve and one `Delta rho` measurement;
- `2` refinements: three complete solves, three `Delta rho` measurements, two
  density mixes, and two potential updates;
- fixed metric/kinetic/nonlocal/SIPG fingerprints at every epoch;
- fixed fragment-construction-WF and generalized-complement fingerprints at
  every epoch, with zero additional construction Wannier90 calls;
- with RT output requested, zero complete-system RT-localizer calls on
  nonterminal epochs and exactly one on the final epoch, independent of the
  refinement count;
- a new local/total operator fingerprint after each density update;
- no intermediate publication; and
- final coefficients and eigenpair receipt from the final operator epoch.

In the same RED step, extend the three RT fixtures to distinguish the stored
potential-generating density from the final orbital density.  Require
checkpoint authentication to bind both arrays to their epochs, RT startup to
validate the stored eigenpair before rebuilding the TD Hamiltonian from the
orbital density, and stationarity to enforce strict drift only for reference
continuation while measuring finite drift for divided one-shot/refined modes.
Require non-finite drift and broken density/operator fingerprints to fail in
every mode.

**Step 2: Run RED**

Run:

```text
python3 tests/dg/run_dg_hybrid_lcfo_refinement_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/run_rt_dg_hybrid_initialization_mpi.py
python3 tests/dg/run_rt_dg_hybrid_stationarity_mpi.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
```

Expected: compile/source-contract failure because the input and controller are
absent.

**Step 3: Implement the explicit state machine**

```fortran
do epoch=0,refine_steps
  final_epoch=(epoch==refine_steps)
  call solve_once(epoch,prepare_rt_checkpoint=&
    final_epoch.and.rt_checkpoint_requested)
  call measure_shared_dc_density_difference(rho_potential,rho_lcfo,delta_rho)
  if(final_epoch)exit
  call mix_density(rho_potential,rho_lcfo,rho_next)
  call update_potential(rho_next)
  call assemble_local_and_complete_rows(rho_next,epoch=epoch+1)
enddo
call publish_final_epoch_only()
```

Every refinement reuses the already frozen Hybrid basis.  It is a small,
explicit number of complete LCFO density corrections; it does not rerun the
fragment density SCF, expand fragment state caches, or invoke either
construction or certified-RT Wannier90 between epochs.

Record route code, requested/completed refinement count, final `Delta rho`,
potential-density fingerprint, and final operator epoch in the existing
version-3 continuation receipt/fingerprint fields.  Use the two existing v3
density arrays without changing the serialized layout:

- `payload%density`: the density that generated the final stored operator;
- `payload%rt_space%density`: the density reconstructed from the final LCFO
  state and used to initialize RT.

Update authentication so the reconstructed certified density must equal only
`rt_space%density`; authenticate the potential density through the overall
payload fingerprint and epoch receipt.  Existing continuation v3 files, where
the two densities are equal, remain valid.  Do not silently add a refinement
based on `Delta rho`.

At RT startup, first validate the stored eigenpair against the stored operator
and its potential-density epoch.  Then reconstruct the orbital density, verify
`rt_space%density`, update the TD Hamiltonian from that physical density, and
emit the initial density/Hamiltonian change.  Split zero-field stationarity
into measurement and enforcement:

- continuation-reference checkpoints retain the existing strict tolerance
  gate;
- divided one-shot/refined checkpoints require finite authenticated receipts
  and report drift without a magnitude gate; and
- non-finite drift or any broken fingerprint/epoch remains fatal in every
  mode.

Do not implement measurement mode by passing artificially large tolerances.

**Step 4: Run GREEN**

Run:

```text
python3 tests/dg/run_dg_hybrid_lcfo_refinement_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_finalization_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/run_rt_dg_hybrid_initialization_mpi.py
python3 tests/dg/run_rt_dg_hybrid_stationarity_mpi.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
cmake --build build-hybrid-commit -j2
```

Expected: all PASS and the default route still counts one solve.

**Step 5: Commit**

```text
git add src/gs/dc/dg_hybrid_lcfo_refinement.f90 tests/dg/test_dg_hybrid_lcfo_refinement_mpi.f90 tests/dg/run_dg_hybrid_lcfo_refinement_mpi.py
git add -p src/gs/dc/CMakeLists.txt src/io/salmon_global.f90 src/io/inputoutput.f90 src/gs/main_dft.f90 src/rt/dg/rt_dg_hybrid_checkpoint.f90 src/rt/dg/rt_dg_hybrid_initialization.f90 src/rt/dg/rt_dg_hybrid_stationarity.f90 src/rt/main_tddft.f90 tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90 tests/dg/test_rt_dg_hybrid_initialization_mpi.f90 tests/dg/test_rt_dg_hybrid_stationarity_mpi.f90 tests/dg/check_dg_hybrid_localization_first_inputs.py tests/dg/check_dg_hybrid_divided_lcfo_route.py
git diff --cached --check
git diff --cached
git commit -m "feat(dg): add explicit LCFO refinement epochs"
```

### Task 11: Validate Si64 one-shot GS and separate-directory zero-field RT

**Files:**

- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in`
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_dg_zero_field_rt.in`
- Modify existing untracked file: `tests/dg/run_dg_hybrid_si64_divided_lcfo.py`
- Modify: `tests/dg/run_dg_hybrid_si64_continuation_rt.py` only to share parser helpers when duplication would otherwise occur

**Step 1: Extend parser-only acceptance RED**

The runner must require:

- a reused conventional DC seed and no `DC #SCF =` lines in the reused GS;
- the same eight ranks and exact rank--fragment mapping;
- one unconstrained construction-WF localization per fragment communicator and
  basis generation, no complete-cell construction localization, and no
  localization call inside divided SCF;
- a dynamic, non-hard-coded retained rank and fragment state counts, including
  logged whole-shell extension events;
- `dg_hybrid_fragment_cg_steps=3` and no fragment iteration receipt above that
  cap;
- divided-SCF convergence, common chemical potential, and per-iteration
  electron receipts;
- complete SIPG/nonlocal final operator receipts;
- exactly one complete eigensolve when refinement count is zero;
- requested/certified/proof energy-window receipts;
- physical occupied/window/density symmetry receipts;
- a complete-v3 checkpoint and matching GS/RT fingerprints;
- `Delta rho` recorded as a diagnostic; and
- finite zero-field stationarity measurements for every requested RT step,
  without applying the continuation-reference magnitude gate to divided mode.

Add negative parser fixtures for a hidden second solve, stale density/operator
epoch, construction localization on the total communicator, a repeated
fragment-localization call during SCF, a fixed 384-state assumption,
basis-closure used as an acceptance gate, rank/mapping mismatch, and an old
version-2 occupied checkpoint.

**Step 2: Run parser RED**

Run: `python3 tests/dg/run_dg_hybrid_si64_divided_lcfo.py --parser-only`

Expected: FAIL until the new receipts are parsed and required.

**Step 3: Implement the runner and fresh-directory policy**

Reuse the already generated final canonical DC seed; do not repeat
conventional DC.  Create separate fresh directories for one-shot GS and RT.
Hash the seed before and after each run and preserve all logs.  Never overwrite
the current continuation-oracle directory.

Treat identical total MPI rank count and identical rank--fragment mapping as a
permanent production compatibility requirement, not a temporary validation
restriction.  A mismatch must stop before consuming any saved `rwf`, density,
potential, spectrum, or occupation payload.

Record the one-shot zero-field drift even when it exceeds the strict
continuation-reference tolerance.  Run a second fresh GS with a small explicit
refinement count and quantify the improvement.  Non-finite drift remains a
failure.  Do not turn the observation into automatic fallback logic.

**Step 4: Run focused GREEN before Si64**

Run:

```text
python3 tests/dg/run_dc_scf_convergence_mpi.py
python3 tests/dg/run_dc_fragment_occupation_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_wannier_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_wannier_lcfo_mpi.py
python3 tests/dg/run_dg_hybrid_wannier_complement_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_subspace_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/run_dg_hybrid_divided_operator_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_finalization_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_refinement_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/run_rt_dg_hybrid_initialization_mpi.py
python3 tests/dg/run_rt_dg_hybrid_stationarity_mpi.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
cmake --build build-hybrid-commit -j2
```

Expected: all PASS.

**Step 5: Run fresh eight-rank one-shot and RT**

Run only after every focused test and the release-style build above pass:

```text
OMP_NUM_THREADS=1 python3 tests/dg/run_dg_hybrid_si64_divided_lcfo.py --mpi-ranks 8 --binary build-hybrid-commit/salmon --result-dir verification-si64-divided-one-shot-lcfo-20260903
```

Expected: unchanged reused DC seed, one construction Wannier90 operation for
each of the eight fragments, converged bounded-update divided SCF, one final
LCFO solve, complete version-3 checkpoint, and a separate zero-field RT
result.  Preserve failure evidence if any gate fails.

**Step 6: Compare with the reference oracle**

Record energy, band gap, occupied-projector difference, density difference,
certified rank, per-fragment construction localization count/time, fragment
state-count history, local-step residuals, symmetry defects, wall time, and
peak RSS.  Compare construction time with the preserved 384-WF complete-cell
Wannier90 reference log, but do not encode 384 into production acceptance.
The comparison is validation evidence, not a post-LCFO density convergence
gate.

**Step 7: Commit the runner contract**

```text
git add -N tests/dg/run_dg_hybrid_si64_divided_lcfo.py
git add -p tests/dg/run_dg_hybrid_si64_divided_lcfo.py tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_dg_zero_field_rt.in tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_divided_lcfo.in tests/dg/run_dg_hybrid_si64_continuation_rt.py
git diff --cached --check
git diff --cached
git commit -m "test(dg): certify divided one-shot LCFO and RT"
```

### Task 12: Run protected regressions, review, and finish the branch

**Files:**

- Review all files changed in Tasks 1-11
- Preserve all untracked verification output

**Step 1: Run all focused and protected verification fresh**

Run:

```text
git diff --check
git diff --check 0bc6f9a38a97dc6b8c204079440230f180d5367d..HEAD
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
python3 tests/dg/run_dc_scf_convergence_mpi.py
python3 tests/dg/run_dc_fragment_occupation_mpi.py
python3 tests/dg/run_dg_dc_seed_checkpoint_mpi.py
python3 tests/dg/run_dg_dc_seed_state_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_wannier_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_wannier_lcfo_mpi.py
python3 tests/dg/run_dg_hybrid_wannier_complement_mpi.py
python3 tests/dg/run_dg_hybrid_projected_fragment_pipeline_mpi.py
python3 tests/dg/run_dg_hybrid_production_fragment_basis_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_solver_mpi.py
python3 tests/dg/run_dg_hybrid_block_cg_mpi.py
python3 tests/dg/run_dg_hybrid_fragment_subspace_mpi.py
python3 tests/dg/run_dg_hybrid_divided_scf_mpi.py
python3 tests/dg/run_dg_hybrid_divided_operator_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_support_redistribution_mpi.py
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/run_dg_hybrid_occupation_policy_mpi.py
python3 tests/dg/run_dg_hybrid_low_energy_symmetry_mpi.py
python3 tests/dg/run_dg_hybrid_certified_rt_basis_mpi.py
python3 tests/dg/run_dg_hybrid_continuation_controller_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_finalization_mpi.py
python3 tests/dg/run_dg_hybrid_lcfo_refinement_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/run_rt_dg_hybrid_initialization_mpi.py
python3 tests/dg/run_rt_dg_hybrid_stationarity_mpi.py
python3 tests/dg/run_rt_dg_hybrid_length_gauge_mpi.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/check_dg_hybrid_divided_lcfo_route.py
python3 tests/dg/check_dg_hybrid_fragment_wannier_route.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/run_dg_hybrid_si64_divided_lcfo.py --parser-only
cmake --build build-hybrid-commit -j2
ctest --test-dir build-hybrid-commit --output-on-failure
```

Expected: all standalone checks PASS, the build completes, and CTest has no
failures.

**Step 2: Request code review**

Invoke `@superpowers:requesting-code-review`.  Review especially:

- exact conventional-DC convergence equivalence;
- common chemical potential and electron count;
- one unconstrained construction-WF call per fragment/basis generation and no
  total-communicator construction localization;
- preservation of every local WF column and gauge invariance of the
  generalized fragment-union PW projection;
- dynamic occupied-plus-guard state extension with no fixed 384-state path;
- SIPG/nonlocal terms exactly once;
- bounded warm-started `[X,R,P]` updates, safe rollback, and no dense or
  complete solve inside fragment SCF;
- no preliminary complete LCFO before divided density convergence;
- zero-refinement solve count exactly one;
- physical-space rather than individual-WF symmetry acceptance;
- complete-system certified RT localization only when checkpoint output needs
  it;
- version-3 density/operator epoch provenance; and
- exact rank-count/rank--fragment DC-seed reuse.

**Step 3: Resolve findings rigorously**

Use `@superpowers:receiving-code-review`.  For each Critical or Important
finding, reproduce it with a RED test, apply the minimal fix, and rerun every
affected check.  Commit review fixes separately.

**Step 4: Verify before any completion claim**

Invoke `@superpowers:verification-before-completion` and rerun affected commands
fresh.  A preserved prior log is evidence for comparison, not a substitute for
fresh verification of changed code.

**Step 5: Finish without automatic merge**

Invoke `@superpowers:finishing-a-development-branch`.  Present the integration
options and do not merge, delete the worktree, clean output, or push unless the
user explicitly requests it.
