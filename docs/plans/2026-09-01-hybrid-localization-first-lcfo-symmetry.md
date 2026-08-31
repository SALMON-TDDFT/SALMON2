# Hybrid Localization-First LCFO Symmetry Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build Hybrid DG functions for locality first, certify the smallest symmetry-closed full-system LCFO energy space requested by the user, propagate only in a localized unitary gauge of that certified space, and reuse an exactly compatible conventional DC result.

**Architecture:** Keep the first unconstrained WF+PW basis as a construction space.  Solve its complete LCFO spectrum once, reapply the authoritative SALMON occupation policy, extend the requested HOMO-relative energy window to the first degeneracy-complete cluster that passes both subspace-closure and energy-covariance gates, and localize only inside that certified span for RT.  A separate rank-sharded DC seed format may skip the conventional DC solve only when MPI rank count, rank-to-fragment mapping, array bounds, and ownership maps are identical.

**Tech Stack:** Fortran 2008, MPI, ScaLAPACK `PZHEEVD`, Wannier90 library interface, Python contract/MPI runners, CMake.

---

## Execution constraints

- Work only in `/Users/otobetoshihito/SALMON-dev/SALMON2_RTDG/.worktrees/wpw-s-orthogonal-complement`.
- Do not create another worktree.
- Preserve every pre-existing dirty change and every generated verification log.
- Use `git add -p` for pre-existing dirty files and inspect `git diff --cached --name-status` plus `git diff --cached --check` before every commit.
- Never stage `build-hybrid-commit/`, replay artifacts, or `verification-*` directories.
- Invoke `@superpowers:test-driven-development` before implementation tasks.
- **Stop immediately after Task 1 is committed and report that checkpoint before starting Task 2.**

### Task 1: Define the user input contracts

**Files:**

- Create: `tests/dg/check_dg_hybrid_localization_first_inputs.py`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_dg_continuation.in`

**Step 1: Write the failing contract test**

Create a source-level input contract test that requires all three new globals,
their `&dc` namelist membership, defaults, broadcasts, unit handling, logs, and
validation.  It must require these declarations:

```fortran
real(8)        :: dg_hybrid_symmetry_energy_window
character(16)  :: dg_dc_seed_mode
character(256) :: dg_dc_seed_directory
```

The checker must also require the Si64 continuation input to contain an
explicit nonnegative window, initially:

```text
dg_hybrid_symmetry_energy_window=0.2d0
```

Require source patterns proving that the four cutoffs remain distinct:
`energy_cut`, `lambda_cut`, `wannier_pw_cutoff`, and
`dg_hybrid_symmetry_energy_window`.

**Step 2: Run RED**

Run:

```text
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
```

Expected: FAIL because the new input variables and contracts do not exist.

**Step 3: Implement the minimal input plumbing**

Add the three variables to `salmon_global`, the `&dc` namelist, defaults,
global broadcasts, and `variables.log` output.  Use exactly:

```fortran
dg_hybrid_symmetry_energy_window = -1d0
dg_dc_seed_mode = 'off'
dg_dc_seed_directory = ''
```

Validate the raw energy-window value before conversion:

- exactly `-1d0` is the only legacy-rank sentinel;
- a finite value greater than or equal to zero selects energy-window mode;
- every other negative value and every nonfinite value is fatal.

Preserve the sentinel and convert only nonnegative values:

```fortran
if(dg_hybrid_symmetry_energy_window>=0d0) &
  dg_hybrid_symmetry_energy_window = &
    dg_hybrid_symmetry_energy_window*uenergy_to_au
```

Accept only trimmed seed modes `off`, `write`, `read`, and `auto`.  Require a
nonblank `dg_dc_seed_directory` for every mode except `off`.  On the Hybrid
continuation route independently require finite `energy_cut`, positive finite
`lambda_cut`, and positive finite `wannier_pw_cutoff`; do not use one setting
as a fallback for another.

**Step 4: Run GREEN and protected input checks**

Run:

```text
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
python3 tests/dg/check_dg_hybrid_divided_dc_controls.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
python3 tests/dg/check_sawf_input_and_build.py
git diff --check
```

Expected: PASS.

**Step 5: Commit only Task 1 hunks**

Stage the new checker and fixture line normally.  Stage only Task 1 hunks from
the two already-dirty Fortran files.  Commit as:

```text
feat(dg): define localization-first input contracts
```

**Step 6: Stop at the first checkpoint**

Report the exact tests, commit, and remaining dirty status.  Do not begin Task
2 until the checkpoint has been reviewed.

### Task 2: Implement the strict rank-sharded DC seed format

**Files:**

- Create: `src/gs/dc/dg_dc_seed_checkpoint.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_dc_seed_checkpoint_mpi.f90`
- Create: `tests/dg/run_dg_dc_seed_checkpoint_mpi.py`

**Step 1: Write RED MPI round-trip and rejection tests**

Define tests on 1, 2, 4, and 8 ranks for valid round trip, missing manifest,
missing/truncated/corrupt shard, interrupted publication, rank-count mismatch,
rank-to-fragment mismatch, local-bound mismatch, and ownership-map mismatch.
Also require downstream-only control changes to leave compatibility unchanged.

Run:

```text
python3 tests/dg/run_dg_dc_seed_checkpoint_mpi.py
```

Expected: compile failure because the module is absent.

**Step 2: Add the versioned manifest and payload API**

Expose a small API:

```fortran
type,public :: s_dg_dc_seed_contract
  integer :: version=1,mpi_size=0,rank=0,fragment_id=0
  integer :: rwf_bounds(14)=0,rho_bounds(6)=0,vloc_bounds(6)=0
  integer(8) :: immutable_fingerprint=0_8,ownership_fingerprint=0_8
end type

type,public :: s_dg_dc_seed_payload
  real(8),allocatable :: rwf(:,:,:,:,:,:,:)
  real(8),allocatable :: rho_tot(:,:,:),vloc_tot(:,:,:)
  real(8),allocatable :: esp(:,:,:),rocc(:,:,:)
  real(8) :: mu=0d0,residual=huge(0d0)
  integer :: iteration=0
end type

public :: write_dg_dc_seed,read_dg_dc_seed,probe_dg_dc_seed
```

The immutable fingerprint must cover cell/grid, atoms/species, fragment and
buffer topology, electron/state/spin/Gamma-real configuration, occupation/XC,
stencil and pseudopotential content, MPI size, rank-to-fragment mapping, exact
array bounds, and exact global-to-local ownership maps.  Exclude localization,
Wannier90, PW, LCFO, and symmetry-window controls.

Write one temporary shard per rank, validate its digest collectively, rename
all shards to their final names, and atomically publish the manifest last.  A
reader accepts only the complete committed set with matching publication ID,
ordered global digest, shapes, bounds, finite values, electron count, and a
stored residual below the current convergence threshold.  Do not implement
redistribution.

**Step 3: Run GREEN and commit**

Run:

```text
python3 tests/dg/run_dg_dc_seed_checkpoint_mpi.py
git diff --check
```

Commit as:

```text
feat(dg): add strict distributed DC seed checkpoints
```

### Task 3: Integrate DC seed read/write around conventional SCF

**Files:**

- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dcdft.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/check_dg_dc_seed_route.py`
- Create: `tests/dg/test_dg_dc_seed_state_mpi.f90`
- Create: `tests/dg/run_dg_dc_seed_state_mpi.py`

**Step 1: Write RED route and restored-state tests**

Require seed probing after `initialization2_dft` and before
`scf_iteration_dft`.  A valid `read` or `auto` seed must produce no
conventional `DC #SCF =` line.  `read` with no committed seed must fail;
`auto` with no seed must run conventional SCF then write; a present invalid
seed must fail in both `read` and `auto`.  `probe_dg_dc_seed` must distinguish
an absent seed (no recognized seed artifacts) from an invalid interrupted seed
(shards or temporary/manifest artifacts exist without one complete commit).

Run:

```text
python3 tests/dg/check_dg_dc_seed_route.py
python3 tests/dg/run_dg_dc_seed_state_mpi.py
```

**Step 2: Implement the route state machine**

After normal allocation/initialization, build the immutable contract and
apply:

```text
off   -> conventional DC only
write -> conventional DC, then seed publish
read  -> exact seed load or fatal
auto  -> exact load when committed; otherwise conventional DC then publish
```

On load restore exact `spsi%rwf` bounds and values, owned
`dc%rho_tot_s(1)%f`, owned `dc%vloc_tot(1)%f`, `energy%esp`, `system%rocc`,
`system%mu`, residual, and iteration.  Rebuild fragment density/potential and
derived work arrays with existing DC routines; do not deserialize mixing,
communicators, grids, or Hamiltonian workspaces.  On write publish only after
`sum1 < threshold` and before any overlapping-Wannier/Hybrid dispatch.

Emit one receipt:

```text
[DG-DC-SEED] mode=... publication_id=... scf_skipped=... mpi_size=... mapping_fingerprint=...
```

**Step 3: Run GREEN, build, and commit**

Run:

```text
python3 tests/dg/check_dg_dc_seed_route.py
python3 tests/dg/run_dg_dc_seed_state_mpi.py
python3 tests/dg/run_dg_dc_seed_checkpoint_mpi.py
cmake --build build-hybrid-commit -j2
```

Commit as:

```text
feat(dg): reuse exact compatible conventional DC states
```

### Task 4: Make Wannier90 symmetry mode explicit

**Files:**

- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_w90_mpi.py`

**Step 1: Write RED mode-policy tests**

Require public constants and an explicit argument:

```fortran
integer,parameter,public :: DG_W90_CONSTRAINED=1
integer,parameter,public :: DG_W90_UNCONSTRAINED=2
```

Test that constrained mode requires `.dmn` and writes
`site_symmetry=.true.`, while unconstrained mode removes a stale `.dmn`, never
consumes one, and writes `site_symmetry=.false.`.  Invalid or rank-disagreeing
modes fail collectively.

Run:

```text
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
```

**Step 2: Implement without changing legacy behavior**

Add a mandatory `symmetry_mode` to `setup_dg_w90_gamma_library` and its
artifact policy.  Pass `DG_W90_CONSTRAINED` at every existing call site in
this task.  Do not route continuation to unconstrained mode yet.

**Step 3: Run GREEN, build, and commit**

Run:

```text
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
cmake --build build-hybrid-commit -j2
```

Commit as:

```text
refactor(w90): make symmetry constraint mode explicit
```

### Task 5: Make Wannier90 replay bundles mode-aware

**Files:**

- Modify: `src/gs/dc/dg_overlapping_wannier_w90.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_w90_mpi.f90`
- Modify: `tests/dg/replay_dg_wannier90_bundle.py`
- Modify: `tests/dg/test_replay_dg_wannier90_bundle.py`

**Step 1: Write RED export/replay tests**

Require `export_dg_w90_replay_bundle(...,symmetry_mode,...)`.  Constrained
bundles include and hash `.dmn`; unconstrained bundles neither open, copy, nor
hash `.dmn` and reject a stale target `.dmn`.

Add the required CLI choice:

```text
--symmetry-mode constrained|unconstrained
```

Run:

```text
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
python3 tests/dg/test_replay_dg_wannier90_bundle.py
```

**Step 2: Implement and verify**

Normalize unconstrained replay `.win` files to `site_symmetry=false` and remove
`symmetrize_eps`; keep optimizer overrides.  Record `symmetry_mode` in the
replay receipt.

Run the two commands above again and commit as:

```text
feat(w90): make replay bundles symmetry-mode aware
```

### Task 6: Add the localization-first seed and receipt contract

**Files:**

- Create: `src/gs/dc/dg_hybrid_localization_first.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_hybrid_localization_first_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_localization_first_mpi.py`

**Step 1: Write RED MPI tests**

Test a deliberately nonclosed raw occupied+s+p span.  Require metric
orthonormalization to preserve its weighted projector and full rank without
group averaging, characters, orbit completion, or per-WF pruning.  Large but
finite spreads pass; nonfinite centers/spreads, a nonunitary transform,
nonconvergence, or rank loss fail.  Require rank-invariant receipts on 1, 2, 4,
and 8 ranks.

Run:

```text
python3 tests/dg/run_dg_hybrid_localization_first_mpi.py
```

**Step 2: Implement the narrow module**

Expose:

```fortran
type,public :: s_dg_hybrid_localization_receipt
  logical :: valid=.false.,symmetry_constrained=.false.,converged=.false.
  integer :: raw_rank=0,retained_rank=0,iterations=0
  real(8) :: spread_min=0d0,spread_max=0d0
  real(8) :: spread_mean=0d0,spread_total=0d0
  real(8) :: transform_unitarity_defect=huge(0d0)
  integer(8) :: seed_fingerprint=0_8,transform_fingerprint=0_8
end type
```

`prepare_dg_hybrid_localization_first_seed` may call only the existing metric
rank/orthonormalization path.  `build_dg_hybrid_localization_receipt` validates
the fixed-rank W90 result and never applies a material-independent spread
cutoff.

**Step 3: Run GREEN and commit**

```text
python3 tests/dg/run_dg_hybrid_localization_first_mpi.py
```

Commit as:

```text
feat(dg): add localization-first basis contracts
```

### Task 7: Route Hybrid continuation through unconstrained localization

**Files:**

- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_hybrid_localization_first_mpi.f90`
- Modify: `tests/dg/check_dg_hybrid_continuation_route.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write RED route assertions**

Require continuation to retain crystallographic map construction/integrity,
but bypass fixed-center groups, occupied averaging, spectral basin/orbit
adaptation, translation-character reconstruction, `.dmn` generation, and
pre/post-Wannier affine acceptance gates.  Require
`DG_W90_UNCONSTRAINED`; legacy routes must still be constrained.

**Step 2: Run RED**

```text
python3 tests/dg/check_dg_hybrid_continuation_route.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

**Step 3: Implement the explicit arm**

Branch once after raw occupied+s+p seed assembly.  The localization-first arm
calls Task 6 metric preparation, ordinary W90 localization, and retains every
localized column.  Never pass its projected symmetry action to helpers whose
contract assumes a closed representation.  Keep full crystal maps for later
LCFO and operator diagnostics.

Emit:

```text
[HYBRID-WF-LOCALIZATION] symmetry_constraint=off raw_rank=... retained_rank=... iterations=... spread_min=... spread_max=... spread_mean=... spread_total=... unitarity=...
```

Individual-WF, center, and complete construction-basis symmetry values remain
finite diagnostics only.

**Step 4: Run GREEN, build, and commit**

```text
python3 tests/dg/run_dg_hybrid_localization_first_mpi.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
cmake --build build-hybrid-commit -j2
```

Commit as:

```text
feat(dg): localize Hybrid construction functions unconstrained
```

### Task 8: Preserve authoritative reciprocal PW symmetry and cutoff shells

**Files:**

- Modify: `src/common/dg_hybrid_windowed_pw_types.f90`
- Modify: `src/common/dg_hybrid_windowed_pw_basis.f90`
- Modify: `src/common/dg_hybrid_production_pw_basis.f90`
- Modify: `src/common/dg_hybrid_reciprocal_catalog.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_hybrid_production_pw_basis_mpi.f90`
- Modify: `tests/dg/test_dg_hybrid_reciprocal_catalog_mpi.f90`

**Step 1: Write RED shell/orbit tests**

Use nonidentity rotations that split fragment rows but close reciprocal
vectors.  Require identity-only window bookkeeping and the complete physical
reciprocal catalog simultaneously.  Add cutoff-boundary shell, orbit partner,
missing identity, nonclosed catalog, and post-completion `wannier_pw_max`
capacity cases.

Run:

```text
python3 tests/dg/run_dg_hybrid_reciprocal_catalog_mpi.py
python3 tests/dg/run_dg_hybrid_production_pw_basis_mpi.py
```

**Step 2: Implement exact cutoff completion**

Separate `window_operation_count` from the physical reciprocal
`operation_count`.  Never replace a nonidentity crystal catalog by identity or
a fragment-compatible subgroup.  Select with:

```text
tau_pw = 64*epsilon(1d0)*max(1d0,abs(cutoff),abs(E_g))
```

Include `E_g <= cutoff + tau_pw`, complete the boundary kinetic shell, then
complete every authoritative reciprocal orbit.  `wannier_pw_max` is a
post-completion capacity guard and may not clip the result.

Emit requested/effective cutoff, `shell_added`, and `orbit_added`.

**Step 3: Run GREEN and commit**

```text
python3 tests/dg/run_dg_hybrid_reciprocal_catalog_mpi.py
python3 tests/dg/run_dg_hybrid_production_pw_basis_mpi.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
```

Commit as:

```text
fix(dg): preserve reciprocal symmetry in cutoff-complete PW bases
```

### Task 9: Share the authoritative occupation and HOMO policy

**Files:**

- Create: `src/gs/dc/dg_hybrid_occupation_policy.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Modify: `src/gs/dc/dg_hybrid_ground_state_types.f90`
- Create: `tests/dg/test_dg_hybrid_occupation_policy_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_occupation_policy_mpi.py`
- Modify: `tests/dg/test_dg_hybrid_ground_state_types_mpi.f90`

**Step 1: Write RED policy tests**

Cover zero/finite temperature, Fermi degeneracy, reordered input states,
capacity failure, omitted-tail electron count, occupation threshold, and HOMO
selection.

Run:

```text
python3 tests/dg/run_dg_hybrid_occupation_policy_mpi.py
python3 tests/dg/run_dg_hybrid_ground_state_types_mpi.py
```

**Step 2: Implement the final-spectrum kernel**

Factor the existing SALMON occupation/chemical-potential policy through a
Hybrid-facing wrapper rather than inventing a new Fermi rule.  Apply it to the
ascending complete LCFO spectrum.  Define:

```fortran
noccupied = count(occupations > 64d0*epsilon(1d0))
e_homo = eigenvalues(noccupied)
```

Reject insufficient spectrum capacity and an omitted occupation tail larger
than `dg_dc_gs_electron_count_tolerance`.  Do not copy occupations by old
state identity.

**Step 3: Run GREEN and commit**

Commit after the two commands pass:

```text
feat(dg): derive occupations and HOMO from the final LCFO spectrum
```

### Task 10: Select the first symmetry-closed spectral extension

**Files:**

- Modify: `src/gs/dc/dg_hybrid_generalized_eigensystem.f90`
- Modify: `src/gs/dc/dg_hybrid_low_energy_symmetry.f90`
- Modify: `src/gs/dc/dg_hybrid_ground_state_types.f90`
- Modify: `tests/dg/test_dg_hybrid_generalized_eigensystem_mpi.f90`
- Modify: `tests/dg/test_dg_hybrid_low_energy_symmetry_mpi.f90`

**Step 1: Write RED spectrum and cluster tests**

Cover zero window, a cutoff between levels, a cutoff on a level, degenerate and
numerically split clusters, different spectra selecting different ranks, the
first passing cluster versus a later passing cluster, no proof state, and no
passing cluster.  Verify that occupied-projector or density failure cannot be
repaired by adding empty states.

Run:

```text
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/run_dg_hybrid_low_energy_symmetry_mpi.py
```

**Step 2: Publish one complete solve**

The current ScaLAPACK backend must expose all construction-basis eigenpairs
from one `PZHEEVD` call per continuation candidate.  Remove repeated prefix
solves; geometric growth remains documented only for a future true partial
backend.

**Step 3: Implement adaptive complete-cluster search**

Start with `E_cut=E_HOMO+window`, include every state satisfying
`E_i<=E_cut+tau_deg`, and complete adjacent gaps with:

```text
tau_deg = max(dg_dc_gs_final_orbital_tolerance,64*epsilon(1d0))
          * max(1d0,abs(E_i),abs(E_j),abs(E_cut))
```

At successive complete cluster boundaries evaluate both target-subspace
closure and target-energy covariance.  Select their first simultaneous pass.
Require one state above the final certified cluster as a proof state.  The
basis ceiling or absence of a passing cluster is a capacity failure.

Expose a named result containing requested/certified cutoffs and ranks,
extension energy/states, proof energy, defects, and worst operation.

**Step 4: Run GREEN and commit**

```text
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/run_dg_hybrid_low_energy_symmetry_mpi.py
```

Commit as:

```text
feat(dg): certify the first symmetry-closed LCFO energy space
```

### Task 11: Build a localized certified RT basis

**Files:**

- Create: `src/gs/dc/dg_hybrid_certified_rt_basis.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_hybrid_certified_rt_basis_mpi.f90`
- Create: `tests/dg/run_dg_hybrid_certified_rt_basis_mpi.py`

**Step 1: Write RED projector/gauge/operator tests**

Construct a nonclosed construction space whose certified LCFO prefix is
closed.  Require a second unconstrained localization to converge with finite
before/after spread receipts while leaving the certified projector and
symmetry defects invariant.  Record spread improvement without imposing a
material-independent improvement threshold.  Test nonunitary transforms, rank
loss, scalar covariance, vector/tensor transformation laws, and rejection of
independent element pruning.

Run:

```text
python3 tests/dg/run_dg_hybrid_certified_rt_basis_mpi.py
```

**Step 2: Implement the certified gauge**

Expose a result holding row-owned `C_cert`, `U_rt`, `B_rt`, certified
eigenvalues, initial occupied amplitudes, projected exact operators, spread
receipts, and fingerprints.  Enforce:

```text
B_rt = C_cert U_rt
P_cert = C_cert C_cert^H S = B_rt B_rt^H S
A_occ(0) = U_rt^H(:,1:noccupied)
S_rt = I
H_rt(0) = U_rt^H diag(epsilon_cert) U_rt
```

Use ordinary unconstrained localization with fixed certified rank.  Store all
projected rows; do not drop small matrix elements unless a later feature drops
complete symmetry orbits under its own error bound.

**Step 3: Run GREEN and commit**

Commit after the MPI runner passes:

```text
feat(dg): localize the certified RT symmetry space
```

### Task 12: Integrate physical LCFO acceptance and RT-basis publication

**Files:**

- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_hybrid_continuation_controller.f90`
- Modify: `tests/dg/test_dg_hybrid_continuation_controller_mpi.f90`
- Modify: `tests/dg/check_dg_hybrid_continuation_route.py`
- Modify: `tests/dg/run_dg_hybrid_si64_continuation_rt.py`

**Step 1: Write RED route/controller assertions**

Require one full LCFO solve, Task 9 occupations, unconditional occupied and
density gates, Task 10 first-cluster certification, Task 11 final basis, and no
propagation publication of construction-only directions.

**Step 2: Implement final candidate acceptance**

After each converged continuation candidate, perform the steps in that order.
Do not symmetrize Hamiltonian, coefficients, occupations, or density.  Emit:

```text
[HYBRID-LCFO-WINDOW] delta_e=... homo=... requested_cutoff=... requested_rank=... certified_cutoff=... certified_rank=... extension_states=... proof_energy=...
[HYBRID-LCFO-SYMMETRY] occupied=... target=... energy=... density=... worst_operation=...
[HYBRID-RT-BASIS] certified_rank=... localization_spread=... embedding_fingerprint=... operator_covariance=...
```

The final RT rank must equal the certified rank, not the construction rank.

**Step 3: Run GREEN, build, and commit**

```text
python3 tests/dg/run_dg_hybrid_continuation_controller_mpi.py
python3 tests/dg/run_dg_hybrid_certified_rt_basis_mpi.py
python3 tests/dg/check_dg_hybrid_continuation_route.py
cmake --build build-hybrid-commit -j2
```

Commit as:

```text
feat(dg): hand off only the certified LCFO symmetry space
```

### Task 13: Upgrade the Hybrid GS-to-RT checkpoint to version 3

**Files:**

- Modify: `src/rt/dg/rt_dg_hybrid_checkpoint.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_checkpoint_mpi.f90`
- Modify: `tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py`

**Step 1: Write RED v3 serialization tests**

Require round trip and digest coverage for construction provenance,
row-owned `C_cert`, `U_rt`, final `B_rt`, certified eigenvalues, proof energy,
occupations, initial occupied amplitudes, final RT basis values/operators, and
all named certification receipts.  Tamper with each field.  Require new-route
v2 rejection and rank-redistributed v3 equality.

**Step 2: Implement named v3 payloads**

Do not extend the old positional real receipt.  Add named derived types and
update validation, hashing, write/read, coalescing, and redistribution paths
together.  The new continuation route requires v3; legacy routes remain
confined to their existing formats.  This checkpoint may redistribute ranks;
the strict rank/mapping rule belongs only to the DC seed.

**Step 3: Run GREEN and commit**

```text
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
```

Commit as:

```text
feat(rt): checkpoint certified Hybrid RT bases as version 3
```

### Task 14: Propagate only inside the certified RT basis

**Files:**

- Modify: `src/rt/dg/rt_dg_hybrid_initialization.f90`
- Modify: `src/rt/dg/rt_dg_hybrid_density_update.f90`
- Modify: `src/rt/dg/rt_dg_hybrid_length_gauge.f90`
- Modify: `src/rt/dg/rt_dg_hybrid_stationarity.f90`
- Modify: `src/rt/main_tddft.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_initialization_mpi.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_length_gauge_mpi.f90`
- Modify: `tests/dg/test_rt_dg_hybrid_stationarity_mpi.f90`

**Step 1: Write RED startup and evolution tests**

Require RT startup to recompute `H*C_cert-S*C_cert*epsilon`,
`C_cert^H*S*C_cert-I`, embedding/unitarity, certified target closure and energy
covariance, electron count, reconstructed density, occupied projector,
projected fixed-operator covariance, and every fingerprint.  Require all
coefficient/state/operator extents to equal certified rank.

Test zero-field stationarity, a field preserving a subgroup, and a field that
physically lowers equilibrium symmetry while maintaining the vector/tensor
covariance relation between symmetry-related fields.

**Step 2: Implement certified-space RT**

Initialize occupied amplitudes from stored `U_rt`, reconstruct density with
the final localized RT basis values, project time-dependent local potentials
into that same basis, and evolve only:

```text
i S_rt dA(t)/dt = H_rt(t) A(t)
```

Never allocate or populate construction-only coefficient directions.  Full
construction-basis covariance remains diagnostic only.

**Step 3: Run GREEN and commit**

```text
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/run_rt_dg_hybrid_initialization_mpi.py
python3 tests/dg/run_rt_dg_hybrid_stationarity_mpi.py
python3 tests/dg/run_rt_dg_hybrid_length_gauge_mpi.py
```

Commit as:

```text
feat(rt): evolve Hybrid states in the certified symmetry space
```

### Task 15: Run end-to-end reuse, symmetry, and RT verification

**Files:**

- Modify: `tests/dg/run_dg_hybrid_si64_continuation_rt.py`
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_hybrid_dg_continuation.in`
- Modify only when a fresh failing test proves a defect elsewhere.

**Step 1: Extend the Si64 runner**

Run a cold `auto` calculation that publishes a DC seed, then a second run with
the same ranks/mapping that reuses it and emits no conventional `DC #SCF =`
iterations.  Change localization, PW cutoff, and symmetry window separately
and require the same seed to remain valid.  Change MPI size/mapping and require
clear rejection.  Preserve every output under a new untracked verification
directory.

Require requested/certified/proof receipts, nonidentity reciprocal action,
unconstrained WF localization, a v3 payload, matching GS/RT fingerprints,
certified-rank RT extents, and zero-field stationarity.

**Step 2: Run focused and protected verification**

Run:

```text
python3 tests/dg/check_dg_hybrid_localization_first_inputs.py
python3 tests/dg/run_dg_dc_seed_checkpoint_mpi.py
python3 tests/dg/run_dg_dc_seed_state_mpi.py
python3 tests/dg/run_dg_hybrid_localization_first_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_w90_mpi.py
python3 tests/dg/test_replay_dg_wannier90_bundle.py
python3 tests/dg/run_dg_hybrid_reciprocal_catalog_mpi.py
python3 tests/dg/run_dg_hybrid_production_pw_basis_mpi.py
python3 tests/dg/run_dg_hybrid_occupation_policy_mpi.py
python3 tests/dg/run_dg_hybrid_low_energy_symmetry_mpi.py
python3 tests/dg/run_dg_hybrid_generalized_eigensystem_mpi.py
python3 tests/dg/run_dg_hybrid_certified_rt_basis_mpi.py
python3 tests/dg/run_dg_hybrid_continuation_controller_mpi.py
python3 tests/dg/run_rt_dg_hybrid_checkpoint_mpi.py
python3 tests/dg/run_rt_dg_hybrid_initialization_mpi.py
python3 tests/dg/run_rt_dg_hybrid_stationarity_mpi.py
python3 tests/dg/run_rt_dg_hybrid_length_gauge_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
cmake --build build-hybrid-commit -j2
ctest --test-dir build-hybrid-commit --output-on-failure
```

Then run the fresh eight-rank Si64 route in a new output directory:

```text
python3 tests/dg/run_dg_hybrid_si64_continuation_rt.py build-hybrid-commit/salmon verification-si64-localization-first-lcfo-20260901/si64-run --ranks 8
```

**Step 3: Review and finish**

Invoke `@superpowers:requesting-code-review`.  Resolve all Critical and
Important findings with focused RED/GREEN commits.  Then invoke
`@superpowers:verification-before-completion` and rerun every affected command
fresh before making any completion claim.  Finally invoke
`@superpowers:finishing-a-development-branch`; present integration options and
do not merge automatically.

Commit the final test-contract changes as:

```text
test(dg): protect certified localization-first Hybrid RT
```
