# SAWF Symmetry-Closed LCFO Basis Implementation Plan

> **Historical/removed:** This document describes an obsolete experimental DG route
> removed on 2026-07-31. It is retained only as an implementation record and is
> not executable guidance.

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build Wannier90 SAWF seeds from a dynamically sized LCFO basis that is closed under the validated crystal symmetry group while preserving the normal DC-SCF and DC-LCFO outputs.

**Architecture:** Add a tested basis-closure module that orthonormalizes complete symmetry orbits, then add fragment-owner communication that maps source bases onto target fragments. Refactor the Flux-LCFO matrix assembly to accept a dynamic basis descriptor and run it a second time only for the SAWF seed path. All Wannier seed artifacts are generated from that second eigensystem with separate provenance.

**Tech Stack:** Fortran 2003/2008, MPI, BLAS/LAPACK, EigenExa, existing DC-LCFO Flux assembly, Wannier90 SAWF `.dmn`, Python focused test drivers, CMake.

---

### Task 1: Implement rank-revealing symmetry-orbit closure

**Files:**
- Create: `src/gs/dc/lcfo_wannier_sawf_basis.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/check_sawf_closed_basis.py`

**Step 1: Write the failing synthetic tests**

Compile a small Fortran driver that supplies weighted real-space columns and
their symmetry images. Require:

- identity leaves an orthonormal basis unchanged;
- inversion of a nonsymmetric two-column basis adds the missing image columns;
- the returned basis is orthonormal in the supplied `hvol` metric;
- applying inversion to the result has closure leakage below `1d-12`;
- linearly dependent orbit images are removed deterministically;
- insufficient capacity fails without returning a partial basis;
- non-finite candidates fail before LAPACK.

**Step 2: Run the test and verify RED**

```bash
python3 tests/dg/check_sawf_closed_basis.py
```

Expected: FAIL because `lcfo_wannier_sawf_basis.f90` is absent.

**Step 3: Implement the minimal pure numerical module**

Provide a public routine with explicit dimensions:

```fortran
subroutine close_sawf_candidate_basis(candidate,npoint,ncandidate,hvol, &
    rank_tolerance,max_basis,basis,nbasis,singular_values,ok,message)
```

Form the weighted Gram matrix, diagonalize it with `eigen_dsyev`, retain
eigenvalues above `rank_tolerance**2`, construct normalized columns in
descending-eigenvalue order, and run a final modified Gram-Schmidt cleanup.
Allocate outputs transactionally and publish them only on success.

**Step 4: Verify GREEN and build**

```bash
python3 tests/dg/check_sawf_closed_basis.py
cmake --build build-mpi-eigenexa-wannier-lib -j 2
```

Expected: PASS and successful build.

**Step 5: Commit**

```bash
git add src/gs/dc/lcfo_wannier_sawf_basis.f90 src/gs/dc/CMakeLists.txt tests/dg/check_sawf_closed_basis.py
git commit -m "feat: close SAWF candidate basis by symmetry orbit"
```

### Task 2: Build fragment-local candidate orbits with MPI ownership

**Files:**
- Modify: `src/gs/dc/lcfo_wannier_sawf_basis.f90`
- Modify: `src/gs/dc/lcfo_flux.f90`
- Create: `tests/dg/check_sawf_closed_basis_mpi.py`

**Step 1: Write failing fragment-map tests**

Use a two-rank driver with one fragment owner per rank. Cover identity,
self-mapped inversion, pair-exchanging inversion, and a four-operation group.
Check that each target receives every and only mapped source block, independent
of rank numbering. Inject a NaN on rank 1 and require its rank-local diagnostic
before collective failure.

**Step 2: Verify RED**

```bash
python3 tests/dg/check_sawf_closed_basis_mpi.py
```

Expected: FAIL because orbit exchange is absent.

**Step 3: Implement mapped block exchange**

Reuse normalized operations and `validate_sawf_fragment_symmetry_map`. For each
operation, determine the unique source/target owner, transform core-grid
columns with `map_sawf_periodic_grid_point`, exchange only that block, and add
it to the target candidate set. Sort by `(operation,source_fragment,column)`.

**Step 4: Close and validate each target basis**

Call Task 1's numerical routine, then use
`diagnose_sawf_fragment_basis_closure` for every operation. Reject dimensions,
MPI counts, capacity overflow, non-finite data, or closure above tolerance.

**Step 5: Verify and commit**

```bash
python3 tests/dg/check_sawf_closed_basis_mpi.py
python3 tests/dg/check_sawf_fragment_symmetry_map.py
cmake --build build-mpi-eigenexa-wannier-lib -j 2
git add src/gs/dc/lcfo_wannier_sawf_basis.f90 src/gs/dc/lcfo_flux.f90 tests/dg/check_sawf_closed_basis_mpi.py
git commit -m "feat: exchange SAWF fragment symmetry orbits"
```

### Task 3: Introduce a dynamic Flux-LCFO basis descriptor

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Create: `tests/dg/check_sawf_dynamic_flux_basis.py`

**Step 1: Write failing source and numerical regression tests**

Require a basis descriptor containing per-fragment dimensions, capacity,
core/buffer values, and provenance. A small identity-only case must reproduce
the legacy Hamiltonian and eigenvalues to `1d-12`. A synthetic augmented case
must accept `nbasis > nstate_frag` without truncation.

**Step 2: Verify RED**

```bash
python3 tests/dg/check_sawf_dynamic_flux_basis.py
```

**Step 3: Refactor matrix assembly without changing formulas**

Replace fixed `dc%nstate_frag` loop bounds in the basis/Hamiltonian/halo matrix
path with descriptor capacity and actual `n_basis`. Keep normal-path allocation
at the legacy capacity. Add dynamic integer-overflow and allocation checks.

**Step 4: Reapply the physical Hamiltonian**

For an augmented descriptor, call the existing `hpsi` and unchanged Flux
volume/surface/nonlocal assembly on every retained function. Do not construct a
group-averaged Hamiltonian and do not polar-repair symmetry matrices.

**Step 5: Verify and commit**

```bash
python3 tests/dg/check_sawf_dynamic_flux_basis.py
python3 tests/dg/check_sawf_closed_basis.py
cmake --build build-mpi-eigenexa-wannier-lib -j 2
git add src/gs/dc/lcfo_flux.f90 tests/dg/check_sawf_dynamic_flux_basis.py
git commit -m "refactor: support dynamic Flux-LCFO seed bases"
```

### Task 4: Add the isolated SAWF seed eigensystem

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/io/salmon_global.f90`
- Modify: `src/io/inputoutput.f90`
- Modify: `tests/dg/check_sawf_input_and_build.py`
- Create: `tests/dg/check_sawf_seed_provenance.py`

**Step 1: Write failing isolation tests**

Require `wannier_site_symmetry='off'` to skip every closure allocation. In SAWF
mode, require normal outputs to be completed before the second seed pass and
require distinct seed filenames/provenance. Stale normal or SAWF seed files
must be rejected.

**Step 2: Add bounded namelist controls**

Add positive controls for the SAWF basis rank threshold and hard maximum basis
per fragment. Defaults must preserve off-mode behavior and fail clearly when
the symmetry orbit does not fit.

**Step 3: Run the second eigensystem only in SAWF mode**

After normal LCFO publication, construct the closed descriptor, rebuild the
Flux Hamiltonian, run EigenExa, and publish `wavefunctions_sawf_seed.bin` plus
matching basis/grid/Hamiltonian metadata. Restore or deallocate all seed-only
state before returning.

**Step 4: Verify and commit**

```bash
python3 tests/dg/check_sawf_input_and_build.py
python3 tests/dg/check_sawf_seed_provenance.py
python3 tests/dg/check_sawf_dynamic_flux_basis.py
cmake --build build-mpi-eigenexa-wannier-lib -j 2
git add src/io/salmon_global.f90 src/io/inputoutput.f90 src/gs/dc/lcfo_flux.f90 tests/dg/check_sawf_input_and_build.py tests/dg/check_sawf_seed_provenance.py
git commit -m "feat: build isolated symmetry-closed SAWF seed"
```

### Task 5: Bind all Wannier artifacts to the closed seed

**Files:**
- Modify: `src/gs/dc/lcfo_flux.f90`
- Modify: `src/gs/dc/lcfo_wannier_sawf_band.f90`
- Modify: `tests/dg/check_sawf_band_representation.py`
- Modify: `tests/dg/check_sawf_integration_log.py`

**Step 1: Write failing mixed-provenance tests**

Attempt to combine normal `.eig`/coefficients with closed-basis AMN/MMN/DMN
and require a fatal provenance mismatch. Require all artifacts to report the
same seed token, dimensions, basis rank threshold, and operation digest.

**Step 2: Route seed generation consistently**

Make `.eig`, `.amn`, `.mmn`, `D_band`, `.dmn`, exported real-space Wannier
basis, and DG import read only the closed seed when SAWF is enabled.

**Step 3: Verify representation gates**

Require basis closure, `D_band` unitarity, group products, Hamiltonian
covariance, and AMN covariance before activating `site_symmetry` in `.win`.

**Step 4: Verify and commit**

```bash
python3 tests/dg/check_sawf_band_representation.py
python3 tests/dg/check_sawf_integration_log.py --help
cmake --build build-mpi-eigenexa-wannier-lib -j 2
git add src/gs/dc/lcfo_flux.f90 src/gs/dc/lcfo_wannier_sawf_band.f90 tests/dg/check_sawf_band_representation.py tests/dg/check_sawf_integration_log.py
git commit -m "feat: bind Wannier artifacts to closed SAWF seed"
```

### Task 6: Run C64 correctness gates and clean the branch

**Files:**
- Modify: `samples/exercise_dg_fragment_rt/inputfile_gs_w90_pseudo_sawf_aligned_nw576_nb664`
- Modify: `docs/plans/2026-07-10-sawf-dmn.md`

**Step 1: Run focused and full SAWF tests**

```bash
python3 tests/dg/check_sawf_closed_basis.py
python3 tests/dg/check_sawf_closed_basis_mpi.py
python3 tests/dg/check_sawf_dynamic_flux_basis.py
python3 tests/dg/check_sawf_seed_provenance.py
python3 tests/dg/check_sawf_band_representation.py
python3 tests/dg/check_sawf_win_integration.py
cmake --build build-mpi-eigenexa-wannier-lib -j 2
```

**Step 2: Run aligned C64 with MPI 8 and OMP 2**

Use a positive, converged `energy_cut`; record the normal and closed basis
dimensions, closure leakage, EigenExa dimensions, representation residuals,
and Wannier90 result. Do not weaken the closure tolerance to pass the gate.

**Step 3: Run physics gates**

Run field-off propagation first. Only after stationarity passes, run the agreed
long-pulse `dt=2` calculation and compare `Pz` and HHG odd/even structure. Ignore
the known-unreliable DG `Jz` route.

**Step 4: Review and cleanup**

Review the complete diff for fixed-dimension assumptions, mixed provenance,
rank-local error loss, stale experimental environment toggles, and obsolete
diagnostics. Remove only routes superseded by the verified closed-seed path.

**Step 5: Commit the verified integration**

```bash
git add samples/exercise_dg_fragment_rt/inputfile_gs_w90_pseudo_sawf_aligned_nw576_nb664 docs/plans/2026-07-10-sawf-dmn.md
git commit -m "test: validate symmetry-closed C64 SAWF route"
```
