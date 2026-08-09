# Memory-Bounded Global Wannier Symmetry Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace replicated global-LCFO symmetry workspaces with disk-free orbital-by-spatial MPI distribution while preserving exact full-system affine covariance, fragment redistribution, and V3 publication.

**Architecture:** Factor the full affine catalog into translation and point-cogroup data, distribute orbital rows independently of core-grid ownership, and stream one point representative at a time through sparse `MPI_Alltoallv` point exchange.  Store dense metric and representation data in distributed blocks, reduce residual receipts without full image/residual arrays, then redistribute completed Wanniers by periodic center with buffer-only tail replication.

**Tech Stack:** Fortran 2008, MPI, ScaLAPACK, EigenExa, BLAS/LAPACK, spglib, Python contract tests, CMake clean overlays.

---

## Mandatory execution discipline

For every task below:

1. use `@test-driven-development` and record a genuine RED before production changes;
2. run the listed focused verification on MPI 1/2/4/8 where applicable;
3. perform a specification review and a code-quality review;
4. resolve every Critical and Important finding and repeat affected tests;
5. build from a clean `git archive HEAD` overlay plus only the current diff and documented parent prerequisites; and
6. commit only the reviewed task files, preserving unrelated Si64/HHG worktree changes.

### Task 1: Define the two-dimensional ownership and memory contract

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_types.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Modify: `tests/dg/run_dg_overlapping_wannier_construction_mpi.py`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write the failing ownership tests**

Add cases that create independent orbital and spatial ownership maps, require every global orbital and core point to have one owner, and reject invalid process grids or duplicate/missing ownership.  Add allocation accounting assertions forbidding a production workspace shaped as full `Norb x Nglobal_grid` or `Nsym x Norb x Norb`.

**Step 2: Run RED**

Run:

```bash
python3 tests/dg/run_dg_overlapping_wannier_construction_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: FAIL because distributed ownership and memory counters do not exist.

**Step 3: Implement minimal ownership metadata**

Add a distributed layout type containing communicator handles, orbital block bounds, spatial owner metadata, distributed dense-block descriptors, and current/peak workspace bytes.  Provide overflow-checked allocation accounting and collective validation.  Replicate only `O(Norb)` integer metadata and the small symmetry catalog.

**Step 4: Focused verification**

Run the construction fixture on MPI 1/2/4/8 and the route checker.  Require identical canonical ownership maps and decreasing per-rank large-array bytes as ranks increase.

**Step 5: Review and commit**

After both reviews and a clean-first prerequisite overlay build:

```bash
git add src/gs/dc/dg_overlapping_wannier_types.f90 \
  src/gs/dc/dg_overlapping_wannier_construction.f90 \
  tests/dg/test_dg_overlapping_wannier_construction_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_construction_mpi.py \
  tests/dg/check_dg_overlapping_wannier_route.py
git commit -m "feat(dg): distribute global Wannier ownership"
```

### Task 2: Exchange symmetry images without full-basis broadcasts

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Write the failing permutation tests**

Construct operations whose targets cross every spatial rank.  Compare the expected point-permuted orbital blocks with the distributed result and require malformed send counts, missing target IDs, and duplicate target IDs to fail collectively.  Add a source contract forbidding `MPI_Bcast(owner_basis,nstate*nlocal,...)` in the production symmetry path.

**Step 2: Run RED**

Run the construction fixture on 1/2/4/8 ranks.  Expected: FAIL on the missing sparse exchange API and forbidden broadcast.

**Step 3: Implement streamed `MPI_Alltoallv` exchange**

Precompute integer send/receive plans from physical core IDs for one representative.  Pack only locally owned orbital rows and requested point values, exchange them with `MPI_Alltoallv`, validate counts collectively, and unpack into a bounded point tile.  Reuse the buffers for the next tile and operation.

**Step 4: Focused verification**

Require exact rank-count identity for point images, allocation bounds, and collective rejection cases on MPI 1/2/4/8.

**Step 5: Review and commit**

Run both reviews and the clean-first overlay, then commit the two files as `feat(dg): stream distributed symmetry images`.

### Task 3: Stream metric representations and residual receipts

Execute this task in two reviewed commits.  Task 3A removes unused all-operation
storage and tiles operation-local images/residuals.  Task 3B then gives metric,
inverse-square-root, overlap, and representation matrices genuine block-cyclic
ownership from their first allocation; wrappers that first replicate a dense
input do not satisfy Task 3B.  Do not begin Task 4 until both commits pass.

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_metric.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_metric_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Write dense-reference RED tests**

For small deterministic bases, compare distributed block metric, representation, total/boundary/interior residuals, metric-unitarity defects, and cocycle products against the existing dense calculation.  Require operation workspace to return to baseline after each representative and reject non-finite blocks or singular metrics.

**Step 2: Run RED**

Run metric and construction fixtures on 1/2/4/8 ranks.  Expected: FAIL because all representations and full image/residual arrays are still resident.

**Step 3: Implement distributed dense blocks**

For Task 3A, process one point representative and one orbital tile at a time;
update squared residual norms directly, release the image tile before the next
tile, and do not retain the unused occupied-LCFO representations.  Publish an
honest byte receipt for every explicit workspace allocation.

For Task 3B, use ScaLAPACK descriptors for metric and overlap blocks from their
first allocation.  Do not call a wrapper that accepts a replicated dense input.
Validate each streamed representation, update only the cocycle fingerprint and
maximum receipts, and discard the operation block after its required products
are checked.  Represent pure translations through point permutations and the
Gamma phase rule rather than stored dense matrices.

**Step 4: Focused verification**

Run metric, construction, symmetry, solver, and SCF fixtures on MPI 1/2/4/8.  Require dense-reference agreement within existing tolerances and no `Nsym*Norb**2` allocation.

**Step 5: Review and commit**

Resolve all findings, repeat the clean-first overlay build, and commit as `feat(dg): stream distributed affine covariance receipts`.

### Task 4: Distribute localization and center-based fragment redistribution

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_localization.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_localization_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_solver_mpi.f90`

**Step 1: Write failing distributed-localization tests**

Require a distributed orbital transform to reproduce the small dense reference, preserve the occupied projector and occupation vector, close center orbits even when the symmetry center lies outside a fragment, and redistribute each core exactly once with only covered buffer tails duplicated.

**Step 2: Run RED**

Run localization and solver fixtures on MPI 1/2/4/8.  Expected: FAIL because localization still assumes full orbital storage.

**Step 3: Implement bounded localization and redistribution**

Partition transform rows by Wannier index, reduce centers/spreads as `O(Norb)` metadata, and retain orbital values in distributed blocks.  After deterministic periodic-center ownership, use bounded `MPI_Alltoallv` transfers for core blocks and generate only required buffer tails.  Extend the redistribution fingerprint with process-grid-independent canonical ownership and physical IDs.

**Step 4: Focused verification**

Run localization, solver/density, operator, SCF, checkpoint, and coefficient-RT fixtures on MPI 1/2/4/8.  Require byte-identical canonical V3 payloads across rank counts.

**Step 5: Review and commit**

After reviews and clean-first build, commit as `feat(dg): localize and redistribute Wanniers in bounded memory`.

### Task 5: Enforce V3 memory and provenance acceptance

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_checkpoint.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_checkpoint_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write failing publication tests**

Require V3 to record layout-independent orbital/spatial ownership fingerprints, symmetry workspace peak bytes, and redistribution workspace peak bytes.  Reject publication when a forbidden full-size allocation is observed, workspace accounting differs across ranks unexpectedly, or any provenance is zero/noncanonical.

**Step 2: Run RED**

Run checkpoint fixture and route checker.  Expected: FAIL on missing memory receipts.

**Step 3: Implement memory receipts**

Extend manifest serialization, digest, size validation, broadcast, replicated-payload validation, and publication rejection.  Preserve legacy rejection and the existing generalized-eigenvalue Exp RT reader boundary.

**Step 4: Focused verification**

Run checkpoint, metric, construction, solver/density, SCF, operator, and coefficient-RT fixtures on MPI 1/2/4/8 plus both retained-route contracts.

**Step 5: Review and commit**

After reviews and clean-first build, commit as `feat(dg): require bounded-memory V3 provenance`.

### Task 6: Genuine ideal-Si64 memory, GS, LR, and HHG acceptance

**Files:**
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_gs.in`
- Modify: `tests/dg/run_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/check_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/test_si64_harmonic_morphology.py`
- Modify: `docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md`

**Step 1: Record physical and memory RED**

Require ideal undisplaced Si64, 384 retained/128 occupied states, full-system inversion, all V3 GS receipts, no forbidden allocation, and a documented per-rank/module peak below the replicated baseline.  Retain polarization as the primary LR/HHG source, current only as a secondary `dP/dt` check, and explicit H2/H4 peak/dip/slope classification.

**Step 2: Build a fresh acceptance overlay**

Create it from committed HEAD plus the current task diff and parent prerequisites.  Configure Release with MPI, ScaLAPACK, EigenExa, and spglib enabled, Wannier90 disabled.  Build EigenExa with `-j1`, then run `cmake --build <overlay> --clean-first -j4`.

**Step 3: Run genuine physics evidence**

On eight MPI ranks and one OpenMP thread, run normal DC-SCF, LCFO+EigenExa, distributed global Wannier GS, and V3.  Require memory receipts, 384/128 ranks, affine closure, inversion, operator covariance, stationarity, and checkpoint publication.  Then run field-off, impulse LR, and long-pulse Exp coefficient RT.  Produce polarization-derived spectra and an absolute-path semi-log HHG figure.

**Step 4: Final verification and reviews**

Run every retained metadata, metric, construction, symmetry, localization, operator, solver/density, SCF, checkpoint, and coefficient-RT fixture on MPI 1/2/4/8; the route/removal contracts; genuine Si64 checker; morphology checker; `git diff --check`; and one final clean-first committed-HEAD overlay.  Perform final specification and code-quality reviews and resolve all Critical/Important findings.

**Step 5: Commit and dual push**

```bash
git add src/io/inputoutput.f90 tests/dg \
  docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md
git commit -m "test(dg): validate bounded-memory Wannier polarization HHG"
git push origin codex/wpw-s-orthogonal-complement
git push upstream HEAD:codex/wpw-s-orthogonal-complement
```

Verify local HEAD, origin branch, and upstream branch resolve to the identical commit before reporting the inversion, memory, LR, H2/H4 morphology, and figure evidence.
