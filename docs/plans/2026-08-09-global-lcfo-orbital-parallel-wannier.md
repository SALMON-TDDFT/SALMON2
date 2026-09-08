# Global LCFO Orbital-Parallel Wannier Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Construct every retained overlapping-Wannier orbital from one full-system DC-LCFO eigenspace, localize without changing the one-particle density matrix, redistribute by periodic center, and publish only an accepted V3 for generalized-eigenvalue Exp RT.

**Architecture:** Keep LCFO basis data spatially distributed by fragment, form the localization gauge in the bounded retained coefficient space, and transpose bounded real-space orbital batches from fragment contraction ownership to Wannier-index ownership and finally to center-fragment ownership. Derive symmetry only from the full atomic configuration, represent translations separately from the point co-group, and accept the GS only after projector, density, covariance, redistribution, and stationarity gates.

**Tech Stack:** Fortran 2008, MPI `Allreduce`/`Alltoallv`, OpenMP, LAPACK/ScaLAPACK, EigenExa, spglib, CMake, Python contract/evidence tests.

---

Every task below follows the same non-optional completion protocol:

1. record a behavioral RED that fails for the intended missing behavior;
2. implement only the task scope and run its focused verification on 1/2/4/8 MPI ranks where applicable;
3. perform a specification review against `docs/plans/2026-08-09-global-lcfo-orbital-parallel-wannier-design.md`;
4. perform a separate code-quality review;
5. resolve every Critical and Important finding and repeat affected tests;
6. run `git diff --check` and a fresh clean-first parent-prerequisite overlay build;
7. commit only after all preceding gates pass.

The overlay is made from `git archive HEAD`, with only the current task diff applied. Configure Release with MPI, ScaLAPACK, EigenExa, and spglib enabled and Wannier90 disabled. Build the EigenExa prerequisite with `-j1`, then run `cmake --build <overlay-build> --clean-first -j4`. Never use a prior in-tree object as acceptance evidence.

### Task 1: Retain the coherent LCFO eigenspace and occupation invariant

**Files:**
- Modify: `src/gs/dc/lcfo.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Create: `tests/dg/check_global_lcfo_wannier_contract.py`

**Step 1: Write and run the RED**

Add focused cases requiring an explicit retained rank independent of `dc%nstate_tot`, all requested EigenExa columns, the matching occupation vector, and reconstruction from the single LCFO coefficient matrix. Require rejection when the LCFO metric rank is too small, a coefficient is non-finite, or a Gamma-real coefficient has a non-negligible imaginary part. For Si64 require retained rank 384 and the existing normal DC LCFO output to remain unchanged.

Run:

```bash
python3 tests/dg/check_global_lcfo_wannier_contract.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: FAIL because the current optional `occupied_count` path neither carries the full retained-rank/occupation contract nor proves all 384 columns have one coherent provenance.

**Step 2: Implement the bounded coefficient interface**

Replace the occupied-only request with an explicit retained-column request and returned occupation metadata. Keep each fragment root's LCFO-basis row block and requested coefficient columns; do not gather fragment basis functions or full real-space orbitals. Validate dimensions, metric rank, finiteness, and real-Gamma scope before publication. Preserve the unchanged normal `lcfo.f90` write path and EigenExa route.

**Step 3: Verify, review, overlay-build, and commit**

Run the two contract checks, the construction fixture on 1/2/4/8 ranks, normal LCFO source checks, and the mandatory completion protocol. Commit:

```bash
git add src/gs/dc/lcfo.f90 src/gs/main_dft.f90 \
  tests/dg/check_global_lcfo_wannier_contract.py \
  tests/dg/check_dg_overlapping_wannier_route.py \
  tests/dg/test_dg_overlapping_wannier_construction_mpi.f90
git commit -m "feat(dg): retain coherent global LCFO Wannier space"
```

### Task 2: Localize only inside equal-occupation blocks

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_localization.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_localization_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Write and run the RED**

Add integer-occupation, fractional equal-occupation, unequal-occupation rejection, singleton-block, nonidentity metric, degenerate rotation, deterministic sign, and incomplete atomic-projector cases. Compare `C f C^T` before and after localization and require metric-orthogonality and the complete periodic-spread gradient gate.

Run the localization and construction MPI fixtures on one rank. Expected: FAIL because the current path treats the retained space as one symmetry/localization block and does not expose a density-matrix-preservation receipt.

**Step 2: Implement block-metric localization**

Partition occupations using a tolerance fixed before localization. In each block, measure the LCFO metric, transform to an orthonormal frame, use atomic projections only for initialization, minimize the complete periodic spread with real antisymmetric generators, and map the gauge back with metric bookkeeping. Reject complex, cross-occupation, non-monotone, or unconverged results. Order and sign results deterministically using periodic centers and global grid identifiers.

**Step 3: Verify, review, overlay-build, and commit**

Run localization and construction on 1/2/4/8 ranks, finite-difference gradient checks, route checks, and the mandatory completion protocol. Commit:

```bash
git add src/gs/dc/dg_overlapping_wannier_localization.f90 \
  src/gs/dc/dg_overlapping_wannier_construction.f90 src/gs/main_dft.f90 \
  tests/dg/test_dg_overlapping_wannier_localization_mpi.f90 \
  tests/dg/test_dg_overlapping_wannier_construction_mpi.f90
git commit -m "feat(dg): preserve occupations during global localization"
```

### Task 3: Materialize by Wannier index and redistribute by periodic center

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_types.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`
- Create: `tests/dg/test_dg_overlapping_wannier_redistribution_mpi.f90`
- Create: `tests/dg/run_dg_overlapping_wannier_redistribution_mpi.py`

**Step 1: Write and run the RED**

Cover balanced ownership when the orbital count is and is not divisible by MPI size, bounded batches, `Alltoallv` count overflow, periodic centers crossing cell faces, exact face/edge/corner tie-breaking, centers outside the construction rank, norm/checksum conservation, and complete destination core+buffer coverage. Include a memory assertion that no rank holds all 384 full-system orbitals.

Run the new fixture on 1/2/4/8 ranks. Expected: FAIL because no orbital/spatial transpose and center-fragment redistribution API exists.

**Step 2: Implement the two transposes**

Contract each bounded transformed-coefficient batch on fragment-resident LCFO basis tails, integrate only unique cores, and use `MPI_Alltoallv` to assemble complete orbitals on balanced Wannier owners. Compute periodic moments and deterministic owners, then transpose each orbital's destination core+periodic buffer to its center fragment. Check all integer counts before allocation/communication and release full-system batches after receipts pass. Use OpenMP over grid contractions without changing deterministic reductions.

**Step 3: Verify, review, overlay-build, and commit**

Run redistribution and construction fixtures on 1/2/4/8 ranks with floating-point traps and memory diagnostics, then the mandatory completion protocol. Commit:

```bash
git add src/gs/dc/dg_overlapping_wannier_construction.f90 \
  src/gs/dc/dg_overlapping_wannier_types.f90 src/gs/main_dft.f90 \
  tests/dg/test_dg_overlapping_wannier_construction_mpi.f90 \
  tests/dg/test_dg_overlapping_wannier_redistribution_mpi.f90 \
  tests/dg/run_dg_overlapping_wannier_redistribution_mpi.py
git commit -m "feat(dg): redistribute global Wanniers by periodic center"
```

### Task 4: Represent general full-system crystallographic symmetry

**Files:**
- Modify: `src/gs/dc/dg_overlapping_wannier_symmetry.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_construction.f90`
- Modify: `tests/dg/check_dg_crystallographic_point_groups.py`
- Modify: `tests/dg/test_dg_overlapping_wannier_fragment_symmetry_mpi.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_construction_mpi.f90`

**Step 1: Write and run the RED**

Test `C1`, inversion, noncentrosymmetric rotations, nonsymmorphic screw/glide, centers whose symmetry center lies outside their fragment, and deliberately displaced lower-symmetry structures. Require translation-subgroup and point-co-group closure including the translation cocycle, plus measured `H` and density-projector commutators. Separate interior covariance failure from boundary-layer LCFO stitching error.

Run crystallographic and MPI symmetry fixtures. Expected: FAIL until the full affine operation data, cocycle, commutator gates, and boundary calibration are connected to the global LCFO space.

**Step 2: Implement general affine covariance**

Build valid operations only from the full instantaneous atomic configuration. Factor pure translations from point representatives while retaining fractional translations and product cocycles. Apply real-space permutations at Gamma, measure representations in the LCFO metric, and require center-orbit closure and final-operator covariance without demanding individual Wannier invariance or a common fixed point. Do not introduce displaced-structure acceptance for an undisplaced target.

**Step 3: Verify, review, overlay-build, and commit**

Run crystallographic checks and symmetry/construction fixtures on 1/2/4/8 ranks, then the mandatory completion protocol. Commit:

```bash
git add src/gs/dc/dg_overlapping_wannier_symmetry.f90 \
  src/gs/dc/dg_overlapping_wannier_construction.f90 \
  tests/dg/check_dg_crystallographic_point_groups.py \
  tests/dg/test_dg_overlapping_wannier_fragment_symmetry_mpi.f90 \
  tests/dg/test_dg_overlapping_wannier_construction_mpi.f90
git commit -m "feat(dg): enforce full-system affine covariance"
```

### Task 5: Gate reconstructed GS and mandatory V3 provenance

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_checkpoint.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_density.f90`
- Modify: `src/gs/dc/dg_overlapping_wannier_operators.f90`
- Modify: `tests/dg/test_dg_overlapping_wannier_checkpoint_mpi.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`
- Modify: `tests/dg/check_si64_overlapping_wannier_gate.py`

**Step 1: Write and run the RED**

Require exact rank, occupation spectrum, electron count, LCFO-projector/density/occupied-energy agreement, inversion-odd density where applicable, center ownership, buffer coverage, all operator covariance, and field-off stationarity receipts. Require all new fingerprints in V3 and reject legacy V3 files missing them.

Run route, checkpoint, and synthetic GS gates. Expected: FAIL because the current V3 schema lacks the complete global-LCFO, occupation-block, affine-cocycle, and redistribution provenance.

**Step 2: Integrate the accepted GS route**

Assemble metric and Hamiltonian from center-owned core+buffer orbitals, solve the fixed-rank generalized eigenproblem, compare against the authoritative LCFO density matrix, and publish V3 only after every gate passes. Keep normal DC LCFO+EigenExa unchanged and retain only V3-backed generalized-eigenvalue Exp coefficient RT; do not add a fallback route.

**Step 3: Verify, review, overlay-build, and commit**

Run checkpoint, metric, density, operator, solver, SCF, RT, route, and obsolete-route fixtures on their supported 1/2/4/8 ranks, then the mandatory completion protocol. Commit:

```bash
git add src/gs/main_dft.f90 src/gs/dc/dg_overlapping_wannier_checkpoint.f90 \
  src/gs/dc/dg_overlapping_wannier_density.f90 \
  src/gs/dc/dg_overlapping_wannier_operators.f90 tests/dg
git commit -m "feat(dg): gate global LCFO Wannier V3 publication"
```

### Task 6: Genuine ideal-Si64 GS, LR, and polarization HHG evidence

**Files:**
- Modify: `tests/dg/data/si64_overlapping_wannier_rt/input_gs.in`
- Modify: `tests/dg/run_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/check_si64_overlapping_wannier_response_hhg.py`
- Modify: `tests/dg/test_si64_harmonic_morphology.py`
- Modify: `docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md`

**Step 1: Record the RED evidence contract**

Require strict undisplaced ideal Si64, 384 retained/128 occupied states, successful global-LCFO and V3 receipts, field-off stationarity, polarization-derived LR, a longer laser pulse with an Exp-compatible large time step, polarization-derived semi-log HHG, and explicit even-order peak/dip/slope classification. Reject displaced inputs and current-derived primary spectra.

Run the morphology and evidence checkers. Expected: FAIL until fresh genuine artifacts satisfy the new provenance and spectral contract.

**Step 2: Run genuine physics evidence**

From a successful clean-first overlay binary, run DC-SCF, LCFO+EigenExa, global Wannier GS, and V3 on eight ranks. Then run field-off, LR, and long-pulse laser Exp coefficient RT. Compute spectra from polarization; retain current only as a secondary cross-check. Generate the semi-log HHG figure and report whether H2/H4 are peaks or dips, their local slopes, odd/even suppression, inversion-odd polarization, and field-off drift. Treat the small system as qualitative evidence, not quantitative material prediction.

**Step 3: Final reviews, verification, commit, and dual push**

Repeat all focused contracts and MPI fixtures, perform the final clean-first committed-HEAD overlay build, repeat specification and code-quality reviews, and resolve every Critical/Important finding. Commit the evidence and verify that local, origin branch, and upstream branch resolve to the identical commit before reporting completion.

```bash
git add tests/dg docs/plans/2026-08-01-overlapping-wannier-polarization-hhg-results.md
git commit -m "test(dg): validate global LCFO Wannier polarization HHG"
git push origin codex/wpw-s-orthogonal-complement
git push upstream HEAD:codex/wpw-s-orthogonal-complement
```
