# Nonlocal Projector Range and Symmetry Diagnostic Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Measure remote-fragment Wannier/projector coupling and operation-5 nonlocal symmetry directly from the total-system projector representation.

**Architecture:** Add a bounded-tile MPI diagnostic that uses the same `dc%ppg_tot` data as production `hpsi`, accumulates projector overlaps collectively, and publishes only distance-class and symmetry-pair receipts. Run it once after Wannier centers are known and remove the normalization-incompatible fragment reference from the SCF Hamiltonian builder.

**Tech Stack:** Fortran 2008, MPI, SALMON `s_pp_grid`, Python route checks.

---

### Task 1: Specify distance-class and symmetry-pair receipts

**Files:**
- Create: `tests/dg/test_dg_nonlocal_projector_range_mpi.f90`
- Create: `tests/dg/run_dg_nonlocal_projector_range_mpi.py`
- Modify: `src/gs/dc/CMakeLists.txt`

**Step 1: Write the failing synthetic MPI test**

Construct a periodic two-fragment grid with four atoms, short-range synthetic
projectors, and Wannier functions whose projector norms are known separately
for local, adjacent, and remote atoms.  Require a tile-width-independent result.

**Step 2: Add operation pairing cases**

Require zero pair defect for an exact twofold rotation, and collective rejection
for a missing same-species atom or projector partner.

**Step 3: Run RED**

Run:

```bash
python3 tests/dg/run_dg_nonlocal_projector_range_mpi.py
```

Expected: compilation fails because `dg_nonlocal_projector_range_diagnostic`
does not exist.

**Step 4: Commit the failing test**

```bash
git add tests/dg/test_dg_nonlocal_projector_range_mpi.f90 tests/dg/run_dg_nonlocal_projector_range_mpi.py src/gs/dc/CMakeLists.txt
git commit -m "test: specify nonlocal projector range diagnostic"
```

### Task 2: Implement bounded total-system projector analysis

**Files:**
- Create: `src/gs/dc/dg_nonlocal_projector_range.f90`
- Test: `tests/dg/test_dg_nonlocal_projector_range_mpi.f90`

**Step 1: Define the result type and contract**

Store local/adjacent/remote contribution totals, maximum remote fraction,
operation pair defect, unmatched-channel count, and peak workspace bytes.
Validate all dimensions, finite inputs, unique physical grid ownership, and the
one-to-one same-species atom map collectively.  Accept explicit dense Wannier
and projector representations so angular-channel rotations are never inferred
from atom indices alone.

**Step 2: Accumulate overlap tiles**

For a tile of at most 16 Wannier functions, accumulate local
`conjg(projector_value)*wannier_value`, perform one `MPI_Allreduce` for all
channels in the tile, and immediately reduce contributions into distance
classes and symmetry-pair defects.

**Step 3: Enforce bounded storage**

Use `int64` extent checks before allocation.  Do not retain a
`Nwannier*Nprojector` array outside the current tile and do not issue per-channel
collectives.

**Step 4: Run GREEN**

Run:

```bash
python3 tests/dg/run_dg_nonlocal_projector_range_mpi.py
```

Expected: PASS for MPI 1, 2, 4, and 8.

**Step 5: Commit**

```bash
git add src/gs/dc/dg_nonlocal_projector_range.f90 tests/dg/test_dg_nonlocal_projector_range_mpi.f90 tests/dg/run_dg_nonlocal_projector_range_mpi.py src/gs/dc/CMakeLists.txt
git commit -m "feat: diagnose distributed nonlocal projector range"
```

### Task 3: Connect the one-shot production diagnostic

**Files:**
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/check_dg_overlapping_wannier_route.py`

**Step 1: Write failing route assertions**

Require the diagnostic exactly once after periodic Wannier centers and final
row-owned values exist, before Hybrid-SCF.  Reject calls from
`ow_build_hamiltonian`.  Require removal of
`full_cell/reference_nonlocal_difference` and its fragment reference assembly.

**Step 2: Run RED**

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
```

Expected: failure because the temporary reference diagnostic is still inside
every Hamiltonian build.

**Step 3: Implement the production adapter**

Map `ow_core_values` to `dc%mg_tot` in tiles using the existing cached
redistribution schedule, call the new diagnostic with `dc%ppg_tot`, physical
atomic positions/species, Wannier centers, and affine operation 5, and print one
rank-zero receipt.  Construct complete radial `(2l+1)` projector blocks in
SALMON's real-harmonic convention and combine them with the exact periodic atom
permutation.  Do not change `hrows`, `ow_srows`, `ow_rhorows`, density, or
any acceptance tolerance.

**Step 4: Verify focused tests and build**

```bash
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/run_dg_nonlocal_projector_range_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_full_cell_mpi.py
cmake --build /Users/otobetoshihito/SALMON-dev/verification/20260818-onepass-production-build -j 4
git diff --check
```

Expected: all PASS and `Built target salmon`.

**Step 5: Commit**

```bash
git add src/gs/main_dft.f90 tests/dg/check_dg_overlapping_wannier_route.py
git commit -m "diag: report remote nonlocal Wannier coupling"
```

### Task 4: Verify Si64 without changing physics

**Files:**
- Runtime artifact only: `verification/20260823-si64-nonlocal-range-symmetry-mpi8-omp1/`

**Step 1: Run Si64**

Run the existing Si64 Hybrid-SCF input with MPI 8,
`OMP_NUM_THREADS=1`, standard input redirection, and no time cutoff.

**Step 2: Record receipts**

Require finite local/adjacent/remote totals, complete projector pairing, the
operation-5 pair defect, peak workspace, wall time, and maximum RSS per rank.

**Step 3: Interpret before fixing**

If remote coupling is material, reject the fragment reference as an oracle.  If
operation-5 pairing is complete but defective, trace the total-system p-channel
action.  Do not symmetrize or relax tolerances in either case.

**Step 4: Run regressions**

```bash
python3 tests/dg/run_dg_nonlocal_projector_range_mpi.py
python3 tests/dg/run_dg_overlapping_wannier_full_cell_mpi.py
python3 tests/dg/check_dg_overlapping_wannier_route.py
python3 tests/dg/check_dg_hybrid_self_consistent_route.py
git diff --check
```

Expected: all PASS.
