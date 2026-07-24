# DG-DC Local-Basis SIPG Ground-State Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Produce the Task 5 DG ground state by retaining the existing DC
fragment-local orbital machinery and adding only physical-face SIPG coupling,
without concatenating local candidates into global real-space states.

**Architecture:** DC continues to generate and update
`natom_fragment * dg_dc_candidate_orbitals_per_atom` local orbitals per
fragment.  A new default-off local-basis adapter assembles the global
Hermitian Hamiltonian from existing DC volume blocks plus SIPG neighbor blocks,
solves the requested total-system `nstate` bands, and reconstructs density for
the existing DC Hartree/XC/`vlocal`/mixing path.  The adaptive continuation and
transactional checkpoint consume this representation.

**Tech Stack:** Fortran 2008, MPI, BLAS/LAPACK, existing SALMON DC structures,
existing physical-face SIPG helpers, Python contract runners.

---

### Task 1: Separate fragment-local basis size from global band count

**Files:**
- Create: `src/gs/dc/dg_dc_local_basis_ground_state.f90`
- Modify: `src/gs/dc/CMakeLists.txt`
- Create: `tests/dg/test_dg_dc_local_basis_layout_mpi.f90`
- Create: `tests/dg/run_dg_dc_local_basis_layout_mpi.py`
- Modify: `docs/plans/2026-07-24-dg-dc-local-basis-sipg-implementation.md`

**Step 1: Write the failing layout test**

Construct two MPI fragments with unequal local basis sizes.  Request three
global bands and verify the production layout reports:

```text
local_basis_count(rank) = supplied local count
global_basis_count = sum(local counts)
global_band_count = requested nstate
```

Require `global_band_count` to remain three rather than
`global_basis_count`.  Reject zero basis counts, rank-disagreeing band counts,
and duplicate global basis offsets.

**Step 2: Verify RED**

Run:

```bash
python3 tests/dg/run_dg_dc_local_basis_layout_mpi.py
```

Expected: FAIL because the local-basis layout module does not exist.

**Step 3: Implement the minimal distributed layout**

Add `s_dg_dc_local_basis_layout` containing local/global basis counts, global
band count, rank offsets, fragment IDs, and fingerprints.  Build offsets with
MPI collectives and validate them collectively in O(number of fragments).

Do not integrate the layout into the production solver yet: unequal local row
counts are unsafe until the distributed coefficient representation exists.
Task 3 performs the atomic production switch and removes
`expand_dg_dc_global_candidate_axis` from production.

**Step 4: Verify GREEN**

Run:

```bash
python3 tests/dg/run_dg_dc_local_basis_layout_mpi.py
git diff --check
```

Expected: PASS.

**Step 5: Run specification and code-quality reviews**

Confirm the new layout API keeps local basis rows and global band columns
distinct and that no partial production integration was introduced.  Resolve
every Critical/Important finding before continuing.

**Step 6: Commit**

```bash
git add src/gs/dc/dg_dc_local_basis_ground_state.f90 \
  src/gs/dc/CMakeLists.txt \
  tests/dg/test_dg_dc_local_basis_layout_mpi.f90 \
  tests/dg/run_dg_dc_local_basis_layout_mpi.py \
  docs/plans/2026-07-24-dg-dc-local-basis-sipg-implementation.md
git commit -m "refactor(dg): separate local basis from global bands"
```

### Task 2: Assemble DC volume and SIPG interface matrices

**Files:**
- Modify: `src/gs/dc/dg_dc_local_basis_ground_state.f90`
- Modify: `src/gs/dc/dg_dc_ground_state_adapter.f90`
- Create: `tests/dg/test_dg_dc_local_basis_sipg_mpi.f90`
- Create: `tests/dg/run_dg_dc_local_basis_sipg_mpi.py`

**Step 1: Write the failing SIPG assembly test**

Use two fragments with analytic local basis functions and one physical shared
face.  Supply existing DC volume blocks and verify:

- `H(lambda=0)` equals the block-diagonal DC volume Hamiltonian;
- `H(lambda=1)` equals volume plus the analytic consistency, symmetry, and
  penalty blocks;
- opposite-face contributions occur once under canonical ownership;
- the global matrix is Hermitian;
- an auxiliary buffer boundary contributes nothing;
- periodic wrapping selects the physical neighbor and image exactly once.

**Step 2: Verify RED**

Run:

```bash
python3 tests/dg/run_dg_dc_local_basis_sipg_mpi.py
```

Expected: FAIL because local-basis SIPG assembly is absent.

**Step 3: Implement minimal assembly**

Reuse the existing DC fragment orbital arrays and existing volume Hamiltonian
action.  Exchange only boundary values and normal derivatives needed by
`evaluate_dg_dc_local_sipg_face`.  Contract the face action with the local
basis functions to form neighbor matrix blocks, then assemble

```text
H(lambda) = H_DC_volume + lambda * H_SIPG_interface
```

Do not call LCFO, WPW, checkpoint consumers, or RT code.

**Step 4: Verify GREEN**

Run the new MPI fixture plus:

```bash
python3 tests/dg/run_dg_dc_ground_state_adapter_mpi.py
python3 tests/dg/check_dg_dc_local_periodic_route.py
git diff --check
```

Expected: PASS.

**Step 5: Review and commit**

Run specification and quality reviews.  Resolve Critical/Important findings.

```bash
git add src/gs/dc/dg_dc_local_basis_ground_state.f90 \
  src/gs/dc/dg_dc_ground_state_adapter.f90 \
  tests/dg/test_dg_dc_local_basis_sipg_mpi.f90 \
  tests/dg/run_dg_dc_local_basis_sipg_mpi.py
git commit -m "feat(dg): add SIPG coupling to DC local basis"
```

### Task 3: Solve global bands and reuse the DC SCF path

**Files:**
- Modify: `src/gs/dc/dg_dc_local_basis_ground_state.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `src/gs/scf_iteration_dft.f90`
- Create: `tests/dg/test_dg_dc_local_basis_scf_mpi.f90`
- Create: `tests/dg/run_dg_dc_local_basis_scf_mpi.py`
- Modify: `tests/dg/check_dg_dc_local_periodic_route.py`

**Step 1: Write the failing SCF test**

Use a small two-fragment Hermitian model with known eigenpairs.  Verify that:

- the solver returns the requested global band count, not the local basis sum;
- empty bands are retained;
- eigenvectors are orthonormal in the local-basis metric;
- density reconstructed from coefficients and fragment-local orbitals matches
  the analytic density;
- the existing DC density, Hartree, XC, `vlocal`, and mixing callbacks are
  invoked once per SCF iteration;
- lambda zero reproduces the DC result and lambda one includes SIPG coupling.

**Step 2: Verify RED**

Run:

```bash
python3 tests/dg/run_dg_dc_local_basis_scf_mpi.py
```

Expected: FAIL because the global local-basis solve is absent.

**Step 3: Implement the minimal solve**

Use the existing DC CG and BLAS Gram--Schmidt to update fragment-local
orbitals.  Assemble the local-basis overlap and Hamiltonian matrices and solve
the Hermitian generalized eigenproblem for total-system `nstate` bands with the
configured distributed eigensolver.  Reconstruct core density from the
eigenvectors and occupations, then call the existing DC potential update and
mixing path.

Remove the production calls to:

- `expand_dg_dc_global_candidate_axis`;
- `solve_nodal_ground_state_cg_mpi`;
- nodal Cholesky orthogonalization.

Keep those neutral RT/common modules available only for their existing
non-production tests and consumers.

**Step 4: Verify GREEN**

Run:

```bash
python3 tests/dg/run_dg_dc_local_basis_scf_mpi.py
python3 tests/dg/run_dg_dc_ground_state_mpi.py
python3 tests/dg/check_dg_dc_local_periodic_route.py
git diff --check
```

Expected: PASS.

**Step 5: Review and commit**

Run specification and quality reviews.  Resolve Critical/Important findings.

```bash
git add src/gs/dc/dg_dc_local_basis_ground_state.f90 \
  src/gs/main_dft.f90 src/gs/scf_iteration_dft.f90 \
  tests/dg/test_dg_dc_local_basis_scf_mpi.f90 \
  tests/dg/run_dg_dc_local_basis_scf_mpi.py \
  tests/dg/check_dg_dc_local_periodic_route.py
git commit -m "feat(dg): run SIPG ground state through DC SCF"
```

### Task 4: Adapt continuation and checkpoint to the local-basis state

**Files:**
- Modify: `src/gs/dc/dg_dc_ground_state.f90`
- Modify: `src/common/dg_ground_state_checkpoint.f90`
- Modify: `src/gs/main_dft.f90`
- Modify: `tests/dg/test_dg_dc_ground_state_mpi.f90`
- Modify: `tests/dg/test_dg_ground_state_checkpoint_mpi.f90`
- Modify: `tests/dg/run_dg_ground_state_checkpoint_mpi.py`

**Step 1: Write RED checkpoint and continuation tests**

Require the checkpoint to store local-basis ownership, offsets, overlap and
Hamiltonian fingerprints, requested global band count, eigenvalues,
occupations, coefficient shards, density/potentials, topology, continuation
history, and every acceptance diagnostic.

Reject concatenated-nodal labeling, incomplete coefficient shards, stale basis
fingerprints, non-Hermitian matrix diagnostics, unaccepted lambda, mixed final
verification, corruption, and rank disagreement.

**Step 2: Verify RED**

Run:

```bash
python3 tests/dg/run_dg_ground_state_checkpoint_mpi.py
python3 tests/dg/run_dg_dc_ground_state_mpi.py
```

Expected: FAIL on the new local-basis schema requirements.

**Step 3: Implement the schema and production publication**

Keep transactional generation-specific payloads and atomic manifest
publication.  Replace global nodal-orbital payload meaning with explicit
fragment-local basis and global coefficient payloads.  Publish only after
lambda one, accepted SCF diagnostics, and an unmixed fixed-point verification.

**Step 4: Verify GREEN**

Run both MPI fixtures, the route contract, and `git diff --check`.

**Step 5: Review and commit**

Run specification and quality reviews.  Resolve Critical/Important findings.

```bash
git add src/gs/dc/dg_dc_ground_state.f90 \
  src/common/dg_ground_state_checkpoint.f90 src/gs/main_dft.f90 \
  tests/dg/test_dg_dc_ground_state_mpi.f90 \
  tests/dg/test_dg_ground_state_checkpoint_mpi.f90 \
  tests/dg/run_dg_ground_state_checkpoint_mpi.py
git commit -m "feat(dg): checkpoint local-basis SIPG ground state"
```

### Task 5: Parent-prerequisite overlay and Si64 gate

**Files:**
- Modify: `docs/plans/2026-07-24-dg-dc-local-periodic-wannier.md`
- Modify: `docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md`

**Step 1: Run the complete focused verification**

Run every new local-basis MPI fixture, existing Task 1--5 focused fixtures, and
`git diff --check`.  Run final specification and code-quality reviews and
resolve every Critical/Important finding.

**Step 2: Build a fresh parent-prerequisite overlay**

Record branch, HEAD, clean status, and parent porcelain hash.  Build a new
read-only-parent `/tmp` overlay containing committed branch changes plus only
the required parent prerequisites:

```bash
cmake -S <overlay-source> -B <overlay-build> \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_Fortran_COMPILER=/opt/homebrew/bin/mpifort \
  -DCMAKE_Fortran_FLAGS=-fallow-invalid-boz \
  -DUSE_MPI=ON -DUSE_SCALAPACK=ON -DUSE_EIGENEXA=ON -DUSE_WANNIER90=OFF
cmake --build <overlay-build> --target salmon -j1
```

Expected: `[100%] Built target salmon`.

**Step 3: Run the required Si64 matrix**

Use fresh directories, eight MPI ranks, and `OMP_NUM_THREADS=1`.  Compare
handoff tolerances `1E-2`, `1E-3`, `1E-4`, at least two buffer widths, and at
least two equivalent fragment decompositions.  Retain 40 local orbitals per
atom and the requested total-system global band count.

Record handoff state, local/global basis dimensions, global band count, lambda
history, rollbacks, SCF/CG counts, residuals, charge, energy, orthogonality,
Hermiticity, face balance, fixed-point result, runtime, and checkpoint.

**Step 4: Enforce the stop gate**

If no configuration passes every ground-state gate, document the failures,
stop here, and do not begin Wannier Tasks 6--8.  Do not promote normal DC,
LCFO, WPW, checkpoint consumers, or RT.

**Step 5: Document and commit**

Update both handoff documents with commands, raw-log-derived values, ratios,
overlay provenance, and the gate decision.

```bash
git add docs/plans/2026-07-24-dg-dc-local-periodic-wannier.md \
  docs/plans/2026-07-19-wpw-density-carrying-seed-session-handoff.md
git commit -m "docs(dg): record local-basis Si64 ground-state gate"
```
