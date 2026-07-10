# SAWF Fragment Alignment Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add an explicit, validated periodic translation workflow that makes supplied crystal symmetry operations act as whole-fragment permutations before constructing the fragment-local SAWF representation.

**Architecture:** A standalone Python preparation tool searches periodic translations, translates the atom coordinates and symmetry translations together, and writes dedicated output files only after atom, integer-grid, group, and fragment-map checks pass. SALMON reads the selected symmetry file through a namelist variable and independently validates that each normalized operation induces an integer grid permutation whose image of every fragment core is exactly one complete fragment core. For even meshes, inversion centers are expected at half-grid positions; correctness is expressed through the integer index map, never through a requirement that the center coincide with a grid point.

**Tech Stack:** Python 3 standard library and NumPy-free exact/rational geometry where practical, Fortran 2003/2008, SALMON namelist and MPI helpers, existing SAWF symmetry types, focused Python/Fortran driver tests, CMake MPI+EigenExa+Wannier90 build.

---

### Task 1: Add an explicit SAWF symmetry-file input

**Files:**
- Modify: `src/io/salmon_global.f90:529`
- Modify: `src/io/inputoutput.f90:665-681`
- Modify: `src/io/inputoutput.f90:1184-1200`
- Modify: `src/io/inputoutput.f90:1913-1934`
- Modify: `src/io/inputoutput.f90:2982-3000`
- Modify: `src/gs/dc/lcfo_flux.f90:6180-6270`
- Modify: `tests/dg/check_sawf_input_and_build.py`

**Step 1: Write the failing source/input tests**

Extend `check_sawf_input_and_build.py` to require:

```python
required_global = r"character\s*\([^)]*\)\s*::\s*wannier_symmetry_file"
required_default = r"wannier_symmetry_file\s*=\s*['\"]sym\.dat['\"]"
required_bcast = r"call\s+comm_bcast\s*\(\s*wannier_symmetry_file"
required_log = r"['\"]wannier_symmetry_file['\"]"
```

Also require the file-mode loader to resolve `wannier_symmetry_file`, not a hardcoded `sym.dat`.

**Step 2: Run the test and verify it fails**

Run:

```bash
python3 tests/dg/check_sawf_input_and_build.py
```

Expected: FAIL because `wannier_symmetry_file` is absent.

**Step 3: Implement the namelist variable**

Add:

```fortran
character(256) :: wannier_symmetry_file
```

to `salmon_global`, include it in `namelist/dc/`, default it to `sym.dat`, broadcast it, and write it to `variables.log`. In `build_sawf_dmn_import`, resolve relative paths against `import_run_root_dir()` while preserving absolute paths. Reject an empty path when `wannier_site_symmetry='file'`.

**Step 4: Run focused tests and build**

Run:

```bash
python3 tests/dg/check_sawf_input_and_build.py
python3 tests/dg/check_sawf_win_integration.py
cmake --build build-mpi-eigenexa-wannier-lib -j 2
```

Expected: all tests PASS and `salmon` builds.

**Step 5: Commit the input slice**

```bash
git add src/io/salmon_global.f90 src/io/inputoutput.f90 src/gs/dc/lcfo_flux.f90 tests/dg/check_sawf_input_and_build.py
git commit -m "feat: configure SAWF symmetry input file"
```

### Task 2: Implement exact grid and fragment compatibility checks in the preparation tool

**Files:**
- Create: `tools/align_periodic_structure_to_fragments.py`
- Create: `tests/dg/check_sawf_fragment_alignment_tool.py`

**Step 1: Write failing unit tests for the discrete geometry**

The tests import the tool module and cover these cases:

```python
def test_even_grid_inversion_uses_half_grid_center():
    # 32 points, two 16-point fragments, i' = 15-i.
    op = SymOp(W=neg_identity(), tau=(15/32, 15/32, 15/32))
    assert grid_shift(op, mesh=(32, 32, 32)) == (15, 15, 15)
    assert inversion_center_grid_units(op) == (7.5, 7.5, 7.5)
    assert fragment_targets(op, (32, 32, 32), (2, 2, 2)) == {
        0: {0}, 1: {1}, 2: {2}, 3: {3},
        4: {4}, 5: {5}, 6: {6}, 7: {7},
    }

def test_unaligned_inversion_splits_fragment():
    op = SymOp(W=neg_identity(), tau=(1/8, 1/8, 1/8))
    assert max(map(len, fragment_targets(op, (32, 32, 32), (2, 2, 2)).values())) == 8

def test_even_grid_does_not_require_center_on_grid_point():
    # A validator that requires integer center coordinates must fail this test.
    assert is_integer_grid_permutation(aligned_inversion, (32, 32, 32))
```

Also test identity, non-integral grid shifts, mixed species, atom bijection,
invalid fragment divisibility, and an axis-swapping operation on unequal mesh
counts that must be rejected when it is not a discrete-grid bijection.

**Step 2: Run the test and verify it fails**

```bash
python3 tests/dg/check_sawf_fragment_alignment_tool.py
```

Expected: FAIL because the tool module does not exist.

**Step 3: Implement the geometry core**

Represent fractional coordinates with `fractions.Fraction` while parsing decimal input. For a symmetry operation `(W,tau)` and mesh `N`, derive the index map through fractional coordinates:

```text
j_alpha = sum_beta W(alpha,beta) i_beta N_alpha/N_beta
          + N_alpha tau_alpha                         (mod N_alpha)
```

after accounting for the repository's real-space grid origin convention.  The
simple `j=W i+q` form is valid only when the operation mixes axes with equal
mesh counts, as in cubic C64.  Accept an operation only when the exact rational
expression is an integer and a bijection for every grid index.  Do not
calculate an inversion center as a grid-point index. For `W=-I`, report the
physical center as `q/2`, which is half-integer when `q` is odd.

Enumerate every source-core grid index for the small preprocessing validation and collect target fragment IDs. Accept an operation only when each source set has cardinality one and its image contains exactly the target core's number of points.

**Step 4: Implement atom and operation validation**

For translation `a`, compute:

```text
r'_atom = r_atom + a                 (mod 1)
tau'_g  = tau_g + (I - W_g) a       (mod 1)
```

Verify a same-species periodic atom bijection for every transformed operation, unique identity, inverse existence, and multiplication closure within the declared tolerance.

**Step 5: Run the unit tests**

```bash
python3 tests/dg/check_sawf_fragment_alignment_tool.py
```

Expected: PASS, including the `7.5` grid-spacing inversion center.

**Step 6: Commit the geometry core**

```bash
git add tools/align_periodic_structure_to_fragments.py tests/dg/check_sawf_fragment_alignment_tool.py
git commit -m "feat: validate symmetry-compatible fragment maps"
```

### Task 3: Add deterministic periodic-translation search and safe file output

**Files:**
- Modify: `tools/align_periodic_structure_to_fragments.py`
- Modify: `tests/dg/check_sawf_fragment_alignment_tool.py`

**Step 1: Add failing CLI tests**

Use a temporary C64-like atom file and `sym.dat`. Require the CLI to:

- find `a=(11/64,11/64,11/64)` for the current 32-grid, 2-fragment inversion case;
- produce `tau'=(15/32,15/32,15/32)`;
- leave input files byte-for-byte unchanged;
- refuse to overwrite outputs without `--force`;
- write no partial output when no compatible translation exists;
- preserve atom species and non-coordinate columns.
- accept explicit nonnegative buffer widths and reject an axis-exchanging
  operation when the exchanged axes have unequal buffer widths.

**Step 2: Run the CLI tests and verify they fail**

```bash
python3 tests/dg/check_sawf_fragment_alignment_tool.py
```

Expected: FAIL because search/output behavior is absent.

**Step 3: Implement deterministic search**

Search a bounded rational translation lattice derived from the mesh and
operation denominators.  Rank valid candidates first by the number of
source fragments that every operation maps onto themselves, then by wrapped
Euclidean displacement, then lexicographically.  This deliberately selects
the C64 `a=(11/64)^3`, `q=(15,15,15)` placement, where each 16-point core is
closed about its own half-grid center, over the shorter
`a=(-5/64)^3`, `q=(31,31,31)` placement, which merely exchanges the two cores
on each axis.  Evaluate the full supplied operation set before ranking so
repeated runs choose the same globally compatible translation.

If no candidate exists, exit nonzero with a diagnostic naming the first failing gate and do not create output files.

**Step 4: Implement validated pair publication with rollback**

Write temporary files beside the destination and validate them by re-reading
before publishing either file.  Two independent filesystem paths cannot be
replaced atomically as one operation, so publish them sequentially and roll
back both destinations on every catchable failure or interruption.  Inputs
must tell users not to consume the pair while the tool is running; an
uncatchable process or machine failure remains outside this guarantee. Include
provenance with input paths, mesh, fragment counts, buffer widths, tolerance,
and translation.  Do not put a leading comment in the atom file because SALMON
reads a fixed number of atom rows directly; put provenance in the symmetry file
comments and stdout.

**Step 5: Run tests and commit**

```bash
python3 tests/dg/check_sawf_fragment_alignment_tool.py
git add tools/align_periodic_structure_to_fragments.py tests/dg/check_sawf_fragment_alignment_tool.py
git commit -m "feat: align periodic structures to fragment partitions"
```

### Task 4: Add SALMON-side whole-fragment runtime validation

**Files:**
- Modify: `src/gs/dc/lcfo_wannier_sawf_band.f90`
- Modify: `src/gs/dc/lcfo_flux.f90:6200-6330`
- Create: `tests/dg/check_sawf_fragment_symmetry_map.py`

**Step 1: Write a failing Fortran-driver test**

Compile a small driver against `lcfo_wannier_sawf_band.f90`. It must check:

```fortran
! identity: accepted
! W=-I, q=(15,15,15): accepted on 32^3 / 2x2x2
! W=-I, q=(4,4,4): rejected because every source core splits
! W=-I, q=(31,31,31): accepted as a one-to-one fragment permutation
! q=(0.5,0,0): rejected because it is not an integer index shift
```

The accepted inversion must report center coordinates `(7.5d0,7.5d0,7.5d0)` for diagnostics without requiring them to be integers.

**Step 2: Run the test and verify it fails**

```bash
python3 tests/dg/check_sawf_fragment_symmetry_map.py
```

Expected: FAIL because the validator is absent.

**Step 3: Implement the numerical helper**

Add a LAPACK-free helper that converts the normalized fractional operation to an integer index map within tolerance, enumerates each regular fragment core, validates the actual fragment buffer widths, and returns:

```fortran
logical :: grid_map_ok, fragment_map_ok
integer :: max_targets_per_source
integer :: source_to_target(nfrag)
real(8) :: max_grid_residual
```

Keep it independent of C/Si and inversion. A whole-fragment permutation is
valid even when source and target fragment numbers differ; reject splitting,
duplicate target ownership, unequal source/target core sizes, and buffer widths
that are incompatible with an exchanged axis. For non-inversion operations,
do not manufacture a center.

**Step 4: Call it before `D_band` construction**

After symmetry normalization and identity ordering, validate every operation before preparing the fragment state cache. Rank 0 logs:

```text
[DC-LCFO-SAWF-ALIGN] operation=2 grid_map_ok=T fragment_map_ok=T max_targets_per_source=1 max_grid_residual=...
```

Broadcast failure and stop before `.dmn` publication if any operation splits a core.

**Step 5: Run focused tests and build**

```bash
python3 tests/dg/check_sawf_fragment_symmetry_map.py
python3 tests/dg/check_sawf_band_representation.py
python3 tests/dg/check_sawf_mpi_band_accumulation.py
cmake --build build-mpi-eigenexa-wannier-lib -j 2
```

Expected: all PASS and `salmon` builds.

**Step 6: Commit the runtime guard**

```bash
git add src/gs/dc/lcfo_wannier_sawf_band.f90 src/gs/dc/lcfo_flux.f90 tests/dg/check_sawf_fragment_symmetry_map.py
git commit -m "feat: reject symmetry-split SAWF fragments"
```

### Task 5: Generate and lock the aligned C64 sample

**Files:**
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/atom_sawf_aligned.dat`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/sym_sawf_aligned.dat`
- Create: `samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_gs_w90_pseudo_sawf_aligned_nw576_nb664`
- Modify: `tests/dg/check_sawf_fragment_alignment_tool.py`
- Modify: `tests/dg/check_sawf_integration_log.py`

**Step 1: Add a failing repository-data regression**

Require the committed files to reproduce from `atom.dat` and `sym.dat`, with translation `(11/64)^3`, aligned inversion translation `(15/32)^3`, 64 same-species atom matches, and eight single-target fragment images.

**Step 2: Generate the dedicated files**

Run the tool with the real C64 cell, `32 32 32` mesh, and `2 2 2` fragments. Update the new input only:

```fortran
file_atom_coor = 'atom_sawf_aligned.dat'
wannier_symmetry_file = 'sym_sawf_aligned.dat'
wannier_site_symmetry = 'file'
wannier_symmetry_tolerance = 1.0d-6
```

Do not modify the historical `atom.dat`, `sym.dat`, or prior input.

**Step 3: Make the integration checker provenance-aware**

Teach `check_sawf_integration_log.py` to accept the symmetry path from the input/provenance and validate the aligned inversion rather than hardcoding run-root `sym.dat` and `tau=(1/8)^3`.

**Step 4: Run sample tests and commit**

```bash
python3 tests/dg/check_sawf_fragment_alignment_tool.py
python3 tests/dg/check_sawf_integration_log.py --self-test
git add samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/atom_sawf_aligned.dat \
  samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/sym_sawf_aligned.dat \
  samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac/inputfile_gs_w90_pseudo_sawf_aligned_nw576_nb664 \
  tests/dg/check_sawf_fragment_alignment_tool.py tests/dg/check_sawf_integration_log.py
git commit -m "test: add fragment-aligned C64 SAWF case"
```

### Task 6: Run the formal C64 SAWF integration gate

**Files:**
- Modify only if evidence requires it: `docs/plans/2026-07-10-sawf-dmn.md`
- Produce untracked run output under a new dedicated run directory.

**Step 1: Rebuild the production executable**

```bash
cmake --build build-mpi-eigenexa-wannier-lib -j 2
```

Expected: `[100%] Built target salmon`.

**Step 2: Launch with the established rank/thread layout**

Use MPI 8 and:

```bash
export OMP_NUM_THREADS=2
```

Run the aligned `576/664` input in a fresh directory. Do not reuse an earlier `.dmn`, `.amn`, checkpoint, or Wannier90 output.

**Step 3: Check the early structural gates**

Require before waiting for Wannier90:

- SCF reaches the established convergence criterion;
- identity and inversion log `fragment_map_ok=T` and `max_targets_per_source=1`;
- inversion does not show the prior `closure_residual ~= 0.99375` failure;
- no polar repair or post-hoc inversion symmetrization is used.

**Step 4: Check SAWF completion**

Run:

```bash
python3 tests/dg/check_sawf_integration_log.py RUN_LOG --expected-operations 2 --expected-bands 664
```

Expected: PASS, `.dmn` accepted by Wannier90, SAWF localization completed, current-run seed provenance present, and SALMON ends normally.

**Step 5: Record evidence and commit documentation**

Record exact SCF count, energy residual, fragment-map lines, singular-value/closure metrics, W90 completion, and output paths in `docs/plans/2026-07-10-sawf-dmn.md`.

```bash
git add docs/plans/2026-07-10-sawf-dmn.md
git commit -m "docs: record aligned C64 SAWF integration"
```

### Task 7: Verify field-off stationarity and HHG symmetry

**Files:**
- Create aligned RT inputs beside the C64 sample only after Task 6 passes.
- Reuse existing polarization and HHG analysis scripts after confirming their timestep normalization.

**Step 1: Run a short field-off DG propagation**

Use the aligned current-run SAWF seed. Check electron number, polarization drift, and DG eigenstate residual before applying a field. A nonstationary seed blocks the laser run.

**Step 2: Run the long laser pulse**

Use the established `dt=2` setup, MPI 8, and `OMP_NUM_THREADS=2`. Run through the end of the pulse and retain `Pz`; ignore the known-unreliable DG `Jz` route.

**Step 3: Compute and plot HHG**

Use polarization-derived response with the repository's verified FFT normalization. Plot `Pz`, the field/vector potential used by the input, and harmonic intensities. Quantify even-order peaks relative to adjacent odd orders rather than relying only on visual inspection.

**Step 4: Review before claiming the physical bug fixed**

Invoke `superpowers:requesting-code-review` and review symmetry-file provenance, half-grid index conventions, fragment mapping, absence of forced inversion projection, field-off residuals, and FFT normalization.

**Step 5: Run final verification**

Invoke `superpowers:verification-before-completion`, rerun all focused SAWF tests, the production build, the formal integration checker, and the field-off/HHG analysis commands. Only then report whether even-order HHG is physically suppressed.
